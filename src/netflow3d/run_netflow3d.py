"""
Find heat-diffusion subnetworks from prepared gene scores and network inputs.
"""
from __future__ import annotations

import argparse
import os
import shutil
from argparse import Namespace
from collections.abc import Sequence
from os import PathLike
from pathlib import Path

import networkx as nx
import pandas as pd
from diffusion import (
    BACKGROUND_NULLS, HEAT_SHUFFLE_MODES, background_settings, check_background_settings,
    estimate_bg_delta, identify_hot_modules_sparse, log, round_delta, select_background_delta,
)

def args_parser() -> Namespace:
    """Parse network inputs, analysis settings, and the output destination."""
    parser = argparse.ArgumentParser(description="NetFlow3D", allow_abbrev=False,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-j', '--job-name', required=True, help="Job directory and subnetwork filename prefix")
    parser.add_argument('-g', '--gene-score', required=True, help="Gene/Score TSV")
    parser.add_argument('-o', '--output-dir', required=True, help="Parent directory; results go into <output-dir>/<job-name>")
    parser.add_argument('-i', '--disrupted-interaction', help="P1/P2 TSV")
    reference_dir = Path(__file__).resolve().parents[2] / "resources/netflow3d"
    parser.add_argument('-b', '--binary-interactome',
                        default=str(reference_dir / "HomoSapiens_binary_HINThq.txt"))
    parser.add_argument('-p', '--id-mapping',
                        default=str(reference_dir / "HUMAN_9606_idmapping.dat.gz"))
    parser.add_argument('-X', '--expressed-genes', help="Expressed-gene list; omit to disable")
    parser.add_argument('-R', '--resolution', choices=['low', 'high'], default='low',
                        help="Minimum (low) or median (high) of finite background deltas")
    parser.add_argument('-r', '--restart-prob', type=float, default=0.5)
    parser.add_argument('-w', '--weighted-edge', type=float, default=5.0)
    parser.add_argument('-x', '--max-subnetwork-size', type=int, default=5,
                        help="Module-size cutoff for background calibration")
    parser.add_argument('-D', '--delta-trial', type=int, default=1000, help="Number of background trials")
    parser.add_argument('--swaps-per-edge', type=int, default=10,
                        help="Degree-preserving swaps per edge in each background trial")
    parser.add_argument('--background-null', choices=BACKGROUND_NULLS, default='both',
                        help="Randomize edges keeping heat (rewire), shuffle heat across "
                             "proteins keeping edges (heat), or do both (both)")
    parser.add_argument('--background-weights', choices=['unit', 'shuffle', 'keep'], default='unit',
                        help="Background edge weights: all one (unit); rewire or both null: "
                             "original weights shuffled across rewired edges (shuffle); "
                             "heat null: original weights left in place (keep)")
    parser.add_argument('--heat-shuffle', choices=HEAT_SHUFFLE_MODES, default='uniform',
                        help="Heat or both null: shuffle within degree bins (binned) or "
                             "across all proteins (uniform)")
    parser.add_argument('-S', '--random-seed', type=int, default=4001, help="First background trial seed")
    parser.add_argument('-t', '--n-workers', type=int, default=8, help="Background worker processes")
    parser.add_argument('--background-delta-file', help="Reuse trials instead of recalculating")
    parser.add_argument('--delta-decimals', type=int, default=3,
                        help="Round the selected delta to this many decimals (half up) and keep "
                             "transfers that round to at least it")
    return parser.parse_args()


def load_id_mapping(id_mapping_file: str | PathLike[str]) -> pd.DataFrame:
    """Load a UniProt-to-gene-tag table used to map ``Gene/Score`` inputs."""
    df_annot = pd.read_csv(id_mapping_file, sep = "\t", dtype = str, header = None)
    df_annot[0] = df_annot[0].apply(lambda x: x.split("-")[0])
    df_annot[2] = df_annot[2].apply(lambda x: x.split(".")[0])
    df_annot = df_annot[[0, 2]].drop_duplicates().rename(columns = {0: "UniProt", 2: "tag"})
    return df_annot


def build_interactome(binary_interactome_file: str | PathLike[str]) -> nx.Graph:
    """Build an initialized undirected protein-protein interaction graph."""
    df_interactome = pd.read_csv(binary_interactome_file, sep="\t", header=None)

    edges = set()
    for _, row in df_interactome.iterrows():
        pro1 = str(row[0]).split("-")[0]
        pro2 = str(row[1]).split("-")[0]
        edges.add(tuple(sorted([pro1, pro2])))

    G = nx.Graph(sorted(edges))

    for node in G.nodes():
        G.nodes[node]["heat_score"] = 0.0

    for u, v in G.edges():
        G[u][v]["weight"] = 1.0

    return G


def filter_expressed_genes(
    G: nx.Graph,
    df_annot: pd.DataFrame,
    expressed_genes_file: str | PathLike[str] | None,
) -> nx.Graph:
    """Filter a protein graph to proteins mapped from an expressed-gene list."""
    if expressed_genes_file is None:
        return G

    expr_genes = set(
        pd.read_csv(expressed_genes_file, sep="\t", header=None)[0]
        .astype(str)
        .str.split(".")
        .str[0]
    )

    expr_uniprots = set(df_annot[df_annot["tag"].isin(expr_genes)]["UniProt"])
    return nx.Graph(G.subgraph(expr_uniprots))


def assign_heat_scores(
    G: nx.Graph,
    df_annot: pd.DataFrame,
    gene_score_file: str | PathLike[str],
) -> nx.Graph:
    """Map ``Gene/Score`` rows to proteins and assign their maximum scores."""
    gene_score = pd.read_csv(gene_score_file, sep="\t", header=0)
    gene_score_annot = pd.merge(gene_score, df_annot, left_on="Gene", right_on="tag", how="inner")

    for uniprot, df_one in gene_score_annot.groupby("UniProt"):
        if uniprot in G:
            G.nodes[uniprot]["heat_score"] = float(df_one["Score"].max())

    return G


def assign_disrupted_edge_weights(
    G: nx.Graph,
    disrupted_interaction_file: str | PathLike[str] | None,
    weighted_edge: float,
) -> nx.Graph:
    """Upweight eligible disrupted interactions and record diagnostic counts."""
    if disrupted_interaction_file is None:
        G.graph["disrupted_edge_stats"] = {
            "input_unique": 0,
            "in_interactome": 0,
            "zero_heat_skipped": 0,
            "weighted": 0,
        }
        return G

    df_disppi = pd.read_csv(disrupted_interaction_file, sep="\t", header=0)
    required_columns = {"P1", "P2"}
    if not required_columns.issubset(df_disppi.columns):
        raise ValueError("Disrupted interaction file must contain P1 and P2 columns")

    edges = set()
    for p1, p2 in df_disppi[["P1", "P2"]].dropna().itertuples(index=False, name=None):
        p1 = str(p1).split("-")[0]
        p2 = str(p2).split("-")[0]
        edges.add(tuple(sorted((p1, p2))))

    in_interactome = 0
    zero_heat_skipped = 0
    weighted = 0
    for p1, p2 in sorted(edges):
        if not G.has_edge(p1, p2):
            continue
        in_interactome += 1
        if G.nodes[p1]["heat_score"] <= 0 or G.nodes[p2]["heat_score"] <= 0:
            zero_heat_skipped += 1
            continue
        G[p1][p2]["weight"] = weighted_edge
        weighted += 1

    G.graph["disrupted_edge_stats"] = {
        "input_unique": len(edges),
        "in_interactome": in_interactome,
        "zero_heat_skipped": zero_heat_skipped,
        "weighted": weighted,
    }

    return G


def write_initial_graphml(
    G: nx.Graph,
    graph_output: str | PathLike[str],
) -> None:
    """Write GraphML with disrupted-edge counts temporarily scalarized."""
    stats = G.graph.pop("disrupted_edge_stats")
    graphml_stats = {
        f"disrupted_edges_{key}": value
        for key, value in stats.items()
    }
    G.graph.update(graphml_stats)
    try:
        nx.write_graphml(G, graph_output)
    finally:
        for key in graphml_stats:
            G.graph.pop(key)
        G.graph["disrupted_edge_stats"] = stats


def save_subnetworks(
    all_subnetworks: Sequence[set[str]],
    output_dir: str | PathLike[str],
    job_name: str,
) -> str:
    """Write subnetworks in descending size order."""
    out_file = os.path.join(output_dir, f"{job_name}_subnetworks.txt")
    with open(out_file, "w") as f:
        f.write("\t".join(["Subnetwork_genes", "Subnetwork_size"]) + "\n")
        for subnet in sorted(all_subnetworks, key=lambda z: len(z), reverse=True):
            f.write(",".join(sorted(subnet)) + "\t" + str(len(subnet)) + "\n")
    return out_file


def main(args: Namespace) -> None:
    """Run graph preparation, threshold selection, and module extraction."""
    background = background_settings(
        args.background_null, args.swaps_per_edge, args.background_weights, args.heat_shuffle
    )
    if args.delta_decimals is not None and args.delta_decimals < 0:
        raise ValueError("--delta-decimals must be non-negative")
    if args.job_name in {"", ".", ".."} or os.path.basename(args.job_name) != args.job_name:
        raise ValueError("--job-name must be a single directory name")
    output_dir = os.path.join(args.output_dir, args.job_name)
    os.makedirs(output_dir, exist_ok=True)
    assert os.path.isfile(args.binary_interactome), f"No such file or directory: {args.binary_interactome}"
    if args.expressed_genes != None:
        assert os.path.isfile(args.expressed_genes), f"No such file or directory: {args.expressed_genes}"

    df_annot = load_id_mapping(args.id_mapping)
    G = build_interactome(args.binary_interactome)
    G = filter_expressed_genes(G, df_annot, args.expressed_genes)
    G = assign_heat_scores(G, df_annot, args.gene_score)
    G = assign_disrupted_edge_weights(G, args.disrupted_interaction, args.weighted_edge)
    stats = G.graph["disrupted_edge_stats"]
    log(
        "Disrupted edges: "
        f"input_unique={stats['input_unique']}, "
        f"in_interactome={stats['in_interactome']}, "
        f"zero_heat_skipped={stats['zero_heat_skipped']}, "
        f"weighted={stats['weighted']}"
    )
    if args.disrupted_interaction is not None and stats["weighted"] == 0:
        log("WARNING: no disrupted interaction received an increased edge weight.")

    initial_graph_file = os.path.join(output_dir, "initial_state.graphml.gz")
    write_initial_graphml(G, initial_graph_file)

    if args.background_delta_file is None:
        if G.number_of_edges() < args.max_subnetwork_size:
            previous_background = os.path.join(output_dir, "bg_delta.txt")
            if os.path.isfile(previous_background):
                os.remove(previous_background)
            log("Background trials skipped: graph is smaller than the calibration cutoff; no bg_delta.txt is written.")
        delta = estimate_bg_delta(
            G=G,
            output_dir=output_dir,
            restart_prob=args.restart_prob,
            max_subnetwork_size=args.max_subnetwork_size,
            delta_trial=args.delta_trial,
            random_seed=args.random_seed,
            n_workers=args.n_workers,
            resolution=args.resolution,
            swaps_per_edge=args.swaps_per_edge,
            background_weights=args.background_weights,
            background_null=args.background_null,
            heat_shuffle=args.heat_shuffle,
        )
    else:
        assert os.path.isfile(args.background_delta_file), \
            f"No such file or directory: {args.background_delta_file}"
        check_background_settings(args.background_delta_file, background)
        delta = select_background_delta(args.background_delta_file, args.resolution)
        reused_delta_output = os.path.join(output_dir, "bg_delta.txt")
        if os.path.abspath(args.background_delta_file) != os.path.abspath(reused_delta_output):
            shutil.copyfile(args.background_delta_file, reused_delta_output)
        log(
            f"Reusing background delta table {args.background_delta_file} "
            f"with resolution={args.resolution}; background trials skipped."
        )
    log(f"The chosen delta is {delta} (resolution={args.resolution}, trials={args.delta_trial}).")
    if args.delta_decimals is not None:
        delta = round_delta(delta, args.delta_decimals)
        half_step = 0.5 * 10 ** -args.delta_decimals
        log(
            f"Delta rounded to {args.delta_decimals} decimals with ties passing: "
            f"{delta + half_step:.{args.delta_decimals}f}; edges need transfer >= "
            f"{delta:.{args.delta_decimals + 1}f}."
        )

    final_graph_file = os.path.join(output_dir, "final_state.graphml.gz")
    all_subnetworks = identify_hot_modules_sparse(
        G,
        final_graph_file,
        beta=args.restart_prob,
        delta=delta,
    )

    subnet_file = save_subnetworks(all_subnetworks, output_dir, args.job_name)

    if len(all_subnetworks) > 0:
        sizes = [len(x) for x in all_subnetworks]
        log(
            f"{len(all_subnetworks)} subnetworks are identified, "
            f"whose sizes range from {min(sizes)} to {max(sizes)}."
        )
    else:
        log("No subnetworks were identified.")

    log(f"Subnetwork results saved to {subnet_file}")
    log("Job finished!")

if __name__ == '__main__':
    main(args_parser())
