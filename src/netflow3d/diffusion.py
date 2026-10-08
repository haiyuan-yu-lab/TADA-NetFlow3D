"""
Heat diffusion, randomized background calibration, and module extraction.

Heat transfer is E = beta * (I - (1-beta)W)^-1 * diag(heat).
Columns are sources; rows are targets. Only transfers strictly above delta
form edges. The runner selects background settings and low/high resolution;
round_delta implements half-up threshold rounding with ties passing.
"""
from __future__ import annotations

import os
from collections.abc import Iterator, Sequence
from datetime import datetime
from decimal import ROUND_HALF_UP, Decimal
from multiprocessing import Pool
from os import PathLike
from typing import Any

import networkx as nx
import numpy as np
import pandas as pd
from scipy import linalg, sparse
from scipy.sparse import csgraph
from scipy.sparse.linalg import splu

DeltaResult = dict[str, Any]
SparseComponent = tuple[sparse.csc_matrix, list[str], np.ndarray]
DENSE_SOLVE_MAX_NODES = 20_000
BACKGROUND_NULLS = ("rewire", "heat", "both")
BACKGROUND_WEIGHT_MODES = {"rewire": ("unit", "shuffle"), "heat": ("unit", "keep"), "both": ("unit", "shuffle")}
HEAT_SHUFFLE_MODES = ("binned", "uniform")
HEAT_SHUFFLE_BINS = 10
_WORKER_GRAPH = None
_WORKER_RESTART_PROB = None
_WORKER_MAX_SUBNETWORK_SIZE = None
_WORKER_BACKGROUND = None

def log(msg: object, tag: str = "[NetFlow3D]", flush: bool = True) -> None:
    """Print a timestamped NetFlow3D log message."""
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{now}] [{tag}] {msg}", flush=flush)


def get_sparse_diffusion_sources(
    G: nx.Graph,
    beta: float,
) -> tuple[sparse.csc_matrix, list[str], np.ndarray]:
    """Compute sparse diffusion columns for positive-heat sources."""
    nodes = list(G.nodes())
    heat_scores = np.array([G.nodes[node]["heat_score"] for node in nodes], dtype=float)
    source_positions = np.flatnonzero(heat_scores > 0)

    if len(nodes) == 0 or len(source_positions) == 0:
        return sparse.csc_matrix((len(nodes), 0), dtype=float), nodes, source_positions

    adjacency = nx.to_scipy_sparse_array(G, nodelist=nodes, dtype=float, format="csc")
    col_sums = np.asarray(adjacency.sum(axis=0)).ravel()
    if np.any(col_sums == 0):
        raise ValueError("Sparse diffusion received a graph with isolated nodes.")

    W = adjacency @ sparse.diags(1.0 / col_sums, format="csc")
    system = sparse.eye(len(nodes), format="csc") - (1.0 - beta) * W

    rhs = sparse.csc_matrix(
        (
            beta * heat_scores[source_positions],
            (source_positions, np.arange(len(source_positions))),
        ),
        shape=(len(nodes), len(source_positions)),
    )
    if len(nodes) <= DENSE_SOLVE_MAX_NODES:
        # Use dense LAPACK for components up to the configured size threshold.
        solution = linalg.solve(
            system.toarray(), rhs.toarray(), overwrite_a=True, overwrite_b=True, check_finite=False
        )
    else:
        solution = splu(system).solve(rhs.toarray())
    return sparse.csc_matrix(solution), nodes, source_positions


def sparse_diffusion_to_adjacency(
    E_sources: sparse.spmatrix,
    source_positions: np.ndarray,
    delta: float | None = None,
) -> sparse.csr_matrix:
    """Convert source-column diffusion output to directed sparse adjacency."""
    if E_sources.shape[1] == 0:
        return sparse.csr_matrix((E_sources.shape[0], E_sources.shape[0]), dtype=bool)

    E = E_sources.tocsc()
    row_chunks = []
    col_chunks = []

    for source_col, source_node_idx in enumerate(source_positions):
        start = E.indptr[source_col]
        end = E.indptr[source_col + 1]
        values = E.data[start:end]
        targets = E.indices[start:end]

        if delta is not None:
            keep = values > delta
        else:
            keep = values != 0

        if not np.any(keep):
            continue

        kept_targets = targets[keep]
        row_chunks.append(np.full(len(kept_targets), source_node_idx, dtype=np.int32))
        col_chunks.append(kept_targets.astype(np.int32, copy=False))

    if len(row_chunks) == 0:
        return sparse.csr_matrix((E_sources.shape[0], E_sources.shape[0]), dtype=bool)

    rows = np.concatenate(row_chunks)
    cols = np.concatenate(col_chunks)
    data = np.ones(len(rows), dtype=bool)
    return sparse.csr_matrix((data, (rows, cols)), shape=(E_sources.shape[0], E_sources.shape[0]))


def iter_sparse_threshold_edges(
    E_sources: sparse.spmatrix,
    source_positions: np.ndarray,
    nodes: Sequence[str],
    delta: float | None = None,
) -> Iterator[tuple[str, str, dict[str, float]]]:
    """Yield thresholded source-to-target edges from sparse diffusion columns."""
    E = E_sources.tocsc()

    for source_col, source_node_idx in enumerate(source_positions):
        start = E.indptr[source_col]
        end = E.indptr[source_col + 1]
        values = E.data[start:end]
        targets = E.indices[start:end]

        if delta is not None:
            keep = values > delta
        else:
            keep = values != 0

        for target_idx, weight in zip(targets[keep], values[keep]):
            yield nodes[source_node_idx], nodes[target_idx], {"weight": float(weight)}


def sparse_strong_component_sets(
    adj: sparse.spmatrix,
    nodes: Sequence[str],
) -> list[set[str]]:
    """Return strongly connected node sets from source-row adjacency."""
    n_components, labels = csgraph.connected_components(
        adj,
        directed=True,
        connection="strong",
        return_labels=True,
    )
    components = []
    for label in range(n_components):
        idx = np.flatnonzero(labels == label)
        components.append({nodes[i] for i in idx})
    return components


def sparse_strong_component_sizes(adj: sparse.spmatrix) -> list[int]:
    """Return strongly connected component sizes from directed adjacency."""
    _, labels = csgraph.connected_components(
        adj,
        directed=True,
        connection="strong",
        return_labels=True,
    )
    return np.bincount(labels).astype(int).tolist()


def compute_sparse_diffusion_components(
    G: nx.Graph,
    beta: float,
) -> tuple[list[SparseComponent], np.ndarray]:
    """Compute sparse diffusion data for eligible connected components."""
    component_data = []
    heat_edge_chunks = []

    for component in nx.connected_components(G):
        component = list(component)
        heat_scores = [G.nodes[node]["heat_score"] for node in component]

        if len(component) <= 1 or max(heat_scores) <= 0:
            continue

        G_sub = G.subgraph(component)
        E_sources, nodes, source_positions = get_sparse_diffusion_sources(G_sub, beta=beta)
        if E_sources.nnz == 0:
            continue

        component_data.append((E_sources, nodes, source_positions))
        heat_edge_chunks.append(E_sources.data[E_sources.data != 0])

    if len(heat_edge_chunks) == 0:
        heat_edges = np.array([], dtype=float)
    else:
        heat_edges = np.sort(np.concatenate(heat_edge_chunks))[::-1]
    return component_data, heat_edges


def get_module_sizes_from_sparse_components(
    component_data: Sequence[SparseComponent],
    delta: float | None = None,
) -> list[int]:
    """Collect SCC sizes after thresholding sparse source-column components."""
    module_sizes = []

    for E_sources, nodes, source_positions in component_data:
        adj = sparse_diffusion_to_adjacency(E_sources, source_positions, delta=delta)
        module_sizes.extend(sparse_strong_component_sizes(adj))

    return module_sizes


def find_suitable_delta_from_sparse_components(
    component_data: Sequence[SparseComponent],
    heat_edges: Sequence[float],
    size_cutoff: int,
    upper_bound: float = 1.0,
    lower_bound: float = 0.0,
) -> DeltaResult:
    """Select a sparse-backend threshold by binary search over heat transfers.

    heat_edges is sorted from largest to smallest, and delta = heat_edges[k]
    keeps the k transfers strictly above it. The search ends with lower_bound
    at the largest tested fraction whose modules are all at most size_cutoff,
    and the returned delta is that tested value. Module sizes are recomputed
    at the returned threshold.
    """
    if len(component_data) == 0 or len(heat_edges) == 0:
        return {
            "delta": -np.inf,
            "top_edge_cutoff": np.nan,
            "subnetwork_sizes": [],
        }

    module_sizes = get_module_sizes_from_sparse_components(component_data, delta=None)

    if len(module_sizes) == 0:
        log("no modules identified without thresholding; returning -inf", tag="delta", flush=True)
        return {
            "delta": -np.inf,
            "top_edge_cutoff": np.nan,
            "subnetwork_sizes": [],
        }

    if max(module_sizes) <= size_cutoff:
        log(f"no thresholding needed; max module size={max(module_sizes)} <= cutoff={size_cutoff}", tag="delta", flush=True)
        return {
            "delta": -1.0,
            "top_edge_cutoff": 1.0,
            "subnetwork_sizes": module_sizes,
        }

    diff = max(2 / len(heat_edges), 1e-5)
    log(f"starting binary search with diff threshold={diff}", tag="delta", flush=True)

    while upper_bound - lower_bound > diff:
        midpoint = (upper_bound + lower_bound) / 2
        idx = min(int(midpoint * len(heat_edges)), len(heat_edges) - 1)
        delta = heat_edges[idx]

        module_sizes = get_module_sizes_from_sparse_components(component_data, delta=delta)
        max_size = max(module_sizes) if len(module_sizes) > 0 else 0

        if len(module_sizes) == 0 or max_size <= size_cutoff:
            lower_bound = midpoint
        else:
            upper_bound = midpoint

    final_midpoint = lower_bound
    final_idx = min(int(final_midpoint * len(heat_edges)), len(heat_edges) - 1)
    final_delta = heat_edges[final_idx]
    final_sizes = get_module_sizes_from_sparse_components(component_data, delta=final_delta)

    return {
        "delta": final_delta,
        "top_edge_cutoff": final_midpoint,
        "subnetwork_sizes": final_sizes,
    }


def get_suitable_delta_sparse(
    G: nx.Graph,
    beta: float,
    size_cutoff: int,
) -> DeltaResult:
    """Compute sparse diffusion components and select their delta threshold."""
    component_data, heat_edges = compute_sparse_diffusion_components(G, beta=beta)
    return find_suitable_delta_from_sparse_components(
        component_data=component_data,
        heat_edges=heat_edges,
        size_cutoff=size_cutoff,
    )


def background_settings(
    background_null: str = "rewire",
    swaps_per_edge: int = 10,
    background_weights: str = "unit",
    heat_shuffle: str = "binned",
) -> dict[str, Any]:
    """Validate background settings and return the values recorded per trial.

    Settings that do not apply to the chosen null are recorded as
    swaps_per_edge=0 (heat null) or heat_shuffle="none" (rewire null).
    """
    if background_null not in BACKGROUND_NULLS:
        raise ValueError(f"background_null must be one of {BACKGROUND_NULLS}")
    allowed = BACKGROUND_WEIGHT_MODES[background_null]
    if background_weights not in allowed:
        raise ValueError(f"background_weights for the {background_null} null must be one of {allowed}")
    if background_null != "heat" and swaps_per_edge < 1:
        raise ValueError("swaps_per_edge must be at least 1")
    if background_null != "rewire" and heat_shuffle not in HEAT_SHUFFLE_MODES:
        raise ValueError(f"heat_shuffle must be one of {HEAT_SHUFFLE_MODES}")
    return {
        "background_null": background_null,
        "swaps_per_edge": 0 if background_null == "heat" else int(swaps_per_edge),
        "background_weights": background_weights,
        "heat_shuffle": "none" if background_null == "rewire" else heat_shuffle,
    }


def shuffle_background_heat(
    G_original: nx.Graph,
    seed: int | Sequence[int],
    heat_shuffle: str,
    background_weights: str = "keep",
) -> nx.Graph:
    """Return a copy with heat scores shuffled across proteins and edges unchanged.

    "binned" ranks proteins by number of distinct partners (ties in node order)
    and shuffles heat within HEAT_SHUFFLE_BINS equal-size bins, so hubs keep
    hub-like heat; "uniform" shuffles across all proteins. Zero heat scores are
    shuffled too. background_weights="keep" leaves every edge weight in place;
    "unit" sets them to one.
    """
    if heat_shuffle not in HEAT_SHUFFLE_MODES:
        raise ValueError(f"heat_shuffle must be one of {HEAT_SHUFFLE_MODES}")
    if background_weights not in BACKGROUND_WEIGHT_MODES["heat"]:
        raise ValueError(f"background_weights must be one of {BACKGROUND_WEIGHT_MODES['heat']}")
    G = G_original.copy()
    nodes = list(G.nodes())
    heat = np.array([G.nodes[node]["heat_score"] for node in nodes], dtype=float)
    if heat_shuffle == "uniform":
        groups = [np.arange(len(nodes))]
    else:
        partners = np.array([len(set(G[node]) - {node}) for node in nodes])
        groups = np.array_split(np.argsort(partners, kind="stable"), min(HEAT_SHUFFLE_BINS, len(nodes)))

    rng = np.random.default_rng(seed)
    shuffled = heat.copy()
    for idx in groups:
        shuffled[idx] = heat[rng.permutation(idx)]
    for node, value in zip(nodes, shuffled):
        G.nodes[node]["heat_score"] = float(value)

    if background_weights == "unit":
        nx.set_edge_attributes(G, 1.0, "weight")
    return G


def make_background_graph(
    G_original: nx.Graph,
    seed: int,
    settings: dict[str, Any],
) -> nx.Graph:
    """Build one background graph from settings returned by background_settings.

    The combined null shuffles heat with its own random stream ([seed, 1]), so
    the heat order is independent of the edge-weight order drawn from seed.
    """
    if settings["background_null"] == "heat":
        return shuffle_background_heat(
            G_original, seed, settings["heat_shuffle"], settings["background_weights"]
        )
    G = rewire_background_graph(
        G_original, seed, settings["swaps_per_edge"], settings["background_weights"]
    )
    if settings["background_null"] == "both":
        G = shuffle_background_heat(G, [seed, 1], settings["heat_shuffle"], "keep")
    return G


def rewire_background_graph(
    G_original: nx.Graph,
    seed: int,
    swaps_per_edge: int,
    background_weights: str = "unit",
) -> nx.Graph:
    """Return a degree-preserving rewired copy of the graph.

    Self-loops (homodimers) are held out of the swaps and restored afterwards.
    networkx would otherwise turn a self-loop into two partner edges, or merge
    two self-loops into one edge and lower both degrees. Only edges between
    different proteins are rewired, and nswap counts only those edges.

    background_weights="unit" sets every edge weight to one. With "shuffle",
    rewired edges receive the original between-protein weights in random order, and
    self-loops keep their own weights, so the number of upweighted edges is
    unchanged.

    networkx counts only attempts that reach its parallel-edge check, which a
    sparse graph rarely rejects, so ten tries per swap leaves ample headroom.
    """
    if background_weights not in BACKGROUND_WEIGHT_MODES["rewire"]:
        raise ValueError(f"background_weights must be one of {BACKGROUND_WEIGHT_MODES['rewire']}")
    G = G_original.copy()
    self_loops = list(nx.selfloop_edges(G, data=True))
    G.remove_edges_from(self_loops)
    weights = [w for _, _, w in G.edges(data="weight", default=1.0)]
    nswap = swaps_per_edge * G.number_of_edges()

    nx.double_edge_swap(
        G,
        nswap=nswap,
        max_tries=10 * nswap,
        seed=seed,
    )

    if background_weights == "shuffle":
        order = np.random.default_rng(seed).permutation(len(weights))
        for (u, v), i in zip(G.edges(), order):
            G[u][v]["weight"] = weights[i]
        G.add_edges_from(self_loops)
    else:
        G.add_edges_from(self_loops)
        nx.set_edge_attributes(G, 1.0, "weight")
    return G


def estimate_one_bg_delta(
    args: tuple[nx.Graph, float, int, dict[str, Any], int],
) -> DeltaResult:
    """Estimate one delta on a randomized background graph."""
    G_original, restart_prob, max_subnetwork_size, settings, seed = args
    log(f"seed={seed}: starting background delta estimation", tag="worker", flush=True)

    G = make_background_graph(G_original, seed, settings)

    result = get_suitable_delta_sparse(
        G,
        beta=restart_prob,
        size_cutoff=max_subnetwork_size,
    )
    result["seed"] = seed
    log(
        f"seed={seed}: done, delta={result['delta']}, "
        f"max_module_size={max(result['subnetwork_sizes']) if result['subnetwork_sizes'] else 'NA'}",
        tag="worker",
        flush=True,
    )
    return result


def init_bg_worker(
    G: nx.Graph,
    restart_prob: float,
    max_subnetwork_size: int,
    settings: dict[str, Any],
) -> None:
    """Initialize per-process globals shared by background-trial workers."""
    global _WORKER_GRAPH
    global _WORKER_RESTART_PROB
    global _WORKER_MAX_SUBNETWORK_SIZE
    global _WORKER_BACKGROUND

    _WORKER_GRAPH = G
    _WORKER_RESTART_PROB = restart_prob
    _WORKER_MAX_SUBNETWORK_SIZE = max_subnetwork_size
    _WORKER_BACKGROUND = settings


def estimate_one_bg_delta_from_seed(seed: int) -> DeltaResult:
    """Run one background trial using the initialized worker globals."""
    return estimate_one_bg_delta((
        _WORKER_GRAPH,
        _WORKER_RESTART_PROB,
        _WORKER_MAX_SUBNETWORK_SIZE,
        _WORKER_BACKGROUND,
        seed,
    ))


def estimate_bg_delta(
    G: nx.Graph,
    output_dir: str | PathLike[str],
    restart_prob: float,
    max_subnetwork_size: int,
    delta_trial: int,
    random_seed: int,
    n_workers: int,
    resolution: str,
    resolution_quantile: float | None = None,
    swaps_per_edge: int = 10,
    background_weights: str = "unit",
    background_null: str = "rewire",
    heat_shuffle: str = "binned",
) -> float:
    """Estimate a background delta from randomized background trials."""
    settings = background_settings(background_null, swaps_per_edge, background_weights, heat_shuffle)
    log("Estimating background delta")
    bg_delta_file = os.path.join(output_dir, "bg_delta.txt")

    if G.number_of_edges() < max_subnetwork_size:
        log(
            f"Graph has only {G.number_of_edges()} edges, smaller than "
            f"max_subnetwork_size={max_subnetwork_size}; returning -inf"
        )
        return -np.inf

    seeds = [random_seed + i for i in range(delta_trial)]
    log(
        f"Launching {delta_trial} background trials with {n_workers} workers using sparse backend "
        f"({', '.join(f'{k}={v}' for k, v in settings.items())})"
    )

    results = []
    if n_workers <= 1:
        for i, seed in enumerate(seeds, start=1):
            result = estimate_one_bg_delta((G, restart_prob, max_subnetwork_size, settings, seed))
            results.append(result)
            log(f"Background trial {i}/{delta_trial} completed", tag="delta")
    else:
        with Pool(
            n_workers,
            initializer=init_bg_worker,
            initargs=(G, restart_prob, max_subnetwork_size, settings),
        ) as pool:
            for i, result in enumerate(pool.imap_unordered(estimate_one_bg_delta_from_seed, seeds), start=1):
                results.append(result)
                log(f"Background trial {i}/{delta_trial} completed", tag="delta")

    log("All background delta trials completed")
    delta_df = pd.DataFrame({
        "seed": [r["seed"] for r in results],
        "delta": [r["delta"] for r in results],
        "top_edge_cutoff": [r["top_edge_cutoff"] for r in results],
        "subnetwork_sizes": [
            ",".join(map(str, r["subnetwork_sizes"])) for r in results
        ],
        **settings,
    })
    delta_df.to_csv(bg_delta_file, sep="\t", index=False)

    return select_background_delta(bg_delta_file, resolution, resolution_quantile)


def check_background_settings(
    background_delta_file: str | PathLike[str],
    settings: dict[str, Any],
) -> None:
    """Refuse to reuse a background-delta table built with other settings.

    settings comes from background_settings. A table missing one of its
    columns predates that setting (built with the rewire null, one swap per
    edge, unit weights and no heat shuffling) and is accepted with a log message.
    """
    delta_df = pd.read_csv(background_delta_file, sep="\t")
    legacy_values = {"background_null": "rewire", "swaps_per_edge": 1,
                     "background_weights": "unit", "heat_shuffle": "none"}
    for column, requested in settings.items():
        legacy = legacy_values[column]
        if column not in delta_df.columns:
            log(
                f"{background_delta_file} has no {column} column; it predates the setting "
                f"and was built with {column}={legacy}"
            )
            continue
        recorded = set(delta_df[column].astype(type(requested)))
        if recorded != {requested}:
            raise ValueError(
                f"{background_delta_file} was built with {column} {sorted(recorded)}, "
                f"but {column}={requested} was requested"
            )


def select_background_delta(
    background_delta_file: str | PathLike[str],
    resolution: str,
    resolution_quantile: float | None = None,
) -> float:
    """Load a background-delta table and select its requested resolution.

    low: minimum trial delta. high: median. quantile: the trial delta at
    'resolution_quantile' with lower interpolation, so the result is always an
    observed trial delta (quantile 0 equals low).
    """
    delta_df = pd.read_csv(background_delta_file, sep="\t")
    if "delta" not in delta_df.columns:
        raise ValueError("Background delta table must contain a 'delta' column")

    valid_deltas = delta_df["delta"].replace([np.inf, -np.inf], np.nan).dropna()
    if len(valid_deltas) == 0:
        return -np.inf

    if resolution == "low":
        return valid_deltas.min()

    if resolution == "high":
        return valid_deltas.median()

    if resolution == "quantile":
        if resolution_quantile is None or not 0.0 <= resolution_quantile <= 1.0:
            raise ValueError("resolution='quantile' needs a resolution_quantile in [0, 1]")
        return float(np.quantile(valid_deltas.to_numpy(), resolution_quantile, method="lower"))

    raise ValueError(f"Unknown resolution: {resolution}")


def round_delta(delta: float, decimals: int) -> float:
    """Return the edge cutoff for delta rounded to `decimals` places, ties passing.

    delta is rounded half up, and an edge is kept when its transfer, rounded the
    same way, is at least the rounded delta, i.e. transfer >= rounded delta -
    0.5 * 10**-decimals. Edges are kept when strictly above the returned value, so
    it sits one floating-point step below that boundary. Non-finite deltas are
    returned unchanged.
    """
    if decimals < 0:
        raise ValueError("decimals must be non-negative")
    if not np.isfinite(delta):
        return float(delta)
    step = Decimal(1).scaleb(-decimals)
    rounded = Decimal(repr(float(delta))).quantize(step, rounding=ROUND_HALF_UP)
    return float(np.nextafter(float(rounded - step / 2), -np.inf))


def identify_hot_modules_sparse(
    G: nx.Graph,
    graph_output: str | PathLike[str],
    beta: float,
    delta: float,
) -> list[set[str]]:
    """Extract sparse-backend hot SCCs and write their directed diffusion graph."""
    all_modules = []
    G_final = nx.DiGraph()

    for component in nx.connected_components(G):
        component = list(component)
        heat_scores = [G.nodes[item]["heat_score"] for item in component]

        if max(heat_scores) <= 0:
            continue

        if len(component) == 1:
            all_modules.append(set(component))
            continue

        G_sub = G.subgraph(component)
        E_sources, nodes, source_positions = get_sparse_diffusion_sources(G_sub, beta=beta)
        adj = sparse_diffusion_to_adjacency(E_sources, source_positions, delta=delta)

        G_final.add_nodes_from(nodes)
        G_final.add_edges_from(iter_sparse_threshold_edges(E_sources, source_positions, nodes, delta=delta))

        for module in sparse_strong_component_sets(adj, nodes):
            if len(module) > 1 or G.nodes[next(iter(module))]["heat_score"] > 0:
                all_modules.append(module)

    nx.write_graphml(G_final, graph_output)
    return all_modules
