# Input formats

All tables are tab-separated with a header unless stated otherwise. Gene IDs
must use one consistent namespace and omit version suffixes where possible.

## TADA gene information

Required metadata: `gene_id`, `LOEUF`, `pLI`, and a gene-symbol column.
The runner resolves `Hugo_Symbol`, then `gene`, then `gene_gencodeV33`.
Noncoding-only runs use this table as the output gene backbone.

## TADA coding table

One row per gene, with gene metadata. Fields for selected components cover
mutation rates, priors, de novo case and
sibling counts, and transmitted/untransmitted counts for `ptv`, `misb`, and
`misa`. The canonical names are demonstrated in
`examples/input/tada/SSC_KOR_STAR/coding_counts.tsv`.

## TADA noncoding table

One row per regulatory element. The mutation-rate column is selected explicitly with `--mut-col` (`mut.em`
in the bundled examples). Component columns include `gene_id`,
`dn.em`, `dn.em.sib`, `in.t.em`, and `in.u.em`. Multiple target genes are
comma-separated in `gene_id`.

## Priors and sample sizes

Two-column tables named `variable` and `value`. Required variables are shown in
the example files. Sample sizes count analyzed individuals, not variants.

## NetFlow3D gene scores

Columns: `Gene` and `Score`. `Gene` may contain a gene symbol or supported gene
identifier present in the ID mapping. Supply finite numeric scores; positive
scores initiate diffusion. If several input rows map to the same protein,
NetFlow3D uses the largest score. Network proteins without a matching score
start with zero heat and remain available to receive heat from other proteins.

## Binary interactome

Headerless two-column table of UniProt accessions, one interaction per row.

## UniProt ID mapping

Headerless three-column UniProt mapping table: UniProt accession, identifier
type, and identifier. The runner reads columns one and three, removes UniProt
isoform suffixes from the first column and version suffixes from the third.
A gzip-compressed mapping table is also accepted.

## Disrupted interactions

Optional table with `P1` and `P2` headers containing UniProt accessions.
Isoform suffixes are removed and pairs deduplicated. Only existing network edges
whose endpoints both have positive heat receive the specified edge weight.

## Expressed genes

Optional headerless, one-column list. Values may be gene symbols, Ensembl gene
IDs, or other identifiers included in the ID mapping.
