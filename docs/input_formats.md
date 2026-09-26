# Input formats

All tables are tab-separated with a header unless stated otherwise. Gene IDs
must use one consistent namespace and omit version suffixes where possible.

## TADA gene information

Required columns: `gene_id`, `gene`, `LOEUF`, and `pLI`.

## TADA coding table

One row per gene. Required fields cover mutation rates, priors, de novo case and
sibling counts, and transmitted/untransmitted counts for `ptv`, `misb`, and
`misa`. The canonical names are demonstrated in
`examples/input/coding_counts.tsv`.

## TADA noncoding table

One row per regulatory element. Required columns are `gene_id`, `mut.em`,
`dn.em`, `dn.em.sib`, `in.t.em`, and `in.u.em`. Multiple target genes are
comma-separated in `gene_id`.

## Priors and sample sizes

Two-column tables named `variable` and `value`. Required variables are shown in
the example files. Sample sizes count analyzed individuals, not variants.

## NetFlow3D gene scores

Columns: `Gene` and `Score`. `Gene` may contain a gene symbol or supported gene
identifier present in the ID mapping. `Score` must be finite and numeric.

## Binary interactome

Headerless two-column table of UniProt accessions, one interaction per row.

## UniProt ID mapping

Headerless three-column UniProt mapping table: UniProt accession, identifier
type, and identifier. The provisional implementation reads columns one and
three and strips isoform/version suffixes.

## Disrupted interactions

Optional table with `P1` and `P2` headers containing UniProt accessions.

## Expressed genes

Optional headerless, one-column list. Values may be gene symbols, Ensembl gene
IDs, or other identifiers included in the ID mapping.
