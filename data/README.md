# Data files (not in git — download from GEO)

Place these files in this `data/` directory before knitting the notebooks:

| File | Source |
|------|--------|
| `GSE118189_ATAC_counts.txt.gz` | [GSE118189](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118189) — consensus peak count matrix |
| `GSE118165_raw_counts_GRCh38.p13_NCBI.tsv.gz` (or `.tsv`) | [GSE118165](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118165) — raw gene counts (NCBI GeneID) |
| `Supplementary_tables.xlsx` | Paper supplementary tables (sheet **ATAC-seq samle QC** used for sample metadata and TSS enrichment) |

Already included in the repository:

- `GSE118165_GSM_to_sample_id.tsv` — maps GEO GSM column names to `sample_id` labels used in the ATAC matrix.

Optional (for full Stanford atlas matrices): [immune atlas data archive](https://web.stanford.edu/group/pritchardlab/dataArchive/immune_atlas_web/).
