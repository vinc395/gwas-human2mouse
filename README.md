# gwas-human2mouse

Pipeline for filtering selected GWAS Catalog studies, preserving local provenance, standardizing human coordinates, and generating conservative human-to-mouse mm39 liftOver outputs.

## Included in this repository
- `scripts/` pipeline and helper scripts
- `manifests/` selected study manifests
- `reports/` compact run reports and summaries
- `tests/` regression test for large TSV CSV field handling

## Excluded from version control
Large raw and generated artifacts are intentionally excluded:
- `downloads/` raw GCST downloads and liftOver assets
- `tmp/` temporary files
- `parsed/` large parsed association tables
- `human_bed/` generated human BED/detail outputs
- `mm39_bed/` generated mm39 BED/detail outputs
- `reports/unmapped_human_features.tsv` large generated report

## Notes
This repository tracks the workflow and compact reports, not the full downloaded GWAS payloads. The mapping outputs represent orthologous mouse genomic intervals carrying human trait metadata; they are not direct mouse disease-SNP claims.
