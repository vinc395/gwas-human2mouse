#!/usr/bin/env python3
import csv
import json
from pathlib import Path


def main() -> int:
    output_root = Path('/home/vince395/gwas_human2mouse')
    input_path = output_root / 'downloads' / 'unpublished_studies_v1.0.3.1.tsv'
    reports_dir = output_root / 'reports'
    reports_dir.mkdir(parents=True, exist_ok=True)

    if not input_path.exists():
        raise FileNotFoundError(f'Missing input manifest: {input_path}')

    with input_path.open('r', encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        columns = reader.fieldnames or []
        first5 = []
        row_count = 0
        full_summary_yes = 0
        non_placeholder_assoc_count = 0
        for row in reader:
            row_count += 1
            if len(first5) < 5:
                first5.append(row)
            if (row.get('FULL SUMMARY STATISTICS') or '').strip().lower() in {'yes', 'true', '1'}:
                full_summary_yes += 1
            assoc_count = (row.get('ASSOCIATION COUNT') or '').strip().lower()
            if assoc_count and assoc_count not in {'not yet curated', 'na', 'n/a'}:
                non_placeholder_assoc_count += 1

    summary_cols = [c for c in columns if 'SUMMARY' in c.upper() or 'FTP' in c.upper() or 'LOCATION' in c.upper()]
    trait_cols = [c for c in columns if 'TRAIT' in c.upper()]
    association_like_cols = [
        c for c in columns if any(token in c.upper() for token in ['P-VALUE', 'PVALUE', 'CHR', 'CHROM', 'POSITION', 'BP', 'RSID', 'SNP'])
    ]
    study_cols = [c for c in columns if 'STUDY' in c.upper() or 'ACCESSION' in c.upper()]

    if association_like_cols:
        content_type = 'association-level rows'
    elif summary_cols and study_cols:
        content_type = 'study metadata / study manifest with summary-statistics locations'
    else:
        content_type = 'study metadata only'

    report_txt = reports_dir / 'inspection_report.txt'
    report_json = reports_dir / 'inspection_summary.json'

    text_lines = [
        f'Input file: {input_path}',
        f'Number of data rows: {row_count}',
        'Column names:',
        *[f'- {col}' for col in columns],
        '',
        f'Content assessment: {content_type}',
        f'Summary statistics location columns: {summary_cols}',
        f'Trait columns: {trait_cols}',
        f'Association-like columns: {association_like_cols}',
        f'Rows with FULL SUMMARY STATISTICS=yes: {full_summary_yes}',
        f'Rows with non-placeholder ASSOCIATION COUNT: {non_placeholder_assoc_count}',
        '',
        'First 5 rows (JSON lines):',
        *[json.dumps(row, ensure_ascii=False) for row in first5],
    ]
    report_txt.write_text('\n'.join(text_lines) + '\n', encoding='utf-8')
    report_json.write_text(
        json.dumps(
            {
                'input_file': str(input_path),
                'row_count': row_count,
                'columns': columns,
                'content_assessment': content_type,
                'summary_location_columns': summary_cols,
                'trait_columns': trait_cols,
                'association_like_columns': association_like_cols,
                'rows_with_full_summary_statistics_yes': full_summary_yes,
                'rows_with_non_placeholder_association_count': non_placeholder_assoc_count,
                'first5': first5,
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding='utf-8',
    )

    print(f'Rows: {row_count}')
    print('Columns:')
    for col in columns:
        print(f'- {col}')
    print('First 5 rows:')
    for row in first5:
        print(json.dumps(row, ensure_ascii=False))
    print(f'Content assessment: {content_type}')
    print(f'Reports written: {report_txt} and {report_json}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
