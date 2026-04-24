#!/usr/bin/env python3
import csv
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

OUTPUT_ROOT = Path('/home/vince395/gwas_human2mouse')
INPUT_MANIFEST = OUTPUT_ROOT / 'downloads' / 'unpublished_studies_v1.0.3.1.tsv'
MANIFEST_DIR = OUTPUT_ROOT / 'manifests'
REPORT_DIR = OUTPUT_ROOT / 'reports'

TRAIT_FIELDS = ['DISEASE/TRAIT', 'MAPPED_TRAIT']
CONTEXT_FIELDS = ['STUDY', 'BACKGROUND TRAIT', 'MAPPED BACKGROUND TRAIT']
ALL_MATCH_FIELDS = TRAIT_FIELDS + CONTEXT_FIELDS

KIDNEY_PATTERNS: Sequence[Tuple[str, str]] = [
    ('chronic_kidney_disease', r'\bchronic kidney disease\b'),
    ('ckd', r'\bckd\b'),
    ('renal_insufficiency', r'\brenal insufficien(?:cy|t)?\b'),
    ('renal_function', r'\brenal function\b'),
    ('kidney_function', r'\bkidney function\b'),
    ('egfr', r'\begfr(?:crea|cys)?\b|estimated glomerular filtration rate'),
    ('serum_creatinine', r'serum creatinine|(?<!urine )\bcreatinine\b'),
    ('cystatin_c', r'\bcystatin c\b'),
    ('bun', r'\bblood urea nitrogen\b|\bbun\b'),
    ('albuminuria', r'albuminuria|microalbuminuria'),
    ('uacr', r'\buacr\b|urine albumin creatinine ratio|urinary albumin creatinine ratio|albumin creatinine ratio'),
    ('proteinuria', r'proteinuria|urine protein\s*[:/]\s*creatin'),
    ('diabetic_kidney_disease', r'diabetic kidney disease|diabetic nephropathy'),
    ('nephropathy', r'\bnephropathy\b'),
    ('esrd', r'\besrd\b|\beskd\b|end[- ]stage kidney disease|end[- ]stage renal disease'),
    ('kidney_failure', r'kidney failure|renal failure|acute renal failure'),
    ('rapid_decline', r'rapid kidney function decline'),
    ('nephrolithiasis', r'nephrolithiasis|kidney stones?'),
]

COMPARATOR_PATTERNS: Sequence[Tuple[str, str]] = [
    ('hypertension', r'\bhypertension\b'),
    ('systolic_blood_pressure', r'systolic blood pressure|\bcsbp\b|central systolic blood pressure'),
    ('diastolic_blood_pressure', r'diastolic blood pressure|\bdbp\b'),
    ('pulse_pressure', r'\bpulse pressure\b|\bpp\b'),
    ('type2_diabetes', r'\btype 2 diabetes\b|\bt2d\b|type ii diabetes'),
    ('fasting_glucose', r'fasting glucose'),
    ('hba1c', r'\bhba1c\b|glycated h[ae]moglobin|hemoglobin a1c|haemoglobin a1c'),
    ('bmi', r'\bbmi\b|body mass index'),
    ('obesity', r'\bobesity\b'),
    ('waist_to_hip_ratio', r'waist[- ]to[- ]hip ratio|waist hip ratio'),
    ('urate', r'\burate\b|uric acid'),
    ('gout', r'\bgout\b'),
    ('cad_cvd', r'coronary artery disease|cardiovascular disease'),
]

COMPARATOR_EXCLUDE_PATTERNS: Sequence[re.Pattern[str]] = [
    re.compile(pattern, flags=re.IGNORECASE)
    for pattern in [
        r'portal hypertension',
        r'pulmonary hypertension',
        r'gestational',
        r'gestational hypertension',
        r'pregnancy-induced hypertension',
        r'preeclampsia',
        r'eclampsia',
        r'maternal hypertension',
        r'oesophageal varices',
        r'childbirth',
        r'puerperium',
    ]
]

EXCLUDE_IF_ONLY_CONTEXT_MATCH = True


def compile_patterns(patterns: Sequence[Tuple[str, str]]):
    return [(label, re.compile(pattern, flags=re.IGNORECASE)) for label, pattern in patterns]


def normalize_value(value: Optional[str]) -> str:
    if value is None:
        return ''
    stripped = value.strip()
    if stripped.lower() == 'not yet curated':
        return ''
    return stripped


def matches_in_fields(row: Dict[str, str], field_names: Sequence[str], compiled_patterns) -> List[Tuple[str, str, str]]:
    hits: List[Tuple[str, str, str]] = []
    for field_name in field_names:
        value = normalize_value(row.get(field_name, ''))
        if not value:
            continue
        for label, pattern in compiled_patterns:
            if pattern.search(value):
                hits.append((label, field_name, value))
    return hits


def choose_category(kidney_hits: List[Tuple[str, str, str]], comparator_hits: List[Tuple[str, str, str]]) -> Optional[str]:
    if kidney_hits:
        return 'kidney'
    if comparator_hits:
        return 'comparator'
    return None


def selection_reason(category: str, hits: List[Tuple[str, str, str]]) -> str:
    parts = []
    seen = set()
    for label, field_name, value in hits:
        key = (label, field_name, value)
        if key in seen:
            continue
        seen.add(key)
        parts.append(f'{label} matched in {field_name}: {value}')
    return '; '.join(parts)


def clean_sample_size(value: str) -> str:
    text = normalize_value(value)
    return text


def main() -> int:
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    if not INPUT_MANIFEST.exists():
        raise FileNotFoundError(f'Missing input manifest: {INPUT_MANIFEST}')

    compiled_kidney = compile_patterns(KIDNEY_PATTERNS)
    compiled_comparator = compile_patterns(COMPARATOR_PATTERNS)

    with INPUT_MANIFEST.open('r', encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        source_columns = reader.fieldnames or []
        base_columns = source_columns + [
            'category',
            'selection_reason',
            'kidney_trait_hits',
            'comparator_trait_hits',
            'summary_stats_available',
            'summary_stats_location',
            'sample_size_combined',
        ]
        all_candidates: List[Dict[str, str]] = []
        selected_rows: List[Dict[str, str]] = []

        for row in reader:
            kidney_trait_hits = matches_in_fields(row, TRAIT_FIELDS, compiled_kidney)
            comparator_trait_hits = matches_in_fields(row, TRAIT_FIELDS, compiled_comparator)
            kidney_context_hits = matches_in_fields(row, CONTEXT_FIELDS, compiled_kidney)
            comparator_context_hits = matches_in_fields(row, CONTEXT_FIELDS, compiled_comparator)

            kidney_hits = list(kidney_trait_hits)
            comparator_hits = list(comparator_trait_hits)
            if not EXCLUDE_IF_ONLY_CONTEXT_MATCH:
                kidney_hits.extend(kidney_context_hits)
                comparator_hits.extend(comparator_context_hits)

            if not kidney_hits and kidney_context_hits:
                row['_context_only_kidney_reason'] = selection_reason('kidney', kidney_context_hits)
            if not comparator_hits and comparator_context_hits:
                row['_context_only_comparator_reason'] = selection_reason('comparator', comparator_context_hits)

            category = choose_category(kidney_hits, comparator_hits)
            if category is None:
                continue

            trait_text_for_exclusion = ' | '.join(
                normalize_value(row.get(field_name, '')) for field_name in TRAIT_FIELDS
            )
            if category == 'comparator' and any(pattern.search(trait_text_for_exclusion) for pattern in COMPARATOR_EXCLUDE_PATTERNS):
                continue

            hits = kidney_hits if category == 'kidney' else comparator_hits
            enriched = dict(row)
            enriched['category'] = category
            enriched['selection_reason'] = selection_reason(category, hits)
            enriched['kidney_trait_hits'] = '; '.join(sorted({label for label, _, _ in kidney_hits}))
            enriched['comparator_trait_hits'] = '; '.join(sorted({label for label, _, _ in comparator_hits}))
            enriched['summary_stats_available'] = normalize_value(row.get('FULL SUMMARY STATISTICS', '')).lower() in {'yes', 'true', '1'} and 'yes' or 'no'
            enriched['summary_stats_location'] = normalize_value(row.get('SUMMARY STATS LOCATION', ''))
            initial_size = clean_sample_size(row.get('INITIAL SAMPLE SIZE', ''))
            repl_size = clean_sample_size(row.get('REPLICATION SAMPLE SIZE', ''))
            enriched['sample_size_combined'] = ' | '.join([part for part in [initial_size, repl_size] if part])
            all_candidates.append(enriched)
            selected_rows.append(enriched)

    selected_rows.sort(key=lambda r: (r['category'], r.get('DISEASE/TRAIT', ''), r.get('STUDY ACCESSION', '')))
    all_candidates.sort(key=lambda r: (r['category'], r.get('DISEASE/TRAIT', ''), r.get('STUDY ACCESSION', '')))

    all_candidates_tsv = MANIFEST_DIR / 'all_candidate_kidney_and_comparator_studies.tsv'
    selected_tsv = MANIFEST_DIR / 'selected_kidney_and_comparator_studies.tsv'
    selected_csv = MANIFEST_DIR / 'selected_kidney_and_comparator_studies.csv'
    manifest_report = REPORT_DIR / 'manifest_filter_report.txt'

    def write_table(path: Path, rows: List[Dict[str, str]], delimiter: str):
        with path.open('w', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=base_columns, delimiter=delimiter, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(rows)

    write_table(all_candidates_tsv, all_candidates, '\t')
    write_table(selected_tsv, selected_rows, '\t')
    write_table(selected_csv, selected_rows, ',')

    counts = {'kidney': 0, 'comparator': 0}
    full_summary_counts = {'kidney': 0, 'comparator': 0}
    for row in selected_rows:
        counts[row['category']] += 1
        if row['summary_stats_available'] == 'yes':
            full_summary_counts[row['category']] += 1

    context_excluded = []
    with INPUT_MANIFEST.open('r', encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            kidney_trait_hits = matches_in_fields(row, TRAIT_FIELDS, compiled_kidney)
            comparator_trait_hits = matches_in_fields(row, TRAIT_FIELDS, compiled_comparator)
            if kidney_trait_hits or comparator_trait_hits:
                continue
            kidney_context_hits = matches_in_fields(row, CONTEXT_FIELDS, compiled_kidney)
            comparator_context_hits = matches_in_fields(row, CONTEXT_FIELDS, compiled_comparator)
            if kidney_context_hits or comparator_context_hits:
                context_excluded.append(
                    {
                        'study_accession': row.get('STUDY ACCESSION', ''),
                        'disease_trait': row.get('DISEASE/TRAIT', ''),
                        'study': row.get('STUDY', ''),
                        'kidney_context_reason': selection_reason('kidney', kidney_context_hits) if kidney_context_hits else '',
                        'comparator_context_reason': selection_reason('comparator', comparator_context_hits) if comparator_context_hits else '',
                    }
                )

    lines = [
        f'Input manifest: {INPUT_MANIFEST}',
        f'Selected kidney studies: {counts["kidney"]}',
        f'Selected comparator studies: {counts["comparator"]}',
        f'Kidney studies with FULL SUMMARY STATISTICS=yes: {full_summary_counts["kidney"]}',
        f'Comparator studies with FULL SUMMARY STATISTICS=yes: {full_summary_counts["comparator"]}',
        f'Context-only rows excluded to avoid trait leakage: {len(context_excluded)}',
        '',
        'First 20 context-only excluded rows:',
    ]
    for row in context_excluded[:20]:
        lines.append(str(row))
    manifest_report.write_text('\n'.join(lines) + '\n', encoding='utf-8')

    print(f'Wrote {all_candidates_tsv}')
    print(f'Wrote {selected_tsv}')
    print(f'Wrote {selected_csv}')
    print(f'Kidney selected: {counts["kidney"]}')
    print(f'Comparator selected: {counts["comparator"]}')
    print(f'Context-only excluded: {len(context_excluded)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
