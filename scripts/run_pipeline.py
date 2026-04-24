#!/usr/bin/env python3
import csv
import gzip
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

import requests
import yaml


def configure_csv_field_limit() -> int:
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return limit
        except OverflowError:
            limit //= 10


CSV_FIELD_LIMIT = configure_csv_field_limit()

OUTPUT_ROOT = Path('/home/vince395/gwas_human2mouse')
SELECTED_MANIFEST = OUTPUT_ROOT / 'manifests' / 'selected_kidney_and_comparator_studies.tsv'
INVENTORY_TSV = OUTPUT_ROOT / 'reports' / 'selected_study_inventory.tsv'
DOWNLOADS_DIR = OUTPUT_ROOT / 'downloads'
PARSED_DIR = OUTPUT_ROOT / 'parsed'
REPORTS_DIR = OUTPUT_ROOT / 'reports'
LOGS_DIR = OUTPUT_ROOT / 'logs'
HUMAN_BED_DIR = OUTPUT_ROOT / 'human_bed'
MM39_BED_DIR = OUTPUT_ROOT / 'mm39_bed'
TMP_DIR = OUTPUT_ROOT / 'tmp'
SCRIPTS_DIR = OUTPUT_ROOT / 'scripts'

P_PRIMARY = 5e-8
P_SECONDARY = 1e-6
LOCUS_CLUSTER_DISTANCE = 500_000
FALLBACK_WINDOW = 50_000
SENSITIVITY_WINDOW = 500_000
BED_SCORE_CAP = 1000

TARGET_HUMAN_BUILD = 'GRCh38'
TARGET_HUMAN_UCSC = 'hg38'
TARGET_MOUSE_BUILD = 'mm39'

CHAIN_URL = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToMm39.over.chain.gz'
CHAIN_PATH = DOWNLOADS_DIR / 'hg38ToMm39.over.chain.gz'
REVERSE_CHAIN_URL = 'https://hgdownload.soe.ucsc.edu/goldenPath/mm39/liftOver/mm39ToHg38.over.chain.gz'
REVERSE_CHAIN_PATH = DOWNLOADS_DIR / 'mm39ToHg38.over.chain.gz'
HG19_TO_HG38_CHAIN_URL = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz'
HG19_TO_HG38_CHAIN_PATH = DOWNLOADS_DIR / 'hg19ToHg38.over.chain.gz'
LIFTOVER_URL = 'https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver'
LIFTOVER_BIN = DOWNLOADS_DIR / 'liftOver'

CHROM_ALLOWED = {str(i) for i in range(1, 23)} | {'X', 'Y', 'MT', 'M'}
LOG_PREFIX = '[gwas-human2mouse]'


def log(msg: str) -> None:
    print(f'{LOG_PREFIX} {msg}', flush=True)


def ensure_dirs() -> None:
    for d in [DOWNLOADS_DIR, PARSED_DIR, REPORTS_DIR, LOGS_DIR, HUMAN_BED_DIR, MM39_BED_DIR, TMP_DIR, SCRIPTS_DIR]:
        d.mkdir(parents=True, exist_ok=True)


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open('r', encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: '' if v is None else v for k, v in row.items()})


def write_csv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: '' if v is None else v for k, v in row.items()})


def normalize_text(value: Optional[str]) -> str:
    if value is None:
        return ''
    text = str(value).strip()
    if text.lower() == 'not yet curated':
        return ''
    return text


def chrom_to_ucsc(chrom_value: object) -> Optional[str]:
    if chrom_value is None:
        return None
    chrom = str(chrom_value).strip()
    chrom = chrom.removeprefix('chr').upper()
    if chrom == 'M':
        chrom = 'MT'
    if chrom not in CHROM_ALLOWED:
        return None
    if chrom == 'MT':
        return 'chrM'
    return f'chr{chrom}'


def parse_float(value: object) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text in {'NA', 'N/A', '#NA', 'null', 'None'}:
        return None
    try:
        return float(text)
    except Exception:
        return None


def parse_int(value: object) -> Optional[int]:
    f = parse_float(value)
    if f is None:
        return None
    return int(round(f))


def p_from_row(row: Dict[str, object]) -> Optional[float]:
    p = parse_float(row.get('p_value'))
    if p is not None and p > 0:
        return p
    neglog = parse_float(row.get('neg_log_10_p_value') or row.get('neg_log10_p_value'))
    if neglog is not None:
        return 10 ** (-neglog)
    mlog = parse_float(row.get('minus_log10_pvalue'))
    if mlog is not None:
        return 10 ** (-mlog)
    return None


def bed_score_from_p(p: Optional[float]) -> int:
    if p is None or p <= 0:
        return 0
    score = min(BED_SCORE_CAP, int(round(min(1000.0, -math.log10(p) * 100.0))))
    return max(0, score)


def detect_delimiter(header_line: str) -> str:
    if '\t' in header_line:
        return '\t'
    if ',' in header_line:
        return ','
    return '\t'


def open_table_text(path: Path):
    if path.suffix == '.gz':
        return gzip.open(path, 'rt', encoding='utf-8-sig', newline='')
    return path.open('r', encoding='utf-8-sig', newline='')


def iter_nonempty_lines(handle: Iterable[str]) -> Iterator[str]:
    for line in handle:
        if not line.strip():
            continue
        yield line


def sniff_table_format(path: Path) -> Tuple[str, str, bool]:
    sample_chunks: List[str] = []
    with open_table_text(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            sample_chunks.append(line)
            if sum(len(chunk) for chunk in sample_chunks) >= 65536 or len(sample_chunks) >= 100:
                break
    sample = ''.join(sample_chunks)
    if not sample:
        return '\t', '"', False
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters='\t,;')
        return dialect.delimiter, getattr(dialect, 'quotechar', '"') or '"', getattr(dialect, 'skipinitialspace', False)
    except csv.Error:
        return detect_delimiter(sample_chunks[0] if sample_chunks else sample), '"', False


def make_unique_fieldnames(fieldnames: List[object]) -> List[str]:
    seen: Dict[str, int] = {}
    cleaned: List[str] = []
    for idx, field in enumerate(fieldnames):
        base = str(field or '').strip().strip('"').strip("'")
        if not base:
            base = f'column_{idx + 1}'
        seen[base] = seen.get(base, 0) + 1
        if seen[base] > 1:
            cleaned.append(f'{base}__{seen[base]}')
        else:
            cleaned.append(base)
    return cleaned


def iter_local_table(path: Path) -> Iterator[Dict[str, str]]:
    delimiter, quotechar, skipinitialspace = sniff_table_format(path)
    with open_table_text(path) as handle:
        reader = csv.reader(
            iter_nonempty_lines(handle),
            delimiter=delimiter,
            quotechar=quotechar,
            skipinitialspace=skipinitialspace,
        )
        headers_raw = next(reader, None)
        if headers_raw is None:
            return
        headers = make_unique_fieldnames(list(headers_raw))
        for parts in reader:
            if not parts:
                continue
            row_values = list(parts)
            extras: List[str] = []
            if len(row_values) > len(headers):
                extras = row_values[len(headers):]
                row_values = row_values[:len(headers)]
            if len(row_values) < len(headers):
                row_values.extend([''] * (len(headers) - len(row_values)))
            row = {headers[idx]: row_values[idx] for idx in range(len(headers))}
            if extras:
                row['__extra_columns'] = delimiter.join(extras)
            yield row


def normalize_coordinate_system(value: object) -> str:
    text = normalize_text(value).lower().replace(' ', '')
    if text in {'0-based', '0based'}:
        return '0-based'
    if text in {'1-based', '1based'}:
        return '1-based'
    return ''


def meta_candidates_for_selected_file(data_file_name: str, inventory_row: Dict[str, str]) -> List[str]:
    meta_files = [m.strip() for m in inventory_row.get('meta_files', '').split(';') if m.strip().endswith('.yaml')]
    expected_exact = f'{data_file_name}-meta.yaml' if data_file_name else ''
    expected_basename = f'{Path(data_file_name).name}-meta.yaml' if data_file_name else ''
    expected_parent = str(Path(data_file_name).parent) if data_file_name else ''

    def sort_key(meta_file: str) -> Tuple[int, int, str]:
        meta_parent = str(Path(meta_file).parent)
        meta_basename = Path(meta_file).name
        if expected_exact and meta_file == expected_exact:
            rank = 0
        elif expected_basename and meta_basename == expected_basename and meta_parent == expected_parent:
            rank = 1
        elif expected_basename and meta_basename == expected_basename:
            rank = 2
        else:
            rank = 3
        return rank, len(meta_file), meta_file

    return sorted(meta_files, key=sort_key)


def collect_companion_files(data_file_name: str, inventory_row: Dict[str, str]) -> List[str]:
    companions: List[str] = []
    seen = set()

    def add(relpath: str) -> None:
        relpath = relpath.strip()
        if not relpath or relpath in seen:
            return
        seen.add(relpath)
        companions.append(relpath)

    for entry in [m.strip() for m in inventory_row.get('meta_files', '').split(';') if m.strip()]:
        if entry.endswith('.yaml'):
            continue
        if Path(entry).name.lower() == 'md5sum.txt':
            add(entry)

    for token in [t.strip() for t in inventory_row.get('available_entries', '').split(';') if t.strip()]:
        if token.lower().endswith('readme.txt'):
            add(token)

    default_md5 = 'harmonised/md5sum.txt' if data_file_name.startswith('harmonised/') else 'md5sum.txt'
    add(default_md5)
    return companions


def load_yaml_file(path: Path) -> Dict[str, object]:
    try:
        return yaml.safe_load(path.read_text(encoding='utf-8')) or {}
    except Exception:
        return {}


def download_study_artifacts(study_row: Dict[str, str], inventory_row: Dict[str, str]) -> Dict[str, object]:
    accession = study_row['STUDY ACCESSION']
    base_url = study_row['summary_stats_location'].rstrip('/').replace('ftp://', 'https://') + '/'
    data_file_name = inventory_row.get('data_file_name', '').strip()
    study_dir = DOWNLOADS_DIR / accession
    study_dir.mkdir(parents=True, exist_ok=True)

    context: Dict[str, object] = {
        'study_accession': accession,
        'base_url': base_url,
        'study_dir': study_dir,
        'data_file_name': data_file_name,
        'local_data_path': None,
        'local_files_downloaded': [],
        'download_notes': [],
        'metadata': {},
        'meta_file_used': '',
        'local_meta_path': None,
    }

    if not data_file_name or data_file_name.lower() == 'readme.txt':
        context['download_notes'].append('no_downloadable_association_table_listed')
        return context

    def attempt_download(relpath: str, required: bool = False) -> Optional[Path]:
        url = base_url + relpath
        local_path = study_dir / relpath
        try:
            download_file(url, local_path)
            context['local_files_downloaded'].append(relpath)
            return local_path
        except Exception as exc:
            note = f'{relpath}: {exc}'
            context['download_notes'].append(note)
            if required:
                raise
            return None

    try:
        context['local_data_path'] = attempt_download(data_file_name, required=True)
    except Exception:
        return context

    meta_candidates = meta_candidates_for_selected_file(data_file_name, inventory_row)
    for meta_relpath in meta_candidates:
        local_meta_path = attempt_download(meta_relpath, required=False)
        if local_meta_path and local_meta_path.exists():
            context['metadata'] = load_yaml_file(local_meta_path)
            context['meta_file_used'] = meta_relpath
            context['local_meta_path'] = local_meta_path
            if context['metadata']:
                break

    for companion in collect_companion_files(data_file_name, inventory_row):
        attempt_download(companion, required=False)

    return context


def build_detection(study_row: Dict[str, str], inventory_row: Dict[str, str], local_context: Dict[str, object]) -> Dict[str, object]:
    meta = local_context.get('metadata') or {}
    meta_file = str(local_context.get('meta_file_used') or '')
    data_file_name = inventory_row.get('data_file_name', '')
    build = normalize_text(meta.get('genome_assembly'))
    coordinate_system = normalize_coordinate_system(meta.get('coordinate_system'))
    notes = []
    if build:
        notes.append(f'meta:{meta_file}')
    else:
        notes.append('no_meta_build')
    if not build:
        if data_file_name.startswith('harmonised/'):
            build = TARGET_HUMAN_BUILD
            notes.append('assumed_from_harmonised_reference_release95')
        else:
            build = 'unknown'
    if coordinate_system:
        notes.append(f'coordinate_system:{coordinate_system}')
    else:
        coordinate_system = '1-based'
        notes.append('coordinate_system_defaulted_to_1-based')
    ancestry = ''
    sample_size = ''
    if isinstance(meta.get('samples'), list) and meta['samples']:
        first_sample = meta['samples'][0]
        ancestry = '; '.join(first_sample.get('sample_ancestry_category', []) if isinstance(first_sample.get('sample_ancestry_category'), list) else [])
        sample_size = str(first_sample.get('sample_size', ''))
    trait_description = ''
    if isinstance(meta.get('trait_description'), list):
        trait_description = '; '.join(str(x) for x in meta.get('trait_description', []))
    return {
        'study_accession': study_row['STUDY ACCESSION'],
        'category': study_row['category'],
        'trait_name': study_row['DISEASE/TRAIT'],
        'detected_build': build,
        'coordinate_system': coordinate_system,
        'meta_file_used': meta_file,
        'detection_notes': '; '.join(notes),
        'ancestry_from_meta': ancestry,
        'sample_size_from_meta': sample_size,
        'trait_description_from_meta': trait_description,
        'source_path': study_row['summary_stats_location'],
        'data_file_name': data_file_name,
        'local_data_file': str(local_context.get('local_data_path') or ''),
        'local_meta_file': str(local_context.get('local_meta_path') or ''),
    }


def choose_source_type(inventory_row: Dict[str, str]) -> str:
    data_name = inventory_row.get('data_file_name', '')
    if data_name.startswith('harmonised/'):
        return 'harmonised_api'
    if data_name.endswith(('.tsv', '.tsv.gz', '.txt', '.txt.gz', '.csv', '.csv.gz')):
        return 'full_summary_stats'
    if data_name == 'readme.txt':
        return 'unavailable'
    return 'unavailable'


def download_file(url: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0:
        return
    log(f'downloading {url} -> {path}')
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        with path.open('wb') as handle:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)


def ensure_liftover_assets() -> None:
    download_file(CHAIN_URL, CHAIN_PATH)
    download_file(REVERSE_CHAIN_URL, REVERSE_CHAIN_PATH)
    download_file(HG19_TO_HG38_CHAIN_URL, HG19_TO_HG38_CHAIN_PATH)
    download_file(LIFTOVER_URL, LIFTOVER_BIN)
    LIFTOVER_BIN.chmod(0o755)


def normalize_human_interval(chrom: str, start: int, end: int, coordinate_system: str) -> Tuple[int, int]:
    if coordinate_system == '0-based':
        normalized_start = start + 1
        normalized_end = max(start + 1, end)
    else:
        normalized_start = start
        normalized_end = max(start, end)
    return normalized_start, normalized_end


def interval_to_bed(start_1based: int, end_1based_inclusive: int) -> Tuple[int, int]:
    bed_start = max(0, start_1based - 1)
    bed_end = max(bed_start + 1, end_1based_inclusive)
    return bed_start, bed_end


def feature_identifier(row: Dict[str, object]) -> str:
    feature_type = str(row.get('feature_type') or 'snp')
    if feature_type == 'interval':
        return normalize_text(row.get('variant_id')) or f"{row['chromosome']}:{row['base_pair_start']}-{row['base_pair_end']}"
    return normalize_text(row.get('rsid')) or normalize_text(row.get('variant_id')) or f"{row['chromosome']}:{row['base_pair_location']}"


def dedup_key(row: Dict[str, object]) -> str:
    feature_type = str(row.get('feature_type') or 'snp')
    if feature_type == 'interval':
        return f"interval:{row['chromosome']}:{row['base_pair_start']}-{row['base_pair_end']}:{normalize_text(row.get('variant_id'))}"
    snp_id = normalize_text(row.get('rsid')) or normalize_text(row.get('variant_id')) or f"{row['chromosome']}:{row['base_pair_location']}"
    return f"snp:{snp_id}"


def write_interval_records_to_bed(records: List[Dict[str, object]], bed_path: Path) -> None:
    with bed_path.open('w', encoding='utf-8') as handle:
        for record in records:
            handle.write(
                f"{record['chromosome']}\t{record['bed_start']}\t{record['bed_end']}\t{record['bed_name']}\t{record['bed_score']}\t.\n"
            )


def liftover_human_records(records: List[Dict[str, object]], chain_path: Path, prefix: str) -> Tuple[List[Dict[str, object]], int, int]:
    if not records:
        return [], 0, 0
    input_bed = TMP_DIR / f'{prefix}.hg_input.bed'
    mapped_bed = TMP_DIR / f'{prefix}.hg_lifted.bed'
    unmapped_bed = TMP_DIR / f'{prefix}.hg_unmapped.bed'
    write_interval_records_to_bed(records, input_bed)
    liftover_file(input_bed, mapped_bed, unmapped_bed, chain_path)
    lifted_rows = read_bed(mapped_bed)
    mapped_by_name: Dict[str, List[Dict[str, object]]] = {}
    for row in lifted_rows:
        mapped_by_name.setdefault(str(row['name']), []).append(row)
    lifted_records: List[Dict[str, object]] = []
    mapped_count = 0
    unmapped_count = 0
    for record in records:
        lifted_list = mapped_by_name.get(str(record['bed_name']), [])
        if len(lifted_list) != 1:
            unmapped_count += 1
            continue
        lifted = lifted_list[0]
        lifted_record = dict(record)
        lifted_record['chromosome'] = lifted['chrom']
        lifted_record['base_pair_start'] = int(lifted['start']) + 1
        lifted_record['base_pair_end'] = int(lifted['end'])
        if str(record.get('feature_type')) == 'snp':
            lifted_record['base_pair_location'] = lifted_record['base_pair_start']
        lifted_record['bed_start'] = int(lifted['start'])
        lifted_record['bed_end'] = int(lifted['end'])
        lifted_record['normalized_build'] = TARGET_HUMAN_BUILD
        lifted_record['human_coordinate_system'] = '1-based'
        lifted_record['human_coordinate_source'] = f"lifted_from_{record.get('genome_build') or 'unknown'}"
        lifted_record['build_conversion_method'] = 'liftOver_hg19ToHg38'
        lifted_record['build_conversion_status'] = 'lifted'
        mapped_count += 1
        lifted_records.append(lifted_record)
    return lifted_records, mapped_count, unmapped_count


def normalize_significant_records(
    rows: List[Dict[str, object]],
    build_row: Dict[str, object],
    study_accession: str,
) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    detected_build = str(build_row.get('detected_build') or 'unknown')
    coordinate_system = str(build_row.get('coordinate_system') or '1-based')
    normalized_records: List[Dict[str, object]] = []
    skipped_no_coords = 0
    for row in rows:
        chrom = chrom_to_ucsc(row.get('chromosome'))
        if chrom is None:
            skipped_no_coords += 1
            continue
        feature_type = str(row.get('feature_type') or 'snp')
        if feature_type == 'interval':
            raw_start = parse_int(row.get('base_pair_start'))
            raw_end = parse_int(row.get('base_pair_end'))
            if raw_start is None or raw_end is None:
                skipped_no_coords += 1
                continue
            normalized_start, normalized_end = normalize_human_interval(chrom, raw_start, raw_end, coordinate_system)
            if normalized_end < normalized_start:
                normalized_start, normalized_end = normalized_end, normalized_start
            row['base_pair_start'] = normalized_start
            row['base_pair_end'] = normalized_end
            row['base_pair_location'] = normalized_start
        else:
            raw_pos = parse_int(row.get('base_pair_location'))
            if raw_pos is None:
                skipped_no_coords += 1
                continue
            normalized_start, normalized_end = normalize_human_interval(chrom, raw_pos, raw_pos, coordinate_system)
            row['base_pair_start'] = normalized_start
            row['base_pair_end'] = normalized_end
            row['base_pair_location'] = normalized_start
        row['chromosome'] = chrom
        bed_start, bed_end = interval_to_bed(int(row['base_pair_start']), int(row['base_pair_end']))
        row['bed_start'] = bed_start
        row['bed_end'] = bed_end
        row['bed_score'] = bed_score_from_p(parse_float(row.get('p_value')))
        row['bed_name'] = bed_name([row['trait_name'], row['study_accession'], feature_identifier(row), row['source_type']])
        row['human_coordinate_system'] = '1-based'
        row['normalized_build'] = detected_build
        row['human_coordinate_source'] = 'metadata_normalized'
        row['build_conversion_method'] = 'none_needed' if detected_build == TARGET_HUMAN_BUILD else ''
        row['build_conversion_status'] = 'retained' if detected_build == TARGET_HUMAN_BUILD else 'pending'
        normalized_records.append(row)

    conversion_summary = {
        'study_accession': study_accession,
        'input_build': detected_build,
        'target_build': TARGET_HUMAN_BUILD,
        'coordinate_system': coordinate_system,
        'conversion_method': 'none_needed',
        'chain_used': '',
        'status': 'retained' if detected_build == TARGET_HUMAN_BUILD else 'pending',
        'input_feature_count': len(normalized_records),
        'output_feature_count': len(normalized_records),
        'unmapped_feature_count': 0,
        'notes': f'skipped_missing_coordinates={skipped_no_coords}',
    }

    if detected_build == TARGET_HUMAN_BUILD:
        return normalized_records, conversion_summary
    if detected_build != 'GRCh37':
        conversion_summary['conversion_method'] = 'unsupported_input_build'
        conversion_summary['status'] = 'excluded_if_no_coordinates'
        conversion_summary['output_feature_count'] = 0
        conversion_summary['notes'] = f"unsupported_input_build={detected_build}; skipped_missing_coordinates={skipped_no_coords}"
        return [], conversion_summary

    lifted_records, lifted_count, unmapped_count = liftover_human_records(normalized_records, HG19_TO_HG38_CHAIN_PATH, f'{study_accession}.hg19_to_hg38')
    conversion_summary['conversion_method'] = 'liftOver_hg19ToHg38'
    conversion_summary['chain_used'] = HG19_TO_HG38_CHAIN_PATH.name
    conversion_summary['status'] = 'lifted' if lifted_records else 'unmapped'
    conversion_summary['output_feature_count'] = lifted_count
    conversion_summary['unmapped_feature_count'] = unmapped_count
    conversion_summary['notes'] = f"skipped_missing_coordinates={skipped_no_coords};lifted_features={lifted_count};unmapped_after_liftover={unmapped_count}"
    return lifted_records, conversion_summary


def extract_significant_rows(
    study_row: Dict[str, str],
    inventory_row: Dict[str, str],
    build_row: Dict[str, object],
    local_context: Dict[str, object],
) -> Tuple[List[Dict[str, object]], Dict[str, object], Dict[str, object]]:
    accession = study_row['STUDY ACCESSION']
    data_file_name = inventory_row.get('data_file_name', '')
    source_type = choose_source_type(inventory_row)
    evidence_level = 'full_summary_stats' if source_type in {'full_summary_stats', 'harmonised_api'} else 'index_only'
    local_data_path = local_context.get('local_data_path')
    local_files_downloaded = '; '.join(str(x) for x in local_context.get('local_files_downloaded', []))
    record_log = {
        'study_accession': accession,
        'category': study_row['category'],
        'trait_name': study_row['DISEASE/TRAIT'],
        'source_path': study_row['summary_stats_location'],
        'data_file_name': data_file_name,
        'download_status': 'available' if data_file_name and data_file_name != 'readme.txt' else 'unavailable',
        'file_format_guess': data_file_name.split('.')[-1] if data_file_name else '',
        'local_files_downloaded': local_files_downloaded,
        'file_size': inventory_row.get('data_file_size_text', ''),
        'source_type': source_type,
        'notes': '',
    }
    empty_conversion = {
        'study_accession': accession,
        'input_build': str(build_row.get('detected_build') or 'unknown'),
        'target_build': TARGET_HUMAN_BUILD,
        'coordinate_system': str(build_row.get('coordinate_system') or ''),
        'conversion_method': 'not_run',
        'chain_used': '',
        'status': 'no_data',
        'input_feature_count': 0,
        'output_feature_count': 0,
        'unmapped_feature_count': 0,
        'notes': '',
    }
    if not data_file_name or data_file_name == 'readme.txt':
        record_log['notes'] = 'No downloadable association table found'
        empty_conversion['status'] = 'unavailable'
        empty_conversion['notes'] = 'No downloadable association table found'
        return [], record_log, empty_conversion
    if not local_data_path or not Path(local_data_path).exists():
        record_log['download_status'] = 'error'
        record_log['notes'] = '; '.join(local_context.get('download_notes', [])) or 'selected data file was not downloaded locally'
        empty_conversion['status'] = 'download_failed'
        empty_conversion['notes'] = record_log['notes']
        return [], record_log, empty_conversion

    rows: List[Dict[str, object]] = []
    total_seen = 0
    total_primary = 0
    total_secondary = 0

    try:
        for raw_row in iter_local_table(Path(local_data_path)):
            total_seen += 1
            p_value = p_from_row(raw_row)
            if p_value is None or p_value > P_SECONDARY:
                continue
            chrom_raw = raw_row.get('chromosome') or raw_row.get('chr') or raw_row.get('Chromosome')
            chrom = chrom_to_ucsc(chrom_raw)
            if chrom is None:
                continue
            raw_start = parse_int(raw_row.get('base_pair_start') or raw_row.get('bp_start') or raw_row.get('chrom_start'))
            raw_end = parse_int(raw_row.get('base_pair_end') or raw_row.get('bp_end') or raw_row.get('chrom_end'))
            raw_pos = parse_int(raw_row.get('base_pair_location') or raw_row.get('position') or raw_row.get('pos'))
            feature_type = 'interval' if raw_start is not None and raw_end is not None else 'snp'
            if feature_type == 'interval':
                base_pair_start = raw_start
                base_pair_end = raw_end
                base_pair_location = raw_start
            else:
                if raw_pos is None:
                    continue
                base_pair_start = raw_pos
                base_pair_end = raw_pos
                base_pair_location = raw_pos
            variant_id = normalize_text(raw_row.get('variant_id') or raw_row.get('cnv_id') or raw_row.get('VariantID'))
            rsid = normalize_text(raw_row.get('rsid'))
            if feature_type == 'snp' and not rsid and variant_id.lower().startswith('rs'):
                rsid = variant_id
            if not variant_id:
                if feature_type == 'interval':
                    variant_id = normalize_text(raw_row.get('cnv_id')) or f'{chrom}:{base_pair_start}-{base_pair_end}'
                else:
                    variant_id = normalize_text(raw_row.get('variant_id')) or rsid or f'{chrom}:{base_pair_location}'
            row_out: Dict[str, object] = {
                'study_accession': accession,
                'trait_name': study_row['DISEASE/TRAIT'],
                'category': study_row['category'],
                'ancestry': normalize_text(build_row.get('ancestry_from_meta') or study_row.get('COHORT', '')),
                'variant_id': variant_id,
                'rsid': rsid,
                'feature_type': feature_type,
                'chromosome': chrom,
                'base_pair_location': base_pair_location,
                'base_pair_start': base_pair_start,
                'base_pair_end': base_pair_end,
                'p_value': p_value,
                'effect_allele': normalize_text(raw_row.get('effect_allele')),
                'other_allele': normalize_text(raw_row.get('other_allele')),
                'effect_allele_frequency': parse_float(raw_row.get('effect_allele_frequency')),
                'beta': parse_float(raw_row.get('beta')),
                'odds_ratio': parse_float(raw_row.get('odds_ratio')),
                'se': parse_float(raw_row.get('standard_error') or raw_row.get('se')),
                'source_file': str(local_data_path),
                'source_type': source_type,
                'genome_build': build_row['detected_build'],
                'human_coordinate_system': build_row.get('coordinate_system', ''),
                'n': parse_float(raw_row.get('n')),
                'evidence_level': evidence_level,
                'signal_tier': 'primary' if p_value <= P_PRIMARY else 'secondary',
                'raw_row_json': json.dumps(raw_row, ensure_ascii=False),
            }
            rows.append(row_out)
            if row_out['signal_tier'] == 'primary':
                total_primary += 1
            else:
                total_secondary += 1
    except Exception as exc:
        record_log['download_status'] = 'error'
        record_log['notes'] = str(exc)
        empty_conversion['status'] = 'parse_failed'
        empty_conversion['notes'] = str(exc)
        return [], record_log, empty_conversion

    normalized_rows, conversion_row = normalize_significant_records(rows, build_row, accession)
    record_notes = [
        f'total_rows_scanned={total_seen}',
        f'primary={total_primary}',
        f'secondary={total_secondary}',
        f'normalized_retained={len(normalized_rows)}',
        f"build_conversion={conversion_row.get('conversion_method')}",
        f"conversion_status={conversion_row.get('status')}",
    ]
    if local_context.get('download_notes'):
        record_notes.append('download_notes=' + '|'.join(str(x) for x in local_context.get('download_notes', [])))
    record_log['notes'] = ';'.join(record_notes)
    return normalized_rows, record_log, conversion_row


def deduplicate_variants(rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    best: Dict[Tuple[str, str], Dict[str, object]] = {}
    for row in rows:
        key = (str(row['study_accession']), dedup_key(row))
        prev = best.get(key)
        if prev is None:
            best[key] = row
            continue
        row_type = str(row.get('feature_type') or 'snp')
        prev_type = str(prev.get('feature_type') or 'snp')
        if row_type == 'snp' and prev_type != 'snp':
            best[key] = row
            continue
        if row_type == prev_type and (row.get('p_value') or 1.0) < (prev.get('p_value') or 1.0):
            best[key] = row
    return list(best.values())


def make_primary_secondary_tables(rows: List[Dict[str, object]]) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    primary = [r for r in rows if r['signal_tier'] == 'primary']
    secondary = [r for r in rows if r['signal_tier'] == 'secondary']
    return deduplicate_variants(primary), deduplicate_variants(secondary)


def group_loci(rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    grouped: Dict[Tuple[str, str, str], List[Dict[str, object]]] = {}
    for row in rows:
        key = (str(row['study_accession']), str(row['trait_name']), str(row['chromosome']))
        grouped.setdefault(key, []).append(row)
    loci: List[Dict[str, object]] = []
    for subset in grouped.values():
        subset = sorted(subset, key=lambda r: (int(r['base_pair_start']), int(r['base_pair_end'])))
        cluster: List[Dict[str, object]] = []
        cluster_end = 0
        for row in subset:
            row_start = int(row['base_pair_start'])
            row_end = int(row['base_pair_end'])
            if not cluster:
                cluster = [row]
                cluster_end = row_end
                continue
            if row_start - cluster_end <= LOCUS_CLUSTER_DISTANCE:
                cluster.append(row)
                cluster_end = max(cluster_end, row_end)
            else:
                loci.append(cluster_to_locus(cluster))
                cluster = [row]
                cluster_end = row_end
        if cluster:
            loci.append(cluster_to_locus(cluster))
    return loci


def cluster_to_locus(cluster: List[Dict[str, object]]) -> Dict[str, object]:
    lead = min(cluster, key=lambda r: r.get('p_value') or 1.0)
    starts = [int(r['base_pair_start']) for r in cluster]
    ends = [int(r['base_pair_end']) for r in cluster]
    cluster_has_interval = any(str(r.get('feature_type')) == 'interval' for r in cluster)
    lead_feature_type = str(lead.get('feature_type') or 'snp')
    if len(cluster) == 1 and not cluster_has_interval and lead_feature_type == 'snp':
        locus_start = max(1, int(lead['base_pair_location']) - FALLBACK_WINDOW)
        locus_end = int(lead['base_pair_location']) + FALLBACK_WINDOW
    else:
        locus_start = min(starts)
        locus_end = max(ends)
    lead_feature_bed_name = bed_name([lead['trait_name'], lead['study_accession'], feature_identifier(lead), lead['source_type']])
    return {
        'study_accession': lead['study_accession'],
        'trait_name': lead['trait_name'],
        'category': lead['category'],
        'ancestry': lead.get('ancestry', ''),
        'feature_type': 'locus',
        'contains_interval': int(cluster_has_interval),
        'chromosome': lead['chromosome'],
        'lead_snp': feature_identifier(lead),
        'lead_variant_id': normalize_text(lead.get('variant_id')),
        'lead_rsid': normalize_text(lead.get('rsid')),
        'lead_feature_type': lead_feature_type,
        'lead_feature_start': int(lead['base_pair_start']),
        'lead_feature_end': int(lead['base_pair_end']),
        'lead_feature_bed_name': lead_feature_bed_name,
        'lead_position': lead['base_pair_location'],
        'lead_p_value': lead['p_value'],
        'locus_start': locus_start,
        'locus_end': locus_end,
        'cluster_size': len(cluster),
        'evidence_level': lead.get('evidence_level', ''),
        'source_type': lead.get('source_type', ''),
        'genome_build': lead.get('genome_build', ''),
        'normalized_build': lead.get('normalized_build', ''),
        'human_coordinate_system': lead.get('human_coordinate_system', ''),
        'human_coordinate_source': lead.get('human_coordinate_source', ''),
        'build_conversion_method': lead.get('build_conversion_method', ''),
        'build_conversion_status': lead.get('build_conversion_status', ''),
        'supporting_variants': ';'.join(feature_identifier(r) for r in cluster),
        'sensitivity_start': max(1, int(lead['base_pair_location']) - SENSITIVITY_WINDOW),
        'sensitivity_end': int(lead['base_pair_location']) + SENSITIVITY_WINDOW,
    }


def save_variant_tables(all_rows: List[Dict[str, object]], primary_rows: List[Dict[str, object]], secondary_rows: List[Dict[str, object]]) -> None:
    fieldnames = [
        'study_accession', 'trait_name', 'category', 'ancestry', 'variant_id', 'rsid',
        'feature_type', 'chromosome', 'base_pair_location', 'base_pair_start', 'base_pair_end',
        'bed_start', 'bed_end', 'p_value', 'effect_allele', 'other_allele',
        'effect_allele_frequency', 'beta', 'odds_ratio', 'se', 'source_file', 'source_type',
        'genome_build', 'normalized_build', 'human_coordinate_system', 'human_coordinate_source',
        'build_conversion_method', 'build_conversion_status', 'n', 'evidence_level', 'signal_tier', 'raw_row_json'
    ]
    write_tsv(PARSED_DIR / 'all_significant_variants.tsv', all_rows, fieldnames)
    write_tsv(PARSED_DIR / 'primary_variants.tsv', primary_rows, fieldnames)
    write_tsv(PARSED_DIR / 'secondary_variants.tsv', secondary_rows, fieldnames)


def save_locus_table(loci: List[Dict[str, object]]) -> None:
    fieldnames = [
        'study_accession', 'trait_name', 'category', 'ancestry', 'feature_type', 'contains_interval',
        'chromosome', 'lead_snp', 'lead_variant_id', 'lead_rsid', 'lead_feature_type',
        'lead_feature_start', 'lead_feature_end', 'lead_feature_bed_name', 'lead_position', 'lead_p_value',
        'locus_start', 'locus_end', 'cluster_size', 'evidence_level', 'source_type', 'genome_build',
        'normalized_build', 'human_coordinate_system', 'human_coordinate_source',
        'build_conversion_method', 'build_conversion_status', 'supporting_variants',
        'sensitivity_start', 'sensitivity_end'
    ]
    write_tsv(PARSED_DIR / 'primary_loci.tsv', loci, fieldnames)


def bed_name(parts: Iterable[object]) -> str:
    return '|'.join(str(p).replace(' ', '_') for p in parts)


def write_bed_and_detail(
    rows: List[Dict[str, object]],
    bed_path: Path,
    detail_path: Path,
    locus_mode: bool,
) -> None:
    detail_rows: List[Dict[str, object]] = []
    with bed_path.open('w', encoding='utf-8') as bed_handle:
        for row in rows:
            chrom = row['chromosome']
            if locus_mode:
                start, end = interval_to_bed(int(row['locus_start']), int(row['locus_end']))
                label = row['lead_snp']
                score = bed_score_from_p(parse_float(row['lead_p_value']))
                identifier = f"{chrom}:{int(row['locus_start'])}-{int(row['locus_end'])}"
            else:
                start = int(row['bed_start'])
                end = int(row['bed_end'])
                label = feature_identifier(row)
                score = bed_score_from_p(parse_float(row['p_value']))
                identifier = label
            name = bed_name([row['trait_name'], row['study_accession'], identifier, row['source_type']])
            bed_handle.write(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t.\n')
            out = dict(row)
            out['feature_name'] = label if not locus_mode else row.get('lead_snp', '')
            out['bed_name'] = name
            out['bed_start'] = start
            out['bed_end'] = end
            out['bed_score'] = score
            detail_rows.append(out)
    fieldnames = sorted({k for row in detail_rows for k in row.keys()}) if detail_rows else ['bed_name']
    write_tsv(detail_path, detail_rows, fieldnames)


def split_categories(primary_rows: List[Dict[str, object]], loci: List[Dict[str, object]]) -> None:
    kidney_snps = [r for r in primary_rows if r['category'] == 'kidney']
    comp_snps = [r for r in primary_rows if r['category'] == 'comparator']
    kidney_loci = [r for r in loci if r['category'] == 'kidney']
    comp_loci = [r for r in loci if r['category'] == 'comparator']
    write_bed_and_detail(kidney_snps, HUMAN_BED_DIR / 'kidney_lead_snps_human.bed', HUMAN_BED_DIR / 'kidney_lead_snps_human.tsv', locus_mode=False)
    write_bed_and_detail(kidney_loci, HUMAN_BED_DIR / 'kidney_loci_human.bed', HUMAN_BED_DIR / 'kidney_loci_human.tsv', locus_mode=True)
    write_bed_and_detail(comp_snps, HUMAN_BED_DIR / 'comparator_lead_snps_human.bed', HUMAN_BED_DIR / 'comparator_lead_snps_human.tsv', locus_mode=False)
    write_bed_and_detail(comp_loci, HUMAN_BED_DIR / 'comparator_loci_human.bed', HUMAN_BED_DIR / 'comparator_loci_human.tsv', locus_mode=True)


def run_cmd(cmd: List[str], cwd: Optional[Path] = None) -> subprocess.CompletedProcess:
    log('running: ' + ' '.join(cmd))
    return subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=True, text=True, capture_output=True)


def liftover_file(input_bed: Path, mapped_bed: Path, unmapped_bed: Path, chain_path: Path) -> None:
    run_cmd([str(LIFTOVER_BIN), str(input_bed), str(chain_path), str(mapped_bed), str(unmapped_bed)])


def read_bed(path: Path) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    if not path.exists():
        return rows
    with path.open('r', encoding='utf-8') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 6:
                continue
            rows.append({
                'chrom': parts[0],
                'start': int(parts[1]),
                'end': int(parts[2]),
                'name': parts[3],
                'score': parts[4],
                'strand': parts[5],
            })
    return rows


def read_detail_table(path: Path) -> Dict[str, Dict[str, str]]:
    if not path.exists():
        return {}
    rows = read_tsv(path)
    return {r['bed_name']: r for r in rows}


def interval_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def merge_fragments_if_close(fragments: List[Dict[str, object]]) -> Tuple[List[Dict[str, object]], str]:
    if not fragments:
        return [], 'unmapped'
    fragments = sorted(fragments, key=lambda x: (x['chrom'], x['start']))
    chroms = {f['chrom'] for f in fragments}
    if len(chroms) > 1:
        return fragments, 'multi_mapped'
    merged = [dict(fragments[0])]
    for frag in fragments[1:]:
        last = merged[-1]
        if frag['start'] - last['end'] <= 10_000:
            last['end'] = max(last['end'], frag['end'])
        else:
            merged.append(dict(frag))
    status = 'unique' if len(merged) == 1 and len(fragments) == 1 else ('merged_nearby' if len(merged) == 1 else 'multi_mapped')
    return merged, status


def compute_confidence(
    original_len: int,
    mapped_len: int,
    reciprocal_overlap_fraction: float,
    mapping_status: str,
    lead_snp_mapped: bool,
    boundary_score: float,
    lead_feature_available: bool,
) -> Tuple[float, str, float, float, float, float]:
    fraction_mapped = min(1.0, mapped_len / original_len) if original_len > 0 else 0.0
    if mapping_status == 'unique':
        unique_mapping_score = 1.0
    elif mapping_status == 'merged_nearby':
        unique_mapping_score = 0.6
    elif mapping_status == 'multi_mapped':
        unique_mapping_score = 0.3
    else:
        unique_mapping_score = 0.0
    if reciprocal_overlap_fraction >= 0.8:
        reciprocal_score = 1.0
    elif reciprocal_overlap_fraction >= 0.2:
        reciprocal_score = 0.5
    else:
        reciprocal_score = 0.0
    if lead_feature_available:
        lead_snp_score = 1.0 if lead_snp_mapped else 0.0
    else:
        lead_snp_score = 1.0 if lead_snp_mapped else (0.5 if mapped_len > 0 else 0.0)
    confidence = (
        0.35 * fraction_mapped +
        0.20 * unique_mapping_score +
        0.20 * reciprocal_score +
        0.15 * lead_snp_score +
        0.10 * boundary_score
    )
    if confidence >= 0.80:
        label = 'high'
    elif confidence >= 0.50:
        label = 'medium'
    elif confidence > 0:
        label = 'low'
    else:
        label = 'unmapped'
    return confidence, label, fraction_mapped, unique_mapping_score, reciprocal_score, lead_snp_score


def annotate_mm39(
    input_bed: Path,
    detail_tsv: Path,
    output_bed: Path,
    output_detail_tsv: Path,
    unmapped_report: List[Dict[str, object]],
    multi_report: List[Dict[str, object]],
    is_locus: bool,
) -> Tuple[int, int]:
    mapped_bed = TMP_DIR / f'{input_bed.stem}.mm39.raw.bed'
    unmapped_bed = TMP_DIR / f'{input_bed.stem}.unmapped.bed'
    reciprocal_bed = TMP_DIR / f'{input_bed.stem}.reciprocal.hg38.bed'
    reciprocal_unmapped = TMP_DIR / f'{input_bed.stem}.reciprocal.unmapped.bed'

    liftover_file(input_bed, mapped_bed, unmapped_bed, CHAIN_PATH)
    liftover_file(mapped_bed, reciprocal_bed, reciprocal_unmapped, REVERSE_CHAIN_PATH)

    original_rows = read_bed(input_bed)
    mapped_rows = read_bed(mapped_bed)
    reciprocal_rows = read_bed(reciprocal_bed)
    details = read_detail_table(detail_tsv)

    mapped_by_name: Dict[str, List[Dict[str, object]]] = {}
    for row in mapped_rows:
        mapped_by_name.setdefault(row['name'], []).append(row)
    reciprocal_by_name: Dict[str, List[Dict[str, object]]] = {}
    for row in reciprocal_rows:
        reciprocal_by_name.setdefault(row['name'], []).append(row)

    detailed_rows: List[Dict[str, object]] = []
    final_bed_rows: List[Dict[str, object]] = []
    mapped_count = 0
    unmapped_count = 0
    lead_feature_lookup: Dict[str, bool] = {}
    if is_locus:
        lead_detail_path: Optional[Path] = None
        if detail_tsv.name.startswith('kidney_loci_'):
            lead_detail_path = MM39_BED_DIR / 'kidney_lead_snps_mm39_detailed.tsv'
        elif detail_tsv.name.startswith('comparator_loci_'):
            lead_detail_path = MM39_BED_DIR / 'comparator_lead_snps_mm39_detailed.tsv'
        if lead_detail_path and lead_detail_path.exists():
            for row in read_tsv(lead_detail_path):
                bed_name_value = row.get('bed_name')
                if not bed_name_value:
                    continue
                mapped_flag = row.get('mapping_status') != 'unmapped'
                lead_feature_lookup[bed_name_value] = lead_feature_lookup.get(bed_name_value, False) or mapped_flag

    for original in original_rows:
        name = original['name']
        meta = details.get(name, {})
        frags = mapped_by_name.get(name, [])
        merged, mapping_status = merge_fragments_if_close(frags)
        original_bp = int(original['end']) - int(original['start'])
        mapped_bp = sum(int(f['end']) - int(f['start']) for f in merged)
        reciprocal_frags = reciprocal_by_name.get(name, [])
        reciprocal_overlap = 0
        for rf in reciprocal_frags:
            reciprocal_overlap += interval_overlap(int(original['start']), int(original['end']), int(rf['start']), int(rf['end']))
        reciprocal_overlap_fraction = min(1.0, reciprocal_overlap / original_bp) if original_bp > 0 else 0.0
        lead_feature_name = meta.get('lead_feature_bed_name') if is_locus else name
        lead_feature_available = bool(lead_feature_name)
        if is_locus and lead_feature_available:
            lead_snp_mapped = bool(lead_feature_lookup.get(str(lead_feature_name), False))
        else:
            lead_snp_mapped = bool(lead_feature_name and mapped_by_name.get(str(lead_feature_name), []))
        if not is_locus and not lead_feature_available:
            lead_feature_available = True
            lead_snp_mapped = len(frags) > 0
        if mapping_status == 'unique':
            boundary_score = 1.0
        elif mapping_status == 'merged_nearby':
            boundary_score = 0.5
        elif mapping_status == 'multi_mapped':
            boundary_score = 0.0
        else:
            boundary_score = 0.0
        confidence, conf_label, fraction_mapped, unique_mapping_score, reciprocal_score, lead_snp_score = compute_confidence(
            original_bp,
            mapped_bp,
            reciprocal_overlap_fraction,
            mapping_status,
            lead_snp_mapped,
            boundary_score,
            lead_feature_available,
        )
        feature_label = meta.get('feature_type') or ('locus' if is_locus else meta.get('feature_type', 'snp'))
        if not merged:
            unmapped_count += 1
            unmapped_report.append({
                'feature_name': name,
                'feature_type': feature_label,
                'category': meta.get('category', ''),
                'study_accession': meta.get('study_accession', ''),
                'trait_name': meta.get('trait_name', ''),
                'source_type': meta.get('source_type', ''),
                'original_chrom': original['chrom'],
                'original_start': original['start'],
                'original_end': original['end'],
                'reason': 'liftOver_unmapped',
            })
            detailed_rows.append({
                **meta,
                'feature_name': name,
                'mapping_status': 'unmapped',
                'mouse_chrom': '',
                'mouse_start': '',
                'mouse_end': '',
                'original_human_bp': original_bp,
                'mapped_mouse_bp': 0,
                'fraction_mapped': 0,
                'reciprocal_overlap_fraction': reciprocal_overlap_fraction,
                'unique_vs_split_mapping': 'unmapped',
                'lead_snp_mapped': int(lead_snp_mapped),
                'lead_feature_available': int(lead_feature_available),
                'boundary_score': boundary_score,
                'confidence_score': 0,
                'confidence_label': 'unmapped',
            })
            continue

        mapped_count += 1
        if mapping_status == 'multi_mapped' or len(frags) > 1:
            multi_report.append({
                'feature_name': name,
                'feature_type': feature_label,
                'category': meta.get('category', ''),
                'study_accession': meta.get('study_accession', ''),
                'trait_name': meta.get('trait_name', ''),
                'fragment_count': len(frags),
                'merged_fragment_count': len(merged),
                'mapping_status': mapping_status,
            })
        for idx, frag in enumerate(merged, start=1):
            detail_row = {
                **meta,
                'feature_name': name,
                'fragment_index': idx,
                'mapping_status': mapping_status,
                'mouse_chrom': frag['chrom'],
                'mouse_start': frag['start'],
                'mouse_end': frag['end'],
                'original_human_bp': original_bp,
                'mapped_mouse_bp': mapped_bp,
                'fraction_mapped': round(fraction_mapped, 6),
                'reciprocal_overlap_fraction': round(reciprocal_overlap_fraction, 6),
                'unique_vs_split_mapping': mapping_status,
                'lead_snp_mapped': int(lead_snp_mapped),
                'lead_feature_available': int(lead_feature_available),
                'boundary_score': boundary_score,
                'unique_mapping_score': unique_mapping_score,
                'reciprocal_score': reciprocal_score,
                'lead_snp_score': lead_snp_score,
                'confidence_score': round(confidence, 6),
                'confidence_label': conf_label,
            }
            detailed_rows.append(detail_row)
        chosen_for_bed = merged if mapping_status in {'unique', 'merged_nearby'} else frags
        for frag in chosen_for_bed:
            final_bed_rows.append({
                'chrom': frag['chrom'],
                'start': frag['start'],
                'end': frag['end'],
                'name': name,
                'score': original['score'],
                'strand': '.',
            })

    with output_bed.open('w', encoding='utf-8') as handle:
        for row in final_bed_rows:
            handle.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['name']}\t{row['score']}\t{row['strand']}\n")
    fieldnames = sorted({k for row in detailed_rows for k in row.keys()}) if detailed_rows else ['feature_name']
    write_tsv(output_detail_tsv, detailed_rows, fieldnames)
    return mapped_count, unmapped_count


def ensure_build_conversion_reports(conversion_rows: List[Dict[str, object]], build_rows: List[Dict[str, object]]) -> None:
    write_tsv(REPORTS_DIR / 'build_detection.tsv', build_rows, list(build_rows[0].keys()) if build_rows else ['study_accession'])
    write_tsv(
        REPORTS_DIR / 'build_conversion.tsv',
        conversion_rows,
        list(conversion_rows[0].keys()) if conversion_rows else ['study_accession'],
    )


def write_download_manifest(download_rows: List[Dict[str, object]]) -> None:
    fieldnames = [
        'study_accession', 'category', 'trait_name', 'source_path', 'data_file_name',
        'download_status', 'file_format_guess', 'local_files_downloaded', 'file_size', 'source_type', 'notes'
    ]
    write_tsv(REPORTS_DIR / 'download_manifest.tsv', download_rows, fieldnames)


def top_confidence_rows(detail_path: Path, category: str, limit: int = 20) -> List[Dict[str, str]]:
    if not detail_path.exists():
        return []
    rows = read_tsv(detail_path)
    rows = [r for r in rows if r.get('category') == category]
    rows.sort(key=lambda r: float(r.get('confidence_score') or 0), reverse=True)
    best = []
    seen = set()
    for row in rows:
        key = row.get('feature_name')
        if key in seen:
            continue
        seen.add(key)
        best.append(row)
        if len(best) >= limit:
            break
    return best


def write_run_summary(summary: Dict[str, object]) -> None:
    lines = []
    lines.append('# GWAS human-to-mouse kidney/comparator pipeline run summary')
    lines.append('')
    lines.append(f"Starting TSV: {summary['starting_tsv']}")
    lines.append(f"Starting TSV rows: {summary['starting_rows']}")
    lines.append(f"Starting TSV interpretation: {summary['starting_interpretation']}")
    lines.append('')
    lines.append(f"Candidate studies retained: {summary['candidate_total']}")
    lines.append(f"Kidney studies retained: {summary['kidney_studies']}")
    lines.append(f"Comparator studies retained: {summary['comparator_studies']}")
    lines.append(f"Studies with full summary statistics: {summary['full_summary_count']}")
    lines.append(f"Index-only studies: {summary['index_only_count']}")
    lines.append('')
    lines.append(f"Primary SNPs generated: {summary['primary_snp_count']}")
    lines.append(f"Primary loci generated: {summary['primary_locus_count']}")
    lines.append(f"Mapped mm39 features: {summary['mapped_features']}")
    lines.append(f"Unmapped human features: {summary['unmapped_features']}")
    lines.append('')
    lines.append('Confidence score distribution:')
    for label, count in summary['confidence_distribution'].items():
        lines.append(f"- {label}: {count}")
    lines.append('')
    lines.append('Top 20 highest-confidence kidney loci:')
    for row in summary['top_kidney_loci']:
        lines.append(f"- {row.get('feature_name')} | confidence={row.get('confidence_score')} | trait={row.get('trait_name')} | study={row.get('study_accession')}")
    lines.append('')
    lines.append('Top 20 highest-confidence comparator loci:')
    for row in summary['top_comparator_loci']:
        lines.append(f"- {row.get('feature_name')} | confidence={row.get('confidence_score')} | trait={row.get('trait_name')} | study={row.get('study_accession')}")
    lines.append('')
    lines.append('Major caveats:')
    lines.append('- selected study files are preserved locally under downloads/<study_accession>/ and parsed from those local copies')
    lines.append(f'- all retained human coordinates are standardized internally to {TARGET_HUMAN_BUILD}/{TARGET_HUMAN_UCSC} before human BED generation')
    lines.append('- GRCh37 selected files without harmonised GRCh38 input are lifted to hg38 feature intervals before downstream BED generation and documented in build_conversion.tsv')
    lines.append('- interval/CNV rows are retained as exact intervals rather than collapsing them to 1 bp SNPs')
    lines.append('- missing or readme-only studies are logged and skipped rather than crashing the run')
    lines.append('- cross-species liftOver maps coordinates, not disease causality')
    (REPORTS_DIR / 'run_summary.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')


def main() -> int:
    ensure_dirs()
    ensure_liftover_assets()

    selected_rows = read_tsv(SELECTED_MANIFEST)
    inventory_rows = {r['study_accession']: r for r in read_tsv(INVENTORY_TSV)}

    build_rows: List[Dict[str, object]] = []
    build_conversion_rows: List[Dict[str, object]] = []
    all_significant_rows: List[Dict[str, object]] = []
    download_manifest_rows: List[Dict[str, object]] = []

    for study_row in selected_rows:
        accession = study_row['STUDY ACCESSION']
        inventory_row = inventory_rows.get(accession, {})
        local_context = download_study_artifacts(study_row, inventory_row)
        build_row = build_detection(study_row, inventory_row, local_context)
        build_rows.append(build_row)
        significant_rows, download_row, conversion_row = extract_significant_rows(study_row, inventory_row, build_row, local_context)
        download_manifest_rows.append(download_row)
        build_conversion_rows.append({
            'study_accession': accession,
            'category': study_row['category'],
            'trait_name': study_row['DISEASE/TRAIT'],
            **conversion_row,
        })
        all_significant_rows.extend(significant_rows)
        log(f"{accession}: kept {len(significant_rows)} significant rows")

    ensure_build_conversion_reports(build_conversion_rows, build_rows)
    write_download_manifest(download_manifest_rows)

    primary_rows, secondary_rows = make_primary_secondary_tables(all_significant_rows)
    save_variant_tables(all_significant_rows, primary_rows, secondary_rows)

    loci = group_loci(primary_rows)
    save_locus_table(loci)
    split_categories(primary_rows, loci)

    unmapped_report: List[Dict[str, object]] = []
    multi_report: List[Dict[str, object]] = []

    files_to_map = [
        (HUMAN_BED_DIR / 'kidney_lead_snps_human.bed', HUMAN_BED_DIR / 'kidney_lead_snps_human.tsv', MM39_BED_DIR / 'kidney_lead_snps_mm39.bed', MM39_BED_DIR / 'kidney_lead_snps_mm39_detailed.tsv', False),
        (HUMAN_BED_DIR / 'kidney_loci_human.bed', HUMAN_BED_DIR / 'kidney_loci_human.tsv', MM39_BED_DIR / 'kidney_loci_mm39.bed', MM39_BED_DIR / 'kidney_loci_mm39_detailed.tsv', True),
        (HUMAN_BED_DIR / 'comparator_lead_snps_human.bed', HUMAN_BED_DIR / 'comparator_lead_snps_human.tsv', MM39_BED_DIR / 'comparator_lead_snps_mm39.bed', MM39_BED_DIR / 'comparator_lead_snps_mm39_detailed.tsv', False),
        (HUMAN_BED_DIR / 'comparator_loci_human.bed', HUMAN_BED_DIR / 'comparator_loci_human.tsv', MM39_BED_DIR / 'comparator_loci_mm39.bed', MM39_BED_DIR / 'comparator_loci_mm39_detailed.tsv', True),
    ]

    mapped_total = 0
    unmapped_total = 0
    for input_bed, detail_tsv, output_bed, output_detail_tsv, is_locus in files_to_map:
        mapped, unmapped = annotate_mm39(input_bed, detail_tsv, output_bed, output_detail_tsv, unmapped_report, multi_report, is_locus)
        mapped_total += mapped
        unmapped_total += unmapped
        log(f'lifted {input_bed.name}: mapped={mapped} unmapped={unmapped}')

    write_tsv(REPORTS_DIR / 'unmapped_human_features.tsv', unmapped_report, list(unmapped_report[0].keys()) if unmapped_report else ['feature_name'])
    write_tsv(REPORTS_DIR / 'multi_mapped_features.tsv', multi_report, list(multi_report[0].keys()) if multi_report else ['feature_name'])

    confidence_distribution = {'high': 0, 'medium': 0, 'low': 0, 'unmapped': 0}
    for detail_path in [
        MM39_BED_DIR / 'kidney_loci_mm39_detailed.tsv',
        MM39_BED_DIR / 'comparator_loci_mm39_detailed.tsv',
        MM39_BED_DIR / 'kidney_lead_snps_mm39_detailed.tsv',
        MM39_BED_DIR / 'comparator_lead_snps_mm39_detailed.tsv',
    ]:
        if detail_path.exists():
            for row in read_tsv(detail_path):
                label = row.get('confidence_label', 'unmapped')
                if label in confidence_distribution:
                    confidence_distribution[label] += 1

    inspection_path = REPORTS_DIR / 'inspection_summary.json'
    if inspection_path.exists():
        inspection_summary = json.loads(inspection_path.read_text(encoding='utf-8'))
        starting_rows = inspection_summary.get('row_count', 0)
    else:
        starting_rows = 0
    summary = {
        'starting_tsv': str(OUTPUT_ROOT / 'downloads' / 'unpublished_studies_v1.0.3.1.tsv'),
        'starting_rows': starting_rows,
        'starting_interpretation': 'study metadata / study manifest with summary-statistics locations',
        'candidate_total': len(selected_rows),
        'kidney_studies': sum(1 for r in selected_rows if r['category'] == 'kidney'),
        'comparator_studies': sum(1 for r in selected_rows if r['category'] == 'comparator'),
        'full_summary_count': sum(1 for r in selected_rows if r.get('summary_stats_available') == 'yes'),
        'index_only_count': 0,
        'primary_snp_count': len(primary_rows),
        'primary_locus_count': len(loci),
        'mapped_features': mapped_total,
        'unmapped_features': unmapped_total,
        'confidence_distribution': confidence_distribution,
        'top_kidney_loci': top_confidence_rows(MM39_BED_DIR / 'kidney_loci_mm39_detailed.tsv', 'kidney', 20),
        'top_comparator_loci': top_confidence_rows(MM39_BED_DIR / 'comparator_loci_mm39_detailed.tsv', 'comparator', 20),
    }
    write_run_summary(summary)
    log('pipeline complete')
    return 0


if __name__ == '__main__':
    sys.exit(main())
