#!/usr/bin/env python3
import csv
import re
import sys
import time
from html.parser import HTMLParser
from pathlib import Path
from typing import Dict, List, Optional
from urllib.parse import urljoin

import requests

OUTPUT_ROOT = Path('/home/vince395/gwas_human2mouse')
SELECTED_MANIFEST = OUTPUT_ROOT / 'manifests' / 'selected_kidney_and_comparator_studies.tsv'
REPORTS_DIR = OUTPUT_ROOT / 'reports'
LOGS_DIR = OUTPUT_ROOT / 'logs'

SIZE_RE = re.compile(r'^(?P<num>\d+(?:\.\d+)?)(?P<unit>[KMGTP]?)$')


def size_to_bytes(size_text: str) -> Optional[int]:
    size_text = size_text.strip()
    m = SIZE_RE.match(size_text)
    if not m:
        return None
    num = float(m.group('num'))
    unit = m.group('unit')
    mult = {
        '': 1,
        'K': 1024,
        'M': 1024**2,
        'G': 1024**3,
        'T': 1024**4,
        'P': 1024**5,
    }[unit]
    return int(num * mult)


class ApacheIndexParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.in_link = False
        self.current_href: Optional[str] = None
        self.current_text_parts: List[str] = []
        self.entries: List[Dict[str, str]] = []
        self.capture_size_next = False
        self.last_href_text = None
        self._td_depth = 0
        self._current_row: List[str] = []
        self._in_row = False
        self._in_td = False

    def handle_starttag(self, tag, attrs):
        if tag == 'tr':
            self._in_row = True
            self._current_row = []
        elif tag == 'td' and self._in_row:
            self._in_td = True
            self._td_depth += 1
            self.current_text_parts = []
            self.current_href = None
        elif tag == 'a' and self._in_td:
            self.in_link = True
            self.current_text_parts = []
            attrs_dict = dict(attrs)
            self.current_href = attrs_dict.get('href')

    def handle_endtag(self, tag):
        if tag == 'a' and self.in_link:
            self.in_link = False
        elif tag == 'td' and self._in_td:
            text = ''.join(self.current_text_parts).strip()
            self._current_row.append(text)
            self._in_td = False
            self._td_depth = max(0, self._td_depth - 1)
            self.current_text_parts = []
        elif tag == 'tr' and self._in_row:
            if len(self._current_row) >= 4:
                # Apache index rows: icon, name, modified, size, description
                name = self._current_row[1] if len(self._current_row) > 1 else ''
                modified = self._current_row[2] if len(self._current_row) > 2 else ''
                size = self._current_row[3] if len(self._current_row) > 3 else ''
                href_match = re.search(r'href="([^"]+)"', ''.join(self.get_starttag_text() or ''))
                self.entries.append({'name': name, 'modified': modified, 'size': size})
            self._in_row = False

    def handle_data(self, data):
        if self._in_td or self.in_link:
            self.current_text_parts.append(data)


def parse_directory_listing(html: str) -> List[Dict[str, str]]:
    entries: List[Dict[str, str]] = []
    row_re = re.compile(r'<tr>(.*?)</tr>', flags=re.S | re.I)
    cell_re = re.compile(r'<t[dh][^>]*>(.*?)</t[dh]>', flags=re.S | re.I)
    link_re = re.compile(r'href="([^"]+)"', flags=re.I)
    tag_re = re.compile(r'<[^>]+>')
    for row_html in row_re.findall(html):
        cells = cell_re.findall(row_html)
        if len(cells) < 4:
            continue
        name_cell = cells[1]
        href_match = link_re.search(name_cell)
        href = href_match.group(1) if href_match else ''
        text_cells = [re.sub(r'\s+', ' ', tag_re.sub('', c)).strip() for c in cells]
        name = text_cells[1] if len(text_cells) > 1 else ''
        if name in {'Parent Directory', ''}:
            continue
        entries.append({
            'href': href,
            'name': name,
            'modified': text_cells[2] if len(text_cells) > 2 else '',
            'size': text_cells[3] if len(text_cells) > 3 else '',
        })
    return entries


def load_selected_rows() -> List[Dict[str, str]]:
    with SELECTED_MANIFEST.open('r', encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def choose_data_file(entries: List[Dict[str, str]]) -> Dict[str, str]:
    harmonised = [e for e in entries if e['name'].endswith('.h.tsv.gz')]
    if harmonised:
        # Prefer harmonised compressed file.
        return dict(harmonised[0], chosen_type='harmonised')
    compressed = [e for e in entries if re.search(r'\.(tsv|txt|csv)\.gz$', e['name'])]
    if compressed:
        return dict(compressed[0], chosen_type='compressed_table')
    plain = [e for e in entries if re.search(r'\.(tsv|txt|csv)$', e['name'])]
    if plain:
        return dict(plain[0], chosen_type='plain_table')
    return {'name': '', 'href': '', 'size': '', 'modified': '', 'chosen_type': 'missing'}


def main() -> int:
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    rows = load_selected_rows()
    out_path = REPORTS_DIR / 'selected_study_inventory.tsv'
    log_path = LOGS_DIR / 'selected_study_inventory.log'

    session = requests.Session()
    session.headers['User-Agent'] = 'gwas-human2mouse-pipeline/1.0'

    fieldnames = [
        'study_accession', 'category', 'trait_name', 'source_path', 'root_status',
        'data_file_name', 'data_file_type', 'data_file_size_text', 'data_file_size_bytes',
        'harmonised_dir_present', 'meta_files', 'available_entries', 'notes'
    ]
    records: List[Dict[str, str]] = []
    lines: List[str] = []

    for idx, row in enumerate(rows, start=1):
        accession = row['STUDY ACCESSION']
        source_path = row['summary_stats_location']
        directory_url = source_path.rstrip('/') + '/'
        record = {
            'study_accession': accession,
            'category': row['category'],
            'trait_name': row['DISEASE/TRAIT'],
            'source_path': source_path,
            'root_status': '',
            'data_file_name': '',
            'data_file_type': '',
            'data_file_size_text': '',
            'data_file_size_bytes': '',
            'harmonised_dir_present': 'no',
            'meta_files': '',
            'available_entries': '',
            'notes': '',
        }
        try:
            response = session.get(directory_url.replace('ftp://', 'https://'), timeout=60)
            record['root_status'] = str(response.status_code)
            response.raise_for_status()
            entries = parse_directory_listing(response.text)
            record['available_entries'] = '; '.join(e['name'] for e in entries)
            record['harmonised_dir_present'] = 'yes' if any(e['name'] == 'harmonised/' for e in entries) else 'no'
            meta_files = [e['name'] for e in entries if e['name'].endswith('.yaml') or e['name'] == 'md5sum.txt']
            if record['harmonised_dir_present'] == 'yes':
                harm_url = directory_url.replace('ftp://', 'https://') + 'harmonised/'
                harm_resp = session.get(harm_url, timeout=60)
                harm_resp.raise_for_status()
                harm_entries = parse_directory_listing(harm_resp.text)
                entries.extend([dict(e, name='harmonised/' + e['name']) for e in harm_entries])
                meta_files.extend('harmonised/' + e['name'] for e in harm_entries if e['name'].endswith('.yaml') or e['name'] == 'md5sum.txt')
                chosen = choose_data_file([dict(e, href=urljoin(harm_url, e['href'])) for e in harm_entries])
                chosen_name = 'harmonised/' + chosen['name'] if chosen.get('name') else ''
            else:
                chosen = choose_data_file([dict(e, href=urljoin(directory_url.replace('ftp://', 'https://'), e['href'])) for e in entries])
                chosen_name = chosen.get('name', '')
            record['meta_files'] = '; '.join(meta_files)
            record['data_file_name'] = chosen_name
            record['data_file_type'] = chosen.get('chosen_type', '')
            record['data_file_size_text'] = chosen.get('size', '')
            size_bytes = size_to_bytes(chosen.get('size', ''))
            record['data_file_size_bytes'] = '' if size_bytes is None else str(size_bytes)
            if not chosen_name:
                record['notes'] = 'no parsable data file found'
        except Exception as exc:
            record['notes'] = f'ERROR: {exc}'
        records.append(record)
        lines.append(f'[{idx}/{len(rows)}] {accession} -> {record["data_file_name"] or "NONE"} ({record["data_file_size_text"]}) status={record["root_status"]} notes={record["notes"]}')
        print(lines[-1])
        time.sleep(0.1)

    with out_path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(records)
    log_path.write_text('\n'.join(lines) + '\n', encoding='utf-8')

    total_bytes = sum(int(r['data_file_size_bytes']) for r in records if r['data_file_size_bytes'])
    print(f'Estimated chosen data total bytes: {total_bytes}')
    print(f'Estimated chosen data total GB: {total_bytes / (1024**3):.2f}')
    print(f'Wrote {out_path}')
    print(f'Wrote {log_path}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
