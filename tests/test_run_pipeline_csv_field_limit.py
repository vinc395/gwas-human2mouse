import importlib.util
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve().parents[1] / 'scripts' / 'run_pipeline.py'
spec = importlib.util.spec_from_file_location('run_pipeline', SCRIPT_PATH)
run_pipeline = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(run_pipeline)


def test_read_tsv_handles_large_fields(tmp_path):
    large_value = 'x' * 200_000
    input_tsv = tmp_path / 'large.tsv'
    input_tsv.write_text(
        'bed_name\tpayload\n'
        f'example\t{large_value}\n',
        encoding='utf-8',
    )

    rows = run_pipeline.read_tsv(input_tsv)

    assert rows == [{'bed_name': 'example', 'payload': large_value}]
