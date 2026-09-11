"""Keep a small record of the inputs and settings behind one pipeline run."""

import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

# Only these generated files are replaced. Raw inputs and other files are left alone.
RUN_FILES = (
    'ranked_candidates.csv', 'ranked_candidates.json', 'pipeline_summary.png',
    'validation_report.json', 'run_manifest.json', 'run_manifest.json.tmp',
)


def file_sha256(path: Path) -> str:
    """Fingerprint a file without loading the whole count table into memory."""
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _write_record(record: dict, output_dir: Path) -> None:
    # Replace the record only after the new JSON has been written successfully.
    temporary = output_dir / 'run_manifest.json.tmp'
    temporary.write_text(json.dumps(record, indent=2) + '\n')
    temporary.replace(output_dir / 'run_manifest.json')


def start_run_record(config: dict, output_dir: Path) -> dict:
    """Record this run and remove old generated results before analysis starts."""
    counts_path = Path(config['selex']['counts_file']).resolve()
    input_paths = [counts_path]
    if config['target']['input_type'] == 'fasta':
        input_paths.append(Path(config['target']['input_value']).resolve())
    output_paths = [output_dir / name for name in RUN_FILES]
    reserved_paths = {path.resolve() for path in output_paths}
    if any(path in reserved_paths for path in input_paths):
        raise ValueError('An input uses a reserved output filename in this folder. Tip: choose a different output folder.')

    inputs = [{'path': str(path)} for path in input_paths]
    source_dir = Path(__file__).resolve().parent
    source_hashes = {path.name: file_sha256(path) for path in sorted(source_dir.glob('*.py'))}
    revision, dirty = None, None
    try:
        revision = subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=source_dir.parent, stderr=subprocess.DEVNULL, text=True, timeout=5,
        ).strip()
        dirty = bool(subprocess.check_output(
            ['git', 'status', '--porcelain', '--untracked-files=no'], cwd=source_dir.parent,
            stderr=subprocess.DEVNULL, text=True, timeout=5,
        ).strip())
    except (OSError, subprocess.SubprocessError):
        pass  # A downloaded ZIP has no Git history; source fingerprints still identify the code.

    record = {
        'status': 'started',
        'started_at_utc': datetime.now(timezone.utc).isoformat(),
        'config': config,
        'inputs': inputs,
        'python_version': sys.version,
        'git_revision': revision,
        'git_dirty': dirty,
        'source_sha256': source_hashes,
    }
    for path in output_paths:
        path.unlink(missing_ok=True)
    _write_record(record, output_dir)
    for entry in inputs:
        entry['sha256'] = file_sha256(Path(entry['path']))
    _write_record(record, output_dir)
    return record


def finish_run_record(record: dict, output_dir: Path, rounds: list[str], candidate_count: int) -> None:
    """Mark completion only after the requested exports and plots have finished."""
    record.update(
        status='complete',
        completed_at_utc=datetime.now(timezone.utc).isoformat(),
        rounds=rounds,
        shortlisted_candidates=candidate_count,
        outputs=[name for name in RUN_FILES if name != 'run_manifest.json' and (output_dir / name).exists()],
    )
    _write_record(record, output_dir)
