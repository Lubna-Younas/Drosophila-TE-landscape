#!/usr/bin/env python3
"""One sbatch call per library. Default is a metadata-only dry run."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import uuid
from datetime import datetime, timezone
from mapping_worker import metadata, reference_check

HERE = Path(__file__).resolve().parent


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--submit", action="store_true", help="Actually call sbatch; otherwise only inspect/print")
    p.add_argument("--sample", action="append", help="Limit to this exact sample ID (repeatable)")
    p.add_argument("--layout", choices=["PE", "SE"], help="Optional restriction to one read layout")
    args = p.parse_args()
    c = json.loads((HERE / "config.json").read_text())
    with (HERE / "samples.tsv").open() as f:
        samples = list(csv.DictReader(f, delimiter="\t"))
    known = {s["sample_id"] for s in samples}
    if args.sample and not set(args.sample).issubset(known):
        raise ValueError("Unknown sample ID: " + str(set(args.sample) - known))
    samples = [s for s in samples if (not args.sample or s["sample_id"] in args.sample)
               and (not args.layout or s["layout"] == args.layout)]
    if not samples:
        raise ValueError("No selected samples")
    seen = set()
    for row in samples:
        sample = json.loads((HERE / "samples" / (row["sample_id"] + ".json")).read_text())
        for item in sample["input_metadata"]:
            current = metadata(item["path"])
            if current != item:
                raise ValueError("Input changed since inventory: " + item["path"])
            if current["path"] in seen:
                raise ValueError("Input assigned to more than one selected library")
            seen.add(current["path"])
    reference_ready = True
    try:
        reference_check(c)
    except (OSError, ValueError, KeyError) as exc:
        reference_ready = False
        print("Reference not ready:", exc)
    print("Selected libraries:", len(samples), "| FASTQ files:", len(seen))
    if not args.submit:
        for row in samples:
            print("sbatch", HERE / "jobs" / (row["sample_id"] + ".sbatch"))
        print("DRY RUN: no jobs submitted.")
        if not reference_ready:
            print("Build the reference first: sbatch prepare_reference.sbatch")
        return 0
    if not reference_ready:
        raise ValueError("Reference must finish and validate before independent sample jobs can be submitted")
    if not shutil.which("sbatch"):
        raise ValueError("sbatch is unavailable")
    # Freeze a copy of worker/config/environment so later edits cannot alter queued jobs.
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S") + "_" + uuid.uuid4().hex[:8]
    frozen = HERE / "submissions" / stamp
    frozen.mkdir(parents=True)
    for name in ("mapping_worker.py", "environment.sh", "config.json"):
        shutil.copy2(HERE / name, frozen / name)
    source_hashes = {name: hashlib.sha256((frozen / name).read_bytes()).hexdigest()
                     for name in ("mapping_worker.py", "environment.sh", "config.json")}
    (frozen / "code_sha256.json").write_text(json.dumps(source_hashes, indent=2) + "\n")
    reservations = HERE / "submission_records"
    reservations.mkdir(exist_ok=True)
    for row in samples:
        sid = row["sample_id"]
        if (Path(c["results_dir"]) / sid).exists():
            print("SKIP existing result:", sid)
            continue
        record = reservations / (sid + ".json")
        try:
            fd = os.open(record, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        except FileExistsError:
            print("SKIP already submitted/reserved:", sid)
            continue
        with os.fdopen(fd, "w") as f:
            json.dump({"state": "SUBMISSION_PENDING", "sample_id": sid, "snapshot": str(frozen)}, f)
        sample_file = frozen / (sid + ".json")
        shutil.copy2(HERE / "samples" / (sid + ".json"), sample_file)
        script = (HERE / "jobs" / (sid + ".sbatch")).read_text()
        script = script.replace(str(HERE / "environment.sh"), str(frozen / "environment.sh"))
        script = script.replace(str(HERE / "mapping_worker.py"), str(frozen / "mapping_worker.py"))
        script = script.replace(str(HERE / "config.json"), str(frozen / "config.json"))
        script = script.replace(str(HERE / "samples" / (sid + ".json")), str(sample_file))
        job_file = frozen / (sid + ".sbatch")
        job_file.write_text(script)
        result = subprocess.run(["sbatch", "--parsable", str(job_file)], text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        response = result.stdout.strip()
        if result.returncode != 0 or not re.fullmatch(r"\d+(;[A-Za-z0-9_.-]+)?", response):
            record.write_text(json.dumps({"state": "SUBMISSION_FAILED_OR_UNCERTAIN", "sample_id": sid,
                              "stdout": response, "stderr": result.stderr, "snapshot": str(frozen)}, indent=2))
            raise RuntimeError("Submission failed or is uncertain; inspect " + str(record))
        job_id = response.split(";")[0]
        record.write_text(json.dumps({"state": "SUBMITTED", "sample_id": sid, "job_id": job_id,
                          "snapshot": str(frozen), "script": str(job_file)}, indent=2) + "\n")
        print(sid, job_id, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
