#!/usr/bin/env python3
"""Read receipts and file metadata only; safe for a login-node status check."""
import csv
import json
from pathlib import Path
import sys

root = Path(__file__).resolve().parent
config = json.loads((root / "config.json").read_text())
with (root / "samples.tsv").open() as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
writer = csv.writer(sys.stdout, delimiter="\t")
writer.writerow(["sample_id", "run_accession", "layout", "state", "job_id", "strand", "result_directory"])
for row in rows:
    directory = Path(config["results_dir"]) / row["sample_id"]
    receipt = directory / "receipt.json"
    state, job, strand = "NO_FINAL_RECEIPT", "", ""
    submission = root / "submission_records" / (row["sample_id"] + ".json")
    if submission.exists():
        sub = json.loads(submission.read_text())
        state, job = sub["state"], sub.get("job_id", "")
    if receipt.exists():
        r = json.loads(receipt.read_text())
        state, job = r["state"], r["slurm_job_id"]
        strand = r.get("strandedness", {}).get("used_for_RSEM", "")
        for item in r["files"]:
            p = directory / item["relative_path"]
            if not p.is_file() or p.stat().st_size != item["bytes"] or p.stat().st_mtime_ns != item["mtime_ns"]:
                state = "OUTPUT_CHANGED_REVALIDATION_REQUIRED"
                break
    writer.writerow([row["sample_id"], row["run_accession"], row["layout"], state, job, strand, str(directory)])
