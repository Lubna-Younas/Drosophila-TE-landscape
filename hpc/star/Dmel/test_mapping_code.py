#!/usr/bin/env python3
"""Local software tests with tiny synthetic inputs; no scientific tools/jobs run."""
import csv
import gzip
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import mapping_worker as m
import submit_samples as submit

HERE = Path(__file__).resolve().parent


class MappingTests(unittest.TestCase):
    def test_manifest_coverage_and_no_pooling(self):
        with (HERE / "samples.tsv").open() as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual(len(rows), 50)
        self.assertEqual(sum(r["layout"] == "PE" for r in rows), 38)
        used = []
        for row in rows:
            sample = json.loads((HERE / "samples" / (row["sample_id"] + ".json")).read_text())
            self.assertEqual(sample["sample_id"].rsplit("_R", 1)[1], row["biological_replicate"])
            self.assertEqual(sample["input_metadata"][0]["path"], row["read1"])
            used.extend(x["path"] for x in sample["input_metadata"])
            if row["layout"] == "PE":
                self.assertEqual(row["read1"].replace(".fastq1.gz", ".fastq2.gz"), row["read2"])
        inventory = json.loads((HERE / "input_inventory.json").read_text())
        expected = {r["path"] for r in inventory["files"] if r["path"].endswith((".fastq.gz", ".fastq1.gz", ".fastq2.gz"))}
        self.assertEqual(len(used), 88)
        self.assertEqual(len(used), len(set(used)))
        self.assertEqual(set(used), expected)

    def test_slurm_scripts_are_independent(self):
        scripts = [p for p in (HERE / "jobs").glob("*.sbatch") if not p.name.startswith(".")]
        self.assertEqual(len(scripts), 50)
        for script in scripts + [HERE / "prepare_reference.sbatch"]:
            text = script.read_text()
            for forbidden in ("--array", "--dependency", "--partition", "--nodelist", "--nodes", "--exclusive", "--constraint"):
                self.assertNotIn(forbidden, text)
            self.assertIn("--cpus-per-task=8", text)
            self.assertIn("--mem=32G", text)
            self.assertIn("--time=", text)
            self.assertIn("--account=YOUR_ACCOUNT", text)
            subprocess.run(["bash", "-n", str(script)], check=True)

    def test_annotation_compatibility_preserves_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            a = '3R\tFlyBase\tmRNA\t1\t40\t.\t.\t.\tgene_id "g"; transcript_id "bad";\n'
            valid = '3R\tFlyBase\texon\t1\t10\t.\t+\t.\tgene_id "g"; transcript_id "good";\n'
            b = valid.replace('"good"', '"bad"')
            c = b.replace('\t+\t', '\t-\t')
            src = root / 'source.gtf'; src.write_text(a + valid + b + c)
            audit = m.prepare_rsem_gtf(src, root / 'genes.gtf', ['bad'])
            self.assertEqual((root / 'genes.gtf').read_text(), valid)
            self.assertEqual(src.read_text(), a + valid + b + c)
            self.assertEqual(audit['retained_exons'], 1)
            self.assertEqual((root / 'excluded_transcript_exons.gtf').read_text(), b + c)
            with self.assertRaises(ValueError):
                m.prepare_rsem_gtf(src, root / 'unexpected.gtf', [])

    def test_strand_classification_and_ambiguity(self):
        self.assertEqual(m.strand_call(19000, 1000), "forward")
        self.assertEqual(m.strand_call(1000, 19000), "reverse")
        self.assertEqual(m.strand_call(10000, 10000), "none")
        self.assertEqual(m.strand_call(15000, 5000), "ambiguous")
        self.assertEqual(m.strand_call(900, 1), "ambiguous")

    def test_spliced_alignment_blocks(self):
        self.assertEqual(list(m.alignment_blocks(100, "5S20M100N10M2D5M1I4M")),
                         [(100, 120), (220, 230), (232, 237), (237, 241)])

    def test_read_pair_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            a, b = Path(tmp) / "1.gz", Path(tmp) / "2.gz"
            for p, mate in ((a, 1), (b, 2)):
                with gzip.open(p, "wt") as f:
                    f.write(f"@read/{mate}\nACGT\n+\nIIII\n")
            self.assertEqual(m.check_fastq([a, b], "PE")["checked_records_per_file"], 1)
            with gzip.open(b, "wt") as f:
                f.write("@different/2\nACGT\n+\nIIII\n")
            with self.assertRaises(ValueError):
                m.check_fastq([a, b], "PE")
            with gzip.open(a, "wt") as f:
                f.write("@read\nACGT\n+\nI\n")
            with self.assertRaises(ValueError):
                m.check_fastq([a], "SE")

    def test_no_login_node_execution(self):
        env = dict(os.environ)
        env.pop("SLURM_JOB_ID", None)
        p = subprocess.run(["python3", str(HERE / "mapping_worker.py"), "reference", "--config", str(HERE / "config.json")], env=env, capture_output=True, text=True)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("compute node", p.stderr)

    def test_publication_checksums_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp, patch.dict(os.environ, {"SLURM_JOB_ID": "123"}):
            root = Path(tmp); source = root / "scratch"; source.mkdir()
            p = source / "count.tsv"; p.write_text("gene\tcount\na\t10\n")
            dest = root / "final"
            m.publish(source, dest, [p], {"state": "TEST"}, "receipt.json")
            r = json.loads((dest / "receipt.json").read_text())
            self.assertEqual(r["files"][0]["sha256"], m.digest(dest / "count.tsv"))
            with self.assertRaises(FileExistsError):
                m.publish(source, dest, [p], {}, "receipt.json")

    def test_genomic_and_transcript_branches(self):
        rows = [json.loads(p.read_text()) for p in (HERE / "samples").glob("*.json") if not p.name.startswith(".")]
        for layout in ("PE", "SE"):
            sample = next(s for s in rows if s["layout"] == layout)
            cmd = m.star_command({"reference_dir": "/ref"}, sample, Path("/scratch"), 8)
            i = cmd.index("--readFilesIn")
            j = cmd.index("--readFilesCommand")
            self.assertEqual(j - i - 1, 2 if layout == "PE" else 1)
            self.assertIn("TranscriptomeSAM", cmd)
            self.assertEqual(cmd[cmd.index("--outFilterMultimapNmax") + 1], "100")
            self.assertEqual(cmd[cmd.index("--outSAMmultNmax") + 1], "100")

    def test_submission_dry_run_and_duplicate_prevention(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in ("samples", "jobs", "logs"):
                (root / name).mkdir()
            raw = root / "raw.fastq.gz"; raw.write_bytes(b"synthetic input")
            c = {"results_dir": str(root / "results")}
            (root / "config.json").write_text(json.dumps(c))
            (root / "mapping_worker.py").write_text("# test worker\n")
            (root / "environment.sh").write_text("# test environment\n")
            row = {"sample_id": "rep1", "layout": "SE"}
            with (root / "samples.tsv").open("w") as f:
                w = csv.DictWriter(f, fieldnames=list(row), delimiter="\t"); w.writeheader(); w.writerow(row)
            (root / "samples/rep1.json").write_text(json.dumps({**row, "input_metadata": [m.metadata(raw)]}))
            (root / "jobs/rep1.sbatch").write_text("#!/bin/bash\npython3 " + str(root / "mapping_worker.py") + "\n")
            fake = subprocess.CompletedProcess([], 0, "987654\n", "")
            with patch.object(submit, "HERE", root), patch.object(submit, "reference_check", return_value={}), patch.object(submit.shutil, "which", return_value="/bin/sbatch"), patch.object(submit.subprocess, "run", return_value=fake) as sbatch:
                with patch.object(submit.sys, "argv", ["submit_samples.py"]):
                    submit.main()
                sbatch.assert_not_called()
                with patch.object(submit.sys, "argv", ["submit_samples.py", "--submit"]):
                    submit.main(); submit.main()
                self.assertEqual(sbatch.call_count, 1)
                self.assertEqual(sbatch.call_args.args[0][:2], ["sbatch", "--parsable"])
                record = json.loads((root / "submission_records/rep1.json").read_text())
                frozen = Path(record["snapshot"])
                (root / "mapping_worker.py").write_text("# later edit\n")
                self.assertEqual((frozen / "mapping_worker.py").read_text(), "# test worker\n")


if __name__ == "__main__":
    unittest.main(verbosity=2)
