#!/usr/bin/env python3
"""Dpse annotation-free STAR two-pass genomic alignment. Slurm workers only."""
import argparse
import collections
import csv
import gzip
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import shutil
import socket
import subprocess
import sys
import tempfile
from datetime import datetime, timezone


def now():
    return datetime.now(timezone.utc).isoformat()


def save(path, data):
    Path(path).write_text(json.dumps(data, indent=2) + "\n")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def metadata(path):
    p = Path(path)
    s = p.stat()
    if not p.is_file() or s.st_size == 0:
        raise ValueError("Missing/empty input: " + str(p))
    return {"path": str(p.resolve()), "bytes": s.st_size, "mtime_ns": s.st_mtime_ns}


def unchanged(items):
    for old in items:
        if metadata(old["path"]) != {k: old[k] for k in ("path", "bytes", "mtime_ns")}:
            raise ValueError("Input metadata changed: " + old["path"])


def run(args, work, log, stdout=None):
    args = list(map(str, args))
    with (work / "commands.jsonl").open("a") as f:
        f.write(json.dumps({"time": now(), "argv": args}) + "\n")
    with (work / log).open("ab") as err:
        if stdout:
            with Path(stdout).open("wb") as out:
                subprocess.run(args, cwd=work, stdout=out, stderr=err, check=True)
        else:
            subprocess.run(args, cwd=work, stdout=err, stderr=err, check=True)


def versions():
    out = {}
    for name in ("STAR", "samtools", "rsem-calculate-expression"):
        exe = shutil.which(name)
        if not exe:
            raise ValueError("Required executable absent: " + name)
        p = subprocess.run([exe, "--version"], text=True, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, timeout=30, check=True)
        out[name] = {"executable": exe, "version": p.stdout.strip()}
    prep = shutil.which("rsem-prepare-reference")
    if not prep or Path(prep).parent != Path(out["rsem-calculate-expression"]["executable"]).parent:
        raise ValueError("RSEM reference builder must come from the same installation")
    # RSEM 1.3.1's reference builder has no --version option.
    out["rsem-prepare-reference"] = {"executable": prep,
        "version": out["rsem-calculate-expression"]["version"],
        "version_source": "rsem-calculate-expression in the same bin directory"}
    return out


def unpack(src, dest):
    opener = gzip.open if str(src).endswith(".gz") else open
    with opener(src, "rb") as inp, Path(dest).open("wb") as out:
        shutil.copyfileobj(inp, out, 8 * 1024 * 1024)


def normalize_fastq_compression(paths, depths, work):
    """Remove only an audited extra gzip wrapper, keeping sequences unchanged."""
    if len(paths) != len(depths):
        raise ValueError("Missing per-input compression audit")
    normalized, evidence = [], []
    for number, (source, expected_depth) in enumerate(zip(paths, depths), 1):
        with gzip.open(source, "rb") as outer:
            observed_depth = 2 if outer.peek(2)[:2] == b"\x1f\x8b" else 1
            if observed_depth != expected_depth or expected_depth not in (1, 2):
                raise ValueError("FASTQ compression differs from audited structure")
            if expected_depth == 1:
                normalized.append(source)
                evidence.append({"source": str(source), "gzip_layers": 1, "transformation": "none"})
                continue
            folder = work / "normalized_inputs"
            folder.mkdir(exist_ok=True)
            target = folder / ("read" + str(number) + ".fastq.gz")
            checksum = hashlib.sha256()
            with target.open("xb") as out:
                for block in iter(lambda: outer.read(8 * 1024 * 1024), b""):
                    out.write(block); checksum.update(block)
            with gzip.open(target, "rb") as inner:
                if not inner.peek(1).startswith(b"@"):
                    raise ValueError("Removing one extra gzip layer did not expose FASTQ")
            normalized.append(str(target))
            evidence.append({"source": str(source), "gzip_layers": 2,
                "transformation": "outer gzip wrapper removed into job-private scratch; sequence data unchanged",
                "derived_gzip_sha256": checksum.hexdigest(), "derived_bytes": target.stat().st_size})
    return normalized, evidence


def check_fastq(paths, layout, limit=None):
    """Streaming FASTQ/gzip/mate-ID validation, always on a compute node.

    Default is the complete input, including gzip CRC and all paired read IDs.
    Tests may pass a finite limit; production calls must not do so.
    """
    def records(path):
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as f:
            count = 0
            while limit is None or count < limit:
                h = f.readline()
                if h == "":
                    break
                h = h.rstrip()
                seq, plus, qual = (f.readline().rstrip() for _ in range(3))
                if not h.startswith("@") or not plus.startswith("+") or not seq or len(seq) != len(qual):
                    raise ValueError("Malformed FASTQ: " + str(path))
                token = re.sub(r"/[12]$", "", h.split()[0][1:])
                count += 1
                yield token
    count = 0
    if layout == "PE":
        for a, b in itertools.zip_longest(records(paths[0]), records(paths[1])):
            if a is None or b is None or a != b:
                raise ValueError("Mate read IDs/counts disagree in FASTQ records")
            count += 1
    else:
        count = sum(1 for _ in records(paths[0]))
    if count == 0:
        raise ValueError("No FASTQ reads")
    return {"checked_records_per_file": count, "layout": layout,
            "complete_file_check": limit is None}


def exon_index(gtf):
    bins = collections.defaultdict(list)
    with Path(gtf).open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            x = line.rstrip().split("\t")
            if len(x) != 9 or x[2] != "exon" or x[6] not in ("+", "-"):
                continue
            a, b = int(x[3]) - 1, int(x[4])
            for k in range(a // 10000, (b - 1) // 10000 + 1):
                bins[(x[0], k)].append((a, b, x[6]))
    if not bins:
        raise ValueError("No annotated exons for strand inference")
    return bins


def alignment_blocks(start, cigar):
    pos = start
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        length = int(length)
        if op in "M=X":
            yield pos, pos + length
        if op in "MDN=X":
            pos += length


def strand_call(same, opposite, minimum=10000):
    n = same + opposite
    if n < minimum:
        return "ambiguous"
    fraction = same / n
    if fraction >= .90:
        return "forward"
    if fraction <= .10:
        return "reverse"
    if .40 <= fraction <= .60:
        return "none"
    return "ambiguous"


def infer_strand(bam, gtf, layout, work):
    """Compare unique primary read-1 exonic alignments with annotation strands.

    These are operational classification bands, not a biological QC acceptance.
    NH=1, MAPQ>=30; opposite-strand overlapping exons are excluded.
    At most 200,000 eligible alignments in coordinate order are examined.
    """
    bins = exon_index(gtf)
    cmd = ["samtools", "view", "-F", "3844", "-q", "30"]
    if layout == "PE":
        cmd += ["-f", "66"]  # proper pair, read 1 only
    cmd += [str(bam)]
    with (work / "strand_scan.log").open("wb") as err:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err, text=True)
        same = opposite = uninformative = eligible = 0
        stopped_early = False
        try:
            for line in proc.stdout:
                x = line.rstrip().split("\t")
                if "NH:i:1" not in x[11:]:
                    continue
                eligible += 1
                strands = set()
                for a, b in alignment_blocks(int(x[3]) - 1, x[5]):
                    for k in range(a // 10000, (b - 1) // 10000 + 1):
                        for s, e, strand in bins.get((x[2], k), ()):
                            if s < b and e > a:
                                strands.add(strand)
                if len(strands) == 1:
                    read_strand = "-" if int(x[1]) & 16 else "+"
                    if read_strand in strands:
                        same += 1
                    else:
                        opposite += 1
                else:
                    uninformative += 1
                if eligible >= 200000:
                    stopped_early = True
                    break
        finally:
            proc.stdout.close()
            if stopped_early and proc.poll() is None:
                proc.terminate()
            code = proc.wait()
        if code != 0 and not stopped_early:
            raise RuntimeError("samtools strand scan failed")
    out = {"method": "unique primary exonic read-1 annotation concordance",
           "same_strand": same, "opposite_strand": opposite,
           "uninformative": uninformative, "eligible_examined": eligible,
           "inferred_strandedness": strand_call(same, opposite),
           "minimum_informative": 10000, "forward_fraction_min": .90,
           "reverse_fraction_max": .10, "unstranded_fraction_range": [.40, .60],
           "sampling": "first up to 200000 eligible alignments in coordinate order",
           "limitation": "An empirical diagnostic; library protocol metadata takes precedence."}
    save(work / "strandedness.json", out)
    return out


def reference_check(c, with_hashes=False):
    root = Path(c["reference_dir"])
    receipt = json.loads((root / "reference_receipt.json").read_text())
    if receipt["state"] != "REFERENCE_READY":
        raise ValueError("Reference is not ready")
    if receipt["assembly"] != c["assembly"]:
        raise ValueError("Reference assembly mismatch")
    for f in receipt["files"]:
        p = root / f["relative_path"]
        st = p.stat()
        if st.st_size != f["bytes"] or st.st_mtime_ns != f["mtime_ns"]:
            if not with_hashes or digest(p) != f["sha256"]:
                raise ValueError("Reference changed: " + str(p))
    for key in ("genome",):
        if receipt["source_sha256"][key] != c[key + "_sha256"]:
            raise ValueError("Reference source identity mismatch: " + key)
        metadata_now = metadata(c[key])
        old = receipt["source_metadata"][key]
        if metadata_now != old:
            if not with_hashes or digest(c[key]) != receipt["source_sha256"][key]:
                raise ValueError("Source reference changed: " + key)
    return receipt


def publish(source, dest, files, receipt, receipt_name):
    """Stage, verify every copied final file, and rename. Never replace a result."""
    if dest.exists():
        raise FileExistsError(str(dest))
    dest.parent.mkdir(parents=True, exist_ok=True)
    staging = dest.parent / ("." + dest.name + ".publishing." + os.environ["SLURM_JOB_ID"])
    staging.mkdir()
    records = []
    for p in sorted(files):
        rel = p.relative_to(source)
        target = staging / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        checksum = digest(p)
        shutil.copy2(p, target)
        with target.open("rb") as fh:
            os.fsync(fh.fileno())
        if digest(target) != checksum:
            raise ValueError("Final copy checksum mismatch: " + str(target))
        st = target.stat()
        records.append({"relative_path": str(rel), "bytes": st.st_size,
                        "mtime_ns": st.st_mtime_ns, "sha256": checksum})
    receipt["files"] = records
    receipt["published_at"] = now()
    save(staging / receipt_name, receipt)
    with (staging / receipt_name).open("rb") as fh:
        os.fsync(fh.fileno())
    if dest.exists():
        raise FileExistsError(str(dest))
    staging.rename(dest)


def prepare_rsem_gtf(original, target, expected_excluded):
    exons = []
    loci = collections.defaultdict(set)
    genes = collections.defaultdict(set)
    features = collections.Counter()
    with Path(original).open() as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            x = line.rstrip("\n").split("\t")
            if len(x) != 9:
                raise ValueError("Malformed GTF row")
            features[x[2]] += 1
            if x[2] != "exon":
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)"', x[8]))
            tid, gid = attrs["transcript_id"], attrs["gene_id"]
            if x[6] not in ("+", "-") or int(x[3]) < 1 or int(x[4]) < int(x[3]):
                raise ValueError("Unsupported exon record: " + tid)
            loci[tid].add((x[0], x[6]))
            genes[tid].add(gid)
            exons.append((tid, line))
    excluded = {tid for tid, loc in loci.items() if len(loc) != 1 or len(genes[tid]) != 1}
    if excluded != set(expected_excluded):
        raise ValueError("Annotation exceptions differ from the reviewed configuration: " + str(sorted(excluded)))
    with Path(target).open("x") as keep, (Path(target).parent / "excluded_transcript_exons.gtf").open("x") as omit:
        for tid, line in exons:
            (omit if tid in excluded else keep).write(line)
    audit = {"input_features": dict(features), "input_exons": len(exons),
             "retained_exons": sum(tid not in excluded for tid, _ in exons),
             "retained_transcripts": len(loci) - len(excluded),
             "excluded_transcripts": [{"transcript_id": tid, "gene_ids": sorted(genes[tid]),
                "chromosome_strands": sorted(loci[tid]), "reason": "RSEM requires one chromosome and strand per transcript"}
                for tid in sorted(excluded)],
             "changes": "Exon-only GTF; retained rows are unchanged; original GTF and excluded exons are preserved",
             "limitation": ("RSEM excludes the explicitly reviewed transcript list" if excluded else "No mixed-locus or mixed-strand transcript exclusions")}
    save(Path(target).parent / "reference_annotation_audit.json", audit)
    return audit


def build_reference(c, work, cpus, tool_versions):
    if c["analysis_mode"] != "genome_mapping_only" or c.get("gtf"):
        raise ValueError("This Dpse worker is explicitly annotation-free")
    source_metadata = {"genome": metadata(c["genome"])}
    source_hashes = {"genome": digest(c["genome"])}
    if source_hashes["genome"] != c["genome_sha256"]:
        raise ValueError("Reference checksum differs from selected UCBerk genome")
    unpack(c["genome"], work / "genome.fa")
    run(["STAR", "--runThreadN", cpus, "--runMode", "genomeGenerate",
         "--genomeDir", work, "--genomeFastaFiles", work / "genome.fa",
         "--genomeSAindexNbases", "12", "--outFileNamePrefix", work / c["reference_prefix"]],
        work, "star_reference_build.log")
    for name in ("Genome", "SA", "SAindex", "genomeParameters.txt"):
        metadata(work / name)
    unchanged(list(source_metadata.values()))
    publish(work, Path(c["reference_dir"]), [p for p in work.iterdir() if p.is_file()],
            {"state": "REFERENCE_READY", "assembly": c["assembly"],
             "analysis_mode": c["analysis_mode"], "annotation_status": c["annotation_status"],
             "source_sha256": source_hashes, "tools": tool_versions,
             "source_metadata": source_metadata,
             "slurm_job_id": os.environ["SLURM_JOB_ID"]}, "reference_receipt.json")


def star_command(c, sample, work, cpus):
    paths = [sample["read1"]] + ([sample["read2"]] if sample["layout"] == "PE" else [])
    return ["STAR", "--runThreadN", str(cpus), "--genomeDir", c["reference_dir"],
            "--genomeLoad", "NoSharedMemory", "--readFilesIn", *paths,
            "--readFilesCommand", "zcat", "--outFileNamePrefix", str(work / "star."),
            "--outSAMtype", "BAM", "Unsorted", "--outSAMattributes", "NH", "HI", "AS", "nM", "NM", "MD",
            "--outSAMunmapped", "Within", "--outFilterMultimapNmax", "100",
            "--winAnchorMultimapNmax", "200", "--outSAMmultNmax", "100",
            "--twopassMode", "Basic", "--sjdbOverhang", "149",
            "--outSAMattrRGline", "ID:" + sample["sample_id"], "SM:" + sample["sample_id"],
            "LB:" + sample["run_accession"], "PL:ILLUMINA", "--runRNGseed", "42"]


def map_sample(c, sample, work, cpus, tool_versions):
    if sample["layout"] not in ("PE", "SE"):
        raise ValueError("Explicit PE or SE layout is required")
    if sample["strandedness"] not in ("auto", "none", "forward", "reverse"):
        raise ValueError("Invalid strandedness setting")
    ref = reference_check(c, with_hashes=True)
    if tool_versions["STAR"]["version"] != ref["tools"]["STAR"]["version"]:
        raise ValueError("STAR version differs from the reference builder")
    paths = [sample["read1"]] + ([sample["read2"]] if sample["layout"] == "PE" else [])
    frozen = [metadata(p) for p in paths]
    expected = sample["input_metadata"]
    if frozen != expected:
        raise ValueError("FASTQ files changed since inventory/submission")
    record = {"sample": sample, "started_at": now(), "tools": tool_versions,
              "assembly": c["assembly"], "slurm_job_id": os.environ["SLURM_JOB_ID"],
              "host": socket.gethostname(), "cpus": cpus,
              "code_sha256": digest(__file__), "scientific_QC": "NOT_YET_REVIEWED"}
    paths, record["input_compression"] = normalize_fastq_compression(paths, sample["gzip_depths"], work)
    record["fastq_integrity_check"] = check_fastq(paths, sample["layout"])
    execution_sample = dict(sample, read1=paths[0], read2=paths[1] if sample["layout"] == "PE" else "")
    run(star_command(c, execution_sample, work, cpus), work, "STAR.stdout.log")
    match = re.search(r"Number of input reads\s*\|\s*(\d+)", (work / "star.Log.final.out").read_text())
    if not match or int(match[1]) != record["fastq_integrity_check"]["checked_records_per_file"]:
        raise ValueError("STAR input count disagrees with validated FASTQ record/pair count")
    unsorted = work / "star.Aligned.out.bam"
    transcript = work / "star.Aligned.toTranscriptome.out.bam"
    run(["samtools", "quickcheck", "-v", unsorted], work, "quickcheck.log")
    bam = work / (sample["sample_id"] + ".genome.bam")
    run(["samtools", "sort", "-@", cpus - 1, "-m", "1G", "-T", work / "sorttmp",
         "-o", bam, unsorted], work, "sort.log")
    run(["samtools", "quickcheck", "-v", bam], work, "quickcheck.log")
    run(["samtools", "index", "-@", cpus - 1, bam], work, "index.log")
    run(["samtools", "flagstat", "-@", cpus - 1, bam], work, "flagstat.log", work / "flagstat.txt")
    run(["samtools", "idxstats", bam], work, "idxstats.log", work / "idxstats.tsv")
    # A full BAM pass with a primary mapped read count; zero alignments is a failure.
    run(["samtools", "view", "-c", "-F", "2308", bam], work, "mapped_count.log", work / "primary_mapped_reads.txt")
    record["primary_mapped_reads"] = int((work / "primary_mapped_reads.txt").read_text())
    if record["primary_mapped_reads"] == 0:
        raise ValueError("No primary mapped reads")
    strand = "not_inferred_without_annotation"
    record["state"] = "GENOME_MAPPED_NO_GENE_QUANTIFICATION"
    record["analysis_mode"] = c["analysis_mode"]
    record["annotation_status"] = c["annotation_status"]
    record["strandedness"] = {"inferred_strandedness": None, "used_for_RSEM": None,
                             "reason": "No matching gene annotation; genome mapping is not strand-restricted"}
    unchanged(frozen)
    save(work / "sample_provenance.json", record)
    keep = []
    for p in work.rglob("*"):
        if not p.is_file() or "rsem_tmp" in p.parts or "normalized_inputs" in p.parts:
            continue
        if p == unsorted or (p == transcript and strand != "ambiguous"):
            continue
        if any(x.startswith("sorttmp") or x.startswith("star._STAR") for x in p.parts):
            continue
        keep.append(p)
    publish(work, Path(c["results_dir"]) / sample["sample_id"], keep, record, "receipt.json")
    if strand == "ambiguous":
        print("Alignment retained; RSEM awaits strandedness review", flush=True)
        return 2
    return 0


def main():
    a = argparse.ArgumentParser()
    a.add_argument("mode", choices=["reference", "sample"])
    a.add_argument("--config", required=True)
    a.add_argument("--sample-json")
    args = a.parse_args()
    if not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_JOB_PARTITION") == "login":
        raise SystemExit("Run this worker through sbatch on a compute node")
    cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", "0"))
    if cpus < 2:
        raise SystemExit("At least two allocated CPUs are required")
    c = json.loads(Path(args.config).read_text())
    sample = json.loads(Path(args.sample_json).read_text()) if args.sample_json else None
    label = sample["sample_id"] if sample else "reference"
    scratch_root = Path(os.environ.get("TMPDIR") or c["scratch_root"])
    if not scratch_root.is_dir() or not os.access(scratch_root, os.W_OK):
        raise SystemExit("Configured compute-node scratch is unavailable: " + str(scratch_root))
    input_bytes = sum(x["bytes"] for x in sample["input_metadata"]) if sample else 0
    # A conservative admission estimate, not a measured peak-disk guarantee.
    required = max(32 * 1024**3, 6 * input_bytes)
    if shutil.disk_usage(scratch_root).free < required:
        raise SystemExit("Insufficient scratch headroom for this sample: need " + str(required))
    dest = Path(c["results_dir"]) / label if sample else Path(c["reference_dir"])
    if dest.exists():
        raise SystemExit("Refusing to replace existing output: " + str(dest))
    # This catches immediate write/fsync and coarse capacity problems. It is
    # not a reservation or a substitute for fileset quota/pool monitoring.
    shared = Path(c["run_dir"])
    if shutil.disk_usage(shared).free < required + 20 * 1024**3:
        raise SystemExit("Insufficient visible shared-filesystem headroom")
    probe_fd, probe_name = tempfile.mkstemp(prefix=".write_probe_", dir=shared)
    try:
        with os.fdopen(probe_fd, "wb") as f:
            f.write(b"0" * (1024 * 1024))
            f.flush()
            os.fsync(f.fileno())
        with open(probe_name, "rb") as f:
            if len(f.read()) != 1024 * 1024:
                raise ValueError("Shared-filesystem read-back failed")
    finally:
        Path(probe_name).unlink(missing_ok=True)
    lock = Path(c["run_dir"]) / "worker_locks" / (label + ".lock")
    lock.parent.mkdir(parents=True, exist_ok=True)
    lock.mkdir()  # Prevent two jobs publishing the same sample/reference.
    work = Path(tempfile.mkdtemp(prefix="rna_mapping_" + os.environ["SLURM_JOB_ID"] + "_", dir=scratch_root))
    save(lock / "owner.json", {"job_id": os.environ["SLURM_JOB_ID"], "scratch": str(work), "started_at": now()})
    success = False
    try:
        tool_versions = versions()
        if sample:
            code = map_sample(c, sample, work, cpus, tool_versions)
        else:
            build_reference(c, work, cpus, tool_versions)
            code = 0
        success = True
        return code
    except Exception as exc:
        failed = Path(c["run_dir"]) / "failures" / (label + ".job" + os.environ["SLURM_JOB_ID"])
        failed.mkdir(parents=True, exist_ok=False)
        save(failed / "failure.json", {"state": "FAILED", "error": repr(exc), "time": now(),
                                      "scratch": str(work), "sample": sample})
        for p in work.glob("*"):
            if p.is_file() and p.stat().st_size < 10 * 1024**2 and p.suffix in (".log", ".json", ".jsonl", ".out", ".txt"):
                shutil.copy2(p, failed / p.name)
        raise
    finally:
        if success:
            shutil.rmtree(work)  # Only this job's private scratch, after verified publication.
            shutil.rmtree(lock)
        # Failure keeps scratch and a lock for inspection; a killed job may leave
        # the same state. Never automatically delete/requeue an uncertain job.


if __name__ == "__main__":
    sys.exit(main())
