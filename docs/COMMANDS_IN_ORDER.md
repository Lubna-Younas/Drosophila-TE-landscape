# Commands in execution order

The STAR/Telescope commands below record the historical execution snapshots. For the subsequently published 3 November genomic-density figures, follow [Circos figure reproduction](CIRCOS_FIGURES.md). That package is self-contained and identifies its newer Dpse reference separately.

This is the ordered command record for the STAR → RepeatMasker → Telescope workflow. Commands were reconstructed from the executed source and saved configuration/receipts; a complete interactive shell-history transcript is not available. Internal tool commands below explain what the workers executed. They are **not additional steps to run after the worker**, which would duplicate the work.

The source snapshots use `/PROJECT_ROOT`, `/SOFTWARE_ROOT`, `/DFAM_FAMDB`, and `YOUR_ACCOUNT` placeholders. Replace them only in a **new configured deployment**. Existing successful controllers and results remain immutable. Examples labelled **adapted recipe** explain portable submission; they are not presented as historical terminal commands.

## 1. Identify samples, read layout, and the exact assembly

Read [samples.tsv](../metadata/samples.tsv), [references.tsv](../metadata/references.tsv), and [exclusions](../metadata/excluded_samples.tsv) first.

| Species | Reference | Current TE libraries | Gene quantification |
| --- | --- | ---: | --- |
| Dmel | FlyBase r6.69 | 49 | STAR transcriptome BAM → RSEM |
| Dmir | GCF_003369915.1, D.miranda_PacBio2.1 | 22 | STAR transcriptome BAM → RSEM |
| Dpse | GCA_004329205.1, UCBerk_Dpse_1.0 | 53 | Not run; genome BAM only |

Keep replicates separate. `R1/R2/R3` in Dmel library names refer to biological replicate labels; `fastq1/fastq2` refer to paired read mates. Preserve the verified single-end Dmel libraries rather than inventing mates. Preserve Dmir canonical-file/alias handling. Dpse source labels alone do not prove independent biological replication.

Each STAR controller needs `config.json`, one `samples/<sample_id>.json` per selected library, and generated job scripts. The archived `config.example.json` and `sample.example.json` show the schema. Input bytes and nanosecond modification times must come from independently verified actual files; zero modification times in examples are placeholders. FASTQ integrity/mate validation and large checksums belong on compute nodes. Public metadata is not a substitute for these execution manifests.

Do not reinstate `dmel_female_head_R2` / `SRR6056856`: it is excluded pending resolution of paired FASTQ validation and explicit reinstatement.

## 2. Load the recorded mapping environment

Inside the Slurm job:

```bash
source "$CONTROLLER/environment.sh"
STAR --version
samtools --version
rsem-calculate-expression --version
```

Recorded versions: STAR 2.7.11b and SAMtools 1.21. The RSEM module is labelled `1.3.3-foss-2024a`, but its executable reports **1.3.1**. Record the executable version, not just the module label. Dpse's preserved environment also loads RSEM even though its genome-only worker does not run RSEM quantification.

## 3. Prepare and validate the mapping reference

The actual worker entry point is the same for all three species:

```bash
python3 "$CONTROLLER/mapping_worker.py" reference \
  --config "$CONTROLLER/config.json"
```

Run through a configured Slurm script, never as intensive work on a login node. The historical Dmel command was `sbatch prepare_reference.sbatch`; that script also auto-launches sample jobs. The Dmir/Dpse historical `submit_batch.py` submitted the reference and dependent sample jobs together. Preserve those scripts as history; use the separate submission recipe below when adapting to the current independent-job preference.

For **Dmel/Dmir**, the worker verifies and unpacks reference sources, prepares the RSEM-compatible exon GTF, and runs:

```bash
rsem-prepare-reference --gtf "$GENES_GTF" -p "$CPUS" \
  "$GENOME_FASTA" "$RSEM_REFERENCE_PREFIX"

STAR --runThreadN "$CPUS" --runMode genomeGenerate \
  --genomeDir "$REFERENCE_DIR" --genomeFastaFiles "$GENOME_FASTA" \
  --sjdbGTFfile "$GENES_GTF" --sjdbOverhang 149 \
  --genomeSAindexNbases 12 --outFileNamePrefix "$REFERENCE_PREFIX"
```

Dmel's RSEM GTF omits nine explicitly audited mixed-strand `mod(mdg4)` transcripts. The original GTF, omitted rows, and audit remain preserved. See [recovery history](RECOVERY_HISTORY.md). This is not a genomic-read exclusion. Dmir has no such transcript exclusions.

For **Dpse**, the executed reference command is genome-only:

```bash
STAR --runThreadN "$CPUS" --runMode genomeGenerate \
  --genomeDir "$REFERENCE_DIR" --genomeFastaFiles "$GENOME_FASTA" \
  --genomeSAindexNbases 12 --outFileNamePrefix "$REFERENCE_PREFIX"
```

Do not substitute the GCF_009870125.1 gene/repeat coordinates for this GCA_004329205.1 mapping assembly. Reference completion requires its receipt, source identity, and index validation, not just a Slurm `COMPLETED` label.

### Adapted recipe: independent reference and sample jobs

Create a new `reference.manual.sbatch` with the worker entry point above and the site's required account/resources/log paths. It must contain **only reference preparation**, without the Dmel automatic sample-launch line. Then:

```bash
sbatch --parsable "$CONTROLLER/reference.manual.sbatch"
# Wait for the reference receipt and index checks to pass before the next step.
```

After the reference passes, submit one configured job per selected sample:

```bash
sbatch --parsable "$CONTROLLER/jobs/$SAMPLE.sbatch"
```

These sample jobs have no array or scheduler dependency. Their command is:

```bash
python3 -u "$CONTROLLER/mapping_worker.py" sample \
  --config "$CONTROLLER/config.json" \
  --sample-json "$CONTROLLER/samples/$SAMPLE.json"
```

For Dmel, the archived launcher also supports an explicit selected library:

```bash
python3 "$CONTROLLER/submit_samples.py" --sample "$SAMPLE"
python3 "$CONTROLLER/submit_samples.py" --sample "$SAMPLE" --submit
```

The first invocation is a metadata-only dry run. Configure `samples.tsv`, sample JSONs, and `jobs/` before using it. Submit only the approved selection, not the entire original inventory.

## 4. STAR genome alignment, retaining multimappers

The worker validates the actual FASTQs, uses local disk scratch, and runs the following common options. For a single-end library, supply only `READ1` to `--readFilesIn`.

```bash
STAR --runThreadN "$CPUS" --genomeDir "$REFERENCE_DIR" \
  --genomeLoad NoSharedMemory \
  --readFilesIn "$READ1" "$READ2" --readFilesCommand zcat \
  --outFileNamePrefix "$WORK/star." --outSAMtype BAM Unsorted \
  --outSAMattributes NH HI AS nM NM MD --outSAMunmapped Within \
  --outFilterMultimapNmax 100 --winAnchorMultimapNmax 200 \
  --outSAMmultNmax 100 \
  --outSAMattrRGline "ID:$SAMPLE" "SM:$SAMPLE" "LB:$ACCESSION" PL:ILLUMINA \
  --runRNGseed 42
```

The **Dmel/Dmir** worker adds:

```text
--quantMode TranscriptomeSAM GeneCounts
--quantTranscriptomeSAMoutput BanSingleEnd_BanIndels_ExtendSoftclip
```

The **Dpse** worker instead adds:

```text
--twopassMode Basic --sjdbOverhang 149
```

Only Dpse is explicitly two-pass in these executed workers. Keep secondary/multiple genomic alignments for Telescope; STAR still caps reported mappings at 100, so these are not unlimited multimapping BAMs.

## 5. Sort, index, and check genome BAMs

Internal worker commands include:

```bash
samtools quickcheck -v "$WORK/star.Aligned.out.bam"
samtools sort -@ "$((CPUS-1))" -m 1G -T "$WORK/sorttmp" \
  -o "$GENOME_BAM" "$WORK/star.Aligned.out.bam"
samtools quickcheck -v "$GENOME_BAM"
samtools index -@ "$((CPUS-1))" "$GENOME_BAM"
samtools flagstat -@ "$((CPUS-1))" "$GENOME_BAM"
samtools idxstats "$GENOME_BAM"
samtools view -c -F 2308 "$GENOME_BAM"
```

The workers also reconcile STAR input counts with fully validated FASTQ records/pairs and record files, hashes, command logs, tools, and configuration. `quickcheck` is one structural check, not proof of biological quality or all-record correctness.

For Dmel/Dmir, the separate **transcriptome BAM** is passed to RSEM:

```bash
rsem-calculate-expression --alignments --no-bam-output \
  --strandedness "$RSEM_STRAND" --seed 42 -p "$CPUS" \
  --temporary-folder "$WORK/rsem_tmp" --paired-end \
  "$TRANSCRIPTOME_BAM" "$RSEM_REFERENCE_PREFIX" "$SAMPLE_OUTPUT_PREFIX"
```

Omit `--paired-end` for SE data. Strand inference is recorded per library; ambiguous strand evidence holds RSEM quantification while preserving genomic mapping. RSEM transcriptome BAMs are not the Telescope genome inputs.

## 6. Configure the TE runtime and freeze its inputs

The current counting source is [shared_20261029](../hpc/telescope/shared_20261029/). It reuses the isolated runtime installed by [initial_20261028](../hpc/telescope/initial_20261028/). A fresh deployment must build and validate that base runtime, rather than inventing its `READY.json` receipt or copying historical hashes to new files.

The base worker installs this exact environment on compute:

```bash
"$CONDA" create --yes --solver=libmamba --override-channels \
  --strict-channel-priority -c conda-forge -c bioconda \
  --prefix "$TE_ENV" python=3.6 \
  telescope=1.0.4.1=py36hab16a5f_0 pandas=1.1.3
```

Project-local Conda/cache paths are set by `isolated_env()`; the user's home and shared Conda configuration are unchanged. `runtime_patch.py` applies the recorded single-end F/R orientation repair to this isolated Telescope installation, and `runtime_smoke.py` plus `numeric_smoke.py` validate strand handling and numerical convergence. Do not apply this version-specific patch to an unrelated newer runtime.

Within a **configured base controller**, the compute-worker action was:

```bash
python3 pipeline.py runtime Dmel
```

The current shared controller then verifies the base runtime receipt and repeats its synthetic checks. Supply `plan.json`, validated `manifests/<species>.json` and `.tsv`, and a frozen `package_checksums.json`. Required sample fields are used explicitly in `pipeline.py`/`submit.py`: sample/accession/layout, total input fragments, strand mode, BAM/BAI size/mtime/SHA256, and mapping receipt identity. Keep public inventory separate from these execution records.

`package_checksums.original.json` is a **historical identity record**. Editing paths/configuration changes identities. Rebuild them in a fresh deployment after reviewing inputs; never rename that historical file and assume it validates edited files. The source depends on real base-runtime receipts, not merely a list of software versions.

## 7. Build one shared library; annotate each mapped genome

The shared-library producer exports curated Dfam 4.0 Drosophila genus/ancestral/descendant models through FamDB 3.0.0:

```bash
"$FAMDB" -i "$DFAM_FAMDB" -e ch,uc,uh families -ad -c \
  -f fasta_acc --include-class-in-name 7215
```

`landscape_annotation.py` filters the export and combines **376 retained Dfam TE/unknown models** with **303 Hill & Betancourt obscura-group consensus records**. The frozen result has **679 reference models**, with source aliases/classifications preserved. Exact sequence/reverse-complement deduplication is within class; it found no duplicate models in this build. Unknown candidates remain explicitly labelled. See [provenance](PROVENANCE.md) for source links and interpretation limits.

In the configured shared controller, the metadata launcher supports dry-run then submit:

```bash
python3 submit.py prepare Dmel
python3 submit.py prepare Dmel --submit
# Verify shared_library/receipt.json is SHARED_LIBRARY_FROZEN before other producers.
python3 submit.py prepare Dmir
python3 submit.py prepare Dmir --submit
python3 submit.py prepare Dpse
python3 submit.py prepare Dpse --submit
```

The first producer must finish publishing the valid shared-library receipt before other preparation jobs start. This avoids the observed initial publication race; the archived producer is not a general concurrent-publication fix. Each preparation job calls `python3 pipeline.py prepare <species>` on compute and internally runs:

```bash
RepeatMasker -engine rmblast -pa 1 -s -a -gff -no_is \
  -lib "$SHARED_LIBRARY/library.fa" -dir "$WORK/annotation" \
  "$WORK/annotation/genome.fa"
```

`-pa 1` uses four RMBlast CPUs in this setup. Preserve the RepeatMasker module's Python environment for its h5py/FamDB dependencies; isolate Telescope separately.

Before counting, require `annotations/<species>/receipt.json` with `ANNOTATION_VALIDATED`, the correct assembly/package identity, matching BAM reference dictionary, and validated annotation files. Locus construction and filtering are explained in [Outputs](OUTPUTS_EXPLAINED.md).

## 8. Submit independent Telescope jobs

After each species annotation and runtime pass:

```bash
python3 submit.py samples Dmel
python3 submit.py samples Dmel --submit
python3 submit.py samples Dmir
python3 submit.py samples Dmir --submit
python3 submit.py samples Dpse
python3 submit.py samples Dpse --submit
```

This generates separate `jobs/<species>/<sample>.sbatch` scripts and durable submission records. No arrays, fixed partitions/node names, scheduler dependencies, or controller admission throttle are added. Existing/uncertain submission records must be reconciled before any retry; never clear reservations just to force resubmission.

Each compute worker calls:

```bash
python3 pipeline.py sample "$SPECIES" "$SAMPLE"
```

It validates the genome BAM/BAI, mapping receipt, annotation identity, and matching reference dictionary. It collates **all** retained alignments:

```bash
samtools collate -@ 3 -o "$WORK/all_alignments.collated.bam" \
  "$GENOME_BAM" "$WORK/collate_tmp"
samtools quickcheck -v "$WORK/all_alignments.collated.bam"
```

Do not add `collate -f`, which would discard secondary alignments. The executed counting command is:

```bash
"$TE_ENV/bin/telescope" assign "$WORK/all_alignments.collated.bam" "$TE_GTF" \
  --attribute locus --ncpu 1 --stranded_mode "$SAMPLE_MODE" \
  --reassign_mode average --overlap_mode threshold --overlap_threshold 0.2 \
  --max_iter 100000 --em_epsilon 1e-7 --pi_prior 0 --theta_prior 200000 \
  --outdir "$OUTPUT" --exp_tag "$SAMPLE" --tempdir "$WORK"
```

Use the per-library `telescope_stranded_mode` from its validated manifest; do not assign one mode to all species by assumption. The recorded cohort has Dmel RF/FR/None modes, and Dmir/Dpse None. An uncertain protocol counted with None remains an explicit uncertainty.

Temporary collated BAMs and working files stay in compute scratch. Validated final reports, count tables, logs, and receipts are published; failure evidence/checkpoints are preserved separately for diagnosis. No FASTQ deletion is part of these commands.

## 9. Validate convergence and collect matrices

For non-zero-overlap samples, require both `EM converged after ...` and `telescope assign complete` in `quantification.log`, plus valid finite counts and a matching `TE_COUNTS_COMPLETE` receipt. A report left at an iteration cap is not accepted merely because Slurm returned zero. Zero-overlap libraries have an explicit alternative validation path in the worker.

After **all expected libraries of one species** validate, the collector is run on compute:

```bash
python3 collect.py Dmel
python3 collect.py Dmir
python3 collect.py Dpse
```

These are separate invocations; submit each in its own configured Slurm script. The collector refuses a partial matrix and rejects mixed annotation identities. It writes `matrices/<species>/model_group_estimated_counts.tsv`, `model_group_CPM_total_input.tsv`, sample metadata, and its receipt. The receipt explicitly records scientific QC as not yet reviewed.

Before declaring current completion, reconcile expected samples, Slurm accounting, final logs, configuration/annotation identities, and current file metadata. Changed large files need checksum revalidation on compute. Matrix presence alone does not establish acceptance.

## 10. Check job state and diagnose only actual failures

Metadata-only commands can run on the login node:

```bash
squeue -u "$USER" -o '%.18i %.30j %.10T %.10M %.8C %R'
sacct -j "$JOB_ID" --format=JobID,JobName,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS
tail -n 80 "$STDOUT_LOG"
tail -n 80 "$STDERR_LOG"
```

Inspect current attempts and their timestamps; old failure logs do not prove a current worker is stopped. Read [recovery history](RECOVERY_HISTORY.md) before retrying. Preserve originals and successful outputs, and retry only reconciled failures with the cause fixed. A changed TE catalogue requires recounting its samples, but does **not** require remapping unchanged genome BAMs.
