# Recovery history

This record explains why jobs or output versions changed. It is not a command to rerun completed work. Every retry must first reconcile scheduler accounting, current attempt logs, receipts, and files. Preserve failed evidence and successful outputs in separate versioned locations.

| Issue | Applied change | Scientific consequence |
| --- | --- | --- |
| Drosophila STAR reference indexing | Explicit `--genomeSAindexNbases 12` in the validated reference worker | Correct small-genome index setup; earlier failed reference retained |
| RSEM rejected mixed-strand Dmel transcript definitions | Prepared exon-only RSEM GTF and audited nine reviewed mixed-strand `mod(mdg4)` transcripts separately | Those isoforms are omitted from RSEM estimates; original annotation and genomic reads preserved |
| Some transferred reads had an extra gzip wrapper | Validated and removed only the extra compression layer in scratch; rejected unexpected depth, truncation, and mate mismatches | Original files and read content preserved; no pooling of replicates |
| Home quota affected environment setup | Put Conda/cache bookkeeping under project runtime | No scientific threshold change |
| RepeatMasker/FamDB could not import h5py under over-broad environment isolation | Kept module-provided Python paths for annotation, isolated Telescope only | No library/coordinate substitution |
| Telescope single-end F/R orientation handling | Version-specific `SE_STRAND_FIX_V1` in isolated Telescope 1.0.4.1; synthetic strand tests | Corrected SE orientation; paired/unstranded modes unchanged |
| EM reached the original 1,000-iteration cap | Raised `max_iter` to 100,000 with `em_epsilon=1e-7` unchanged; checked actual convergence | More iterations, not relaxed convergence acceptance |
| Eight late sample failures were outside the initial recovery selection while originals were still running | Reconciled those late failures and recovered only the failed attempts | Successful/active work was not blindly restarted |
| Dpse custom library contained `?/?` classification labels | Normalized 74 such labels to Unknown, preserving all 303 source IDs/sequences and classification provenance | Explicit uncertainty; models not silently discarded |
| Annotation sources differed across species | Built a 679-model shared library, reannotated the exact assemblies, and recounted unchanged genome BAMs | New feature catalogue and new RNA assignments; old/current counts cannot be mixed |
| Simultaneous shared-library producers raced during first publication | Let one producer publish and validate the library before releasing other preparation jobs | Operational readiness gate; archived code is not claimed to solve arbitrary concurrent publication |

The nine Dmel transcript IDs are:

```text
FBtr0084079 FBtr0084080 FBtr0084081 FBtr0084082 FBtr0084083
FBtr0084084 FBtr0084085 FBtr0307759 FBtr0307760
```

The runtime patch verifies the original Telescope `model.py` SHA256 `5bf65c71c4e84d8e95371819aa869f0bfa1b9e78ff82e08125a35bf0d594847e` and the patched hash `efb4eb672cce427ce9699f005f2b28478eb22dacc09b396fbea45ce369741b8f`. It changes two SE F/R orientation expressions. It must not be applied indiscriminately to other package versions.

SRR6056856 remains excluded for paired FASTQ validation and is not an automatic retry target. Ambiguous strandedness or missing scientific inputs must remain explicit holds until their evidence is resolved. A technical completion receipt is separate from scientific acceptance.

Historical launch details are preserved: Dmel reference completion triggered its library launcher; Dmir/Dpse originally used independent jobs with `afterok` reference dependencies. The current manual submission sequence waits for reference validation, then submits independent jobs without dependencies. Archived Telescope launchers already create separate jobs without arrays/dependencies.

No new mapping/counting jobs, QC-threshold changes, FASTQ deletions, or automated retries were performed while publishing this repository.
