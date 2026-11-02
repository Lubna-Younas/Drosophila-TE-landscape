# Provenance and reproducibility

## Source export

This release was assembled on **2 November 2026** from the executed LiSC controllers. [source_manifest.tsv](../metadata/source_manifest.tsv) records the original controller/file SHA256, the published SHA256, and whether path/account substitutions were required. Original and published hashes can differ because private paths were replaced. Scientific parameters were preserved.

[source_verification.json](../metadata/source_verification.json) records a live comparison of STAR worker source hashes with mapping receipts for all **49 Dmel, 22 Dmir, and 53 Dpse** selected libraries. All matched. This verifies which source produced those records; it is not a rerun or a new biological QC decision.

| Public snapshot | Original controller label | Purpose |
| --- | --- | --- |
| `hpc/star/Dmel` | `mapping_lisc_20261023_r3` | STAR + RSEM |
| `hpc/star/Dmir` | `mapping_Dmir_20261023` | STAR + RSEM |
| `hpc/star/Dpse` | `mapping_Dpse_20261023` | Genome-only STAR |
| `hpc/telescope/initial_20261028` | `te_expression_20261028` | First runtime/counting implementation |
| `hpc/telescope/recovery_20261028_r2` | `te_expression_recovery_20261028_r2` | Numerical/runtime recovery |
| `hpc/telescope/legacy_Dpse_20261028_r3` | `te_expression_recovery_20261028_r3_Dpse` | Earlier Dpse assembly-specific annotation |
| `hpc/telescope/shared_20261029` | `comparative_TE_annotation_20261029` | **Current comparison annotation/counting** |

These snapshots are not all meant to execute in succession. The historical folders explain the origin of earlier output versions. The current branch builds a shared catalogue and recounts against unchanged genome BAMs.

Example plans/configurations preserve original scientific source hashes and illustrate private controller relationships, but public path substitutions mean their receipt identities do not validate a new installation. Original runtime manifests/receipts, full FASTQ validation records, and generated per-sample scripts remain external. No complete terminal-history transcript is claimed. Reproducing elsewhere requires configuring and independently validating these inputs, then creating new immutable identities.

## References and cohort

Exact compressed reference basenames and SHA256 values are in [references.tsv](../metadata/references.tsv). File compression changes can change those source hashes; the worker also records the staged reference it actually uses.

The original inventory has 125 libraries (50 Dmel, 22 Dmir, 53 Dpse); the current selected TE cohort has 124. `SRR6056856` is excluded. [samples.tsv](../metadata/samples.tsv) records inclusion separately from the original inventory and retains accessions and source sample names. Blank biological metadata remains unknown rather than inferred as fact. Sample labels and accessions are public; private storage, login, account, and compute-host details are omitted.

## Software identity

| Component | Recorded identity |
| --- | --- |
| STAR | 2.7.11b, module `STAR/2.7.11b-GCC-13.3.0` |
| SAMtools | 1.21, module `SAMtools/1.21-GCC-13.3.0` |
| RSEM | Executable 1.3.1; module labelled `RSEM/1.3.3-foss-2024a` |
| RepeatMasker | 4.2.6 |
| RMBlast | 2.14.1 |
| FamDB / Dfam | 3.0.0 / 4.0 curated installed mirror |
| Telescope | 1.0.4.1, build `py36hab16a5f_0` |
| Telescope environment | Python 3.6, pandas 1.1.3, isolated version-specific strand repair |

These are recorded execution versions, not a claim that the latest releases behave identically. Original package checksum files preserve runtime code/configuration identities. Installed software, database binaries, and third-party source code are not vendored in this repository.

## TE annotation sources

The earlier workflow combined unlike annotation sources: curated FlyBase TE copies for Dmel, an assembly-associated precomputed RepeatMasker annotation for Dmir, and a new Dpse RepeatMasker run with the 303-sequence obscura-group library. Their catalogue sizes were not directly comparable. The exact library release behind the earlier Dmir precomputed annotation was not established and must not be invented.

The current branch searches **the same frozen library** against all three mapped assemblies:

- Curated Dfam 4.0 export for Drosophila taxon 7215 with ancestors/descendants: 403 exported records, including 394 curated DF accessions and nine auxiliary bacterial entries. Excluding the auxiliary entries and 18 non-TE models leaves 376 TE/unknown models.
- 303 obscura-group consensus records from Hill & Betancourt (2018), Additional file 4. Source SHA256: `9d73b8ad683bcf58e52ce362ee48afc8c859479200089591719a7879f054101c`.
- Combined 679 reference models, with aliases/classification provenance retained; 81 Unknown/novel candidate groups remain explicitly labelled. Exact sequence or reverse-complement deduplication within class found no duplicates in this build. Related models are not automatically merged into biological families.

`landscape_annotation.py` records library export, filtering, sequence/model identities, RepeatMasker commands, locus construction, genomic interval-union coverage, and validation receipts. The native `.out`, `.gff`, `.tbl`, and compressed `.align` outputs are preserved with the annotation externally; they are not uploaded here. There is no de novo RepeatModeler discovery step in this release.

## Primary project/tool references

- [STAR source and manual](https://github.com/alexdobin/STAR)
- [RSEM source and documentation](https://github.com/deweylab/RSEM)
- [SAMtools source](https://github.com/samtools/samtools)
- [RepeatMasker source](https://github.com/Dfam-consortium/RepeatMasker)
- [FamDB source](https://github.com/Dfam-consortium/FamDB)
- [Telescope source](https://github.com/mlbendall/telescope)
- [Hill & Betancourt (2018), obscura-group TE library](https://link.springer.com/article/10.1186/s13100-018-0123-6)
- [FlyBase downloads](https://ftp.flybase.net/releases/)
- [NCBI Dmir assembly GCF_003369915.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_003369915.1/)
- [NCBI Dpse assembly GCA_004329205.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_004329205.1/)

Third-party tools, annotations, and libraries retain their own terms and citations. This initial source release does not choose a new software license on behalf of the project owner.
