# Drosophila TE landscape

HPC workflow record for **Spatiotemporal regulation of transposable elements across Drosophila species**.

This repository contains the STAR RNA-seq mapping and Telescope TE-expression code used for **Drosophila species**. 

**Start with [Commands in execution order](docs/COMMANDS_IN_ORDER.md).** It connects each command to its inputs, outputs, and validation gate.

| Step | Code | Result |
| --- | --- | --- |
| 1. Reference preparation and STAR mapping | [Dmel](hpc/star/Dmel/), [Dmir](hpc/star/Dmir/), [Dpse](hpc/star/Dpse/) | One genome BAM and mapping receipt per library; RSEM gene/isoform estimates for Dmel and Dmir |
| 2. Common TE library and genome annotation | [Current shared workflow](hpc/telescope/shared_20261029/landscape_annotation.py) | The same 679 reference models searched against each exact mapped assembly |
| 3. Telescope quantification | [Current counting worker](hpc/telescope/shared_20261029/pipeline.py) | Estimated fragments per inferred TE locus and reference-model group |
| 4. Validation and matrix collection | [Collector](hpc/telescope/shared_20261029/collect.py) | Species-specific count and total-input CPM matrices |

The current TE cohort has **124 libraries: 49 Dmel, 22 Dmir, and 53 Dpse**. Biological replicates remain separate. [Sample metadata](metadata/samples.tsv) also records the one excluded Dmel library from the original 125-library inventory. The STAR worker source hashes match the mapping receipts for all 124 selected libraries; this checks provenance and does not replace biological QC.

The Dpse branch used genome-only STAR alignment because a matching gene annotation was not supplied. It therefore has no RSEM gene/isoform outputs. Its genome BAMs can still be counted against the assembly-matched TE annotation.

These files preserve the executed workflow and its recovery history. **They are source snapshots, not a configured one-command installation.** Before reuse, supply the exact references, validated sample manifests, local paths, runtime receipts, and newly frozen package identities described in the command guide. Example metadata and historical checksums must not be treated as validation of new inputs. Large sequencing files, indexes, installed software, and scientific result matrices are stored outside GitHub.

Useful guides:

- [Gene and TE density figures: two Circos designs (3 November 2026)](analysis/gene_TE_density_circos_20261103/README.md)

- [Output files and why their row counts differ](docs/OUTPUTS_EXPLAINED.md)
- [Assemblies, software versions, and source provenance](docs/PROVENANCE.md)
- [What changed during recovery](docs/RECOVERY_HISTORY.md)
- [HPC resources and independent jobs](docs/RESOURCES.md)
- [Repository maintenance commands](docs/GITHUB_MAINTENANCE.md)

The **shared 29 October annotation** is the version selected for current cross-species comparisons. Earlier Telescope folders are retained to explain prior results, not to mix their counts with the current catalogue. Reference-model groups still require biological family reconciliation; a shared library alone does not establish comparable biological families or remove study, tissue, and library biases.

## Gene and TE density figures — 3 November 2026

- **A — chromosome composition:** broad solid segments show each TE class's share of a chromosome's TE-covered bases; segment positions are not genomic locations.
- **B — dominant class by position:** each 100 kb window is coloured by its largest TE class; minority classes are omitted.

[Two-page PDF](analysis/gene_TE_density_circos_20261103/figures/both_circos_designs.pdf) · [Figures, editable SVGs and data](analysis/gene_TE_density_circos_20261103/) · [Reproduction commands](docs/CIRCOS_FIGURES.md)

Both figures use shared 0–1 scales: gene density is genes per 100 kb divided by 140; TE density is covered bases divided by assembled window span. They retain all 62 Dmir Y-associated scaffolds. The package includes exact plotting inputs, class catalogues, source hashes, code and validation for all 6,117 displayed windows.

This figure package uses the **new Dpse UCI_MV25_Y composite assembly and its matching annotations**, with the same 679-model TE library. It is distinct from the earlier GCA_004329205.1 Dpse STAR/Telescope execution snapshot documented above. These figures show genomic annotation densities, not RNA expression or conserved insertions. No mapping jobs are launched by figure reproduction.

Repository checks use only Python's standard library and Bash:

```bash
python3 scripts/validate_repository.py
python3 scripts/run_synthetic_checks.py
```

These check published source hashes, syntax, metadata, documentation links, publication hygiene, and 32 existing synthetic cases. They do not run mapping or establish scientific QC acceptance. The synthetic selection excludes historical tests that require private sample inventories or generated jobs.
