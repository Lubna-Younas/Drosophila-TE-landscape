# Reproduce the two gene and TE density Circos figures

Published 3 November 2026 for **Spatiotemporal regulation of transposable elements across Drosophila species**.

[View the two-page PDF](../analysis/gene_TE_density_circos_20261103/figures/both_circos_designs.pdf) or browse the [complete figure package](../analysis/gene_TE_density_circos_20261103/README.md).

## Commands in order

From the repository root, create a local Python environment (Python 3.10 or 3.11 recommended), install the plotting dependencies, then enter the figure package:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r analysis/gene_TE_density_circos_20261103/code/requirements.txt
cd analysis/gene_TE_density_circos_20261103

# Optional: recalculate class coverage from the included frozen TE catalogues.
# This must reproduce the included class-window hashes and union coverage.
python code/prepare_class_windows.py

# Draw A (chromosome composition) and B (dominant class by position).
python code/plot_combined_circos.py

# Combine the individual figures into a two-page PDF and side-by-side preview.
python code/assemble_comparison.py

# Reconcile density values, class coverage, chromosomes, Y scaffolds and files.
python code/validate_package.py
```

Run in a copied checkout if you want to preserve the delivered artifacts byte for byte. The scripts write under this package's `figures/`, `tables/` and `provenance/` directories. They do not contact clusters, submit jobs, read FASTQ/BAM files, or rerun mapping or RepeatMasker.

## Interpretation

- Both designs retain Dmir outside, Dpse in the middle, Dmel inside, and the dated central tree.
- Genes use a shared normalization: genes per 100 kb divided by 140. This is not gene base-pair coverage.
- TE density is union-covered bases divided by assembled window span, including N bases.
- A allocates solid segments in proportion to each class's TE bases per chromosome. These segments are composition summaries, not genomic coordinates.
- B colours each window by the class with the most TE bases. It omits minority classes. White indicates no annotated TE; pale beige indicates tied leaders.
- All 62 Dmir Y scaffolds are included; their length-based display order does not reconstruct Y chromosome order or orientation.
- The new Dpse UCI_MV25_Y genome and matching annotation are identified in this package's source receipt. Do not substitute the older Dpse mapping assembly from the historical HPC scripts.

The package README gives annotation definitions, chromosome/Muller mappings, dating-model source, limitations and table descriptions. The public manifest records transformations and original/published hashes. No scientific workflow was rerun for GitHub publication.
