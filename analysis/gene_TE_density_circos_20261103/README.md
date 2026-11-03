# Combined Circos, version 3: two solid TE-class ring designs

Two alternative figures use exactly the same inputs and density values. Version 2 and all earlier inputs remain unchanged.

## Choose the figure by the question

| Figure | Fourth ring | What it shows | What it does not show |
|---|---|---|---|
| **A — chromosome composition** | Broad, solid coloured segments within each chromosome sector | Each class's proportion of that chromosome's total TE-covered bases | Local positions or physical blocks of that class |
| **B — dominant class by position** | Constant-thickness coloured tiles along the chromosome | The class with the most covered bases in each 100 kb window | Minority classes within the same window |

A is simplest for comparing overall TE composition. B is appropriate for locating changes in the predominant class along a chromosome. Neither is an expression plot, and neither implies evolutionary conversion of one class into another.

The class palette is identical in both figures and all species: LTR, LINE, DNA, RC / Helitron, Unknown and Mixed. Unknown remains explicit. Mixed denotes genomic bases overlapped by annotations from more than one class, not a tied ranking. No separate SINE category exists in the retained catalogue.

In B, **white = no annotated TE** (30 displayed windows). **Pale beige = tied highest class coverage** (3 displayed windows; no arbitrary winner). Exact tied classes and the fraction of TE bases belonging to the dominant class are provided in `tables/plotted_windows_100kb.tsv`. The majority may be less than 50%; 'dominant' means largest coverage, not necessarily an absolute majority.

## Shared 0–1 density scales

All three species use identical scaling constants; no species-specific or chromosome-specific min–max normalization is used.

- **Gene density = (genes per 100 kb) / 140.** A value of 1 means 140 genes per 100 kb, 0.5 means 70. This is a normalized gene count density, not the fraction of bases occupied by genes. The common reference of 140 preserves the previous figure's scale; the largest observed density is 137, so no window is clipped.
- **TE density = union of TE-covered bases / assembled window span.** A value of 1 means 100% coverage, 0.5 means 50%. Overlapping TE annotations are counted once; the denominator includes N bases.

Raw gene counts, original densities, TE bases, coordinates and window widths are preserved alongside new `gene_density_0_1` and `TE_density_0_1` columns. Both displayed density axes run from 0 to 1, with 0, 0.5 and 1 marks beside each species' tracks.

The ring order remains **Dmir outside, Dpse middle, Dmel inside**. Each species has four tracks, outside to inside: chromosome span, gene density, total TE density, TE-class ring. Chromosome sectors remain proportional to their assembled span. In A only, the subdivisions of the class ring are compositional segments and must not be interpreted using the adjacent Mb coordinate ticks.

## Chromosome and Y-scaffold representation

All **62 Dmir Y-associated scaffolds (110,504,817 bp)** are included. Their display concatenation is ordered by descending length, then accession. Small marks identify the scaffold boundaries. The Y-sector coordinates are cumulative display span, not a reconstructed chromosome coordinate. No order or orientation is inferred, and original windows are never merged across scaffold boundaries.

Dpse X is a single assembly sequence containing XL and XR (A+D). Dmir XR is predominantly Muller D with A-derived material and is labelled D*. Dmir neo-X is chromosome 3 / Muller C. Non-Y unplaced scaffolds and mitochondria remain in the source tables but are not plotted. Each species independently spans the same angular range; angular alignment is not synteny.

## Dated central tree

The tree uses horizontal branch lengths in **millions of years**, with contemporaneous tips at zero. This version uses the chronology in [Gao et al. (2007), Figure 4 and its accompanying discussion](https://link.springer.com/article/10.1186/1471-2148-7-87).

| Node | Range in the cited dating model | Position used for drawing |
|---|---:|---:|
| Dmel versus the displayed obscura-group species | 30–35 Ma calibration | 32.5 Ma |
| Dpse versus Dmir | 3.76–4.93 Ma estimated ages reported in the text | 4.345 Ma, labelled 4.35 |

The plotted node positions are arithmetic range midpoints chosen for illustration, not new statistical age estimates. Grey horizontal marks show the reported ranges; these are **not confidence intervals**. Published ages differ across calibrations and clock models. This figure uses one explicitly identified chronology rather than claiming a consensus date. The pruned tree contains only these three taxa.

The tree is `(Dmel:32.5,(Dpse:4.345,Dmir:4.345):28.155);`, with branch lengths in Myr. All three root-to-tip path lengths equal 32.5 Myr. Full provenance is in `source/divergence_ages.json`.

## Data and validation

The existing validated gene and total-TE windows are unchanged. Dmel uses FlyBase r6.69, Dmir uses D. miranda PacBio2.1, and Dpse uses the newly selected `Dpse.UCI_MV25_Y.chr.fa` with its matching annotation. The exact genomes, annotations and SHA256 values are recorded in the frozen export receipt. All three retain the same 679-model comparative TE library.

The preserved class windows were calculated from each species' matching frozen `locus_catalog.tsv`, using its annotated blocks. Catalogue hashes must equal those in the original density-export receipt before calculation. The calculation uses only small local annotation tables; it does not remap reads, run RepeatMasker again or scan FASTQ/BAM files.

The class coverage calculation partitions overlaps into mutually exclusive genomic intervals. All overlapping annotations of one class count once in that class; overlaps involving different classes count once in Mixed. For each source window, the sum of all class bases must exactly equal its previously validated total TE-covered bases. Per-chromosome class totals must also match the original export. Both checks are required before plotting a class ring.

Gene density counts one midpoint per unique gene ID, including all annotated gene biotypes. Nonoverlapping windows are 100 kb except at sequence ends, where density uses the actual window length. Total TE density uses assembled span, including N bases, as its denominator. These are genomic annotation densities, not RNA expression measurements.

## Files and reproducibility

- `figures/A_composition_circos.pdf`, `.svg`, `.png`: full-resolution chromosome-composition design.
- `figures/B_dominant_circos.pdf`, `.svg`, `.png`: full-resolution dominant-class positional design.
- `figures/A_composition_preview.png`, `B_dominant_preview.png`: compact individual previews.
- `figures/compare_both_designs.png`: side-by-side preview.
- `figures/both_circos_designs.pdf`: two-page comparison PDF (A, then B).
- `tables/plotted_windows_100kb.tsv`: all 6,117 displayed windows, raw and normalized density, class bases, dominant class, ties and display offsets.
- `tables/chromosome_TE_class_composition.tsv`: exact class base totals, proportions and composition-only angular segments for A.
- `tables/plotted_chromosomes.tsv`: chromosome totals and angular spans.
- `tables/Dmir_Y_scaffold_display_offsets.tsv`: all 62 Y-scaffold accessions and display offsets.
- `tables/*_TE_classes_100kb.tsv`, `source/`, `provenance/`: frozen source data, inherited validation and revision audit.
- `code/`: reproducible preparation, rendering and independent validation.

Run with Python and `code/requirements.txt`:

```bash
python code/plot_combined_circos.py
python code/assemble_comparison.py
python code/validate_package.py
```

The main figures regenerate directly from frozen inputs. The bundled small density export is verified against the included input manifest; no private path or previous results directory is needed. No FASTQ/BAM scans, mapping, RepeatMasker runs, or expression recounts are performed for this figure revision. New checks independently reconcile normalized values to raw counts, class composition to the earlier superfamily coverage export, dominant classes to each window, and Y-scaffold membership to the assembly report.

## Public-package provenance

The figures, class catalogues and numeric tables are unchanged from the validated local version 3 package. Personal paths, cluster hostnames and job identifiers were removed from public provenance records; original scientific input basenames, assembly identities, library identity and checksums remain. Transfer-session details are omitted. The renderer's input-loading step was adapted to verify the bundled inputs without a private sibling directory; plotting calculations and layout are unchanged. `provenance/publication_manifest.json` records original and published SHA256 values and identifies each transformation. Original execution and public reproduction validation are recorded separately.

The frozen `input_manifest.json` describes scientific inputs. The publication manifest describes this published release, including delivered figures. Re-running the code may change timestamps or generated PDF/SVG metadata; those changed files are a new local output and should not be silently substituted for the published release.
