# What the output files count

**Annotation rows, genomic loci, expressed loci, and RNA fragments are different quantities.** Counting lines with `wc -l` does not make them equivalent.

| File | A row represents | Includes zeros? |
| --- | --- | --- |
| RepeatMasker `.out` | A repeat alignment fragment, with headers and possibly lower-priority overlapping hits | Not an expression table |
| `TE.gtf` | An annotated genomic block; several blocks may share one `locus` attribute | Not an expression table |
| `locus_catalog.tsv` | One inferred genomic locus model | Every retained locus |
| Native `*-telescope_report.tsv` | A locus represented in Telescope's assignment report, plus metadata/special rows | Can contain zero estimates; does not enumerate the whole catalogue |
| `locus_counts.tsv` | One catalogue locus, with estimated RNA fragments | Yes; absent native-report loci are filled with zero |
| `model_group_counts.tsv` | Sum over loci assigned to one reference-model/class group | Yes |
| `model_group_estimated_counts.tsv` | One model/class group across all selected sample columns | Yes |

Telescope uses `--attribute locus`, so two GTF blocks with the same locus ID generate **one locus count**, not two. The complete count table should match **unique GTF locus IDs**, after allowing for its single header. Its row count should not equal the number of GTF blocks.

## Dpse example: the two annotation versions

| Quantity | Earlier obscura-only library | Current common library |
| --- | ---: | ---: |
| Reference models searched | 303 | 679 |
| GTF blocks | 99,474 | 79,828 |
| Unique locus models | 80,262 | 65,286 |
| `locus_counts.tsv` physical lines | 80,263 | 65,287 |
| Represented model/class groups | 231 | 491 |

The current Dpse `.out` has 303,470 numeric rows plus three headers. After excluding 30,454 lower-priority overlapping rows, 173,268 simple-repeat rows, and 19,382 low-complexity rows, 80,366 retained fragments become 79,828 blocks and 65,286 inferred loci. Adjacent/overlapping fragments within a locus are merged; gaps remain unannotated. The grouping key is **contig + RepeatMasker ID + model name + strand**.

For the example library `Dpse.RNA_female_AdultHead_Rep1_DRR055248`, the earlier native report has 18,889 TE rows; the current report has 17,754 TE rows. Both also contain zero estimates. Their complete locus tables have exactly the catalogue IDs, with no missing or extra locus IDs in the 2 November audit. Native report physical line counts add a comment, header, and `__no_feature` row.

A concrete current locus, `Dpse_RM_92a1d7b6c6562e0e527e`, has two blocks on `CM014738.1` (12039–12230 and 21813–22009), but one count-table row, with 5.0 estimated fragments in that example library.

Changing the library changes the search results, overlap competition, retained annotations, and RNA assignment. Therefore the old and current count values are different **analyses**, not duplicate copies of one analysis. Match a result to its receipt and annotation hash rather than inferring the annotation from its parent folder's name.

## Counts, model groups, and biological interpretation

`estimated_fragments` is Telescope's fractional RNA-fragment estimate after multimapping reassignment. It is not the number of annotated copies or necessarily an integer. In paired-end data, fragments refer to read pairs; do not label them individual reads.

`model_unique_fragments` is the model's unique-assignment statistic; it is not interchangeable with STAR's genome-wide uniquely mapped read count.

`CPM_total_input_fragments = estimated_fragments / total_input_fragments × 1,000,000`. This is not TPM and does not adjust for locus length. Model-group values are sums of the full locus estimates. Native, complete-locus, and model-group estimated-fragment sums must agree within numeric precision.

The shared library standardizes the models searched. It does not prove that 491 Dpse reference-model groups are 491 distinct biological TE families, nor that unknown candidates are all confirmed TEs. Model-family reconciliation, mappability, assembly differences, biological replication, protocol, and study effects remain necessary for comparisons. Genomic coverage uses interval unions and non-N assembly bases; inferred loci are not validated intact insertions. TE-overlapping RNA can include host-gene/read-through transcription and does not by itself prove autonomous TE activation.
