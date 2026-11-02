# HPC source map

Use [the ordered command guide](../docs/COMMANDS_IN_ORDER.md) as the entry point.

| Directory | Role |
| --- | --- |
| `star/Dmel` | Executed STAR plus RSEM worker, reference script, independent-library launcher |
| `star/Dmir` | Executed STAR plus RSEM worker and historical batch launcher |
| `star/Dpse` | Executed annotation-free, two-pass STAR worker and historical batch launcher |
| `telescope/shared_20261029` | **Current comparison version:** common library, RepeatMasker, Telescope, matrices |
| `telescope/initial_20261028` | Historical first counting implementation |
| `telescope/recovery_20261028_r2` | Historical numerical/runtime recovery |
| `telescope/legacy_Dpse_20261028_r3` | Historical Dpse-specific library and annotation |

The STAR `slurm_examples` are representative actual job scripts with path/account substitutions. The Dmel reference script automatically invokes its sample launcher; the Dmir/Dpse historical launchers attach `afterok` reference dependencies. The manual sequence in the guide separates reference validation from independent sample submission. Do not run historical launchers against an active/completed controller.

The Telescope launchers generate one Slurm job per sample, without arrays, scheduler dependencies, fixed partitions, or node names. The archived resource and walltime settings remain visible as execution history.
