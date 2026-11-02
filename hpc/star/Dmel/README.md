# Dmel STAR source snapshot

Executed controller: `mapping_lisc_20261023_r3`.

`mapping_worker.py` is the scientific worker copied from the executed workflow.
Its source hash was compared with the mapping receipts for every selected TE
sample. `config.example.json`, `sample.example.json` and `slurm_examples/`
contain path/account placeholders and require local configuration. A zero
`mtime_ns` is a schema placeholder, never a validation result. Create input
metadata only after verifying the intended files independently on compute.

Read [the ordered command record](../../../docs/COMMANDS_IN_ORDER.md) before
using these files. The reference, samples, input metadata and runtime receipts
are external to GitHub. No FASTQs, BAMs, indexes or scientific result matrices
are distributed here.

The historical reference Slurm script automatically calls the sample launcher
on success. Its original inventory had 50 libraries; SRR6056856 is excluded
from the current 49-library TE cohort. Use the documented explicit sample
selection when adapting the workflow rather than automatically releasing the
entire historical inventory.
