# Dmir STAR source snapshot

Executed controller: `mapping_Dmir_20261023`.

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

The preserved `submit_batch.py` is a historical launcher: it attached an
`afterok` reference dependency to each independent sample job. The current
manual sequence submits sample scripts only after reference validation, with
no scheduler dependencies. This distinction is documented rather than silently
rewriting the command history. `config.template.json` and the generated
`sample_manifest.json`, `samples/` and `jobs/` were site-specific inputs.
