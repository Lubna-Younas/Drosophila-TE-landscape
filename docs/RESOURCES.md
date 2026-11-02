# HPC resources and storage

The archived scripts retain the resources actually requested; they have not been silently edited into a different execution history.

| Work | Recorded request | Notes |
| --- | --- | --- |
| STAR reference | 8 CPUs, 32 GB, 4 h | Verify index-build receipt before samples |
| STAR sample | 8 CPUs, 32 GB, 48 h | Separate library jobs, local disk scratch |
| RepeatMasker preparation | 4 CPUs, 32 GB, 48 h | `-pa 1` invokes four RMBlast CPUs in this setup |
| Telescope sample | 4 CPUs, per-sample memory, 48 h | Collation uses `-@ 3`; Telescope EM uses `--ncpu 1` |

The TE plan records that the site rejected a missing walltime and `time=0`, so the historical scripts use 48 h. This is a scheduler-required request, not an instruction to impose a new arbitrary analysis cutoff. Do not add array throttles, partitions, fixed nodes, or named clusters when adapting these jobs.

The resource audit suggested **future** Telescope jobs could request one CPU and use single-thread collation (`samtools collate -@ 0`), with memory based on measured peak use plus roughly 40% headroom. This proposal was not applied to the archived workers or active jobs here. Collation and Telescope phases have different CPU use; asking for four CPUs does not make a single-threaded EM use all four. Use completed-job `sacct`/live `sstat` measurements and sample size to size future jobs, rather than requesting the same large memory for every library.

Mapping, FASTQ scans, RepeatMasker, collation, checksum revalidation, Telescope, and matrix collection run on compute nodes. Login nodes are for editing, metadata inspection, scheduler queries, and submission. Temporary uncompressed reads, unsorted/collated BAMs, sorting files, and tool work directories belong in allocated local disk scratch. Successful final BAM/index/count files and their validation records are copied to persistent storage; failed attempts retain diagnostic evidence separately.

This code release does not delete raw reads or other outputs. Any later storage cleanup requires an explicit file-level inventory tying recoverable raw accessions/checksums to independently validated outputs and current downstream needs.
