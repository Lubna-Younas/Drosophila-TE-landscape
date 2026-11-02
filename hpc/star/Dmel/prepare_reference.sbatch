#!/bin/bash -l
#SBATCH --job-name=RNA_reference
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/PROJECT_ROOT/mapping_lisc_20261023_r3/logs/RNA_reference.%j.out
#SBATCH --error=/PROJECT_ROOT/mapping_lisc_20261023_r3/logs/RNA_reference.%j.err
set -euo pipefail
source /PROJECT_ROOT/mapping_lisc_20261023_r3/environment.sh
scratch="${TMPDIR:-/tmp}"
if [[ "$(findmnt -n -o FSTYPE -T "$scratch")" == "tmpfs" ]]; then
 echo 'Refusing RAM-backed scratch' >&2; exit 1
fi
python3 /PROJECT_ROOT/mapping_lisc_20261023_r3/mapping_worker.py reference --config /PROJECT_ROOT/mapping_lisc_20261023_r3/config.json

# Submit independent libraries only after successful reference validation.
python3 /PROJECT_ROOT/mapping_lisc_20261023_r3/submit_samples.py --submit
