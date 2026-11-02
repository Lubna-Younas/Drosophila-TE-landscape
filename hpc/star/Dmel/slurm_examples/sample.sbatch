#!/bin/bash -l
#SBATCH --job-name=RNA_dmel_emb.4to6h_st9_R1
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/PROJECT_ROOT/mapping_lisc_20261023_r3/logs/RNA_dmel_emb.4to6h_st9_R1.%j.out
#SBATCH --error=/PROJECT_ROOT/mapping_lisc_20261023_r3/logs/RNA_dmel_emb.4to6h_st9_R1.%j.err
set -euo pipefail
source /PROJECT_ROOT/mapping_lisc_20261023_r3/environment.sh
scratch="${TMPDIR:-/tmp}"
if [[ "$(findmnt -n -o FSTYPE -T "$scratch")" == "tmpfs" ]]; then
 echo 'Refusing RAM-backed scratch' >&2; exit 1
fi
python3 /PROJECT_ROOT/mapping_lisc_20261023_r3/mapping_worker.py sample --config /PROJECT_ROOT/mapping_lisc_20261023_r3/config.json --sample-json /PROJECT_ROOT/mapping_lisc_20261023_r3/samples/dmel_emb.4to6h_st9_R1.json
