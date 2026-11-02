#!/bin/bash -l
#SBATCH --job-name=RNA_Dpse.RNA_SexNA_St2_Rep1_SRR1055004
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --chdir=/PROJECT_ROOT/mapping_Dpse_20261023
#SBATCH --output=/PROJECT_ROOT/mapping_Dpse_20261023/logs/Dpse.RNA_SexNA_St2_Rep1_SRR1055004.%j.out
#SBATCH --error=/PROJECT_ROOT/mapping_Dpse_20261023/logs/Dpse.RNA_SexNA_St2_Rep1_SRR1055004.%j.err
set -euo pipefail
source '/PROJECT_ROOT/mapping_Dpse_20261023/environment.sh'
scratch="${TMPDIR:-/tmp}"
if [[ "$(findmnt -n -o FSTYPE -T "$scratch")" == "tmpfs" ]]; then echo 'Refusing RAM-backed scratch' >&2; exit 1; fi
python3 -u '/PROJECT_ROOT/mapping_Dpse_20261023/mapping_worker.py' sample --config '/PROJECT_ROOT/mapping_Dpse_20261023/config.json' --sample-json '/PROJECT_ROOT/mapping_Dpse_20261023/samples/Dpse.RNA_SexNA_St2_Rep1_SRR1055004.json'
