#!/bin/bash -l
#SBATCH --job-name=RNA_Dpse_reference
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --chdir=/PROJECT_ROOT/mapping_Dpse_20261023
#SBATCH --output=/PROJECT_ROOT/mapping_Dpse_20261023/logs/reference.%j.out
#SBATCH --error=/PROJECT_ROOT/mapping_Dpse_20261023/logs/reference.%j.err
set -euo pipefail
source '/PROJECT_ROOT/mapping_Dpse_20261023/environment.sh'
python3 -u '/PROJECT_ROOT/mapping_Dpse_20261023/prepare_reference_sources.py'
python3 -u '/PROJECT_ROOT/mapping_Dpse_20261023/mapping_worker.py' reference --config '/PROJECT_ROOT/mapping_Dpse_20261023/config.json'
