#!/usr/bin/env bash
#SBATCH --job-name=gbm-setup
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=__BASE__/logs/Slurm-%j.out
#SBATCH --error=__BASE__/logs/Slurm-%j.err

# =============================================================================
# Run setup.sh as a SLURM job for clean monitoring + offload from login node.
# Sources config.sh to know BASE; cd's to BASE; runs setup.sh.
# =============================================================================

cd __BASE__

echo "[$(date)] Starting setup.sh on $(hostname)"
bash setup.sh
echo "[$(date)] setup.sh exit code: $?"
