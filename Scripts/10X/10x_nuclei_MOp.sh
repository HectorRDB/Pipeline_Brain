#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --nodes=1

module load hdf5

# Loading and filtering
loc="/pylon5/ib5phhp/hectorrb/10x_nuclei_MOp/"
out1="/pylon5/ib5phhp/hectorrb/out/10x_nuclei_MOp_filt.rds"
Rscript --no-save --verbose  1-load_filter.R -l $loc -o $out1 > 1b.out 2>&1
echo "Step 1 done"

# zinbWave
out_norm="/pylon5/ib5phhp/hectorrb/out/10x_nuclei_MOp_norm.rds"
out_rd="/pylon5/ib5phhp/hectorrb/out/10x_nuclei_MOp_zinbWs"
MEMORYFILE="2b_memoryLogger.txt"

while true; do free -h >> $MEMORYFILE; sleep 15; done & \
Rscript --no-save --verbose  2-zinb.R -n 20 -l $out1 -o $out_norm -r $out_rd > 2b.out 2>&1
echo "Step 2 done"

# sc3
out3="/pylon5/ib5phhp/hectorrb/out/10x_nuclei_MOp_sc3.rds"
Rscript --no-save --verbose  3-sc3.R -n 20 -l $out1 -o $out3 > 3b.out 2>&1
echo "Step 3 done"

# Seurat
out4="/pylon5/ib5phhp/hectorrb/out/10x_nuclei_MOp_seurat.rds"
Rscript --no-save --verbose  4-seurat.R -l $out1 -o $out4 > 4b.out 2>&1

# Rsec
MEMORYFILE="5b_memoryLogger.txt"
out5="/pylon5/ib5phhp/hectorrb/out/10x_nuclei_MOp_RSEC.rds"

while true; do free -h >> $MEMORYFILE; sleep 15; done & \
Rscript --no-save --verbose  5-RSEC.R -n 20 -l $out_rd -o $out5 > 5b.out 2>&1

# Consensus Clustering
MEMORYFILE="6b_memoryLogger.txt"
loc="/pylon5/ib5phhp/hectorrb/out/10x_nuclei_MOp"

while true; do free -h >> $MEMORYFILE; sleep 15; done & \
Rscript --no-save --verbose  6-ConsensusClustering.R -n 20 -l $loc -o $loc > 6b.out 2>&1
