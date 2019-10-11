datasets="SMARTer_cells_MOp SMARTer_nuclei_MOp 10x_cells_MOp 10x_nuclei_MOp"
# datasets="SMARTer_cells_MOp SMARTer_nuclei_MOp Regev"

for dataset in $datasets
do
  echo $dataset
  Rscript --no-save --quiet render.R -d $dataset
done
