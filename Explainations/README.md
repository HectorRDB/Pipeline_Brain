# What script should I run and in which order if I want to reproduce that analysis

The scripts to generate the datasets in 
1. All the scripts in [Scripts/Smart-Seq/](/Scripts/Smart-Seq/) should be run in order (1 through 9) for both the SMARTer_cells_MOp and SMARTer_nuclei_MOp datasets.
2. All the scripts in [Scripts/10X](/Scripts/10X/) should be run in order (1 through 9) for both the 10x_cells_MOp and 10x_nuclei_MOp datasets. Note that pipelines 1/ and 2/ are independent, and that running a dataset through the pipeline can be done in the same time as running another dataset through that pipeline.
3. Run [Script/Viz/01-zinb_tsne.R](/Script/Viz/01-zinb_tsne.R)
4. Run [Script/Dune-Comparison/Dune.R](/Script/Dune-Comparison/Dune.R) for two sets of parameters.
5. Run all scripts in [Scripts/Replicability](/Scripts/Replicability/) in order.
6. Run [Script/Viz/02-tsne_replicability.R](/Script/Viz/02-tsne_replicability.R) and [Script/Viz/03-make_gifs.R](/Script/Viz/03-make_gifs.R) 


You can now generate the reports in the [Report](/Report/) folder. The [dataset_analysis.Rmd](/Report/IndividualDatasets/dataset_analysis.Rmd) file can be run on all four datasets. The [DuneVSHierarchical](/Report/IndividualDatasets/DuneVSHierarchical.Rmd), [singleMethod](/Report/IndividualDatasets/singleMethod.Rmd) and [MiniAtlasFigure](/Report/IndividualDatasets/MiniAtlasFigure.Rmd) files generate the figures.  

More details on workflow can be found in the powerpoint presentation [Workflow](/Explainations/Workflow.pptx) and a list of all the data generated and how to generate it can be found in [generateDatasets](/Explainations/generateDatasets.csv). 
