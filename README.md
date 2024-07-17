Download the repo and unpack it. Install RStudio with R version 4.4.1 (2024-06-14 ucrt)

Double click the file interview-code.Rproj to open the project. In the RStudio Console, run the command 

(1) install.packages("renv")

(2) renv::init()

And select that you want to restore from the lockfile (Option 1). A number of packages and their dependencies will be installed. 

You should now be able to carry out the analysis:

(1) The R script entanglement-assoc-pathways.R can be run to generate the file processed-data/analysis-results.csv

(2) The R script create-figures reads in the data in analysis-results.csv to generate the figures in figures/*
