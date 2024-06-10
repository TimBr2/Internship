# Internship
Welcome to the github repository for my internship at the Pediatric Precision Oncology Lab (PPOL) Ghent. \
The goal of the internship is to do a shiny app development fo single cell RNAseq data of IPSC models.
This repository contains the source code and documentation of the SingleCell DiffTrack ShinyApp for the Wild Type single cell RNAseq data from the publication: https://www.sciencedirect.com/science/article/pii/S2589004223021739

## Exercises
To get familiar with ShinyApp development, working with single cell RNA seq data and publishing the app, i did tutorials on https://shiny.posit.co/, https://satijalab.org/seurat/, https://www.appsilon.com/post/r-shiny-docker-getting-started and https://rpubs.com/LoanRobinson/shiny_python. All different folder for the exercises can be found under the folder Exercises.

## ShinyApp
In the ShinyApp folder, two ShinyApp scripts can be found, one with ALK data (BruggemanTim_ShinyApp.R) and one without ALK data (BruggemanTim_ShinyApp_NoALK.R). The image made from the SupportScripts/Pseudotime_Calclations.R can also be found in this folder. This image is saved as Pseudotime.png. Also the logo of the lab can be found here, the PPOL_Handjes.png.
The actual used datafiles were to big to store in this repository and can be found on the PPOL github: https://github.com/PPOLLabGhent/SingleCell_DiffTrack.

## SupportScripts
In the SupportScripts repository, the scripts that helped make the ShinyApps can be found. The DEGenes.R was used to make the DEGenes.rds for the DEGenes datatable. The Pseudotime_Calculations.R was used to make the Pseudotime.png and the sce_slingshot.rds. The Profiling.R script has the package profvis to identify and optimize the performance of the R script.

## Usage
The published ShinyApp can be found on https://ccgg.ugent.be/shiny/singlecell_difftrack/. For this app, the script BruggemanTim_ShinyApp_NoALK.R was used.
