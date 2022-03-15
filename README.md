# Analysis_of_single_cell_data_for_Usp22_project

The analysis pipeline can be run to reproduce the results of the scRNA-seq analysis in our manuscript.
Data required to run the analysis will be made available as soon as the manuscript is published.

One way to run the pipeline is to create a new R-project in RStudio using the Git option. When creating the R-project, specify the URL of our git repository to create a local copy in your R-project. After cloning the repository successfully, the whole analysis can be executed by running the main.R script:
1. Run the first section of the main.R file until all dependencies are installed and R is restarted.
2. Before running the rest of the main.R file, adjust the paths to the data and - if desired - change the name of the output folder     which will be created by the script.
3. Afterwards run the rest of the main.R script and the complete pipeline should be executed automatically.
Make sure that your computer remains connected to the internet while running the pipeline as some functions connect to databases.

known issue:
Sometimes the Biomart server is not available and therefore the cell annotation function fails. This issue should resolve when running the script at a later time.
