# Publish ShinyApp
This file explains step by step how to publish a ShinyApp on the vesalius.ugent.be server. All code in this file was made to make the ShinyApp for SingleCell_DiffTrack public. For different projects, different names can be used.

## Step 1: Create a directory for the ShinyApp
Make a new directory for the ShinyApp, put all the files in this directory that the ShinyApp uses. This is the actual app (named specifically app.R, nothing else), the data, any other scripts the ShinyApp sources, pictures the ShinyApp uses, ... Be aware, the name of this directory will also be the name of the public ShinyApp.
```
mkdir SingleCell_DiffTrack
```

## Step 2: Create a Shiny environment YAML file
Make a .yml file using your favorite text editor. This file defines the environment specifications for the Shiny App. It includes all the packages used. Here the file is called shiny.yml.
The first line is "channels:". Beneath this line, the channels are defined to install the packages. Beneath these channels there is a line "dependencies:". Under this dependencies line, the packages that the ShinyApp uses can be defined. The versions of the packages can also be defined but I recommend to do this only at the end because this can make it hard for Conda to solve the environment. Only the r-base version should be mentionned at the start. This is the R version of your environment.
An example on how to make this file:
```
channels:
  - conda-forge
  - bioconda
dependencies:
  - r-base=4.2
  - r-shiny
  - r-seurat
  - r-bslib
  - r-shinydashboard
  - r-dt
  - r-dplyr
  - r-ggplot2
  - r-ggpubr
  - bioconductor-slingshot
  - r-viridis
  - bioconductor-destiny
  - bioconductor-escape
  - bioconductor-ucell
  - r-gridextra
  - r-shinyjs
  - r-openxlsx
```

## Step 3: Create a Conda environment from the YAML file
The command uses the shiny.yml file to create a new Conda environment named shinyenv. the -f specifies the YAML file to use and -p specifies the path where the environment should be created. The environment will contain all the packages listed in shiny.yml
```
conda env create -f shiny.yml -p shinyenv
```

## Step 4: Transfer files to the remote server
The folder with the ShinyApp and the other files can now be transferred to the remote server(Vesalius). To transfer files to the remote server, a username and password is needed. This can be given by the host of the server. In this case, the username is sarahlee.
To copy a secure copy is used, this is scp. The -r means "recursive", allowing to copy directories with the contents. The directory SingleCell_DiffTrack is copied to the /home/sarahlee/projects folder on the Vesalius server.
```
scp -r ./SingleCell_DiffTrack sarahlee@vesalius.ugent.be:/home/sarahlee/projects
```
After this command you will be prompted to enter your password.

Now you can go to the remote server with the ssh command.
```
ssh sarahlee@vesalius.ugent.be
```
Then go to the directory with your ShinyApp using the cd command

## Step 5: Create a Conda environment on the remote server
The command uses the shiny.yml file to create a new Conda environment named shinyenv. the -f specifies the YAML file to use and -p specifies the path where the environment should be created. The environment will contain all the packages listed in shiny.yml. In this case a Conda environment was created in the SingleCell_DiffTrack directory on the Vesalius server using the shiny.yml file.
```
conda env create -f shiny.yml -p SingleCell_DiffTrack
```

## Step 6: Activate the Conda environment in RStudio
This command will activate the Conda environment in RStudio. Once it is activated, RStudio will use the packages specified in the shiny.yml file.
```
conda_rstudio activate
```

## Step 7: Access RStudio and test your ShinyApp
Open your browser and navigate to the URL to access RStudio in your directories environment (here SingleCell_DiffTrack). Test your ShinyApp to ensure it works correctly before publicating it.
The URL is made of https://vesalius.ugent.be/rstudio/directoryname_username
In this case the URL is https://vesalius.ugent.be/rstudio/SingleCell_DiffTrack_sarahlee with SingleCell_DiffTrack the directory name and sarahlee the username.

## Step 8: Make the ShinyApp public
When the app works perfectly in RStudio, you can publicate the app with following command:
``` I Don't have this command yet ```

# Troubleshooting
The troubleshooting can maybe help you solve errors. This part consist of errors I had to solve.

