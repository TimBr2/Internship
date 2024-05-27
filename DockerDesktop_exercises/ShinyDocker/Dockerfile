### This is my first ever selfmade dockerfile for a shinyapp ###
################################################################
# a couple of distinct keywords you have to learn:
##################################################

# FROM: A starting point for every Dockerfile. Used to tell Docker from which base image # you'll build your own. A rocker/shiny image already has R and Shiny installed, so      # that's the one we'll use.

# RUN: Runs a shell command. Useful for creating files and folders, installing
# dependencies, and more.

# COPY: Copies contents from your local machine to the Docker container. You first
# specify the path to the local file, followed by the file on the container.

# EXPOSE: Exposes a port of a Docker container so you can access it from your local
# machine.

# CMD: Command that will be used every time you launch the container. We'll use it to
# launch the R Shiny app.

#######################################################################################
# The actual docker script:
###########################
# Base R Shiny image
FROM rocker/shiny

# Make a directory in the container
RUN mkdir /home/shiny-app

# Install R dependencies
RUN R -e "install.packages(c('dplyr', 'ggplot2', 'gapminder'))"

# Copy the Shiny app code
COPY example1.R /home/shiny-app/example1.R

# Expose the application port
EXPOSE 8180

# Run the R Shiny app
CMD Rscript /home/shiny-app/example1.R