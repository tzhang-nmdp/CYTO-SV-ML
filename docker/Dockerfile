###############################################################################
# Dockerfile for CYTO-SV-ML R Shiny Web Portal
#
# Builds a Docker image containing an R Shiny application for interactive
# exploration of SV classification results.
#
# Base image: rocker/shiny:3.6.1 (R 3.6.1 with Shiny Server)
#
# Build:
#   sudo docker build -t cyto-sv-ml-app:<cohort_name> .
#
# Run:
#   sudo docker run -d -p 8000:80 cyto-sv-ml-app:<cohort_name>
#   # Open http://localhost:8000/ in a browser
#
# Disk space: ~2 GB for the built image
###############################################################################

FROM docker.io/rocker/shiny:3.6.1
MAINTAINER "Tao Zhang tzhang@nmdp.org"

# ---- System dependencies ----
# Use archived Debian stretch repos (EOL distribution)
RUN echo "deb http://archive.debian.org/debian stretch main" > /etc/apt/sources.list

# Install Java (required for some R packages) and SSL libraries
RUN apt-get update && apt-get install -y libxtst6
RUN apt-get update && apt-get install -y openjdk-8-jdk libssl-dev && apt-get autoremove

# Configure R to find Java
RUN R CMD javareconf

# ---- R package installation ----
# Install remotes for GitHub package installation
RUN install2.r remotes

# lineupjs: interactive ranking visualization widget
RUN R -e 'install.packages("lineupjs", dependencies=TRUE, repos="https://cran.r-project.org/")'

# Core Shiny app dependencies
RUN install2.r ggforce plotly DT ggplot2 stringr shinydashboard shinyjs shinyWidgets dplyr Hmisc

# ---- Application deployment ----
# Copy Shiny Server configuration
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf

# Copy the Shiny app (includes data/cyto_sv_ml.RData, ui.R, server.R)
COPY /cyto-sv-ml /srv/shiny-server/

# Set ownership for the shiny user
RUN chown -R shiny:shiny /srv/shiny-server/

# Expose port 80 for the Shiny Server
EXPOSE 80

# Copy and set the startup script
COPY shiny-server.sh /usr/bin/shiny-server.sh

CMD ["/usr/bin/shiny-server.sh"]
