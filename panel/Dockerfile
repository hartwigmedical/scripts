# Use official R base image
FROM r-base:4.3.2

# Set working directory
WORKDIR /app

# Install system dependencies required by R packages
RUN apt-get update && apt-get install -y \
    cmake \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    zlib1g-dev \
    libcairo2-dev \
    libgit2-dev \
    libgmp-dev \
    libmpfr-dev \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager first
RUN R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/', dependencies=TRUE)" && \
    R -e "if (!require('BiocManager')) stop('BiocManager installation failed')"

# Install CRAN packages with dependencies and verification
RUN R -e "install.packages('tidyr', repos='http://cran.rstudio.com/', dependencies=TRUE)" && \
    R -e "if (!require('tidyr')) stop('tidyr installation failed')"

RUN R -e "install.packages('dplyr', repos='http://cran.rstudio.com/', dependencies=TRUE)" && \
    R -e "if (!require('dplyr')) stop('dplyr installation failed')"

RUN R -e "install.packages('stringi', repos='http://cran.rstudio.com/', dependencies=TRUE)" && \
    R -e "if (!require('stringi')) stop('stringi installation failed')"

RUN R -e "install.packages('ggplot2', repos='http://cran.rstudio.com/', dependencies=TRUE)" && \
    R -e "if (!require('ggplot2')) stop('ggplot2 installation failed')"

RUN R -e "install.packages('gridExtra', repos='http://cran.rstudio.com/', dependencies=TRUE)" && \
    R -e "if (!require('gridExtra')) stop('gridExtra installation failed')"

RUN R -e "install.packages('cowplot', repos='http://cran.rstudio.com/', dependencies=TRUE)" && \
    R -e "if (!require('cowplot')) stop('cowplot installation failed')"

RUN R -e "install.packages('ggpubr', repos='http://cran.rstudio.com/', dependencies=TRUE)" && \
    R -e "if (!require('ggpubr')) stop('ggpubr installation failed')"

# Install Bioconductor packages
RUN R -e "BiocManager::install('VariantAnnotation', ask=FALSE, update=TRUE)" && \
    R -e "if (!require('VariantAnnotation')) stop('VariantAnnotation installation failed')"

# Copy the R script
COPY createSampleQcReport.R /app/

# Make the script executable
RUN chmod +x /app/createSampleQcReport.R
RUN chmod +x /app/createSampleQcReport_Deamination.R

# Set the entrypoint to the R script
ENTRYPOINT ["/app/createSampleQcReport.R"]

# Default command shows help if no arguments provided
CMD []
