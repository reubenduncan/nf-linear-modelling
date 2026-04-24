FROM rocker/r-ver:4.3.2

LABEL maintainer="Reuben Duncan <reuben.duncan25@outlook.com>"
LABEL description="ecologyflow-lm: CODA-GLMNET and best-subset regression for microbiome data"

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libhdf5-dev \
        zlib1g-dev \
        libgit2-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        libglpk-dev \
        libgmp-dev \
        libmpfr-dev \
        libzstd-dev \
        liblz4-dev \
    && rm -rf /var/lib/apt/lists/*

# Set CRAN mirror
RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org"))' > /etc/R/Rprofile.site

# Install BiocManager
RUN Rscript -e "install.packages('BiocManager', dependencies = TRUE)"

# Install CRAN packages
RUN Rscript -e " \
    pkgs <- c('optparse', 'leaps', 'dplyr', 'caret', 'purrr', 'stringr', 'coda4microbiome'); \
    install.packages(pkgs, dependencies = TRUE) \
"

# Install Bioconductor packages
RUN Rscript -e " \
    BiocManager::install(c('phyloseq', 'mixOmics'), ask = FALSE, update = FALSE) \
"

# Install arrow (pre-built C++ library; LIBARROW_BINARY avoids 30-min source compile)
RUN LIBARROW_BINARY=true Rscript -e " \
    install.packages('arrow', repos = 'https://cloud.r-project.org', \
        Ncpus = max(1L, parallel::detectCores() - 1L)) \
"

# Verify key packages load
RUN Rscript -e " \
    library(optparse); \
    library(leaps); \
    library(dplyr); \
    library(caret); \
    library(purrr); \
    library(stringr); \
    library(coda4microbiome); \
    library(phyloseq); \
    library(mixOmics); \
    library(arrow); \
    message('All packages verified successfully') \
"

# Copy R scripts into container
COPY src/ /opt/ecology-scripts/src/

# Default working directory
WORKDIR /data

CMD ["R"]
