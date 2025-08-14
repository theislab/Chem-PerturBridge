FROM python:3.11

# Install system dependencies for R
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libblas-dev \
    liblapack-dev \
    gfortran \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    pkg-config \
    vim \
    && rm -rf /var/lib/apt/lists/*


# Install Jupyter for Python
RUN pip install --no-cache-dir notebook lamindb scanpy gcsfs rpy2 anndata2ri decoupler pubchempy

# Install R's Jupyter kernel
RUN R -e "install.packages('IRkernel', repos='http://cran.r-project.org')" \
    && R -e "IRkernel::installspec(user = FALSE)" \
    && R -e "install.packages('BiocManager')" \
    && R -e "BiocManager::install('scDblFinder')"

# Set working directory
WORKDIR /home/ubuntu 

# Expose Jupyter port
EXPOSE 4222
