# Use the specified Miniconda base image
FROM continuumio/miniconda3:23.5.2-0

# Set environment variables to non-interactive and define Quarto version
ENV DEBIAN_FRONTEND=noninteractive \
    QUARTO_VERSION=1.4.553 \
    LC_ALL=C \
    PATH=/opt/conda/bin/:$PATH

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    pandoc \
    wget \
    original-awk \
    wget \
    pandoc-citeproc \
    fonts-liberation \
    curl \
    gdebi-core \
    librsvg2-bin \
    texlive-xetex \
    texlive-fonts-recommended \
    texlive-plain-generic \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Determine system architecture and download the appropriate Quarto version
RUN ARCH=$(uname -m) && \
    if [ "$ARCH" = "x86_64" ]; then \
        QUARTO_ARCH="amd64"; \
    elif [ "$ARCH" = "aarch64" ]; then \
        QUARTO_ARCH="arm64"; \
    else \
        echo "Unsupported architecture: $ARCH"; \
        exit 1; \
    fi && \
    curl -o quarto-linux-$QUARTO_ARCH.deb -L https://github.com/quarto-dev/quarto-cli/releases/download/v$QUARTO_VERSION/quarto-$QUARTO_VERSION-linux-$QUARTO_ARCH.deb \
    && dpkg --add-architecture $QUARTO_ARCH \
    && gdebi --non-interactive quarto-linux-$QUARTO_ARCH.deb \
    && rm -f quarto-linux-$QUARTO_ARCH.deb

# Install TinyTeX and install packages with tlmgr
RUN wget -qO- "https://yihui.org/tinytex/install-bin-unix.sh" | sh && \
    ARCH=$(uname -m) && \
    if [ "$ARCH" = "x86_64" ]; then \
        PATH=/root/.TinyTeX/bin/x86_64-linux:$PATH; \
    elif [ "$ARCH" = "aarch64" ]; then \
        PATH=/root/.TinyTeX/bin/aarch64-linux:$PATH; \
    fi && \
    tlmgr install fvextra footnotebackref pagecolor sourcesanspro sourcecodepro titling

USER root

RUN set -e -x && \
    mkdir -p /input

# Configure conda channels
RUN conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge

# Install Mamba into the base environment = faster solver
RUN /opt/conda/bin/conda install mamba -n base -c conda-forge -y

# Install Jupyter and scientific libraries with conda
RUN mamba install jupyterlab==4.0.11 \
       nbconvert \
       r-base==4.3.2 \
       r-reticulate \
       r-shiny \
       r-htmltools \
       r-knitr \
       r-jsonlite \
       r-rmarkdown \
       r-renv \
       r-remotes \
       r-vegan \
       r-codetools \
       r-irkernel \
       r-tidyverse==2.0.0 \
       r-ggpubr \
       r-ggbreak \
       r-patchwork \
       r-egg \
       r-scales \
       r-pheatmap \
       r-ggpubr \
       r-glue \
       r-cnogpro \
       r-corrplot \
       r-ggrepel \
       r-ape \
       r-ggbeeswarm \
       r-viridislite \
       r-lme4 \
       rpy2 \
       itables \
       matplotlib-venn \
       openpyxl \
       plotnine==0.13.0 \
       pandas==2.2.0 \ 
       seaborn==0.13.2 \
       numpy==1.26.3 \
       scipy==1.12.0 \
       matplotlib==3.8.2 \
       numpy==1.26.3 \
       scikit-learn==1.4.0 \
       plotly==5.18.0 \
       statsmodels==0.14.1 \
       scikit-bio==0.6.2 \
       venn \
       Pyarrow \
       biopython \
       -y --quiet

RUN pip3 install venn
# Set a CRAN mirror and install BiocManager
RUN R -e "install.packages(c('BiocManager','viridis', 'CNOGpro', 'glue', 'ggrepel','Hmisc', 'ape','paletteer'), repos='http://cran.rstudio.com/')"

# Install packages using BiocManager
RUN R -e "BiocManager::install(c('edgeR', 'ggtree', 'ggtreeExtra','treeio'))"
# Final cleanup
RUN apt-get autoremove -y \
    && apt-get clean

# Docker file adapted from https://dev.to/cavo789/running-quarto-markdown-in-docker-4igj