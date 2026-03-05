# Use specific R version for reproducibility
FROM rocker/r-ver:4.5.0

# Install system dependencies needed for R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
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
  libglpk-dev \
  libcairo2-dev \
  libmagick++-dev \
  libx11-dev \
  pandoc \
  git \
  curl \
  && rm -rf /var/lib/apt/lists/*

# Install Quarto
ARG QUARTO_VERSION=1.5.57
RUN curl -LO https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb \
  && dpkg -i quarto-${QUARTO_VERSION}-linux-amd64.deb \
  && rm quarto-${QUARTO_VERSION}-linux-amd64.deb

# Set working directory
WORKDIR /project

# Copy renv files first for better caching
COPY renv.lock renv.lock
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json

# Install renv and restore R package environment
RUN R -e "install.packages('renv', repos='https://cloud.r-project.org')" \
  && R -e "renv::restore()"

# Copy project files
COPY . .

# Create non-root user for security
RUN useradd -m -u 1000 rstudio \
  && chown -R rstudio:rstudio /project

# Switch to non-root user
USER rstudio

# Set default command to run the workflow
CMD ["R", "-e", "targets::tar_make()"]
