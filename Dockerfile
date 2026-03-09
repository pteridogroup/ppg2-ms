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
  git \
  curl \
  && rm -rf /var/lib/apt/lists/*

# Install Quarto (with architecture detection)
# Quarto bundles its own pandoc, so no need to install system pandoc
ARG QUARTO_VERSION=1.5.57
ARG TARGETARCH
RUN if [ "$TARGETARCH" = "arm64" ]; then \
  QUARTO_ARCH="arm64"; \
  else \
  QUARTO_ARCH="amd64"; \
  fi && \
  curl -LO https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-${QUARTO_ARCH}.deb && \
  dpkg -i quarto-${QUARTO_VERSION}-linux-${QUARTO_ARCH}.deb && \
  rm quarto-${QUARTO_VERSION}-linux-${QUARTO_ARCH}.deb

# Set working directory
WORKDIR /project

# Set renv paths to persist packages outside /project
# This prevents them from being overwritten when mounting volumes
ENV RENV_PATHS_LIBRARY=/usr/local/lib/R/renv-library
ENV RENV_PATHS_CACHE=/usr/local/lib/R/renv-cache

# Copy renv files first for better caching
COPY renv.lock renv.lock
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json

# Install renv and restore R package environment to persistent location
RUN R -e "install.packages('renv', repos='https://cloud.r-project.org')" \
  && R -e "renv::restore()"

# Copy project files
COPY . .

# Create non-root user for security
RUN useradd -m -u 1000 rstudio \
  && chown -R rstudio:rstudio /project \
  && chmod -R 755 /usr/local/lib/R/renv-library \
  && chmod -R 755 /usr/local/lib/R/renv-cache

# Switch to non-root user
USER rstudio

# Set default command to run the workflow
# The .Rprofile in the mounted directory will activate renv automatically
CMD ["R", "--no-save", "-e", "targets::tar_make()"]

