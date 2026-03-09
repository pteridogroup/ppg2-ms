# Use Rocker image with R and Quarto pre-installed
# This avoids architecture-specific Quarto installation issues
FROM rocker/verse:4.5.0

# Install system dependencies needed for R packages
# rocker/verse includes most dependencies, but add extras needed for this project
RUN apt-get update && apt-get install -y --no-install-recommends \
  libglpk-dev \
  libcairo2-dev \
  libmagick++-dev \
  libx11-dev \
  && rm -rf /var/lib/apt/lists/*

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

# Set permissions for rstudio user (already exists in rocker/verse)
RUN chown -R rstudio:rstudio /project \
  && chmod -R 755 /usr/local/lib/R/renv-library \
  && chmod -R 755 /usr/local/lib/R/renv-cache

# Switch to non-root user
USER rstudio

# Set default command to run the workflow
# The .Rprofile in the mounted directory will activate renv automatically
CMD ["R", "--no-save", "-e", "targets::tar_make()"]

