# PPG II Manuscript

This repository contains the source code and data for the manuscript "An
update of the Pteridophyte Phylogeny Group classification for extant
ferns and lycophytes: PPG II" (submitted).

## Overview

The [Pteridophyte Phylogeny Group (PPG)](https://pteridogroup.github.io/) is an international collaborative
effort to develop a comprehensive and phylogenetically-based
classification of ferns and lycophytes. This repository contains all
code and data needed to reproduce the manuscript describing PPG II, an
updated classification building on [PPG I (2016)](https://doi.org/10.1111/jse.12229).

## Repository Structure

- `ppg2_ms.Qmd`: Main manuscript file (Quarto format)
- `fig_s1.qmd`, `table_s1.qmd`: Supplementary materials
- `R/`: R code for data processing and analysis
  - `functions.R`: Custom functions for taxonomy processing, plotting,
    and phylogenetic tree manipulation
  - `packages.R`: Package dependencies
- `data/`: Input data files
  - `compare_ppg_wf_2025-12-21.csv`: Comparison between PPG II and
    [World Ferns](https://www.worldplants.de/world-ferns/ferns-and-lycophytes-list)
  - `ppg2_author_list.csv`: List of manuscript authors
  - `ppg_comments.csv`: Comments and notes for specific taxa
  - `ppg_issues_edited.csv`: Manually categorized taxonomic proposals
  - `ppgi_taxonomy_original.csv`: Original PPG I classification (2016)
  - `species_count_updates.csv`: Manual corrections to species counts
  - `table_uncertainty.csv`: Descriptions of uncertain phylogenetic
    nodes
  - `wf_taxa_count.csv`: Taxon counts from World Ferns database
- `tests/`: Unit tests for data validation
- `_targets.R`: Main workflow orchestration using the
  [`targets`](https://docs.ropensci.org/targets/) package
- `_targets_pre.R`: Pre-workflow for generating data files from raw data (not all raw data made public due to personally identifiable information)
- `images/`: Image files used in the manuscript (not figures produced by
  analysis)
- `renv.lock`: Package dependency specification for reproducibility

## Reproducibility

This project uses the [`targets`](https://docs.ropensci.org/targets/)
package for reproducible workflow management. The workflow automatically:

1. Downloads the latest PPG taxonomy from the [PPG
repository](https://github.com/pteridogroup/ppg)
2. Fetches and processes taxonomic proposals from GitHub Issues
3. Validates taxonomic changes against the previous classification
4. Generates figures and tables
5. Renders the manuscript

Once the workflow has successfully finished, the following files should be produced:

- `ppg2_ms.docx`: Manuscript file
- `fig_1.pdf`: Figure 1
- `fig_2.pdf`: Figure 2
- `fig_s1.pdf`: Figure S1
- `table_s1.pdf`: Table S1

### Requirements

- R (≥ 4.0.0)
- [Quarto](https://quarto.org/docs/get-started/) (≥ 1.4) - document
  rendering system
- [`renv`](https://rstudio.github.io/renv/) for R package management

All R package dependencies are managed with `renv` and recorded in
`renv.lock`.

### Running the Analysis

To reproduce the analysis and manuscript:

1. Install [R](https://www.r-project.org/) (≥ 4.0.0) and
   [Quarto](https://quarto.org/docs/get-started/) (≥ 1.4) if not
   already installed

2. Set up the R environment:

```r
# Install renv if needed
install.packages("renv")

# Restore the package environment
renv::restore()
```

3. Run the workflow:

```r
# Run the workflow
targets::tar_make()
```

## Contributing

This repository contains the manuscript code. For taxonomic proposals
and discussions, please visit the main [PPG
repository](https://github.com/pteridogroup/ppg).

## License

The code in this repository is licensed under the [MIT License](LICENSE). Please
cite the published manuscript when using this work.

## Authors

This is a collaborative work by the Pteridophyte Phylogeny Group. See
the manuscript for the complete author list.

## Citation

*Citation information will be added upon publication.*

## Contact

For questions about the manuscript or code, please [open an issue](https://github.com/pteridogroup/ppg2-ms/issues/new) in this
repository.
