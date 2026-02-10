source("R/packages.R")
source("R/functions.R")

Sys.setenv(TAR_PROJECT = "main")

tar_plan(
  # Download PPG taxonomy ----
  tar_url(
    ppg_csv_url,
    "https://github.com/pteridogroup/ppg/raw/refs/heads/main/data/ppg.csv"
  ),
  ppg_raw = read_csv(ppg_csv_url),
  # Clean PPG data
  # - remove invalid nomenclatural status and unchecked taxa
  # - remove nothogenera that haven't received votes
  ppg = clean_ppg(ppg_raw),

  # Load PPG I taxonomy ----
  tar_file_read(
    ppg_i,
    "data/ppgi_taxonomy_original.csv",
    read_csv(!!.x)
  ),

  # Load PPG II vs World Ferns comparison
  tar_file_read(
    ppg_ii_vs_wf_raw,
    "data/compare_ppg_wf_2025-12-21.csv",
    read_csv(!!.x)
  ),

  ppg_ii_vs_wf = clean_wf_comp_list(ppg_ii_vs_wf_raw),

  # Load WF taxa count

  tar_file_read(
    wf_taxa_count,
    "data/wf_taxa_count.csv",
    read_csv(!!.x)
  ),

  # Process GitHub issues ----

  # - Download issues
  tar_url(
    ppg_repo_url,
    "https://github.com/pteridogroup/ppg/",
  ),

  tar_target(
    ppg_issues_raw,
    fetch_issues(ppg_repo_url),
    cue = tar_cue("always")
  ),

  # - Remove invalid and TBD issues
  ppg_issues = remove_invalid_issues(ppg_issues_raw),

  # - Get github user names of commenters for each valid issue
  commenters = fetch_commenters(ppg_issues$number),

  # - Get voting results for each valid issue
  voting_results = fetch_voting_results(ppg_issues$number),

  # - Generate initial automatic count of isses by type (sink/split)
  ppg_issue_count_raw = count_issues(ppg_issues),

  # - Load manual issues count by type (sink/split)
  tar_file_read(
    ppg_issues_count,
    "data/ppg_issues_edited.csv",
    read_csv(!!.x)
  ),

  # - Check that all passed issues have been included in issue type count
  issue_type_count_check = check_issue_type_count(
    ppg_issues,
    ppg_issues_count
  ),

  # Make family-level tree ----
  phy_family = make_family_tree(),
  # Get families in 'phylogenetic' order
  families_in_phy_order = get_ladderized_tips(phy_family),

  # Load taxon comments ----
  tar_file_read(
    comments_raw,
    "data/ppg_comments.csv",
    read_csv(!!.x)
  ),
  comments = format_ppg_comments(comments_raw),

  # Format data ----
  # - convert DarwinCore format to dataframe in taxonomic order for printing
  #   (only includes accepted taxa at genus and higher)
  ppg_tl = dwc_to_tl(ppg, families_in_phy_order),
  # - count number of children taxa, while excluding any hybrids
  #   (names with '×')
  children_tally = count_children_ppg(
    ppg,
    families_in_phy_order,
    exclude_hybrids = TRUE
  ),
  # - count number of taxa per rank in PPG II
  ppg_2_taxa_count = count_ppg2_taxa(ppg, exclude_hybrids = TRUE),
  ppg_i_taxa_count = count_ppgi(ppg_i),

  # Compare PPG I and PPG II ----
  # - Convert PPG II to wide hierarchical format
  ppg_ii = ppg_to_wide(ppg),
  # - Check for classification changes and verify they have passed
  #   issues
  genus_check = check_ppg_classification_changes(
    ppg_ii,
    ppg_i,
    ppg_issues,
    ppg
  ),
  # - Summarize classification check results
  classification_summary = summarize_classification_check(
    genus_check
  ),

  # Check higher taxonomic ranks for changes
  higher_taxa_check = check_ppg_higher_tax_changes(
    ppg_ii,
    ppg_i,
    ppg_issues
  ),

  # Generate figures ----

  # - proposal count
  issues_count_plot = make_issues_plot(ppg_issues),

  # - tree figure
  tree_fig = make_tree_figure(phy_family, ppg, ppg_tl, children_tally),

  tree_fig_appendix = make_tree_figure_appendix(
    phy_family,
    ppg,
    ppg_ii,
    ppg_issues_count
  ),

  # Format classification ----
  ppg_classification = format_ppg_classification(
    ppg,
    ppg_tl,
    children_tally,
    comments,
    families_in_phy_order
  ),

  # Check that all PPG I taxa are accounted for ----
  ppg_i_check = check_ppg_i_taxa_accounted(
    ppg_i,
    ppg_classification,
    comments_path = "data/ppg_comments.csv"
  ),

  # Output manuscript ----
  tar_file_read(
    uncertainty_table,
    "data/table_uncertainty.csv",
    read_csv(!!.x)
  ),

  tar_file(
    references,
    "references.yaml"
  ),

  tar_quarto(
    ppg2_ms,
    "ppg2_ms.Qmd",
    quiet = FALSE
  ),

  tar_quarto(
    appendix,
    "appendix.qmd",
    quiet = FALSE
  ),

  tar_quarto(
    appendix_2,
    "appendix_2.qmd",
    quiet = FALSE
  )
)
