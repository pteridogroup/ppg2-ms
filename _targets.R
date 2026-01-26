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
  ppg_issues_count_raw = count_issues(ppg_issues),

  # - Load manual issues count by type (sink/split)
  tar_file_read(
    ppg_issues_count,
    "data/ppg_issues_edited.csv",
    read_csv(!!.x)
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

  # Generate figures ----

  # - proposal count
  issues_count_plot = make_issues_plot(ppg_issues),

  # - tree figure
  tree_fig = make_tree_figure(phy_family, ppg, ppg_tl, children_tally),

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
  tar_quarto(
    ppg2_ms,
    "ppg2_ms.Qmd",
    quiet = FALSE
  )
)
