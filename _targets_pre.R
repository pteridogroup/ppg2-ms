source("R/packages.R")
source("R/functions_pre.R")

Sys.setenv(TAR_PROJECT = "pre")

# Pre-work for main workflow. This generates various data files in data/ used
# in the main workflow. The raw data files include PII and are not made public.

tar_plan(
  # Load emails to exclude from author list (can't identify owner, etc)
  tar_file_read(
    exclude_emails,
    "data_raw/emails_exclude.csv",
    read_csv(!!.x)
  ),
  # Load CSV of URLs for Google ballots
  tar_file_read(
    ballot_url,
    "data_raw/ballot-url.csv",
    read_csv(!!.x)
  ),
  # Load voting data
  votes = load_votes(ballot_url),
  # Load email list
  # always run, to catch up any updates (tar_url doesn't like google links)
  tar_target(
    email_url,
    paste0(
      "https://docs.google.com/spreadsheets/d/",
      "1vxlmf8QPndiE6dIeDcjoFT7GA3ZE4pSc_Z1raf4iuwA/edit?usp=sharing"
    )
  ),
  tar_target(
    ppg_emails,
    load_emails(email_url),
    cue = tar_cue("always")
  ),
  # Tally votes
  vote_tally = tally_votes(votes, ppg_emails, exclude_emails),
  # Format author list
  manual_authors = format_manual_authors(),
  author_list = format_author_list(vote_tally, ppg_emails, manual_authors),
  # Write out author list
  tar_file(
    author_csv,
    tar_write_csv(author_list, "data/ppg2_author_list.csv")
  ),
  # Taxon comments
  tar_target(
    taxon_comments_url,
    paste0(
      "https://docs.google.com/spreadsheets/d/",
      "1LJIQcEQkGmqh3e99Q3IqHJczUDv2OD9rvpRS3RrP4io/edit?usp=sharing"
    )
  ),
  tar_target(
    taxon_comments,
    load_taxon_comments(taxon_comments_url),
    cue = tar_cue("always")
  ),
  tar_file(
    taxon_comments_csv,
    tar_write_csv(taxon_comments, "data/ppg_comments.csv")
  ),
  # Issues
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
  # - Write out CSV files
  tar_file(
    ppg_issues_csv,
    tar_write_csv(ppg_issues, "data/ppg_issues.csv")
  ),
  tar_file(
    commenters_csv,
    tar_write_csv(commenters, "data/commenters.csv")
  ),
  tar_file(
    voting_results_csv,
    tar_write_csv(voting_results, "data/voting_results.csv")
  ),
  tar_file(
    ppg_issue_count_raw_csv,
    tar_write_csv(voting_results, "data_raw/ppg_issues_type_raw.csv")
  )
)
