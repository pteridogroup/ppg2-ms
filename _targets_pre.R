source("R/packages.R")
source("R/functions_pre.R")

Sys.setenv(TAR_PROJECT = "pre")

# Pre-work for main workflow. The main purpose is to generate a list of
# PPG 2 authors and to fetch the taxon comments.

tar_plan(
  # Load emails to exclude from author list (can't identify owner, etc)
  tar_file_read(
    exclude_emails,
    "data/emails_exclude.csv",
    read_csv(!!.x)
  ),
  # Load CSV of URLs for Google ballots
  tar_file_read(
    ballot_url,
    "data/ballot-url.csv",
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
  author_list = format_author_list(vote_tally, ppg_emails),
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
)
