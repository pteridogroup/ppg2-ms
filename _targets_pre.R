source("R/packages.R")
source("R/functions_pre.R")

Sys.setenv(TAR_PROJECT = "pre")

# Pre-work for main workflow. The main purpose is to generate a list of
# PPG 2 authors.

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
  email_url = paste0(
    "https://docs.google.com/spreadsheets/d/",
    "1vxlmf8QPndiE6dIeDcjoFT7GA3ZE4pSc_Z1raf4iuwA/edit?usp=sharing"
  ),
  ppg_emails = load_emails(email_url),
  # Tally votes
  vote_tally = tally_votes(votes, ppg_emails, exclude_emails),
  # Format author list
  author_list = select(vote_tally, -n_votes),
  # Write out author list
  tar_file(
    author_csv,
    tar_write_csv(author_list, "data/ppg2_author_list.csv")
  )
)
