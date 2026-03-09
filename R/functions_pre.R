# Fix malformed emails in voting data
fix_email <- function(data) {
  data |>
    dplyr::mutate(
      # Automatically fix typos in email address
      email = tolower(email) |>
        stringr::str_replace_all("\\.ocm$", ".com") |>
        stringr::str_replace_all("\\.combr$", ".com") |>
        stringr::str_replace_all(
          "dar\\.sanin\\@gmail\\.com$",
          "dav.sanin@gmail.com"
        ) |>
        stringr::str_replace_all(
          "fittmatos\\@gmail\\.com$",
          "fbittmatos@gmail.com"
        ) |>
        stringr::str_replace_all("\\.utexxas\\.", ".utexas.") |>
        stringr::str_replace_all(
          "ralf\\.knap\\@gmail\\.com",
          "ralf.knapp@gmail.com"
        ) |>
        stringr::str_replace_all(
          "zuozhengyu\\@mail\\.kib\\.ac\\.cn",
          "zuozhengyu@sdau.edu.cn"
        )
    )
}

load_votes <- function(ballot_url) {
  ballot_url |>
    mutate(
      res = map(url, googlesheets4::read_sheet)
    ) |>
    select(-url) |>
    unnest(res) |>
    select(
      ballot,
      timestamp = Timestamp,
      email = `Email Address`,
      gh_username = `GitHub username`
    ) |>
    mutate(email = str_to_lower(email))
}

load_taxon_comments <- function(taxon_comments_url) {
  googlesheets4::read_sheet(taxon_comments_url) |>
    janitor::clean_names() |>
    select(taxon_id:comment)
}

load_emails <- function(email_url) {
  ppg_emails <- googlesheets4::read_sheet(email_url) |>
    janitor::clean_names() |>
    rename(
      email = email_our_ppg,
      secondary = email_official,
      tertiary = email_other
    ) |>
    rename_with(
      ~ str_remove_all(., "_s"),
      matches("name")
    ) |>
    select(
      given_name,
      surname,
      ppg2,
      institution,
      country,
      email,
      secondary,
      tertiary
    ) |>
    filter(ppg2 == "Yes") |>
    mutate(name = paste(given_name, surname)) |>
    assert(is_uniq, name) |>
    assert(not_na, name)

  emails_long <-
    ppg_emails |>
    select(name, email, secondary, tertiary) |>
    pivot_longer(
      values_to = "email",
      names_to = "which",
      -name
    ) |>
    filter(!is.na(email)) |>
    select(-which)

  ppg_emails |>
    select(given_name, surname, name, institution, country) |>
    left_join(
      emails_long,
      by = join_by(name),
      relationship = "one-to-many"
    ) |>
    mutate(email = str_to_lower(email)) |>
    select(email, given_name, surname, name, institution, country)
}

tally_votes <- function(votes, ppg_emails, exclude_emails) {
  author_institutions <- ppg_emails |>
    select(email, name, institution, country) |>
    distinct() |>
    group_by(name) |>
    summarize(
      # For each name, take the first non-NA value for each column
      email = first(email[!is.na(email)]),
      institution = first(institution[!is.na(institution)]),
      country = first(country[!is.na(country)]),
      .groups = "drop"
    ) |>
    select(name, institution, country) |>
    distinct() |>
    assert(not_na, name) |>
    assert(is_uniq, name)

  votes |>
    fix_email() |>
    anti_join(
      exclude_emails,
      by = "email"
    ) |>
    count(email) |>
    rename(n_votes = n) |>
    left_join(
      select(ppg_emails, email, name),
      by = "email",
      relationship = "one-to-one"
    ) |>
    assert(not_na, name) |>
    select(-email) |>
    group_by(name) |>
    summarize(n_votes = sum(n_votes), .groups = "drop") |>
    left_join(
      author_institutions,
      by = "name",
      relationship = "one-to-one"
    ) |>
    assert(not_na, name, n_votes) |>
    assert(is_uniq, name) |>
    mutate(
      across(
        c(institution, country),
        ~ tidyr::replace_na(.x, "?")
      )
    )
}

format_manual_authors <- function() {
  tibble(
    name = c(
      "Bertrand Black",
      "Muhammad Irfan",
      "Jean-Yves Dubuisson",
      "Sabine Hennequin",
      "Jovani B. de S. Pereira",
      "Shi-Yong Dong"
    )
  )
}

format_author_list <- function(vote_tally, ppg_emails, manual_authors) {
  # Account for all manually added authors and all authors who voted
  assertthat::assert_that(
    isTRUE(all(manual_authors$name %in% ppg_emails$name)),
    msg = "A manually specified author is missing from the author list"
  )

  assertthat::assert_that(
    isTRUE(all(vote_tally$name %in% ppg_emails$name)),
    msg = "A person who voted is missing from the author list"
  )

  ppg_emails |>
    select(given_name, surname, name, institution, country) |>
    unique() |>
    # Make sure spelling of given_name + surname = name
    mutate(
      name_check = paste(given_name, surname) |> str_squish(),
      name_check = name == name_check
    ) |>
    assert(in_set(TRUE), name_check) |>
    select(-name_check) |>
    assert(is_uniq, name) |>
    assert(not_na, everything()) |>
    arrange(name)
}

tar_write_csv <- function(x, file, ...) {
  write_csv(x = x, file = file, ...)
  file
}


# Extract useful information to dataframe
fetch_issues <- function(repo_url, n_max = 1000) {
  repo <- str_match(
    repo_url,
    "github\\.com\\/(.*)"
  ) |>
    magrittr::extract2(2) |>
    str_remove_all("\\/$")

  issues_json <-
    glue::glue(
      "https://api.github.com/repos/{repo}/issues?state=all&page=1&per_page={n_max}"
    ) |>
    jsonlite::fromJSON()

  # Create initial tibble of issues (may include PRs)
  issues_df <- tibble::tibble(
    number = issues_json$number,
    title = issues_json$title,
    url = issues_json$url,
    created_at = issues_json$created_at,
    user = issues_json$user$login,
    state = issues_json$state,
    body = issues_json$body
  )

  if (nrow(issues_df) == n_max) {
    stop("Maximum number of issues fetched; increase n_max")
  }

  # If any PRs exist, remove them
  if (!is.null(issues_json$draft)) {
    issues_df <-
      issues_df |>
      dplyr::mutate(draft = issues_json$draft) |>
      dplyr::filter(is.na(draft)) |>
      dplyr::select(-draft)
  }

  # Format final data frame
  issues_df |>
    dplyr::mutate(
      url = stringr::str_replace_all(
        url,
        "https://api.github.com/repos/",
        "https://github.com/"
      ),
      name = stringr::str_match(body, "Name of taxon[\r|\n]*(×*.*)") |>
        magrittr::extract(, 2),
      rank = stringr::str_match(body, "Rank of taxon[\r|\n]*(\\w+)[\r|\n]*") |>
        magrittr::extract(, 2),
      no_species = stringr::str_match(
        body,
        "number of species affected[\r|\n]*(.*)"
      ) |>
        magrittr::extract(, 2),
      description = stringr::str_match(
        body,
        "Description of change[\r|\n]*(.*)"
      ) |>
        magrittr::extract(, 2)
    ) |>
    dplyr::select(-body)
}


# Make tibble of commenters on an issue
fetch_commenter_single <- function(issue_num) {
  data <- gh::gh(
    "GET /repos/pteridogroup/ppg/issues/{issue_num}/timeline",
    issue_num = issue_num,
    .limit = Inf
  )

  commenters <- purrr::map(data, ~ purrr::pluck(.x, "user", "login")) |>
    purrr::flatten() |>
    unlist()

  tibble::tibble(
    issue = issue_num,
    commenter = commenters
  )
}

fetch_commenters <- function(issue_num) {
  issue_num |>
    unique() |>
    purrr::map_df(fetch_commenter_single)
}


# Make tibble of voting results on an issue
fetch_voting_results_single <- function(issue_num) {
  data <- gh::gh(
    "GET /repos/pteridogroup/ppg/issues/{issue_num}/timeline",
    issue_num = issue_num,
    .limit = Inf
  )

  results <- purrr::map_df(data, function(x) {
    # Only process comments that have a body field
    if (
      !is.null(x$body) &&
        stringr::str_detect(x$body, "voted on during PPG Ballot")
    ) {
      tibble::tibble(
        total_votes = stringr::str_extract(x$body, "A total of (\\d+) votes") |>
          stringr::str_extract("\\d+") |>
          as.numeric(),
        yes_votes = stringr::str_extract(x$body, "(\\d+) 'Yes' votes") |>
          stringr::str_extract("\\d+") |>
          as.numeric(),
        yes_percent = stringr::str_extract(
          x$body,
          "'Yes' votes \\((\\d+(?:\\.\\d+)?)%\\)"
        ) |>
          stringr::str_extract("\\d+(?:\\.\\d+)?") |>
          as.numeric(),
        no_votes = stringr::str_extract(x$body, "(\\d+) 'No' votes") |>
          stringr::str_extract("\\d+") |>
          as.numeric(),
        no_percent = stringr::str_extract(
          x$body,
          "'No' votes \\((\\d+(?:\\.\\d+)?)%\\)"
        ) |>
          stringr::str_extract("\\d+(?:\\.\\d+)?") |>
          as.numeric(),
        ballot = stringr::str_extract(x$body, "PPG Ballot (\\d+)") |>
          stringr::str_extract("\\d+") |>
          as.numeric(),
        result = stringr::str_extract(
          x$body,
          "The proposal (passes|does not pass)\\."
        ) |>
          stringr::str_remove("The proposal ") |>
          stringr::str_remove("\\.")
      )
    } else {
      NULL
    }
  })

  # If no voting results found, return tibble with issue number and NAs
  if (nrow(results) == 0) {
    results <- tibble::tibble(
      total_votes = NA_real_,
      yes_votes = NA_real_,
      yes_percent = NA_real_,
      no_votes = NA_real_,
      no_percent = NA_real_,
      ballot = NA_real_,
      result = NA_character_
    )
  }

  # Add issue number to results
  results |>
    dplyr::mutate(issue = issue_num, .before = 1)
}

fetch_commenters <- function(issue_num) {
  issue_num |>
    unique() |>
    purrr::map_df(fetch_commenter_single)
}

fetch_voting_results <- function(issue_num) {
  issue_num |>
    unique() |>
    purrr::map_df(fetch_voting_results_single)
}

#' Remove Invalid Issues from Issue Tracker Data
#'
#' Filters a data frame of GitHub issues to retain only those with a valid
#' status designation. Valid statuses are "PASSED" or "NOT PASSED", which
#' should appear in square brackets at the end of the issue title (e.g.,
#' "[PASSED]").
#'
#' Also removes issues that have not been voted on yet ('TBD' issues).
#'
#' @param issues A data frame containing GitHub issue data, expected to have
#'   at least a `title` column containing the issue title and a `number`
#'   column for the issue number.
#'
#' @return A tibble containing only the issues with valid status
#'   designations (PASSED or NOT PASSED in the title).
remove_invalid_issues <- function(issues) {
  issues |>
    mutate(
      status = str_extract(title, "\\[([^\\]]+)\\]$"),
      status = str_remove_all(status, "\\[|\\]"),
      status = replace_na(status, "TBD")
    ) |>
    filter(!status %in% c("NOT VALID", "RETRACTED", "TBD"))
}

count_issues <- function(issues) {
  issues |>
    filter(str_detect(title, "\\[PASSED\\]")) |>
    mutate(
      name = str_remove_all(name, regex("and|,", ignore_case = TRUE)) |>
        str_squish()
    ) |>
    select(number, name, rank, description) |>
    mutate(
      change = case_when(
        str_detect(description, "should be sunk") ~ "sink",
        str_detect(description, "subsumed") ~ "sink",
        str_detect(description, fixed("recogni", ignore_case = TRUE)) ~ "split",
        str_detect(
          description,
          fixed("Resurrect", ignore_case = TRUE)
        ) ~ "split",
        str_detect(description, "to be a distinct group") ~ "split",
        str_detect(description, "accommodate") ~ "split",
        str_detect(description, "segregate") ~ "split",
        str_detect(description, "into synonymy") ~ "sink",
        str_detect(description, "transferred to") ~ "sink",
        str_detect(description, "accommodate") ~ "split",
        str_detect(description, "lump") ~ "sink",
        .default = NA_character_
      )
    ) |>
    separate_rows(name)
}