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
    )
}

load_taxon_comments <- function(taxon_comments_url) {
  googlesheets4::read_sheet(taxon_comments_url) |>
    janitor::clean_names() |>
    select(taxon_id:comment)
}

load_emails <- function(email_url) {
  ppg_emails <- googlesheets4::read_sheet(email_url) |>
    janitor::clean_names() |>
    rename_with(
      ~ str_remove_all(., "_email"),
      matches("secondary|tertiary")
    ) |>
    select(
      name,
      institution,
      country,
      email,
      secondary,
      tertiary,
      status
    )
  
  ppg_emails |>
    dplyr::bind_rows(
      dplyr::select(ppg_emails, name, email = secondary)
    ) |>
    dplyr::bind_rows(
      dplyr::select(ppg_emails, name, email = tertiary)
    ) |>
    dplyr::filter(!is.na(email)) |>
    unique() |>
    fix_email() |>
    assertr::assert(assertr::is_uniq, email) |>
    assertr::assert(assertr::not_na, name) |>
    select(email, name, institution, country, status)
}

tally_votes <- function(votes, ppg_emails, exclude_emails) {
  author_institutions <- ppg_emails |>
    select(email, name, institution, country, status) |>
    distinct() |>
    group_by(name) |>
    summarize(
      # For each name, take the first non-NA value for each column
      email = first(email[!is.na(email)]),
      institution = first(institution[!is.na(institution)]),
      country = first(country[!is.na(country)]),
      status = first(status[!is.na(status)]),
      .groups = "drop"
    ) |>
    select(name, institution, country, status) |>
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
        c(institution, country, status),
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
  # Add manually specified authors
  manual_data <- ppg_emails |>
    inner_join(manual_authors, by = join_by(name))

  vote_tally |>
    select(-n_votes) |>
    bind_rows(manual_data) |>
    unique() |>
    assert(not_na, name, institution) |>
    assert(is_uniq, name) |>
    select(name, institution, country)
}

tar_write_csv <- function(x, file, ...) {
  write_csv(x = x, file = file, ...)
  file
}
