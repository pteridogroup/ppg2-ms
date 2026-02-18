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
    rename_with(
      ~ str_remove_all(., "_email"),
      matches("secondary|tertiary")
    ) |>
    rename_with(
      ~ str_remove_all(., "_s"),
      matches("name")
    ) |>
    select(
      given_name,
      surname,
      name,
      ppg2,
      institution,
      country,
      email,
      secondary,
      tertiary
    ) |>
    filter(ppg2 == "Yes") |>
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
