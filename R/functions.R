#' Extract Last Name from Full Name
#'
#' Extracts the last word from a person's name string, typically the
#' surname or family name.
#'
#' @param name A character string containing a person's full name.
#'
#' @return A character string containing the last name.
extract_last_name <- function(name) {
  # Split by spaces and take the last element
  parts <- stringr::str_split(name, "\\s+")[[1]]
  parts[length(parts)]
}

#' Make a Family-Level Fern Phylogeny
#'
#' Constructs a family-level phylogenetic tree for ferns using the FTOL
#' backbone and taxonomy, ensuring each family is represented by a single
#' exemplar species and that all families are monophyletic or monotypic.
#'
#' @return A ladderized phylogenetic tree object of class phylo with
#'   family names as tip labels.
make_family_tree <- function() {
  require(ftolr)
  # Load tree
  phy <- ftolr::ft_tree(branch_len = "ultra", rooted = TRUE, drop_og = TRUE)
  # Load fern taxonomy
  taxonomy <- ftol_taxonomy |>
    # Subset to only species in tree
    filter(species %in% phy$tip.label) |>
    select(species, family) |>
    # Make corrections for new families
    mutate(
      family = case_when(
        str_detect(species, "Arthropteris") ~ "Arthropteridaceae",
        str_detect(species, "Draconopteris") ~ "Pteridryaceae",
        str_detect(species, "Malaifilix") ~ "Pteridryaceae",
        str_detect(species, "Polydictyum") ~ "Pteridryaceae",
        str_detect(species, "Pteridrys") ~ "Pteridryaceae",
        .default = family
      )
    )

  # Analyze monophyly of each family
  family_mono_test <- MonoPhy::AssessMonophyly(
    phy,
    as.data.frame(taxonomy[, c("species", "family")])
  )

  # Check that all families are monophyletic or monotypic
  family_mono_summary <-
    family_mono_test$family$result |>
    tibble::rownames_to_column("family") |>
    as_tibble() |>
    assert(in_set("Yes", "Monotypic"), Monophyly)

  # Get one exemplar tip (species) per family
  rep_tips <-
    taxonomy |>
    group_by(family) |>
    slice(1) |>
    ungroup()

  # Subset phylogeny to one tip per family
  phy_family <- ape::keep.tip(phy, rep_tips$species)

  # Relabel with family names
  new_tips <-
    tibble(species = phy_family$tip.label) |>
    left_join(rep_tips, by = "species") |>
    pull(family)

  phy_family$tip.label <- new_tips

  ape::ladderize(phy_family)
}

#' Convert Darwin Core Taxonomy to Indented Taxonomic Data Frame
#'
#' Converts a cleaned PPG (Pteridophyte Phylogeny Group) Darwin Core
#' taxonomy data frame to an indented, ordered taxonomic data frame
#' suitable for printing or reporting. The function selects higher
#' taxonomic ranks, orders them according to phylogenetic and
#' classification priorities, and formats the output with markdown-style
#' indentation for each rank.
#'
#' @param ppg A cleaned data frame of PPG taxonomic data, typically the
#'   output of clean_ppg().
#' @param families_in_phy_order A character vector of family names in
#'   phylogenetic order (as returned by make_family_tree()).
#' @param higher_tax_levels Character vector; taxonomic levels to include
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{taxonID}{Unique taxon identifier.}
#'     \item{scientificName}{Scientific name of the taxon.}
#'     \item{scientificNameAuthorship}{Authorship of the scientific name.}
#'     \item{taxonRank}{Taxonomic rank (e.g., family, order).}
#'     \item{indent}{Markdown-style header string for indentation.}
#'   }
dwc_to_tl <- function(
  ppg,
  families_in_phy_order,
  higher_tax_levels = c(
    "class",
    "subclass",
    "order",
    "suborder",
    "family",
    "subfamily",
    "tribe",
    "subtribe",
    "genus"
  ),
  return_taxlist = TRUE
) {
  require(taxlist)

  # Format PPG data for printing out
  ppg_print <-
    ppg |>
    # Only keeping higher, accepted taxa
    filter(taxonRank %in% higher_tax_levels) |>
    filter(taxonomicStatus == "accepted") |>
    # TODO fix these in Rhakhis
    # Remove bad taxa
    filter(
      taxonID != "wfo-1000070090", # Todea Bernh., PPG I has Todea Willd. ex Bernh.
      taxonID != "wfo-0001114160", # Duplicate of Selaginella sanguinolenta (L.) Spring in different publication
      taxonID != "wfo-0001110737", # Todea barbara, unplaced and lacks parent
      taxonID != "wfo-0001118486" # Todea papuana, unplaced and lacks parent
    )

  # Identify higher taxonomic levels actually used
  higher_tax_levels_used <- higher_tax_levels[
    higher_tax_levels %in% ppg_print$taxonRank
  ]

  # Convert to taxonlist format
  res <-
    ppg_print |>
    dplyr::select(
      TaxonConceptID = taxonID,
      TaxonUsageID = taxonID,
      TaxonName = scientificName,
      AuthorName = scientificNameAuthorship,
      Level = taxonRank,
      Parent = parentNameUsageID
    ) |>
    mutate(AcceptedName = TRUE) |>
    as.data.frame() |>
    taxlist::df2taxlist(levels = rev(higher_tax_levels_used))

  if (return_taxlist) {
    return(res)
  }

  ppg_print
}

# vectorized version of count_children_single that accepts a character vector
# of taxa to count
count_children <- function(tax_list, taxon, level) {
  purrr::map_dbl(taxon, ~ count_children_single(tax_list, .x, level))
}

count_children_ppg <- function(ppg, families_in_phy_order) {
  # Filter to accepted names and convert to taxlist
  ppg_for_counting_tl <- dwc_to_tl(
    ppg,
    families_in_phy_order,
    c(
      "class",
      "subclass",
      "order",
      "suborder",
      "family",
      "subfamily",
      "tribe",
      "genus",
      "species"
    )
  )

  # Also need a dataframe at genus and higher
  ppg_for_counting_df <- taxlist::taxlist2df(ppg_for_counting_tl) |>
    as_tibble() |>
    janitor::clean_names() |>
    mutate(level = as.character(level)) |>
    filter(level != "species")

  # Loop through dataframe and count number of children taxa at different
  # levels. Exclude if the level to count equals or exceeds the target taxon.
  ppg_for_counting_df |>
    mutate(
      n_species = count_children(ppg_for_counting_tl, taxon_name, "species"),
      n_genera = case_when(
        level == "genus" ~ NaN,
        TRUE ~ count_children(ppg_for_counting_tl, taxon_name, "genus")
      ),
      n_tribe = case_when(
        level == "genus" ~ NaN,
        level == "tribe" ~ NaN,
        TRUE ~ count_children(ppg_for_counting_tl, taxon_name, "tribe")
      ),
      n_subfamily = case_when(
        level == "genus" ~ NaN,
        level == "tribe" ~ NaN,
        level == "subfamily" ~ NaN,
        TRUE ~ count_children(ppg_for_counting_tl, taxon_name, "subfamily")
      ),
      n_family = case_when(
        level == "genus" ~ NaN,
        level == "tribe" ~ NaN,
        level == "subfamily" ~ NaN,
        level == "family" ~ NaN,
        TRUE ~ count_children(ppg_for_counting_tl, taxon_name, "family")
      ),
      n_suborder = case_when(
        level == "genus" ~ NaN,
        level == "tribe" ~ NaN,
        level == "subfamily" ~ NaN,
        level == "family" ~ NaN,
        level == "suborder" ~ NaN,
        TRUE ~ count_children(ppg_for_counting_tl, taxon_name, "suborder")
      ),
      n_order = case_when(
        level == "genus" ~ NaN,
        level == "tribe" ~ NaN,
        level == "subfamily" ~ NaN,
        level == "family" ~ NaN,
        level == "suborder" ~ NaN,
        level == "order" ~ NaN,
        TRUE ~ count_children(ppg_for_counting_tl, taxon_name, "order")
      ),
      n_subclass = case_when(
        level == "genus" ~ NaN,
        level == "tribe" ~ NaN,
        level == "subfamily" ~ NaN,
        level == "family" ~ NaN,
        level == "suborder" ~ NaN,
        level == "order" ~ NaN,
        level == "subclass" ~ NaN,
        TRUE ~ count_children(ppg_for_counting_tl, taxon_name, "subclass")
      )
    ) |>
    select(taxonID = taxon_concept_id, contains("n_"))
}


# function for counting taxa within a group
count_children_single <- function(tax_list, taxon, level) {
  require(taxlist)

  taxon_select <- taxon
  level_select <- level

  # Do subset without subset()
  tax_filter <- tax_list
  tax_filter@taxonNames <- tax_list@taxonNames[
    tax_list@taxonNames$TaxonName == taxon_select,
  ]
  tax_filter <- taxlist::clean(tax_filter)
  tax_filter <- get_children(tax_list, tax_filter)

  # Now count children
  return(taxlist::count_taxa(tax_filter, level = level_select))
}


#' Get Ladderized Tip Labels from a Phylogenetic Tree
#'
#' Returns the tip labels of a phylogenetic tree in ladderized order.
#'
#' @param tree A phylogenetic tree object of class phylo.
#'
#' @return A character vector of tip labels in ladderized order.
get_ladderized_tips <- function(tree) {
  is_tip <- tree$edge[, 2] <= length(tree$tip.label)
  ordered_tips <- tree$edge[is_tip, 2]
  tree$tip.label[ordered_tips]
}

format_ppg_comments <- function(comments_raw) {
  comments_raw |>
    filter(!is.na(comment)) |>
    select(taxonID, comment)
}

# Helper function for format_ppg_taxa_count()
remove_na_elements <- function(
  lst,
  drop_empty = FALSE,
  coerce_character = FALSE
) {
  out <- lapply(lst, function(x) {
    if (coerce_character) {
      x <- as.character(x)
    }
    x[!is.na(x)]
  })
  if (drop_empty) {
    out <- out[lengths(out) > 0]
  }
  out
}

# Helper function for format_ppg_taxa_count()
make_rank_plural <- function(x) {
  patterns <- c(
    genus = "genera",
    tribe = "tribes",
    subfamily = "subfamilies",
    family = "families",
    order = "orders",
    suborder = "suborders",
    subclass = "subclasses"
  )

  patterns_whole <- setNames(
    object = unname(patterns),
    nm = paste0("\\b", names(patterns), "\\b")
  )

  str_replace_all(
    x,
    patterns_whole
  )
}

format_ppg_taxa_count <- function(children_tally) {
  children_tally |>
    select(taxonID, starts_with("n_")) |>
    pivot_longer(names_to = "level", values_to = "n", starts_with("n_")) |>
    filter(!is.na(n)) |>
    filter(n != 0) |>
    mutate(
      level = str_remove(level, "n_") |>
        str_replace_all("genera", "genus"),
      level_eng = case_when(
        n > 1 ~ make_rank_plural(level),
        .default = level
      ),
      n_eng = case_when(
        n < 10 ~ as.character(english::english(n)),
        .default = scales::number(n, big.mark = ",")
      ),
      text = glue("{n_eng} {level_eng}"),
      text = str_replace_all(text, "one species", "monotypic")
    ) |>
    select(taxonID, level, text) |>
    pivot_wider(
      names_from = level,
      values_from = text
    ) |>
    rowwise() |>
    mutate(
      text_list = list(do.call(
        c,
        list(
          subclass,
          order,
          suborder,
          family,
          subfamily,
          tribe,
          genus,
          species
        )
      ))
    ) |>
    ungroup() |>
    mutate(text_list = remove_na_elements(text_list)) |>
    mutate(
      text_c = map_chr(
        text_list,
        ~ knitr::combine_words(words = .x)
      ) |>
        str_to_sentence() |>
        str_replace_all(
          "One genus and monotypic",
          "Single monotypic genus"
        )
    ) |>
    select(taxonID, taxon_count = text_c)
}

count_ppgi <- function(ppgi) {
  ppgi |>
    filter(notes == "original PPGI 2016") |>
    select(class:genus) |>
    mutate(across(everything(), ~ n_distinct(., na.rm = TRUE))) |>
    unique() |>
    pivot_longer(
      names_to = "rank",
      values_to = "n",
      everything()
    )
}


#' Set Taxon Priority Order for Sorting
#'
#' Determines the priority order of higher taxonomic ranks (class, subclass,
#' order, suborder, family) for sorting and reporting, and checks that all
#' expected taxa are present in the provided data. The function returns a
#' character vector of taxon names in the desired order for sorting.
#'
#' @param ppg A cleaned data frame of PPG taxonomic data, typically the output
#'   of clean_ppg().
#' @param families_in_phy_order A character vector of family names in
#'   phylogenetic order, as returned by make_family_tree().
#'
#' @return A character vector of taxon names in the desired priority order for
#'   sorting.
set_taxon_priority <- function(ppg, families_in_phy_order) {
  # Set priorities for sorting by rank ----
  # Also check that all names are in data

  # - class
  priority_class <- c(
    "Lycopodiopsida",
    "Polypodiopsida"
  )

  class_check <- ppg |>
    filter(taxonRank == "class") |>
    filter(taxonomicStatus == "accepted") |>
    assert(
      in_set(priority_class),
      scientificName,
      success_fun = success_logical
    )

  # - subclass
  priority_subclass <- c(
    "Equisetidae",
    "Ophioglossidae",
    "Marattiidae",
    "Polypodiidae"
  )

  subclass_check <- ppg |>
    filter(taxonRank == "subclass") |>
    filter(taxonomicStatus == "accepted") |>
    assert(
      in_set(priority_subclass),
      scientificName,
      success_fun = success_logical
    )

  # - order
  priority_order <- c(
    # (Lycopodiopsida)
    "Lycopodiales",
    "Isoetales",
    "Selaginellales",
    # (Equisetidae)
    "Equisetales",
    # (Ophioglossidae)
    "Psilotales",
    "Ophioglossales",
    # (Marattiidae)
    "Marattiales",
    # (Polypodiidae)
    "Osmundales",
    "Hymenophyllales",
    "Gleicheniales",
    "Schizaeales",
    "Salviniales",
    "Cyatheales",
    "Polypodiales"
  )

  order_check <- ppg |>
    filter(taxonRank == "order") |>
    filter(taxonomicStatus == "accepted") |>
    assert(
      in_set(priority_order),
      scientificName,
      success_fun = success_logical
    )

  # suborder
  priority_suborder <- c(
    "Saccolomatineae",
    "Lindsaeineae",
    "Pteridineae",
    "Dennstaedtiineae",
    "Aspleniineae",
    "Polypodiineae"
  )

  suborder_check <- ppg |>
    filter(taxonRank == "suborder") |>
    filter(taxonomicStatus == "accepted") |>
    assert(
      in_set(priority_suborder),
      scientificName,
      success_fun = success_logical
    )

  # - family
  priority_family <- c(
    # Lycophytes
    "Lycopodiaceae",
    "Isoetaceae",
    "Selaginellaceae",
    # Ferns are determined from FTOL
    rev(families_in_phy_order)
  )

  family_check <- ppg |>
    filter(taxonRank == "family") |>
    filter(taxonomicStatus == "accepted") |>
    assert(
      in_set(priority_family),
      scientificName,
      success_fun = success_logical
    )

  # Compile all priorities
  c(
    priority_class,
    priority_subclass,
    priority_order,
    priority_suborder,
    priority_family
  )
}


#' Count the number of parent taxa in a DwC-style taxonomic dataset
#'
#' Given a Darwin Core (DwC) style table of taxa with `taxonID` and
#' `parentNameUsageID` columns, this function computes the number of
#' parent taxa (ancestors) for each record by traversing the parent-child
#' relationships as a directed graph.
#'
#' The function uses \pkg{igraph} to build a single graph of all taxa and
#' computes distances from a synthetic root node, making it much faster
#' than recursive row-by-row methods even for large datasets.
#'
#' @param ppg A data frame or tibble containing at least two columns:
#'   - `taxonID`: unique identifier for each taxon (character)
#'   - `parentNameUsageID`: the taxonID of the immediate parent (may be `NA` or
#'     empty for root taxa)
#'
#' @return A tibble identical to `ppg` but with an additional integer column:
#'   - `n_parents`: the number of parent taxa (ancestors) for each record.
#'
#' @details
#' If cycles exist in the parent relationships (e.g., data errors),
#' those taxa will receive `NA` for `n_parents` and a warning will be issued.
#'
#' Empty strings in `parentNameUsageID` are treated as missing (`NA`).
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(igraph)
#'
#' ppg_with_counts <- count_taxon_parents(ppg)
#' }
count_taxon_parents <- function(ppg) {
  # Requires: dplyr, tibble, igraph

  # 1) Build a parent→child edge list and add a synthetic ROOT
  edges <- ppg %>%
    dplyr::transmute(
      parent = dplyr::coalesce(dplyr::na_if(parentNameUsageID, ""), "ROOT"),
      child = taxonID
    )

  # 2) Create vertex table and construct directed graph
  verts <- tibble::tibble(name = unique(c(edges$parent, edges$child)))
  g <- igraph::graph_from_data_frame(edges, directed = TRUE, vertices = verts)

  # 3) Compute distances (levels) from ROOT to all nodes
  if (!igraph::is_dag(g)) {
    warning(
      "Cycle(s) detected in parent links; affected taxa will get NA for n_parents."
    )
  }

  d <- igraph::distances(g, v = "ROOT", to = igraph::V(g), mode = "out")[1, ]
  depth <- as.numeric(d) - 1
  depth[is.infinite(depth)] <- NA_real_

  # 4) Attach results back to original table
  ppg %>%
    dplyr::mutate(n_parents = as.integer(depth[match(taxonID, names(d))]))
}

rep_collapse_single <- function(x, n) {
  rep(x, n) |> paste(collapse = "")
}

rep_collapse <- function(x, n) {
  map_chr(n, ~ rep_collapse_single(x, .))
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
      status = replace_na(status, "TBA")
    ) |>
    filter(status != "NOT VALID")
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

count_ppgi <- function(ppg_i) {
  ppg_i |>
    filter(notes == "original PPGI 2016") |>
    select(class:genus) |>
    mutate(across(everything(), ~ n_distinct(., na.rm = TRUE))) |>
    unique() |>
    pivot_longer(
      names_to = "taxonRank",
      values_to = "accepted",
      everything()
    ) |>
    # manual changes
    rows_upsert(
      tribble(
        ~taxonRank , ~accepted ,
        "subclass" ,         4 ,
        "genus"    ,       337 ,
        "species"  ,     11916
      )
    )
}

count_ppg2_taxa <- function(ppg) {
  ppg |>
    mutate(
      nomenclaturalStatus = tidyr::replace_na(
        nomenclaturalStatus,
        "unknown"
      )
    ) |>
    filter(nomenclaturalStatus %in% c("conserved", "valid", "unknown")) |>
    filter(taxonomicStatus %in% c("synonym", "accepted")) |>
    filter(
      taxonRank %in%
        c(
          "class",
          "subclass",
          "order",
          "suborder",
          "family",
          "subfamily",
          "tribe",
          "genus",
          "species"
        )
    ) |>
    mutate(
      taxonRank = case_when(
        taxonRank == "genus" & str_detect(scientificName, "^× ") ~ "nothogenus",
        .default = taxonRank
      ),
      taxonRank = factor(
        taxonRank,
        levels = c(
          "class",
          "subclass",
          "order",
          "suborder",
          "family",
          "subfamily",
          "tribe",
          "genus",
          "nothogenus",
          "species"
        )
      )
    ) |>
    group_by(taxonRank, taxonomicStatus) |>
    count() |>
    ungroup() |>
    pivot_wider(
      names_from = taxonomicStatus,
      values_from = n,
      values_fill = 0
    )
}


# Helper function for fetch_parent
fetch_parent_single <- function(query, rank, ppg, ppg_tl) {
  query_id <- ppg |>
    filter(scientificName == query) |>
    filter(taxonomicStatus == "accepted") |>
    pull(taxonID)

  assertthat::assert_that(
    isTRUE(length(query_id) == 1),
    msg = "Could not find single exact match for query"
  )

  parent_id <- taxlist::parents(
    taxlist = ppg_tl,
    level = rank,
    concept = query_id
  )

  taxlist::print_name(
    ppg_tl,
    parent_id,
    italics = FALSE,
    include_author = FALSE
  )
}

# Retrieve the parent taxon at a specified taxonomic level
fetch_parent <- function(query, rank, ppg, ppg_tl) {
  purrr::map_chr(
    query,
    ~ fetch_parent_single(
      query = .,
      rank = rank,
      ppg = ppg,
      ppg_tl = ppg_tl
    )
  )
}

# Function to create tip labels for tree
format_tip_labels <- function(rank_select, phy_family, ppg, ppg_tl) {
  tips_to_higher_taxon <-
    ppg |>
    filter(taxonRank == "family", taxonomicStatus == "accepted") |>
    select(tip = scientificName) |>
    mutate(label = fetch_parent(tip, rank_select, ppg, ppg_tl)) |>
    add_count(label) |>
    filter(!is.na(label)) %>%
    assertr::verify(nrow(.) > 0) |>
    mutate(
      label_type = if_else(n > 1, "clade", "tip")
    )

  clade_labels <-
    tips_to_higher_taxon |>
    filter(label_type == "clade") |>
    select(-n)

  tip_labels <-
    tips_to_higher_taxon |>
    filter(label_type == "tip") |>
    select(-n)

  if (nrow(clade_labels) > 0) {
    clade_labels <-
      clade_labels |>
      group_by(label) |>
      summarize(node = getMRCA(phy_family, tip)) |>
      mutate(label_type = "clade")
  }

  bind_rows(
    tip_labels,
    clade_labels
  ) |>
    select(all_of(c("label", "label_type", "tip", "node")))
}

# Convert family-level fern tree to a tree including other major tracheophyte
# lineages
fern_to_tracheo_phy <- function(phy_family) {
  # Start with fern phylogeny
  fern_tree <- phy_family
  fern_tree$edge.length <- NULL
  fern_tree$root.edge <- 1

  # Make a lycophyte tree
  lyco_tree <- ape::read.tree(
    text = "(((Selaginellaceae, Isoetaceae), Lycopodiaceae), Algae);"
  )

  # Bind seed plants to ferns
  megaphyll_tree <- phytools::bind.tip(
    fern_tree,
    "Spermatophytes",
    where = "root"
  )

  # Re-root the tree so the root is between spermatophytes and ferns
  megaphyll_tree <- ape::root(
    megaphyll_tree,
    outgroup = "Spermatophytes",
    resolve.root = TRUE
  )

  megaphyll_tree$root.edge <- 1

  # Bind lycophytes to megaphylls
  tracheo_tree <-
    ape::bind.tree(
      megaphyll_tree,
      lyco_tree,
      where = "root"
    )

  # Fix rooting
  tracheo_tree <- ape::root(
    tracheo_tree,
    outgroup = "Algae",
    resolve.root = TRUE
  ) |>
    ape::drop.tip("Algae") |>
    root(
      outgroup = c("Lycopodiaceae", "Isoetaceae", "Selaginellaceae"),
      resolve.root = TRUE
    )

  tracheo_tree
}

rescale_tree <- function(tree, scale) {
  tree$edge.length <-
    tree$edge.length / max(phytools::nodeHeights(tree)[, 2]) * scale
  return(tree)
}

modify_node_height <- function(tree, tax_set, mod_length) {
  # Identify MRCA of taxon set
  mrca_node <- ape::getMRCA(
    tree,
    tax_set
  )

  # Modify length of branches leading to node
  edge_in <- which(tree$edge[, 2] == mrca_node)

  tree$edge.length[edge_in] <- tree$edge.length[edge_in] - mod_length

  # Get edges going out from the mrca node
  out_edges <- which(tree$edge[, 1] == mrca_node)

  # Adjust them to compensate
  tree$edge.length[out_edges] <- tree$edge.length[out_edges] + mod_length

  assertthat::assert_that(
    !isTRUE(any(tree$edge.length < 0)),
    msg = "Modification resulted in negative branchlengths"
  )

  tree
}

add_root_length <- function(tree, length) {
  tree$root.edge <- length
  tree
}

#' Create Combined Plot of PPG Issues Over Time
#'
#' Generates a combined visualization showing the progress of PPG
#' (Pteridophyte Phylogeny Group) proposals over time. The function
#' creates two plots: a cumulative count of submitted proposals and a
#' stacked bar chart showing voting results (passed, not passed, or TBD).
#' The plots are combined side-by-side with subfigure labels using the
#' patchwork package.
#'
#' @param ppg_issues A data frame containing GitHub issue data for PPG
#'   proposals. Expected to have columns: `title` (containing status
#'   designations like "[PASSED]" or "[NOT PASSED]") and `created_at`
#'   (timestamp of issue creation).
#'
#' @return A patchwork object combining two ggplot2 plots:
#'   \describe{
#'     \item{Plot A}{Line plot showing cumulative number of proposals
#'       submitted over time}
#'     \item{Plot B}{Stacked bar chart showing the count of proposals
#'       by voting result (passed, not passed, TBD) for each time
#'       period}
#'   }
make_issues_plot <- function(ppg_issues) {
  require(patchwork)

  # set color palette
  okabe_ito_cols <- c(
    orange = "#E69F00",
    skyblue = "#56B4E9",
    bluishgreen = "#009E73",
    yellow = "#F0E442",
    blue = "#0072B2",
    vermillion = "#D55E00",
    reddishpurple = "#CC79A7"
  )

  # Format issues (proposals)
  issues <-
    ppg_issues |>
    mutate(passed = str_detect(title, "\\[PASSED\\]")) |>
    mutate(not_passed = str_detect(title, "\\[NOT PASSED\\]")) |>
    mutate(voting = !str_detect(title, "PASSED")) |>
    mutate(
      created_at = str_remove_all(created_at, "T.*Z") |>
        lubridate::ymd()
    ) |>
    mutate(
      month = lubridate::month(created_at) |>
        str_pad(side = "left", pad = "0", width = 2)
    ) |>
    mutate(year = lubridate::year(created_at)) |>
    mutate(year_month = paste(year, month, sep = "-"))

  issues_summary <-
    issues %>%
    group_by(year_month) %>%
    summarize(
      passed = sum(passed),
      not_passed = sum(not_passed),
      under_vote = sum(voting)
    )

  # Make plot of cumulative count of submitted proposals
  cum_issue_plot <-
    issues %>%
    group_by(year_month) %>%
    mutate(year_month = lubridate::ym(year_month)) %>%
    summarize(
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(n_total = cumsum(n)) %>%
    mutate(dummy = "a") %>%
    ggplot(aes(x = year_month, y = n_total, group = dummy)) +
    geom_line() +
    geom_point(size = 1) +
    labs(
      y = "Number of proposals"
    ) +
    scale_x_date(
      date_labels = "%Y-%m",
      date_breaks = "4 months"
    ) +
    theme_gray(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank()
    )

  vote_res_plot <-
    issues_summary %>%
    pivot_longer(names_to = "type", values_to = "count", -year_month) %>%
    mutate(year_month = lubridate::ym(year_month)) %>%
    ggplot(aes(x = year_month, y = count, fill = type)) +
    geom_col(position = "stack") +
    labs(
      fill = "Result"
    ) +
    scale_fill_manual(
      values = c(
        "passed" = okabe_ito_cols[["bluishgreen"]],
        "not_passed" = okabe_ito_cols[["vermillion"]],
        "under_vote" = okabe_ito_cols[["skyblue"]]
      ),
      breaks = c("passed", "not_passed", "under_vote"),
      labels = c("Passed", "Not passed", "TBD")
    ) +
    scale_x_date(
      date_labels = "%Y-%m",
      date_breaks = "4 months"
    ) +
    scale_y_continuous(
      breaks = seq(0, 10, by = 2),
      expand = expansion(mult = c(0, 0.05))
    ) +
    theme_gray(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = c(1, 1),
      legend.justification = c(1, 1)
    )

  combined_plot <- cum_issue_plot +
    vote_res_plot +
    plot_annotation(tag_levels = 'A')

  combined_plot
}

#' Create Phylogenetic Tree Figure for PPG II
#'
#' Generates a comprehensive cladogram showing relationships between
#' families of tracheophytes recognized by PPG II. The tree includes
#' various label types (tips, clades, nodes) and visual elements to
#' represent taxonomic diversity and uncertain relationships.
#'
#' @param phy_family A phylo object containing the family-level fern
#'   phylogeny from the Fern Tree of Life.
#' @param ppg A cleaned data frame of PPG taxonomic data, typically the
#'   output of clean_ppg().
#' @param ppg_tl A taxlist object containing PPG taxonomy.
#' @param children_tally A data frame with taxon counts, typically the
#'   output of count_children_ppg().
#'
#' @return A ggtree object showing the phylogenetic tree with families,
#'   clade labels, genus/species counts, and uncertainty indicators.
make_tree_figure <- function(phy_family, ppg, ppg_tl, children_tally) {
  require(ggtree)
  require(ggimage)
  require(ggrepel)
  require(ape)
  require(phytools)

  # Add other tracheophytes to fern tree
  phy_tracheo <- fern_to_tracheo_phy(phy_family)

  # Format family labels with genus/species counts. Keep num species
  # for tip points
  family_tip_labels <-
    ppg |>
    filter(taxonRank == "family", taxonomicStatus == "accepted") |>
    select(taxonID) |>
    left_join(children_tally, by = join_by(taxonID)) |>
    select(taxon_name, n_genera, n_species) |>
    mutate(
      family_label = glue::glue("{taxon_name} ({n_genera}/{n_species})") |>
        as.character()
    ) |>
    select(tip = taxon_name, family_label, n_species) |>
    bind_rows(
      tibble(
        tip = "Spermatophytes",
        family_label = "SEED PLANTS",
        seed_plant_image = "images/starburst.png"
      )
    )

  # Format clade labels for formal taxa
  subclass_labels <- format_tip_labels(
    "subclass",
    phy_tracheo,
    ppg,
    ppg_tl
  )
  order_labels <- format_tip_labels("order", phy_tracheo, ppg, ppg_tl)
  suborder_labels <- format_tip_labels(
    "suborder",
    phy_tracheo,
    ppg,
    ppg_tl
  )

  # Format labels for other clade names
  other_labels <- tribble(
    ~node                                                       , ~label ,
    getMRCA(phy_tracheo, c("Blechnaceae", "Polypodiaceae"))     ,
    "eupolypods"                                                ,
    getMRCA(phy_tracheo, c("Lycopodiaceae", "Selaginellaceae")) ,
    "Lycopodiopsida"                                            ,
    getMRCA(phy_tracheo, c("Equisetaceae", "Polypodiaceae"))    ,
    "Polypodiopsida"                                            ,
    getMRCA(phy_tracheo, c("Spermatophytes", "Polypodiaceae"))  ,
    "euphyllophytes"                                            ,
    getMRCA(phy_tracheo, c("Lycopodiaceae", "Polypodiaceae"))   ,
    "tracheophytes"                                             ,
    # uncertain relationships
    getMRCA(phy_tracheo, c("Equisetaceae", "Ophioglossaceae"))  , "1"    ,
    getMRCA(phy_tracheo, c("Marattiaceae", "Osmundaceae"))      , "2"    ,
    getMRCA(phy_tracheo, c("Gleicheniaceae", "Dipteridaceae"))  , "3"
  ) |>
    mutate(label_type = "clade")

  # Combine label datasets
  all_labs <- subclass_labels |>
    bind_rows(order_labels) |>
    bind_rows(suborder_labels) |>
    bind_rows(other_labels) |>
    rename(taxon = label)

  # - right side clades
  clade_labs_side <- all_labs |>
    filter(label_type == "clade") |>
    filter(
      taxon %in%
        c(
          "Polypodiineae",
          "Aspleniineae",
          "Lindsaeineae",
          "Cyatheales",
          "Salviniales",
          "Schizaeales"
        )
    ) |>
    # ggtree wants different column names for each 'label' in every
    # dataset
    rename(taxon_clade = taxon)

  # - internal nodes
  clade_labs_node <- all_labs |>
    filter(label_type == "clade") |>
    filter(
      taxon %in%
        c(
          "eupolypods",
          "Polypodiales",
          "Polypodiidae",
          "Polypodiopsida",
          "euphyllophytes",
          "tracheophytes"
        )
    ) |>
    mutate(
      taxon = taxon |>
        str_replace_all(
          "Polypodiidae",
          "Polypodiidae\n(leptosporangiates)"
        ) |>
        str_replace_all(
          "Polypodiopsida",
          "Polypodiopsida (ferns)"
        ) |>
        str_replace_all(
          "Polypodiales",
          "Polypodiales\n(cathetogyrates)"
        )
    ) |>
    rename(taxon_node = taxon)

  # - right side tips labels (order or suborder)
  tip_labs_side <- all_labs |>
    filter(label_type == "tip") |>
    select(tip, taxon_tip = taxon) |>
    left_join(
      ppg |>
        filter(taxonomicStatus == "accepted") |>
        select(
          taxon_tip = scientificName,
          rank = taxonRank
        ),
      by = "taxon_tip"
    ) |>
    filter(rank %in% c("order", "suborder"))

  # - node labels above branches
  clade_labs_branch <- all_labs |>
    filter(label_type == "clade") |>
    filter(
      taxon %in%
        c(
          "Lycopodiopsida",
          "Ophioglossidae"
        )
    ) |>
    mutate(
      taxon = str_replace_all(
        taxon,
        "Lycopodiopsida",
        "Lycopodiopsida (lycophytes)"
      )
    ) |>
    rename(taxon_branch = taxon)

  # - node labels to add to uncertain nodes
  node_labs_nums <- all_labs |>
    filter(label_type == "clade") |>
    filter(taxon %in% c("1", "2", "3")) |>
    rename(node_num = taxon)

  # Since we can only adjust x-nudge for an entire layer at a time,
  # need to define separate datasets for each label with a different
  # x-nudge
  branch_labs_equisetaceae <- tribble(
    ~tip           , ~equisetaceae_lab          ,
    "Equisetaceae" , "Equisetidae (horsetails)"
  )

  branch_labs_marattiaceae <- tribble(
    ~tip           , ~marattiaceae_lab ,
    "Marattiaceae" , "Marattiidae"
  )

  # Uncertain nodes for linetype
  line_types_nodes <- as_tibble(ggtree::fortify(phy_tracheo))

  line_types_nodes_add <- all_labs |>
    filter(label_type == "clade") |>
    filter(taxon %in% c("1", "2", "3")) |>
    select(node) |>
    mutate(uncertain = TRUE)

  line_types_nodes <- left_join(
    line_types_nodes,
    line_types_nodes_add,
    by = "node",
    relationship = "one-to-one"
  ) |>
    mutate(
      uncertain = tidyr::replace_na(uncertain, FALSE)
    )

  # Set up branch lengths
  phy_tracheo_rescale <-
    phy_tracheo |>
    # Add branchlengths evenly
    compute.brlen(method = "grafen", power = 0.9) |>
    # Rescale total height (arbitrarily select 20)
    rescale_tree(20) |>
    # Tweak node for Equisetaceae + Ophioglossales
    modify_node_height(
      tax_set = c("Equisetaceae", "Psilotaceae"),
      mod_length = 16
    ) |>
    add_root_length(1)

  # set font sizes etc
  fig_font_size <- 2.5
  clade_lab_offset <- 9
  branch_lab_vjust <- -0.6

  # generate figure
  tree_fig <-
    # Base tree, with line type by uncertainty
    ggtree(phy_tracheo_rescale, aes(linetype = uncertain)) %<+%
    # Add datasets
    line_types_nodes %<+%
    family_tip_labels %<+%
    clade_labs_node %<+%
    clade_labs_branch %<+%
    branch_labs_equisetaceae %<+%
    branch_labs_marattiaceae %<+%
    tip_labs_side %<+%
    node_labs_nums +
    # Diamonds scaled by species count
    geom_tippoint(
      aes(size = n_species),
      shape = 18,
      color = "black"
    ) +
    geom_tiplab(
      aes(image = seed_plant_image),
      geom = "image",
      shape = 11,
      size = 0.019,
      nudge_x = -0.1
    ) +
    # Clade labels (nodes)
    geom_text_repel(
      aes(label = taxon_node),
      min.segment.length = 0,
      position = position_nudge_repel(x = -4, y = 1),
      size = fig_font_size
    ) +
    # Clade labels (above branch)
    geom_nodelab(
      aes(label = taxon_branch),
      size = fig_font_size,
      vjust = branch_lab_vjust,
      nudge_x = -8
    ) +
    # Node labels (at node, showing doubtful relationships)
    geom_nodelab(
      aes(label = node_num),
      geom = "label",
      size = fig_font_size * 0.7,
      nudge_x = -0.1,
      hjust = 0.5,
      vjust = 0.5,
      fill = "grey"
    ) +
    # Clade labels (right side)
    geom_cladelab(
      data = clade_labs_side,
      mapping = aes(node = node, label = taxon_clade),
      offset = clade_lab_offset,
      fontsize = fig_font_size,
      extend = 0.25
    ) +
    # Paraphyletic group labels (right side)
    geom_strip(
      "Gleicheniaceae",
      "Matoniaceae",
      label = "Gleicheniales",
      offset = clade_lab_offset,
      fontsize = fig_font_size,
      extend = 0.25,
      offset.text = 0.2
    ) +
    # Tip labels (family)
    geom_tiplab(
      aes(label = family_label),
      size = fig_font_size,
      offset = 0.5
    ) +
    # Tip labels (higher taxon on right side)
    geom_tiplab(
      aes(label = taxon_tip),
      size = fig_font_size,
      offset = clade_lab_offset
    ) +
    # Branch labels (monotypic groups, above branch)
    geom_tiplab(
      aes(label = equisetaceae_lab),
      size = fig_font_size,
      offset = -11.5,
      vjust = branch_lab_vjust
    ) +
    geom_tiplab(
      aes(label = marattiaceae_lab),
      size = fig_font_size,
      offset = -10,
      vjust = branch_lab_vjust
    ) +
    # seems geom_rootedge() requires brn lengths
    geom_rootedge() +
    xlim(-4, 33) +
    theme(
      legend.position = "none"
    )

  tree_fig
}

#' Format PPG Classification for Printing
#'
#' Creates a formatted, indented list of the PPG classification ready for
#' printing in the manuscript.
#'
#' @param ppg Data frame containing the PPG taxonomy in Darwin Core format
#' @param ppg_tl taxlist object with PPG taxonomy
#' @param children_tally Data frame with counts of children taxa
#' @param comments Data frame with comments for taxa
#' @param families_in_phy_order Character vector of family names in
#'   phylogenetic order
#'
#' @return Character vector of formatted classification lines
format_ppg_classification <- function(
  ppg,
  ppg_tl,
  children_tally,
  comments,
  families_in_phy_order
) {
  children_taxa_count <- format_ppg_taxa_count(children_tally) |>
    assertr::assert(assertr::not_na, taxon_count) |>
    dplyr::mutate(taxon_count = paste0(taxon_count, "."))

  priority_sort <- set_taxon_priority(ppg, families_in_phy_order)

  parent_taxon_count <-
    ppg |>
    dwc_to_tl(families_in_phy_order, return_taxlist = FALSE) |>
    count_taxon_parents() |>
    dplyr::select(taxonID, n_parents) |>
    assertr::assert(assertr::not_na, dplyr::everything()) |>
    assertr::assert(assertr::is_uniq, taxonID)

  ppg_tl |>
    taxlist::sort_taxa(priority = priority_sort) |>
    taxlist::indented_list(print = FALSE) |>
    tibble::as_tibble() |>
    janitor::clean_names() |>
    dplyr::select(
      taxonID = taxon_concept_id,
      scientificName = taxon_name,
      scientificNameAuthorship = author_name,
      taxonRank = level
    ) |>
    dplyr::left_join(
      children_taxa_count,
      by = "taxonID",
      relationship = "one-to-one"
    ) |>
    dplyr::left_join(
      parent_taxon_count,
      by = "taxonID",
      relationship = "one-to-one"
    ) |>
    dplyr::left_join(
      dplyr::select(
        ppg,
        taxonID,
        namePublishedIn
      ),
      by = "taxonID",
      relationship = "one-to-one"
    ) |>
    dplyr::left_join(
      comments,
      by = "taxonID",
      relationship = "one-to-one"
    ) |>
    dplyr::mutate(
      taxonRank_print = stringr::str_replace_all(taxonRank, "genus", "") |>
        stringr::str_to_sentence(),
      indent = rep_collapse("  ", n_parents),
      indent = paste0(indent, "* "),
      comment = tidyr::replace_na(comment, ""),
      taxon_count = tidyr::replace_na(taxon_count, "")
    ) |>
    dplyr::mutate(
      # genera in bold italics, everything else in bold
      name_print = dplyr::case_when(
        taxonRank == "genus" ~ glue::glue("***{scientificName}***"),
        .default = glue::glue("**{scientificName}**")
      ),
      pretty = glue::glue(
        "{taxonRank_print} {name_print}　{scientificNameAuthorship}. {taxon_count} {comment}"
      ) |>
        as.character() |>
        stringr::str_replace_all("\\.+", ".") |>
        stringr::str_squish(),
      pretty = glue::glue("{indent}{pretty}")
    ) |>
    dplyr::select(pretty) |>
    dplyr::pull(pretty)
}
