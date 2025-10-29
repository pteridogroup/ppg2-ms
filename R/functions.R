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