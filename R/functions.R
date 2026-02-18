#' Clean PPG Taxonomic Data
#'
#' Cleans raw PPG (Pteridophyte Phylogeny Group) taxonomic data by
#' removing duplicates, deleting invalid nothogenera, and filtering to
#' valid nomenclatural and taxonomic statuses.
#'
#' @param ppg_raw A data frame in Darwin Core format containing raw PPG
#'   taxonomic data. Expected to have columns taxonID,
#'   nomenclaturalStatus, taxonomicStatus, and other Darwin Core
#'   standard fields.
#'
#' @return A cleaned data frame with duplicate taxa removed, invalid
#'   nothogenera deleted, and filtered to only include taxa with
#'   nomenclatural status of "conserved", "valid", or "unknown" and
#'   taxonomic status of "accepted" or "synonym".
#'
#' @details
#' The function performs the following cleaning steps:
#' \enumerate{
#'   \item Removes a duplicate taxon (Selaginella sanguinolenta)
#'   \item Deletes three nothogenera that have not passed PPG voting:
#'     × Chrinephrium, × Chrismatopteris, and × Glaphyrocyclosorus
#'   \item Replaces NA values in nomenclaturalStatus with "unknown"
#'   \item Filters to retain only taxa with nomenclatural status of
#'     "conserved", "valid", or "unknown"
#'   \item Filters to retain only taxa with taxonomic status of
#'     "accepted" or "synonym"
#' }
clean_ppg <- function(ppg_raw) {
  require(dwctaxon)
  dct_options(stamp_modified = FALSE)
  ppg_raw |>
    # TODO fix these in Rhakhis
    # Remove bad taxa
    filter(
      # Duplicate of Selaginella sanguinolenta (L.) Spring
      # in different publication
      taxonID != "wfo-0001114160"
    ) |>
    # TODO fix these in Rhakhis
    # Delete nothogenera that are still in Rhakhis but have not passed
    # voting
    delete_taxon("× Chrinephrium", quiet = TRUE) |>
    delete_taxon("× Chrismatopteris", quiet = TRUE) |>
    delete_taxon("× Glaphyrocyclosorus", quiet = TRUE) |>
    mutate(
      nomenclaturalStatus = tidyr::replace_na(
        nomenclaturalStatus,
        "unknown"
      )
    ) |>
    filter(nomenclaturalStatus %in% c("conserved", "valid", "unknown")) |>
    filter(taxonomicStatus %in% c("accepted", "synonym"))
}

#' Delete a Taxon and All Its Descendants and Synonyms
#'
#' Removes a taxon from a Darwin Core taxonomic database along with all
#' its descendant taxa and all synonyms that point to it or its
#' descendants.
#'
#' @param dwc A data frame in Darwin Core format with columns taxonID,
#'   scientificName, parentNameUsageID, and acceptedNameUsageID.
#' @param name Character string; the scientific name of the taxon to
#'   delete.
#' @param taxonomic_status Character vector; taxonomic status values to
#'   consider when matching the name. Default is "accepted".
#' @param verbose Logical; if TRUE, prints the scientific names of all
#'   deleted taxa. Default is FALSE.
#' @param quiet Logical; if TRUE, suppresses all messages and warnings.
#'   Default is FALSE.
#'
#' @return A data frame with the specified taxon, its descendants, and
#'   related synonyms removed.
delete_taxon <- function(
  dwc,
  name,
  taxonomic_status = "accepted",
  verbose = FALSE,
  quiet = FALSE
) {
  require(dplyr)

  # Find the taxonID(s) for the given scientific name
  target_ids <- dwc |>
    filter(scientificName == name) |>
    filter(taxonomicStatus %in% taxonomic_status) |>
    pull(taxonID)

  if (length(target_ids) == 0) {
    if (!quiet) {
      warning(paste("No taxon found with name:", name))
    }
    return(dwc)
  }

  if (length(target_ids) > 1) {
    if (!quiet) {
      warning(paste(
        "Multiple taxa found with name:",
        name,
        "- deleting all matches"
      ))
    }
  }

  # Recursively find all descendant taxonIDs
  all_ids <- target_ids
  current_ids <- target_ids

  while (length(current_ids) > 0) {
    # Find children of current set
    children_ids <- dwc |>
      filter(parentNameUsageID %in% current_ids) |>
      pull(taxonID)

    # Add new children to our collection
    new_ids <- setdiff(children_ids, all_ids)
    all_ids <- c(all_ids, new_ids)
    current_ids <- new_ids
  }

  # Find all synonyms pointing to any of these taxa
  synonym_ids <- dwc |>
    filter(acceptedNameUsageID %in% all_ids) |>
    pull(taxonID)

  # Combine all IDs to delete
  ids_to_delete <- unique(c(all_ids, synonym_ids))

  # Print deleted names if verbose
  if (verbose && !quiet) {
    deleted_names <- dwc |>
      filter(taxonID %in% ids_to_delete) |>
      pull(scientificName)
    message("Deleted taxa:")
    message(paste(deleted_names, collapse = "\n"))
  }

  # Remove these taxa from the database
  result <- dwc |>
    filter(!taxonID %in% ids_to_delete)

  if (!quiet) {
    message(paste(
      "Deleted",
      length(ids_to_delete),
      "taxa:",
      length(all_ids),
      "accepted taxa and their descendants,",
      length(synonym_ids),
      "synonyms"
    ))
  }

  return(result)
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

  # Fix water ferns: they should not be sister to tree ferns
  # Extract water ferns clade
  wf <- ape::keep.tip(phy_family, c("Salviniaceae", "Marsileaceae"))

  wf$root.edge <- 10

  # Remove water ferns from main tree
  phy_family <- ape::drop.tip(phy_family, c("Salviniaceae", "Marsileaceae"))

  # Show node IDs on plot
  # plot(tree)
  # Show node IDs on plot
  # nodelabels()

  # Bind water ferns back at node 55
  phy_family <- ape::bind.tree(phy_family, wf, where = 55) |> multi2di()

  # Fix equisetales
  equi <- ape::keep.tip(phy_family, "Equisetaceae")

  equi$root.edge <- 10

  phy_family <- ape::drop.tip(phy_family, "Equisetaceae")

  phy_family <- ape::bind.tree(phy_family, equi, where = "root")

  phy_family <- ape::root(
    phy_family,
    outgroup = "Equisetaceae",
    resolve.root = TRUE
  )

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
    filter(taxonomicStatus == "accepted")

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
count_children <- function(tax_list, taxon, level, exclude_hybrids = FALSE) {
  purrr::map_dbl(
    taxon,
    ~ count_children_single(
      tax_list,
      .x,
      level,
      exclude_hybrids
    )
  )
}

count_children_ppg <- function(
  ppg,
  families_in_phy_order,
  exclude_hybrids = FALSE
) {
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
      is_nothogenus = if_else(
        level == "genus" & str_detect(taxon_name, "^×"),
        TRUE,
        FALSE
      ),
      n_species = case_when(
        is_nothogenus ~ count_children(
          ppg_for_counting_tl,
          taxon_name,
          "species",
          exclude_hybrids = FALSE
        ),
        .default = count_children(
          ppg_for_counting_tl,
          taxon_name,
          "species",
          exclude_hybrids = TRUE
        )
      ),
      n_genera = case_when(
        level == "genus" ~ NaN,
        TRUE ~ count_children(
          ppg_for_counting_tl,
          taxon_name,
          "genus",
          exclude_hybrids = TRUE
        )
      ),
      n_subfamily = case_when(
        level == "genus" ~ NaN,
        level == "subfamily" ~ NaN,
        TRUE ~ count_children(
          ppg_for_counting_tl,
          taxon_name,
          "subfamily",
          # subfamily and higher don't include hybrids (names with ×) anyways
          exclude_hybrids = FALSE
        )
      ),
      n_family = case_when(
        level == "genus" ~ NaN,
        level == "subfamily" ~ NaN,
        level == "family" ~ NaN,
        TRUE ~ count_children(
          ppg_for_counting_tl,
          taxon_name,
          "family",
          exclude_hybrids = FALSE
        )
      ),
      n_suborder = case_when(
        level == "genus" ~ NaN,
        level == "subfamily" ~ NaN,
        level == "family" ~ NaN,
        level == "suborder" ~ NaN,
        TRUE ~ count_children(
          ppg_for_counting_tl,
          taxon_name,
          "suborder",
          exclude_hybrids = FALSE
        )
      ),
      n_order = case_when(
        level == "genus" ~ NaN,
        level == "subfamily" ~ NaN,
        level == "family" ~ NaN,
        level == "suborder" ~ NaN,
        level == "order" ~ NaN,
        TRUE ~ count_children(
          ppg_for_counting_tl,
          taxon_name,
          "order",
          exclude_hybrids = FALSE
        )
      ),
      n_subclass = case_when(
        level == "genus" ~ NaN,
        level == "subfamily" ~ NaN,
        level == "family" ~ NaN,
        level == "suborder" ~ NaN,
        level == "order" ~ NaN,
        level == "subclass" ~ NaN,
        TRUE ~ count_children(
          ppg_for_counting_tl,
          taxon_name,
          "subclass",
          exclude_hybrids = FALSE
        )
      )
    ) |>
    select(taxonID = taxon_concept_id, contains("n_"))
}

#' Apply Manual Species Count Updates
#'
#' Updates species counts in children_tally with manually curated values from
#' a CSV file. This is useful for correcting automated counts that may not
#' reflect the most current taxonomic knowledge.
#'
#' @param children_tally Data frame with taxon counts, typically the output
#'   of count_children_ppg().
#' @param species_count_updates Data frame with columns 'genus' and 'species'
#'   containing manual species count updates.
#' @param ppg Data frame containing the PPG taxonomy in Darwin Core format.
#'
#' @return Updated children_tally data frame with corrected species counts.
apply_species_count_updates <- function(
  children_tally,
  species_count_updates,
  ppg
) {
  # Get taxonIDs for genera to update
  genera_to_update <- ppg |>
    filter(
      scientificName %in% species_count_updates$genus,
      taxonRank == "genus",
      taxonomicStatus == "accepted"
    ) |>
    select(scientificName, taxonID) |>
    left_join(
      species_count_updates,
      by = join_by(scientificName == genus)
    ) |>
    rename(n_species_manual = species)

  # Update children_tally at genus level
  children_tally_updated <- children_tally |>
    left_join(genera_to_update, by = "taxonID") |>
    mutate(
      # Calculate the difference for propagation
      species_diff = case_when(
        !is.na(n_species_manual) ~ n_species_manual - n_species,
        .default = 0
      ),
      n_species = coalesce(n_species_manual, n_species)
    ) |>
    select(-scientificName, -n_species_manual)

  # Get parent relationships to propagate changes up the hierarchy
  parent_map <- ppg |>
    filter(taxonomicStatus == "accepted") |>
    select(taxonID, parentNameUsageID)

  # Calculate total species difference by parent
  species_diff_by_parent <- children_tally_updated |>
    filter(species_diff != 0) |>
    left_join(parent_map, by = "taxonID") |>
    group_by(parentNameUsageID) |>
    summarize(total_species_diff = sum(species_diff), .groups = "drop") |>
    rename(taxonID = parentNameUsageID)

  # Recursively propagate changes up the hierarchy
  while (nrow(species_diff_by_parent) > 0) {
    children_tally_updated <- children_tally_updated |>
      left_join(species_diff_by_parent, by = "taxonID") |>
      mutate(
        n_species = case_when(
          !is.na(total_species_diff) ~ n_species + total_species_diff,
          .default = n_species
        )
      ) |>
      select(-total_species_diff)

    # Find next level up
    species_diff_by_parent <- species_diff_by_parent |>
      left_join(parent_map, by = "taxonID") |>
      filter(!is.na(parentNameUsageID)) |>
      group_by(parentNameUsageID) |>
      summarize(
        total_species_diff = sum(total_species_diff),
        .groups = "drop"
      ) |>
      rename(taxonID = parentNameUsageID)
  }

  # Remove temporary column
  children_tally_updated |>
    select(-species_diff)
}


# function for counting taxa within a group
count_children_single <- function(
  tax_list,
  taxon,
  level,
  exclude_hybrids = FALSE
) {
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

  # Optionally exclude hybrids from children
  if (exclude_hybrids) {
    tax_filter@taxonNames <- tax_filter@taxonNames[
      !str_detect(tax_filter@taxonNames$TaxonName, "×"),
    ]
    tax_filter <- taxlist::clean(tax_filter)
  }

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
    select(taxonID = taxon_id, comment)
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
      text = str_replace_all(text, "one species", "monospecific")
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
          "One genus and monospecific",
          "Single monospecific genus"
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

count_ppgi <- function(ppg_i) {
  ppg_i |>
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
        "species"  ,     11916
      )
    )
}

count_ppg2_taxa <- function(ppg, exclude_hybrids = FALSE) {
  initial_count <- ppg |>
    assert(in_set(c("conserved", "valid", "unknown")), nomenclaturalStatus) |>
    assert(in_set(c("synonym", "accepted")), taxonomicStatus) |>
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
      ),
      is_nothotaxon = str_detect(scientificName, "×")
    ) |>
    group_by(taxonRank, taxonomicStatus, is_nothotaxon) |>
    count()

  # Optionally exclude hybrid species from count
  if (exclude_hybrids) {
    initial_count <- initial_count |>
      filter(!(taxonRank == "species" & is_nothotaxon == TRUE)) |>
      # Should only be one row that is a nothotaxon (taxonRank == "nothgenus")
      verify(sum(is_nothotaxon) == 1) |>
      verify(all((taxonRank == "nothogenus") == is_nothotaxon))
  }

  # Format as wide table
  initial_count |>
    ungroup() |>
    select(-is_nothotaxon) |>
    group_by(taxonRank, taxonomicStatus) |>
    summarize(n = sum(n)) |>
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
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      panel.grid.minor = element_blank()
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
        "passed" = okabe_ito_cols[["blue"]],
        "not_passed" = okabe_ito_cols[["vermillion"]]
      ),
      breaks = c("passed", "not_passed"),
      labels = c("Passed", "Not passed")
    ) +
    scale_x_date(
      date_labels = "%Y-%m",
      date_breaks = "4 months"
    ) +
    scale_y_continuous(
      breaks = seq(0, 10, by = 2),
      expand = expansion(mult = c(0, 0.05))
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      # Move legend just inside top-right corner of plot area
      legend.position = c(0.99, 0.99),
      legend.justification = c("right", "top"),
      panel.grid.minor = element_blank()
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
    ~node                                                          , ~label ,
    getMRCA(phy_tracheo, c("Blechnaceae", "Polypodiaceae"))        ,
    "eupolypods"                                                   ,
    getMRCA(phy_tracheo, c("Lycopodiaceae", "Selaginellaceae"))    ,
    "Lycopodiopsida"                                               ,
    getMRCA(phy_tracheo, c("Equisetaceae", "Polypodiaceae"))       ,
    "Polypodiopsida"                                               ,
    getMRCA(phy_tracheo, c("Spermatophytes", "Polypodiaceae"))     ,
    "euphyllophytes"                                               ,
    getMRCA(phy_tracheo, c("Lycopodiaceae", "Polypodiaceae"))      ,
    "tracheophytes"                                                ,
    # uncertain relationships
    getMRCA(phy_tracheo, c("Equisetaceae", "Ophioglossaceae"))     , "1"    ,
    getMRCA(phy_tracheo, c("Marattiaceae", "Osmundaceae"))         , "2"    ,
    getMRCA(phy_tracheo, c("Hymenophyllaceae", "Gleicheniaceae"))  , "3"    ,
    getMRCA(phy_tracheo, c("Gleicheniaceae", "Dipteridaceae"))     , "4"    ,
    getMRCA(phy_tracheo, c("Hypodematiaceae", "Dennstaedtiaceae")) , "5"    ,
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
    filter(taxon %in% as.character(1:5)) |>
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
    filter(taxon %in% as.character(1:5)) |>
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

  # Get root node for manually drawing root edge
  # (needed because geom_rootedge() doesn't work with branch.length = "none")
  n_tips <- length(phy_tracheo$tip.label)
  root_node_num <- n_tips + 1

  # Get the actual y-coordinate of the root from the tree layout
  tree_layout <- as_tibble(ggtree::fortify(phy_tracheo))
  root_y <- tree_layout |>
    filter(node == root_node_num) |>
    pull(y)

  # set font sizes etc
  fig_font_size <- 2.5
  clade_lab_offset <- 9
  branch_lab_vjust <- -0.6
  root_edge_length <- 1

  # generate figure
  tree_fig <-
    # Base tree, with line type by uncertainty
    ggtree(phy_tracheo, branch.length = "none") %<+%
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
    # Manual root edge (horizontal line at root)
    geom_segment(
      aes(x = -root_edge_length, xend = 0, y = root_y, yend = root_y),
      linewidth = 0.5
    ) +
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

  classification_data <- ppg_tl |>
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
    # Add original order to preserve hierarchy
    dplyr::mutate(original_order = dplyr::row_number()) |>
    # Join parent info
    dplyr::left_join(
      dplyr::select(ppg, taxonID, parentNameUsageID),
      by = "taxonID",
      relationship = "one-to-one"
    ) |>
    dplyr::mutate(
      is_nothogenus = stringr::str_detect(scientificName, "^×")
    )

  # For genera only: create a sort key that moves nothogenera to end of family
  genera_data <- classification_data |>
    dplyr::filter(taxonRank == "genus") |>
    dplyr::group_by(parentNameUsageID) |>
    dplyr::mutate(
      # Within each family, nothogenera get placed after the last regular genus
      genus_position = dplyr::row_number(),
      max_regular_order = max(original_order[!is_nothogenus], na.rm = TRUE)
    ) |>
    dplyr::ungroup()

  # Compute min position for nothogenera separately to avoid warnings
  nothogenus_min_pos <- genera_data |>
    dplyr::filter(is_nothogenus) |>
    dplyr::group_by(parentNameUsageID) |>
    dplyr::summarize(
      min_nothogenus_pos = min(genus_position),
      .groups = "drop"
    )

  # Combine and compute sort keys
  genera_data <- genera_data |>
    dplyr::left_join(nothogenus_min_pos, by = "parentNameUsageID") |>
    dplyr::mutate(
      sort_key = dplyr::if_else(
        is_nothogenus,
        max_regular_order + (genus_position - min_nothogenus_pos + 1) * 0.1,
        as.numeric(original_order)
      )
    ) |>
    dplyr::select(taxonID, sort_key)

  # Apply the new sort key and re-order
  classification_data |>
    dplyr::left_join(genera_data, by = "taxonID") |>
    dplyr::mutate(
      final_order = dplyr::coalesce(sort_key, original_order)
    ) |>
    dplyr::arrange(final_order) |>
    dplyr::select(
      -parentNameUsageID,
      -is_nothogenus,
      -original_order,
      -sort_key,
      -final_order
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

#' Check That All PPG I Taxa Are Accounted For
#'
#' Verifies that all taxa from PPG I (the original classification) are
#' mentioned either in the manuscript or in the comments file. This
#' ensures that readers know what happened to taxa that were synonymized
#' or otherwise removed in PPG II.
#'
#' @param ppg_i Data frame containing PPG I taxonomy with a 'genus'
#'   column.
#' @param ppg_classification Character string; the PPG II classification.
#' @param comments_path Character string; path to the comments CSV file.
#'   Default is "data/ppg_comments.csv".
#'
#' @return Invisible NULL if all taxa are accounted for.
#'
#' @details
#' The function extracts all unique genus names from PPG I, then checks
#' if each name appears in either the manuscript text or the comments
#' file. If any names are missing, it throws an error listing the
#' unaccounted taxa.
check_ppg_i_taxa_accounted <- function(
  ppg_i,
  ppg_classification,
  comments_path = "data/ppg_comments.csv"
) {
  require(dplyr)
  require(readr)
  require(stringr)

  # Get all unique genus names from PPG I
  ppg_i_taxa <- ppg_i |>
    filter(!is.na(genus)) |>
    select(class:genus) |>
    pivot_longer(values_to = "taxon", names_to = "rank", everything()) |>
    filter(!is.na(taxon)) |>
    pull(taxon) |>
    unique()

  # Read comments
  comments_text <- if (file.exists(comments_path)) {
    readr::read_csv(comments_path, show_col_types = FALSE) |>
      pull(comment) |>
      paste(collapse = " ")
  } else {
    ""
  }

  # Combine text sources
  all_text <- paste(
    comments_text,
    ppg_classification,
    collapse = " "
  )

  # Check which genera are missing
  missing_genera <- ppg_i_taxa[
    !stringr::str_detect(
      all_text,
      stringr::fixed(ppg_i_taxa)
    )
  ]

  if (length(missing_genera) > 0) {
    stop(
      "The following PPG I taxa are not accounted for in the PPG II ",
      "classification or comments:\n",
      paste(sort(missing_genera), collapse = "\n"),
      "\n\nPlease add mentions of these taxa or add comments explaining ",
      "what happened to them."
    )
  }

  invisible(NULL)
}

check_issue_type_count <- function(ppg_issues, ppg_issues_count) {
  ppg_issues |>
    filter(str_detect(status, "^PASSED$")) |>
    anti_join(ppg_issues_count, by = "number") %>%
    verify(nrow(.) == 0)
}

#' Convert PPG from Long Format to Wide Hierarchical Format
#'
#' Transforms PPG II Darwin Core taxonomic data from long format (where
#' hierarchy is encoded via parentNameUsageID) to wide format with
#' separate columns for each taxonomic rank. This enables direct
#' comparison with PPG I data.
#'
#' @param ppg A data frame in Darwin Core format containing PPG
#'   taxonomic data with columns taxonID, scientificName, taxonRank,
#'   parentNameUsageID, and taxonomicStatus. Typically the output of
#'   clean_ppg().
#' @param ranks Character vector of taxonomic ranks to include as
#'   columns in wide format, in hierarchical order from highest to
#'   lowest. Default is c("class", "subclass", "order", "suborder",
#'   "family", "subfamily", "genus").
#' @param accepted_only Logical; if TRUE (default), only include
#'   accepted taxa. If FALSE, include both accepted and synonym taxa.
#'
#' @return A tibble in wide format with one row per taxon and columns
#'   for each specified rank, plus additional columns like notes. Taxa
#'   are filtered to only include those matching the lowest specified
#'   rank.
#'
#' @details The function builds the hierarchical structure by
#'   traversing the parent-child relationships encoded in
#'   parentNameUsageID. It creates a lookup table mapping each taxon
#'   to its ancestors at each taxonomic rank, then pivots this to wide
#'   format. Only taxa at the lowest specified rank (typically genus)
#'   are returned as rows.
#'
#' @examples
#' \dontrun{
#' # Convert PPG to wide format matching PPG I structure
#' ppg_wide <- ppg_to_wide(ppg)
#'
#' # Include species and subspecies ranks
#' ppg_wide_species <- ppg_to_wide(
#'   ppg,
#'   ranks = c("class", "order", "family", "genus", "species")
#' )
#' }
ppg_to_wide <- function(
  ppg,
  ranks = c(
    "class",
    "subclass",
    "order",
    "suborder",
    "family",
    "subfamily",
    "genus"
  ),
  accepted_only = TRUE
) {
  require(dplyr)
  require(tidyr)

  # Filter to accepted taxa if requested
  if (accepted_only) {
    ppg_data <- ppg |>
      filter(taxonomicStatus == "accepted")
  } else {
    ppg_data <- ppg
  }

  # Get the lowest rank (should be "genus")
  lowest_rank <- ranks[length(ranks)]

  # Get all unique genera (taxa at the lowest rank, typically genus)
  # This ensures we only get genus-level taxa, not species
  genera <- ppg_data |>
    filter(taxonRank == lowest_rank) |>
    group_by(scientificName) |>
    slice(1) |>
    ungroup() |>
    select(taxonID, genus = scientificName)

  # For each rank, build a lookup table of taxonID -> name at that rank
  # by walking up parent hierarchy
  rank_lookups <- list()

  for (rank_name in ranks) {
    # Get all taxa at this rank
    taxa_at_rank <- ppg_data |>
      filter(taxonRank == rank_name) |>
      select(taxonID, rank_value = scientificName)

    # Build a mapping: for each taxonID in ppg, what is its ancestor at
    # this rank?
    # Walk up to 10 levels (should be enough for any taxonomy)
    current_ids <- ppg_data |>
      select(original_id = taxonID, current_id = taxonID) |>
      mutate(rank_value = NA_character_)

    for (i in 1:10) {
      # Check if current_id is at the target rank
      current_ids <- current_ids |>
        left_join(
          taxa_at_rank,
          by = c("current_id" = "taxonID"),
          suffix = c("", "_new")
        )

      # If we found a rank value, use it
      current_ids <- current_ids |>
        mutate(
          rank_value = coalesce(rank_value, rank_value_new)
        ) |>
        select(-rank_value_new)

      # Get parent for those that don't have rank_value yet
      current_ids <- current_ids |>
        left_join(
          select(ppg_data, current_id = taxonID, parent_id = parentNameUsageID),
          by = "current_id"
        ) |>
        mutate(
          current_id = if_else(is.na(rank_value), parent_id, current_id)
        ) |>
        select(original_id, current_id, rank_value)

      # If all have rank_value or no parents left, stop
      if (all(is.na(current_ids$current_id) | !is.na(current_ids$rank_value))) {
        break
      }
    }

    # Store the final mapping
    rank_lookups[[rank_name]] <- current_ids |>
      select(taxonID = original_id, !!rank_name := rank_value) |>
      distinct()
  }

  # Join all rank lookups to genera
  result <- genera
  for (rank_name in ranks[ranks != "genus"]) {
    result <- result |>
      left_join(rank_lookups[[rank_name]], by = "taxonID")
  }

  # Clean up and select only rank columns (genus already exists)
  result |>
    select(all_of(ranks)) |>
    distinct()
}

#' Check if PPG I to PPG II Classification Changes Have Passed Issues
#'
#' Compares classification changes between PPG I and PPG II at the genus
#' level and verifies whether each change has a corresponding GitHub
#' issue that passed voting. This function helps ensure that all
#' taxonomic changes are properly documented and approved.
#'
#' @param ppg_ii A data frame containing PPG II data in wide format with
#'   columns for taxonomic ranks (class, subclass, order, suborder,
#'   family, subfamily, genus). Typically the output of ppg_to_wide().
#' @param ppg_i A data frame containing PPG I data with columns for
#'   taxonomic ranks (class, order, suborder, family, subfamily, genus).
#' @param ppg_issues A data frame containing GitHub issue data with
#'   columns: number, title (containing [PASSED] for passed issues),
#'   name (taxon name), and rank (taxonomic rank).
#' @param ranks Character vector of ranks to compare. Default is
#'   c("class", "subclass", "order", "suborder", "family", "subfamily").
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{genus}{Genus name}
#'     \item{classification_changed}{Logical; TRUE if any rank changed}
#'     \item{ranks_changed}{Character; which ranks changed}
#'     \item{has_passed_issue}{Logical; TRUE if genus has a passed
#'       issue}
#'     \item{issue_numbers}{Character; comma-separated issue numbers}
#'     \item{needs_attention}{Logical; TRUE if changed but no passed
#'       issue}
#'     \item{*_ii}{PPG II values for each rank}
#'     \item{*_i}{PPG I values for each rank}
#'   }
#'
#' @details The function performs a full join of PPG I and PPG II data
#'   by genus, then identifies which ranks have changed. It cross-
#'   references these changes with passed GitHub issues to flag any
#'   changes lacking proper documentation. New genera in PPG II are
#'   expected to have passed issues.
#'
#' @examples
#' \dontrun{
#' # Convert PPG II to wide format
#' ppg_wide <- ppg_to_wide(ppg)
#'
#' # Check changes
#' changes <- check_ppg_classification_changes(
#'   ppg_wide,
#'   ppg_i,
#'   ppg_issues
#' )
#'
#' # Find genera that need attention
#' changes |> filter(needs_attention)
#' }
check_ppg_classification_changes <- function(
  ppg_ii,
  ppg_i,
  ppg_issues,
  ppg = NULL,
  ranks = c("class", "subclass", "order", "suborder", "family", "subfamily")
) {
  require(dplyr)
  require(tidyr)
  require(stringr)

  # Prepare PPG I data: add _i suffix to rank columns
  ppg_i_comp <- ppg_i |>
    select(all_of(c(ranks, "genus"))) |>
    rename_with(
      ~ paste0(.x, "_i"),
      .cols = -genus
    )

  # Prepare PPG II data: add _ii suffix to rank columns
  # For genera that are in PPG I but not in ppg_ii (because they're now
  # synonyms), add them with all NA values
  genera_in_i_not_in_ii <- ppg_i_comp |>
    filter(!genus %in% ppg_ii$genus) |>
    pull(genus)

  # Create rows for these synonym genera with NA for all ranks
  synonym_genera <- tibble(genus = genera_in_i_not_in_ii)
  for (rank in ranks) {
    synonym_genera[[rank]] <- NA_character_
  }

  # Combine accepted genera from ppg_ii with synonym genera
  ppg_ii_comp <- bind_rows(ppg_ii, synonym_genera) |>
    select(all_of(c(ranks, "genus"))) |>
    rename_with(
      ~ paste0(.x, "_ii"),
      .cols = -genus
    )

  # Get passed issues for all ranks (not just genus)
  # We'll check if the taxon name matches the genus OR any of the
  # higher ranks that changed
  passed_issues_all <- ppg_issues |>
    filter(str_detect(status, "^PASSED$")) |>
    mutate(
      # Standardize rank names
      rank = str_to_lower(rank),
      # Clean up taxon names: remove 'and', commas, normalize spacing
      name = str_remove_all(name, regex("and|,", ignore_case = TRUE)) |>
        str_squish()
    ) |>
    # Separate multiple taxa into individual rows
    separate_rows(name, sep = "\\s+") |>
    mutate(
      # Remove × or x prefix from nothogenera
      name = str_remove(name, "^[×x]") |>
        str_trim()
    ) |>
    filter(name != "") |>
    select(number, rank, name) |>
    distinct()

  # Join and compare
  result <- full_join(ppg_ii_comp, ppg_i_comp, by = "genus")

  # Check each rank for changes - create a helper function
  check_rank_changed <- function(df, rank) {
    ii_col <- paste0(rank, "_ii")
    i_col <- paste0(rank, "_i")

    # Extract the columns
    ii_val <- df[[ii_col]]
    i_val <- df[[i_col]]

    # Consider it changed if:
    # - One is NA and the other isn't, OR
    # - Both are non-NA but values differ
    changed <- (is.na(ii_val) != is.na(i_val)) |
      (!is.na(ii_val) & !is.na(i_val) & ii_val != i_val)

    changed
  }

  # Create change indicator columns for each rank
  for (rank in ranks) {
    result[[paste0(rank, "_changed")]] <-
      check_rank_changed(result, rank)
  }

  # Identify which ranks changed
  rank_change_cols <- paste0(ranks, "_changed")
  result <- result |>
    mutate(
      classification_changed = rowSums(
        pick(all_of(rank_change_cols)),
        na.rm = TRUE
      ) >
        0,
      ranks_changed = pmap_chr(
        pick(all_of(rank_change_cols)),
        function(...) {
          vals <- list(...)
          changed <- ranks[unlist(vals)]
          if (length(changed) == 0) {
            return(NA_character_)
          }
          paste(changed, collapse = ", ")
        }
      )
    ) |>
    select(-all_of(rank_change_cols))

  # Helper function to find matching issues for a single genus
  find_matching_issues <- function(
    genus_val,
    rank_cols,
    issues_df,
    ppg_data = NULL
  ) {
    # Build list of taxa to check, removing × prefix to match issue
    # names (which have × removed during parsing)
    genus_clean <- str_remove(genus_val, "^×") |> str_trim()
    taxa_to_check <- c(genus_clean)

    # If ppg_data provided, check if this genus is a synonym and get its
    # accepted name
    if (!is.null(ppg_data)) {
      genus_record <- ppg_data |>
        filter(
          scientificName == genus_val,
          taxonRank == "genus"
        ) |>
        slice(1)

      if (
        nrow(genus_record) > 0 &&
          genus_record$taxonomicStatus == "synonym" &&
          !is.na(genus_record$acceptedNameUsageID)
      ) {
        # Look up the accepted genus name
        accepted_record <- ppg_data |>
          filter(taxonID == genus_record$acceptedNameUsageID) |>
          slice(1)

        if (nrow(accepted_record) > 0) {
          accepted_name <- str_remove(accepted_record$scientificName, "^×") |>
            str_trim()
          taxa_to_check <- c(taxa_to_check, accepted_name)
        }
      }
    }

    # Add all rank values (both _ii and _i)
    for (rank in ranks) {
      ii_val <- rank_cols[[paste0(rank, "_ii")]]
      i_val <- rank_cols[[paste0(rank, "_i")]]
      if (!is.na(ii_val)) {
        # Clean rank names too (remove × and extra spaces, remove
        # parentheses)
        clean_val <- str_remove(ii_val, "^×") |>
          str_trim() |>
          str_remove("\\s*\\(.*\\)$")
        taxa_to_check <- c(taxa_to_check, clean_val)
      }
      if (!is.na(i_val)) {
        clean_val <- str_remove(i_val, "^×") |>
          str_trim() |>
          str_remove("\\s*\\(.*\\)$")
        taxa_to_check <- c(taxa_to_check, clean_val)
      }
    }

    taxa_to_check <- unique(taxa_to_check)

    # Find matching issues by exact name match
    matches <- issues_df |>
      filter(name %in% taxa_to_check)

    # If no exact matches found and the genus has all NA ranks
    # (new genera), also try partial matching for the genus in case
    # it was mentioned in a broader proposal
    if (nrow(matches) == 0) {
      # Check if all rank columns are NA
      all_rank_values <- unlist(c(
        rank_cols[paste0(ranks, "_ii")],
        rank_cols[paste0(ranks, "_i")]
      ))
      all_ranks_na <- all(is.na(all_rank_values))

      if (all_ranks_na) {
        # Try matching the genus name anywhere in the issue names
        matches <- issues_df |>
          filter(
            str_detect(
              name,
              regex(paste0("^", genus_clean, "$"), ignore_case = TRUE)
            )
          )
      }
    }

    if (nrow(matches) > 0) {
      paste(unique(matches$number), collapse = ", ")
    } else {
      NA_character_
    }
  }

  # Apply the function to each row
  result <- result |>
    rowwise() |>
    mutate(
      issue_numbers = find_matching_issues(
        genus,
        pick(everything()),
        passed_issues_all,
        ppg
      ),
      has_passed_issue = !is.na(issue_numbers),
      needs_attention = classification_changed & !has_passed_issue
    ) |>
    ungroup()

  # Reorder columns for readability
  result |>
    select(
      genus,
      classification_changed,
      ranks_changed,
      has_passed_issue,
      issue_numbers,
      needs_attention,
      everything()
    ) |>
    arrange(desc(needs_attention), genus)
}

#' Summarize classification check results
#'
#' Create a summary of the PPG I vs PPG II classification check,
#' showing how many genera changed and how many have passed issues.
#'
#' @param classification_check Output from check_ppg_classification_changes()
#' @return List with summary statistics
summarize_classification_check <- function(classification_check) {
  require(dplyr)

  list(
    total_genera = nrow(classification_check),
    genera_changed = sum(classification_check$classification_changed),
    genera_with_issues = sum(classification_check$has_passed_issue),
    genera_needing_attention = sum(classification_check$needs_attention),
    match_rate = round(
      sum(
        classification_check$has_passed_issue &
          classification_check$classification_changed
      ) /
        sum(classification_check$classification_changed) *
        100,
      1
    ),
    sample_matched = classification_check |>
      filter(classification_changed, has_passed_issue) |>
      select(genus, ranks_changed, issue_numbers) |>
      head(10)
  )
}

clean_wf_comp_list <- function(ppg_ii_vs_wf_raw) {
  ppg_ii_vs_wf_raw |>
    filter(taxonRank %in% c("order", "family", "genus")) |>
    filter(!str_detect(scientificName, "Japanobotrychium")) |>
    filter(!str_detect(scientificName, "Japanobotrychum")) |>
    pull(scientificName)
}


check_ppg_higher_tax_changes <- function(ppg_ii, ppg_i, ppg_issues) {
  # Get passed issues for all ranks (not just genus)
  # We'll check if the taxon name matches the genus OR any of the
  # higher ranks that changed
  passed_issues_all <- ppg_issues |>
    filter(str_detect(status, "^PASSED$")) |>
    mutate(
      # Standardize rank names
      rank = str_to_lower(rank),
      # Clean up taxon names: remove 'and', commas, normalize spacing
      name = str_remove_all(name, regex("and|,", ignore_case = TRUE)) |>
        str_squish()
    ) |>
    # Separate multiple taxa into individual rows
    separate_rows(name, sep = "\\s+") |>
    mutate(
      # Remove × or x prefix from nothogenera
      name = str_remove(name, "^[×x]") |>
        str_trim()
    ) |>
    filter(name != "") |>
    select(number, rank, name) |>
    distinct()

  ppg_ii_long <- ppg_ii |>
    pivot_longer(names_to = "rank", values_to = "taxon", everything()) |>
    filter(rank != "genus") |>
    unique() |>
    filter(!is.na(taxon))

  ppg_i_long <- ppg_i |>
    select(class:genus) |>
    pivot_longer(names_to = "rank", values_to = "taxon", everything()) |>
    filter(rank != "genus") |>
    unique() |>
    filter(!is.na(taxon))

  ppg_ii_only <-
    ppg_ii_long |>
    anti_join(ppg_i_long) |>
    mutate(diff = "In PPG II but not in PPG I")

  ppg_i_only <-
    ppg_i_long |>
    anti_join(ppg_ii_long) |>
    mutate(diff = "In PPG I but not in PPG II")

  bind_rows(
    ppg_ii_only,
    ppg_i_only
  ) |>
    arrange(rank, taxon) |>
    left_join(
      select(passed_issues_all, -rank),
      by = join_by(taxon == name)
    ) |>
    # manual fixes
    mutate(
      number = case_when(
        taxon == "Ctenitidoideae" ~ 92L,
        taxon == "Drynarioideae" ~ 52L,
        .default = number
      )
    )
}

# Apendix ----

italicize_subgen_single <- function(x) {
  x <- paste0("*", x, "*")
  x <- str_replace_all(
    x,
    "subg. ",
    "*subg. *"
  )
}

italicize_subgen <- function(x) {
  map_chr(x, italicize_subgen_single)
}


italicize_gen_single <- function(x) {
  if (is.na(x)) {
    return(NA)
  }
  is_hybrid <- str_detect(x, "×")
  if (is_hybrid) {
    res <- str_remove_all(x, "×") |> str_squish()
    res <- paste0("*", res, "*")
    res <- paste("×", res)
  } else {
    res <- paste0("*", x, "*")
  }
  str_squish(res)
}

italicize_gen <- function(x) {
  map_chr(x, italicize_gen_single)
}

count_ppgi_gen <- function(ppg_i, ...) {
  ppg_i |>
    filter(...) |>
    summarize(n = n_distinct(genus)) |>
    pull(n)
}

#' Create Phylogenetic Tree Figure for PPG II
#'
#' Generates a comprehensive cladogram showing relationships between
#' families of tracheophytes recognized by PPG II. The tree shows number of
#' taxonomic changes relative to PPG I.
#'
#' @param phy_family A phylo object containing the family-level fern
#'   phylogeny from the Fern Tree of Life.
#' @param ppg A cleaned data frame of PPG taxonomic data, typically the
#'   output of clean_ppg().
#' @param ppg_ii PPG II data in wide format
#' @param ppg_issues_count A data frame with category of taxonomic change
#' (split vs. sink) for each taxonomic proposal
#'
#' @return A ggtree object showing the phylogenetic tree with families
#'
make_tree_figure_appendix <- function(
  phy_family,
  ppg,
  ppg_ii,
  ppg_issues_count
) {
  require(ggtree)
  require(ggimage)
  require(ggrepel)
  require(ape)
  require(phytools)
  require(ggnewscale)

  ppg_higher <- ppg |>
    filter_out(
      scientificName == "Dennstaedtia",
      scientificNameAuthorship == "Bernh."
    ) |>
    mutate(
      acceptedNameUsageID = case_when(
        taxonomicStatus == "accepted" ~ taxonID,
        .default = acceptedNameUsageID
      )
    )

  ppg_issues_count_for_higher <-
    ppg_issues_count |>
    select(
      name,
      rank,
      change
    ) |>
    left_join(
      select(
        ppg_higher,
        name = scientificName,
        acceptedNameUsageID
      ),
      by = "name"
    ) |>
    filter_out(change == "none") |>
    unique() |>
    assert(not_na, acceptedNameUsageID)

  ppg_issues_count_genus <-
    ppg_issues_count_for_higher |>
    filter(rank %in% c("Genus", "Nothogenus")) |>
    left_join(
      select(
        ppg_higher,
        accepted_name = scientificName,
        acceptedNameUsageID = taxonID
      ),
      by = "acceptedNameUsageID"
    ) |>
    assert(not_na, accepted_name) |>
    left_join(
      ppg_ii,
      by = join_by(accepted_name == genus)
    ) |>
    unique() |>
    assert(not_na, family) |>
    assert(is_uniq, name)

  ppg_issues_count_subfamily <-
    ppg_issues_count_for_higher |>
    filter(rank == "Subfamily") |>
    left_join(
      select(
        ppg_higher,
        accepted_name = scientificName,
        acceptedNameUsageID = taxonID
      ),
      by = "acceptedNameUsageID"
    ) |>
    assert(not_na, accepted_name) |>
    left_join(
      unique(select(ppg_ii, -genus)),
      by = join_by(accepted_name == subfamily)
    ) |>
    assert(not_na, family) |>
    assert(is_uniq, name)

  ppg_issues_count_family <-
    ppg_issues_count_for_higher |>
    filter(rank == "Family") |>
    left_join(
      select(
        ppg_higher,
        accepted_name = scientificName,
        acceptedNameUsageID = taxonID
      ),
      by = "acceptedNameUsageID"
    ) |>
    assert(not_na, accepted_name) |>
    left_join(
      unique(select(ppg_ii, -c(genus, subfamily))),
      by = join_by(accepted_name == family)
    ) |>
    mutate(family = accepted_name) |>
    assert(not_na, family) |>
    assert(is_uniq, name)

  ppg_issues_count_higher <-
    ppg_issues_count_genus |>
    bind_rows(ppg_issues_count_subfamily) |>
    bind_rows(ppg_issues_count_family) |>
    assert(not_na, family, name) |>
    assert(is_uniq, name) |>
    # Group change into fct called change_fct,
    # combine values other than "sink" or "split" into "other"
    mutate(
      change_fct = case_when(
        change %in% c("sink", "split") ~ change,
        TRUE ~ "other"
      ) |>
        factor(
          levels = c("sink", "split", "other")
        )
    )

  ppg_issues_count_by_family <-
    ppg_issues_count_higher |>
    group_by(family) |>
    count(change_fct) |>
    pivot_wider(
      names_from = change_fct,
      values_from = n,
      values_fill = 0
    ) |>
    ungroup()

  # Add other tracheophytes to fern tree
  phy_tracheo <- fern_to_tracheo_phy(phy_family)

  # Format family labels with genus/species counts. Keep num species
  # for tip points
  family_tip_labels <-
    tibble(
      taxon_name = phy_tracheo$tip.label
    ) |>
    left_join(
      ppg_issues_count_by_family,
      by = join_by(taxon_name == family),
      relationship = "one-to-one"
    ) |>
    mutate(across(c(split, sink, other), ~ replace_na(., 0))) |>
    filter_out(taxon_name == "Spermatophytes") |>
    bind_rows(
      tibble(
        taxon_name = "Spermatophytes",
        family_label = "SEED PLANTS",
        seed_plant_image = "images/starburst.png"
      )
    ) |>
    mutate(
      family_label = case_when(
        taxon_name != "Spermatophytes" ~ glue::glue(
          "{taxon_name} ({split}/{sink})"
        ) |>
          as.character(),
        .default = taxon_name
      )
    ) |>
    select(tip = taxon_name, family_label, seed_plant_image, split, sink, other)

  # Get root node for manually drawing root edge
  # (needed because geom_rootedge() doesn't work with branch.length = "none")
  n_tips <- length(phy_tracheo$tip.label)
  root_node_num <- n_tips + 1

  # Get the actual y-coordinate of the root from the tree layout
  tree_layout <- as_tibble(ggtree::fortify(phy_tracheo))
  root_y <- tree_layout |>
    filter(node == root_node_num) |>
    pull(y)

  # set font sizes etc
  fig_font_size <- 2.5
  clade_lab_offset <- 9
  branch_lab_vjust <- -0.6
  root_edge_length <- 1

  # generate figure
  tree_fig <-
    # Base tree, with line type by uncertainty
    ggtree(phy_tracheo, branch.length = "none") %<+%
    # Add datasets
    family_tip_labels +
    # Diamonds scaled by number of splits
    geom_tippoint(
      aes(color = split),
      shape = 18,
      size = 4.5,
      position = position_nudge(x = 0.2)
    ) +
    # Set color for splits and sinks to same maximum value
    scale_color_viridis_c(
      option = "D",
      limits = c(0, max(family_tip_labels$split, na.rm = TRUE)),
      na.value = "transparent",
      name = "Number of changes"
    ) +
    # Circles scaled by number of sinks
    # set color for NA values to transparent
    new_scale_color() +
    geom_tippoint(
      aes(color = sink),
      shape = 16,
      size = 4,
      position = position_nudge(x = 1.1)
    ) +
    scale_color_viridis_c(
      option = "D",
      limits = c(0, max(family_tip_labels$split, na.rm = TRUE)),
      na.value = "transparent",
      guide = "none"
    ) +
    # Tip label (seed plants)
    geom_tiplab(
      aes(image = seed_plant_image),
      geom = "image",
      shape = 11,
      size = 0.019,
      nudge_x = -0.1
    ) +
    # Tip labels (family)
    geom_tiplab(
      aes(label = family_label),
      size = fig_font_size,
      offset = 1.5
    ) +
    # Manual root edge (horizontal line at root)
    geom_segment(
      aes(x = -root_edge_length, xend = 0, y = root_y, yend = root_y),
      linewidth = 0.5
    ) +
    xlim(-4, 33) +
    theme(
      # Set legend in upper left
      legend.position = c(0.05, 0.95),
      legend.justification = c("left", "top")
    )

  tree_fig
}

# Run all test files in a directory
# Error if FALSE, return TRUE if passes
run_tests <- function(dir) {
  testthat::test_dir(
    dir,
    reporter = "progress",
    stop_on_failure = TRUE
  )
  TRUE
}
