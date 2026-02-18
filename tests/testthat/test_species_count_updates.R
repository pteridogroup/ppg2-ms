# Test script for apply_species_count_updates()

library(dplyr)
library(tibble)

test_that("apply_species_count_updates works correctly", {
  # Create test data
  # Simple taxonomy: Order -> Family -> 2 Genera
  test_ppg <- tribble(
    ~taxonID , ~scientificName , ~taxonRank , ~taxonomicStatus , ~parentNameUsageID ,
    "1"      , "Testales"      , "order"    , "accepted"       , NA                 ,
    "2"      , "Testaceae"     , "family"   , "accepted"       , "1"                ,
    "3"      , "Testus"        , "genus"    , "accepted"       , "2"                ,
    "4"      , "Mockus"        , "genus"    , "accepted"       , "2"
  )

  # Original counts: Testus has 10 species, Mockus has 5 species
  # Total for family should be 15
  test_children_tally <- tribble(
    ~taxonID , ~n_species , ~n_genera ,
    "1"      ,         15 ,         2 , # Order: 15 species, 2 genera
    "2"      ,         15 ,         2 , # Family: 15 species, 2 genera
    "3"      ,         10 , NaN       , # Testus: 10 species
    "4"      ,          5 , NaN # Mockus: 5 species
  )

  # Manual update: Change Testus to 15 species (increase of 5)
  test_updates <- tribble(
    ~genus   , ~species ,
    "Testus" ,       15
  )

  # Apply the updates
  result <- apply_species_count_updates(
    test_children_tally,
    test_updates,
    test_ppg
  )

  # Run tests
  expect_equal(
    result |> filter(taxonID == "3") |> pull(n_species),
    15,
    label = "Testus updated from 10 to 15"
  )

  expect_equal(
    result |> filter(taxonID == "4") |> pull(n_species),
    5,
    label = "Mockus unchanged at 5"
  )

  expect_equal(
    result |> filter(taxonID == "2") |> pull(n_species),
    20,
    label = "Family total updated to 20"
  )

  expect_equal(
    result |> filter(taxonID == "1") |> pull(n_species),
    20,
    label = "Order total updated to 20"
  )

  expect_equal(
    result |> filter(!is.na(n_genera)) |> pull(n_genera),
    c(2, 2),
    label = "Genus counts unchanged"
  )
})
