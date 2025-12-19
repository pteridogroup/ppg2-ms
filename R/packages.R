# Load conflicted first to manage function conflicts
library(conflicted)

# Set conflict preferences
conflicted::conflicts_prefer(tidyr::replace_na)

# Load other packages
library(tidyverse)
library(targets)
library(tarchetypes)
library(assertr)
library(ape)
library(gh)
