# Load required packages

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# load entanglement data; code uniprot and entanglement variables as characters
AF2_entangled <- readr::read_delim("raw-data/AF_entanglement.dat",
                                   delim = ",",
                                   skip  = 1,
                                   col_names = c("uniprot", 
                                                 "entanglement"),
                                   col_types = list(uniprot      = col_character(),
                                                    entanglement = col_character()))

# test that all rows correspond to a unique value of uniprot
stopifnot(length(unique(AF2_entangled$uniprot)) == nrow(AF2_entangled))

# load the Reactome database file and filter out entries that do not correspond
# to human proteins
human_reactome <- readr::read_delim("raw-data/UniProt2Reactome.txt", 
                                    delim = "\t", 
                                    col_names = c("uniprot", 
                                                  "pathway", 
                                                  "reactome_link", 
                                                  "pathway_name", 
                                                  "evidence_code", 
                                                  "organism"),
                                    col_types = list(uniprot       = col_character(),
                                                     pathway       = col_character(),
                                                     reactome_link = col_character(),
                                                     pathway_name  = col_character(),
                                                     evidence_code = col_character(),
                                                     organism      = col_character())) %>%
                                    dplyr::filter(organism == "Homo sapiens")

# create a joined tibble with one row (entry) for each combination of uniprot
# and pathway identifiers. Note that I replace NA with "dummy"
joined_tibble <- left_join(AF2_entangled, human_reactome, by = "uniprot") %>%
                 mutate(pathway=tidyr::replace_na(pathway, "dummy"))

# make a list of the unique pathways in joined_tibble and remove "dummy"
unique_pathways <- unique(joined_tibble$pathway)

unique_pathways <- unique_pathways[! unique_pathways %in% c("dummy")]

# setup counter and vectors for calculation
i                           <- 1

N_usable_pathways           <- length(unique_pathways)

a_in_pathway_entangled      <- vector(mode = "integer", 
                                      length = N_usable_pathways)

b_out_pathway_entangled     <- vector(mode = "integer", 
                                      length = N_usable_pathways)

c_in_pathway_not_entangled  <- vector(mode = "integer", 
                                      length = N_usable_pathways)

d_out_pathway_not_entangled <- vector(mode = "integer", 
                                      length = N_usable_pathways)

two_sided_p_values          <- vector(mode = "numeric", 
                                      length = N_usable_pathways)

enrich_p_values             <- vector(mode = "numeric", 
                                      length = N_usable_pathways)

deplete_p_values            <- vector(mode = "numeric", 
                                      length = N_usable_pathways)

odds_ratios                 <- vector(mode = "numeric", 
                                      length = N_usable_pathways)

# note that in this instance, as I am able to dimension the vectors before the "for" loop,
# we do not expect vectorization in R to significantly increase speed. The use of a "for"
# loop here maintains readability by people less familiar with R

# loop over pathway identifiers
for (u in unique_pathways) {
  
  # vector of unique uniprot that have this pathway identifier and are entangled
  aVec <- joined_tibble %>% 
          filter(pathway == u & entanglement == "Yes") %>% 
          .$uniprot %>% 
          unique()
  
  a    <- length(aVec)
  
  # vector of unique uniprot that do not have this pathway identifier and are entangled
  # N.B., there may uniprot in both aVec and bVec due to degenerate pathway names; we need to remove dups from bVec
  bVec <- joined_tibble %>% 
          filter(pathway != u & entanglement == "Yes") %>% 
          .$uniprot %>% 
          unique()
  
  bVec <- bVec[! bVec %in% aVec]
  
  b    <- length(bVec)
  
  # vector of unique uniprot that have this pathway identifier and are NOT entangled
  cVec <- joined_tibble %>% 
          filter(pathway == u & entanglement == "No")  %>% 
          .$uniprot %>% 
          unique()
  
  c    <- length(cVec)
  
  # vector of unique uniprot that do not have this pathway identifier and are not entangled
  # N.B., there may be uniprot in both cVec and dVec due to degenerate pathways; we need to remove dups from dVec
  dVec <- joined_tibble %>% filter(pathway != u & entanglement == "No")  %>% .$uniprot %>% unique()
  dVec <- dVec[! dVec %in% cVec]
  d    <- length(dVec)
  
  # compute the sum of {a, b, c, d} and test that value equals number of unique uniprot names
  Ntotal <- a + b + c + d
  stopifnot(Ntotal == length(AF2_entangled$uniprot))
  
  # construct contingency table as a matrix
  contingency_table <- matrix(c(a, c, b, d), nrow = 2, ncol = 2)
  
  # record values of {a, b, c, d} in vectors
  a_in_pathway_entangled[i]      <- a
  b_out_pathway_entangled[i]     <- b
  c_in_pathway_not_entangled[i]  <- c
  d_out_pathway_not_entangled[i] <- d
  
  # perform two-sided test and record results
  double_tail_result    <- stats::fisher.test(contingency_table, alternative = "two.sided")
  two_sided_p_values[i] <- double_tail_result$p.value
  odds_ratios[i]        <- double_tail_result$estimate
  
  # perform one-sided test (greater) and record result (enrichment)
  enrichment_result     <- stats::fisher.test(contingency_table, alternative = "greater")
  enrich_p_values[i]    <- enrichment_result$p.value
  
  # perform one-sided test (less) and record result (depletion)
  depletion_result      <- stats::fisher.test(contingency_table, alternative = "less")
  deplete_p_values[i]   <- depletion_result$p.value
  
  # increment counter
  i <- i + 1
  
}

# correct p-values for False Discovery Rate
corrected_p_values <- stats::p.adjust(two_sided_p_values, method = "BH")

# construct tibble of results and save to file
results <- dplyr::tibble(unique_pathways, 
                         a_in_pathway_entangled,
                         b_out_pathway_entangled, 
                         c_in_pathway_not_entangled, 
                         d_out_pathway_not_entangled,
                         odds_ratios, 
                         two_sided_p_values,
                         corrected_p_values,
                         enrich_p_values,
                         deplete_p_values)

readr::write_csv(results, "processed-data/analysis_results.csv")