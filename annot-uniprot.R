## site/residue centered representation -----------------------------
library(Biostrings)
library(stringr)
library(tibble)
library(dplyr)

fasta_path <- "data/Sequence/homo_sapiens_all_20160725.fasta"

hs_fasta <- Biostrings::readAAStringSet(fasta_path)
hs_fasta <- as.data.frame(hs_fasta) %>% 
    tbl_df() %>% 
    rownames_to_column(var = "protein") %>% 
    select(protein, sequence = x)

regex_uniprot <- ".*([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}).*"
regex_aa <- "A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y"

hs_fasta <- hs_fasta %>% 
    mutate(uniprot_ac = str_replace(protein, pattern = regex_uniprot, replacement = "\\1"))

# hs_fasta %>% count(uniprot_ac) %>% filter(n > 1)

## define modification residues and symbols
mod_residue <- "K"
mod_symbol <- "\\*"

mod_indices <- lapply(hs_fasta$sequence, function(x) as.vector(str_locate_all(x, pattern = mod_residue)[[1]][, 1]))
mod_residues <- lapply(hs_fasta$sequence, function(x) c(str_match_all(x, pattern = mod_residue)[[1]]))
nb_mods <- str_count(hs_fasta$sequence, mod_residue)

## [TODO]: not an accurate representation as a uniprot_ac may match to >1 isoforms
mod_protein <- data_frame(uniprot_ac = rep(hs_fasta$uniprot_ac, nb_mods), 
                          mod_res = unlist(mod_residues), 
                          mod_idx = unlist(mod_indices))

save.image("output/annot_uniprot.RData")