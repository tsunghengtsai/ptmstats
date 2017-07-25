# Load required libraries -------------------------------------------------
library(Biostrings)
library(stringr)
library(tibble)
library(dplyr)
library(purrr)


# Annotation with FASTA ---------------------------------------------------
# http://www.uniprot.org/help/fasta-headers
fasta_path <- "data/Sequence/homo_sapiens_all_20160725.fasta"

hs_fasta <- Biostrings::readAAStringSet(fasta_path)
hs_fasta <- as.data.frame(hs_fasta) %>% 
    tbl_df() %>% 
    rownames_to_column(var = "header") %>% 
    select(header, sequence = x)

regex_uniprot_ac <- regex("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})")
regex_uniprot_iso <- regex("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})([-]\\d{1,}){0,1}")

regex_aa <- "A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y"

hs_fasta <- hs_fasta %>% 
    mutate(
        uniprot_ac = str_extract(header, pattern = regex_uniprot_ac), 
        uniprot_iso = str_extract(header, pattern = regex_uniprot_iso), 
        p_entry = str_extract(header, pattern = "([^\\s\\|]*)(?=\\s)")
    )

# define modification residues and symbols
mod_residue <- "K"
mod_symbol <- "\\*"

# [TODO]: seems no need to export mod_bdy
hs_fasta <- hs_fasta %>% 
    mutate(
        mod_bdy = str_locate_all(sequence, mod_residue), 
        mod_idx = map(mod_bdy, ~.[, "start"]), 
        mod_res = str_match_all(sequence, mod_residue), 
        mod_res = map(mod_res, ~.[, 1])
    )

# save.image("output/uniprot_inf.RData")

