# Load required libraries -------------------------------------------------
library(Biostrings)
library(stringr)
library(tibble)
library(dplyr)


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

hs_fasta <- hs_fasta %>% 
    mutate(
        uniprot_ac = str_extract(header, pattern = regex_uniprot_ac), 
        uniprot_iso = str_extract(header, pattern = regex_uniprot_iso), 
        entry_name = str_extract(header, pattern = "([^\\s\\|]*)(?=\\s)")
    )

# saveRDS(hs_fasta, "output/hs_fasta_20160725.rds")
