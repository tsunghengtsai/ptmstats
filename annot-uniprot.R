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

# regex_uniprot <- ".*([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}).*"
regex_uniprot_ac <- regex("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})")
regex_uniprot_iso <- regex("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})([-]\\d{1,}){0,1}")
# regex_entry <- regex("\\|.* ")

regex_aa <- "A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y"

# hs_fasta <- hs_fasta %>% 
#     mutate(uniprot_ac = str_replace(header, pattern = regex_uniprot, replacement = "\\1"))
hs_fasta <- hs_fasta %>% 
    mutate(
        uniprot_ac = str_extract(header, pattern = regex_uniprot_ac), 
        uniprot_iso = str_extract(header, pattern = regex_uniprot_iso), 
        p_entry = str_extract(header, pattern = "([^\\s\\|]*)(?=\\s)")
    )

# define modification residues and symbols
mod_residue <- "K"
mod_symbol <- "\\*"

hs_fasta <- hs_fasta %>% 
    mutate(
        mod_bdy = str_locate_all(sequence, mod_residue), 
        mod_idx = map(mod_bdy, ~.[, "start"]), 
        mod_res = str_match_all(sequence, mod_residue), 
        mod_res = map(mod_res, ~.[, 1])
    )

# save.image("output/uniprot_inf.RData")


# Site representation for modified peptides -------------------------------
library(tidyverse)
library(stringr)

load("output/uniprot_inf.RData")
load("output/dta_conflu.RData")

df_mod <- dta_conflu %>% 
    distinct(Reference, peptide_sequence) %>% 
    separate(Reference, into = c("p_entry", "uniprot_iso"), sep = "\\|") %>% 
    mutate(
        unmod_peptide = str_replace_all(peptide_sequence, "\\*", ""), 
        # trim_peptide = str_extract(unmod_peptide, regex("\\.[ACDEFGHIKLMNPQRSTVWY]+")) %>% 
        #     str_replace("\\.", ""), 
        trimmed_pep = str_extract(unmod_peptide, "(?<=\\.)([ACDEFGHIKLMNPQRSTVWY\\*]+)"), 
        trimmed_mod = str_extract(peptide_sequence, "(?<=\\.)([ACDEFGHIKLMNPQRSTVWY\\*]+)"), 
        len_peptide = str_length(trimmed_pep)
    ) %>% 
    arrange(uniprot_iso)

hs_fasta <- hs_fasta %>% 
    select(uniprot_ac, uniprot_iso, p_entry, mod_idx, mod_res, header, sequence) %>% 
    arrange(uniprot_iso)


# Not all IsoID from experimental data covered by the fasta
# length(which(!(df_mod$uniprot_iso %in% hs_fasta$uniprot_iso)))

# IsoID and protein entry are not defined unambiguously
#  --> Use the entries from FASTA for consistency
df_mod %>% anti_join(hs_fasta) %>% distinct(uniprot_iso)
df_mod %>% anti_join(hs_fasta %>% select(-p_entry)) %>% distinct(uniprot_iso)
sum(!(unique(df_mod$uniprot_iso) %in% unique(hs_fasta$uniprot_iso)))

df_mod %>% distinct(uniprot_iso, p_entry) %>% filter(uniprot_iso == "Q9UHQ7")
hs_fasta %>% filter(uniprot_iso == "Q9UHQ7") %>% distinct(uniprot_iso, p_entry)

df_mod %>% distinct(uniprot_iso, p_entry) %>% filter(uniprot_iso == "Q969X6")
hs_fasta %>% filter(uniprot_iso == "Q969X6") %>% distinct(uniprot_iso, p_entry)

df_mod %>% distinct(uniprot_iso, p_entry) %>% filter(uniprot_iso == "Q92616")
hs_fasta %>% filter(uniprot_iso == "Q92616") %>% distinct(uniprot_iso, p_entry)

df_mod %>% distinct(uniprot_iso, p_entry) %>% filter(uniprot_iso == "P63244")
hs_fasta %>% filter(uniprot_iso == "P63244") %>% distinct(uniprot_iso, p_entry)

df_mod %>% distinct(uniprot_iso, p_entry) %>% filter(uniprot_iso == "X00760")

df_mod %>% distinct(uniprot_iso, p_entry) %>% filter(uniprot_iso == "P0CW22")

df_mod %>% distinct(uniprot_iso, p_entry) %>% filter(uniprot_iso == "P08107")


# Exploring with matched subset first (ingoring 3 IsoID)
# df_mod %>% semi_join(hs_fasta)
df_fasmod <- df_mod %>% 
    filter(
        len_peptide >= 6, 
        str_detect(p_entry, "HUMAN")
    ) %>% 
    select(-p_entry) %>% 
    inner_join(hs_fasta %>% select(-header))

df_fasmod <- df_fasmod %>% 
    mutate(
        pep_idx = str_locate_all(sequence, trimmed_pep), 
        nb_mch = map_int(pep_idx, ~nrow(.))
    )

# Peptide with no mapping to sequence? (seems a mistake with the first AA)
df_fasmod %>% filter(nb_mch == 0)

# Peptide mapped to multiple locations (using extended AA only helps a bit)
df_fasmod %>% 
    filter(nb_mch > 1) %>% 
    select(uniprot_iso, peptide_sequence, trimmed_pep, nb_mch)

# df_fasmod <- df_fasmod %>% 
#     mutate(ext_peptide = str_replace_all(unmod_peptide, "\\.", ""))

# Ignoring non-specific matching for now
df_fasmod <- df_fasmod %>% filter(nb_mch == 1)

# Checking the boundary of peptide sequence
df_fasmod %>% filter(str_detect(peptide_sequence, "^\\.|\\.$"))
df_fasmod %>% filter(str_detect(peptide_sequence, "^.\\."))  # Left bound is strictly defined
df_fasmod %>% filter(str_detect(peptide_sequence, "[^\\.]{2}$"))  # Right bound is loosely defined

# Making sure no modifications occur on extended AA
df_fasmod %>% filter(str_detect(peptide_sequence, "^K\\*\\.|\\.K\\*$"))

# Making sure modification residues are fully defined
df_fasmod %>% 
    mutate(
        nb_smod = str_count(peptide_sequence, mod_symbol), 
        nb_mod = str_count(peptide_sequence, str_c(mod_residue, mod_symbol))
    ) %>% 
    filter(nb_mod != nb_smod)

# Modification sites
df_fasmod <- df_fasmod %>% 
    mutate(
        site_all = map2(str_locate_all(trimmed_pep, mod_residue), pep_idx,  
                       ~ .x[, "start"] + .y[, "start"] - 1), 
        site_mod = map2(str_locate_all(trimmed_mod, str_c(mod_residue, mod_symbol)), pep_idx,  
                       ~ .x[, "start"] - seq_along(.x[, "start"]) + .y[, "start"]), 
        site_null = map2(site_all, site_mod, setdiff),
        nb_mod = str_count(trimmed_mod, str_c(mod_residue, mod_symbol))
    )


df_fasmod %>% select(uniprot_iso, trimmed_mod, trimmed_pep, site_all, site_mod, nb_mod)
df_fasmod %>% count(nb_mod)
df_fasmod %>% select(uniprot_iso, trimmed_mod, trimmed_pep, site_all, site_mod, nb_mod) %>%
    filter(nb_mod == 0)
df_fasmod %>% select(uniprot_iso, trimmed_mod, trimmed_pep, site_all, site_mod, nb_mod) %>%
    filter(nb_mod == 1)
df_fasmod %>% select(uniprot_iso, trimmed_mod, trimmed_pep, site_all, site_mod, nb_mod) %>%
    filter(nb_mod == 2)


# Annotate modification site
ant_site <- function(index, index_full, residue) {
    if (is_empty(index)) {
        ant_str <- "UNMOD"
    } else {
        ant_len <- max(str_length(index_full))
        ant_idx <- str_pad(index, width = ant_len, pad = "0")
        ant_str <- str_c(residue, ant_idx, collapse = "-")
    }

    return(ant_str)
}

df_fasmod <- df_fasmod %>% 
    mutate(
        site_str = map2_chr(site_mod, mod_idx, ant_site, mod_residue), 
        # site_str = map_chr(site_mod, ant_site, mod_residue), 
        pep_str = map_chr(pep_idx, ~ str_c(., collapse = "-")), 
        pps_str = str_c(uniprot_iso, pep_str, site_str, sep = "_")
    )

# saveRDS(df_fasmod, "output/fasmod.rds")


# Explore FASTA -----------------------------------------------------------
load("output/annot_uniprot.RData")

nrow(hs_fasta)
n_distinct(hs_fasta$protein)
n_distinct(hs_fasta$sequence)
n_distinct(hs_fasta$uniprot_ac)
n_distinct(hs_fasta$uniprot_iso)

# Same sequence with multiple Uniprot AC
dup_seq <- hs_fasta %>% 
    distinct(sequence, uniprot_ac) %>% 
    group_by(sequence) %>% 
    filter(n() > 1) %>% 
    ungroup()

dup_seq %>% arrange(sequence) %>% select(-sequence)

# Same Uniprot AC to multiple sequences 
# Example: http://www.uniprot.org/uniprot/A0A096LP49
dup_uniprot <- hs_fasta %>% 
    distinct(sequence, uniprot_ac) %>% 
    group_by(uniprot_ac) %>% 
    filter(n() > 1) %>% 
    ungroup()

dup_uniprot %>% arrange(uniprot_ac) %>% select(uniprot_ac)

# Use isoform ID to resolve ambiguity due to alternative splicing
hs_fasta %>% 
    distinct(sequence, uniprot_iso) %>% 
    group_by(uniprot_iso) %>% 
    filter(n() > 1) %>% 
    ungroup()

# Still, there are same sequences mapped to multiple isoform ID
dup_seq <- hs_fasta %>% 
    distinct(sequence, uniprot_iso) %>% 
    group_by(sequence) %>% 
    filter(n() > 1) %>% 
    ungroup()

dup_seq %>% arrange(sequence) %>% select(-sequence)


