## load required packages -------------------------------------------
library(readr)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)


## set up working directory and load datasets -----------------------
# work_dir <- "/Users/thtsai/Projects/NEU/ptm/dev/"
# data_dir <- "/Users/thtsai/Projects/NEU/ptm/data/CCCP_confluency/"
# setwd(work_dir)


## vista outputs
df_design <- read.table("data/CCCP_confluency/CellDensity_tft_forTH.txt", 
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)

file_name <- paste0("vistaga_", df_design$run_id, ".txt")

sub_cols <- c("Reference", "run_id", "peptide_id", "vista_peak_area_light", 
              "vista_confidence_score", "peptide_m_z", "peptide_sequence", 
              "peptide_trypticity", "peptide_miscleavages", "peptide_validity", 
              "label1", "ms2_rt", "ms2_charge", "is_xq")

dta_conflu <- vector("list", length(file_name))
for (i in seq_along(file_name)) {
    file_path <- paste0("data/CCCP_confluency/", file_name[i])
    dta_conflu[[i]] <- read.table(file_path, header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE) %>% 
        tbl_df() %>% select_(.dots = sub_cols)
}
dta_conflu <- bind_rows(dta_conflu)
rm(list = c("i", "file_name", "file_path"))

# save.image("output/dta_conflu.RData")
# save.image("~/Projects/NEU/ptm/dev/dta_conflu.RData")


## starting workflow ------------------------------------------------
## load dataset
# work_dir <- "/Users/thtsai/Projects/NEU/ptm/dev/"
# setwd(work_dir)
load("output/dta_conflu.RData")

dta_conflu <- dta_conflu %>% 
    filter(vista_confidence_score >= 83)

# df_design %>% arrange(cond_treatment, cond_confluency, id_subject, id_injectionset) %>% 
#     select(-ann_label)

df_design <- tbl_df(df_design)

dta_conflu <- dta_conflu %>% 
    select(-vista_confidence_score, -peptide_trypticity, -peptide_miscleavages, 
           -peptide_validity)

## OBJ# can be viewed as a batch identity
uniq_subj <- unique(df_design$id_subject)
key_subj <- paste0("S", sprintf("%02d", 1:length(uniq_subj)))
df_design <- df_design %>% 
    mutate(subj = plyr::mapvalues(id_subject, from = uniq_subj, to = key_subj), 
           biorep = ifelse(ann_objective_id == 37122, "B1", "B2"), 
           techrep = paste0("T", id_injectionset), 
           run_bt = paste0(biorep, techrep), 
           run_cbt = paste(id_subject, run_bt, sep = "-"))


## data annotation and filtering ------------------------------------
## [TODO]: need a closer look at `Reference` & `label1`
dta_conflu <- dta_conflu %>% 
    mutate(is_mod = grepl("\\*", peptide_sequence), 
           feature = paste(peptide_sequence, ms2_charge, sep = "_"), 
           log2inty = log2(vista_peak_area_light), 
           cond = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$id_subject), 
           biorep = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$biorep), 
           techrep = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$techrep), 
           run_cbt = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$run_cbt))

## features with no charge state or retention time inforamtion
# dta_conflu %>% filter(is.na(ms2_charge) | is.na(ms2_rt)) %>%
#     select(feature, peptide_sequence, ms2_charge, ms2_rt)

## short matches
# dta_conflu %>% filter(grepl("-\\..$", peptide_sequence)) %>%
#     select(feature, peptide_sequence)

## not belonging to expected uniprot ID
regex_uniprot <- ".*([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}).*"
# dta_conflu %>% filter(!grepl(pattern = regex_uniprot, Reference))

## preliminary filtering
dta_conflu <- dta_conflu %>%
    filter(!is.na(ms2_charge), !is.na(ms2_rt), !grepl("-\\..$", peptide_sequence), grepl(regex_uniprot, Reference))


## extract uniprot accession number
uniq_ref <- dta_conflu %>% distinct(Reference) %>% .$Reference
uniq_unip <- unlist(lapply(uniq_ref, function(x) str_extract(unlist(str_split(x, '\\|'))[2], regex_uniprot)))
dta_conflu <- dta_conflu %>% 
    mutate(uniprot_ac = plyr::mapvalues(Reference, from = uniq_ref, to = uniq_unip))


## data exploration -------------------------------------------------
## note the proportion of modified peptide in different biological replicates
table(dta_conflu$run_cbt, dta_conflu$is_mod)

## check for duplicate features 
dta_conflu %>% 
    count(run_id, feature) %>% ungroup() %>% 
    filter(n > 1)

## 2158 features from 16 runs have duplicates (4984 features in total)
dta_conflu %>% 
    count(run_id, feature) %>% ungroup() %>% 
    filter(n > 1) %>% 
    distinct(feature)

## 814 with additional criterion for max RT diff > 1.5 min
dta_conflu %>% 
    group_by(run_id, feature) %>% 
    summarise(n = n(), d_rt = max(ms2_rt, na.rm = TRUE) - min(ms2_rt, na.rm = TRUE)) %>% 
    ungroup() %>% 
    filter(n > 1, d_rt > 1.5)

## check for mismatched features (different RT in different runs for the same peptide)
## to be carried out after duplicates within run are handled
# dta_conflu %>% 
#     group_by(feature) %>% 
#     summarise(d_rt = max(ms2_rt, na.rm = TRUE) - min(ms2_rt, na.rm = TRUE)) %>% 
#     ungroup() %>% 
#     filter(d_rt > 5)


## protein & protein group
length(unique(dta_conflu$Reference))


## protein (Reference)
## [TODO]: check if uniprot # is a better proxy
dta_conflu %>% select(Reference, run_id, peptide_sequence) %>% distinct() %>% 
    group_by(run_id, peptide_sequence) %>% mutate(nb_prot = n()) %>% ungroup() %>% 
    filter(nb_prot > 1) %>% arrange(peptide_sequence)

dta_conflu %>% select(Reference, peptide_sequence) %>% distinct() %>% 
    group_by(peptide_sequence) %>% mutate(nb_prot = n()) %>% ungroup() %>% 
    filter(nb_prot > 1) %>% arrange(peptide_sequence)


## filtering before checking intensity profile ----------------------
## use the max for duplicate features (per run)
## [TODO]: conisder is_xq, ms2_rt
dta_conflu <- dta_conflu %>% 
    group_by(Reference, uniprot_ac, peptide_sequence, is_mod, feature, ms2_charge, cond, biorep, run_cbt, run_id) %>% 
    summarise(log2inty = max(log2inty)) %>% 
    ungroup()


## data manipulation and transformation -----------------------------
## remove peptides with ambiguous matches (to >1 uniprot_ac) for now
nb_pep2uniprot <- dta_conflu %>% 
    distinct(peptide_sequence, uniprot_ac) %>% 
    count(peptide_sequence)

dta_conflu <- dta_conflu %>% 
    semi_join(filter(nb_pep2uniprot, n == 1))


## normalization (equalizing median)
med_run <- dta_conflu %>% 
    group_by(run_id) %>% 
    summarise(log2inty_med = median(log2inty, na.rm = TRUE)) %>% 
    mutate(log2inty_ref = median(log2inty_med, na.rm = TRUE))

dta_conflu <- dta_conflu %>% 
    left_join(med_run) %>% 
    mutate(log2inty_eqmed = log2inty - log2inty_med + log2inty_ref) %>% 
    select(-log2inty_med, -log2inty_ref)


## visualization ----------------------------------------------------
## QC plot (boxplot)
ggplot(dta_conflu, aes(run_cbt, log2inty)) + 
    geom_boxplot(aes(fill = biorep)) + 
    facet_wrap(~ is_mod) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(dta_conflu, aes(run_cbt, log2inty_eqmed)) + 
    geom_boxplot(aes(fill = biorep)) + 
    facet_wrap(~ is_mod) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


## proportion of modified peptide features
dta_conflu %>% 
    group_by(cond, biorep, run_cbt) %>% 
    summarise(p_mod = sum(is_mod) / n()) %>% 
    ungroup() %>% 
    ggplot(aes(cond, p_mod)) + 
    geom_jitter(aes(colour = biorep), width = 0.1, size = 4) + 
    labs(x = "Condition", y = "Proportion of modified features")

## number of features
dta_conflu %>% 
    group_by(cond, biorep, run_cbt) %>% 
    summarise(nb_feature = n()) %>% 
    ungroup() %>% 
    ggplot(aes(cond, nb_feature)) + 
    geom_jitter(aes(colour = biorep), width = 0.1, size = 4) + 
    labs(x = "Condition", y = "Number of features")

## number of modified/unmodified features
dta_conflu %>% 
    group_by(cond, biorep, run_cbt) %>% 
    summarise(modified = sum(is_mod), unmodified = sum(!is_mod)) %>% 
    ungroup() %>% 
    gather(mod, nb_feature, modified:unmodified) %>% 
    ggplot(aes(cond, nb_feature)) + 
    geom_jitter(aes(colour = biorep, shape = mod), width = 0.1, size = 4) + 
    labs(x = "Condition", y = "Number of features")

# ## number of features for modified peptides
# dta_conflu %>% 
#     group_by(cond, biorep, run_cbt) %>% 
#     summarise(nb_feature = sum(is_mod)) %>% 
#     ungroup() %>% 
#     ggplot(aes(cond, nb_feature)) + 
#     geom_jitter(aes(colour = biorep), width = 0.1, size = 4)
# 
# ## number of features for unmodified peptides
# dta_conflu %>% 
#     group_by(cond, biorep, run_cbt) %>% 
#     summarise(nb_feature = sum(!is_mod)) %>% 
#     ungroup() %>% 
#     ggplot(aes(cond, nb_feature)) + 
#     geom_jitter(aes(colour = biorep), width = 0.1, size = 4)


## functions for profile plots --------------------------------------
plot_profile_eqmed <- function(df, protein, runs) {
    df %>% filter(grepl(pp, Reference)) %>% 
        select(Reference, feature, run_cbt, log2inty_eqmed) %>% 
        complete(Reference, run_cbt = levels_run, feature) %>% 
        ggplot(aes(run_cbt, log2inty_eqmed, group = feature, colour = feature)) + 
        geom_point(size = 4) +
        geom_line() + 
        coord_cartesian(ylim = c(10,35)) + 
        labs(x = "runs", title = protein) + 
        theme(legend.position = "bottom") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

plot_profile_raw <- function(df, protein, runs) {
    df %>% filter(grepl(pp, Reference)) %>% 
        select(Reference, feature, run_cbt, log2inty) %>% 
        complete(Reference, run_cbt = levels_run, feature) %>% 
        ggplot(aes(run_cbt, log2inty, group = feature, colour = feature)) + 
        geom_point(size = 4) +
        geom_line() + 
        coord_cartesian(ylim = c(10,35)) + 
        labs(x = "runs", title = protein) + 
        theme(legend.position = "bottom") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

plot_profiles <- function(df, protein, runs) {
    df %>% filter(grepl(pp, Reference)) %>% 
        select(Reference, feature, run_cbt, log2inty, log2inty_eqmed) %>% 
        complete(Reference, run_cbt = levels_run, feature) %>% 
        gather("norm", "log2inty", log2inty:log2inty_eqmed) %>% 
        ggplot(aes(run_cbt, log2inty, group = feature, colour = feature)) + 
        geom_point(size = 4) +
        geom_line() + 
        facet_wrap(~ norm) + 
        coord_cartesian(ylim = c(10,35)) + 
        labs(x = "runs", title = protein) + 
        theme(legend.position = "bottom") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

plot_profiles_pdf <- function(df, protein, runs) {
    df %>% filter(grepl(pp, Reference)) %>% 
        select(Reference, feature, run_cbt, log2inty, log2inty_eqmed) %>% 
        complete(Reference, run_cbt = levels_run, feature) %>% 
        gather("norm", "log2inty", log2inty:log2inty_eqmed) %>% 
        ggplot(aes(run_cbt, log2inty, group = feature, colour = feature)) + 
        geom_point(size = 2) +
        geom_line() + 
        facet_wrap(~ norm) + 
        coord_cartesian(ylim = c(10,35)) + 
        labs(x = "runs", title = protein) + 
        theme(legend.position = "none")
}

plot_mprofiles <- function(df, protein, runs) {
    df %>% filter(grepl(pp, Reference)) %>% 
        select(Reference, feature, run_cbt, log2inty, log2inty_eqmed) %>% 
        complete(Reference, run_cbt = levels_run, feature) %>% 
        gather("norm", "log2inty", log2inty:log2inty_eqmed) %>% 
        mutate(is_mod = ifelse(grepl("\\*", feature), "modified", "unmodified")) %>% 
        ggplot(aes(run_cbt, log2inty, group = feature, colour = feature)) + 
        geom_point(size = 4) +
        geom_line() + 
        facet_grid(is_mod ~ norm) + 
        coord_cartesian(ylim = c(10,35)) + 
        labs(x = "runs", title = protein) + 
        # theme(legend.position = "bottom") + 
        theme(legend.position = "none") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


## some profile plots -----------------------------------------------
levels_run <- unique(dta_conflu$run_cbt)

## full sets
dta_conflu_full <- dta_conflu %>% 
    group_by(feature) %>% mutate(nb_run = n()) %>% ungroup() %>% 
    filter(nb_run == 16)

pp <- "UBP30_HUMAN"  # USP30
pp <- "Q70CQ3"
plot_profiles(dta_conflu, pp, levels_run)

## not detected...
pp <- "OPA1_HUMAN"  # OPA 1
plot_profiles(dta_conflu, pp, levels_run)

pp <- "MFN1_HUMAN"  # MFN1
pp <- "Q8IWA4"
plot_profiles(dta_conflu, pp, levels_run)

pp <- "PRKN2_HUMAN"  # Parkin
plot_profiles(dta_conflu, pp, levels_run)
plot_profiles(dta_conflu_full, pp, levels_run)
plot_mprofiles(dta_conflu, pp, levels_run)

pp <- "TOM20_HUMAN"  # Tomm20
pp <- "Q15388"
plot_profiles(dta_conflu, pp, levels_run)

pp <- "VDAC3_HUMAN"  # VDAC
plot_profiles(dta_conflu, pp, levels_run)
plot_profiles(dta_conflu_full, pp, levels_run)
plot_mprofiles(dta_conflu, pp, levels_run)

pp <- "VDAC1_HUMAN"  # VDAC
plot_profiles(dta_conflu, pp, levels_run)
plot_mprofiles(dta_conflu, pp, levels_run)

pp <- "VDAC2_HUMAN"  # VDAC
plot_profiles(dta_conflu, pp, levels_run)
plot_mprofiles(dta_conflu, pp, levels_run)

pp <- "G3P_HUMAN"  # GAPDH
plot_profiles(dta_conflu, pp, levels_run)
plot_profiles(dta_conflu_full, pp, levels_run)
plot_mprofiles(dta_conflu, pp, levels_run)

pp <- "ACTB_HUMAN"  # Actin
plot_profiles(dta_conflu, pp, levels_run)
plot_mprofiles(dta_conflu, pp, levels_run)

pp <- "RL40_HUMAN"  # ubiquitin
plot_profiles(dta_conflu, pp, levels_run)
plot_mprofiles(dta_conflu, pp, levels_run)


# plot_profile_eqmed(dta_conflu, pp, levels_run)
# plot_profiles(dta_conflu, pp, levels_run)


## site/residue centered representation -----------------------------

# mod_idx <- lapply(ref_table$sequence, function(x) as.vector(str_locate_all(x, pattern = mod_residue)[[1]][, 1]))
# mod_res <- lapply(ref_table$sequence, function(x) c(str_match_all(x, pattern = mod_residue)[[1]]))
# nb_site <- str_count(ref_table$sequence, mod_residue)
# protein_indices <- data_frame(uniprot_ac = rep(ref_table$uniprot_ac, nb_mods), 
#                               ptm_site = unlist(mod_residues), 
#                               res_index = unlist(mod_indices))

load("~/Projects/NEU/ptm/dev/annot_uniprot.RData")

pep_n_prot <- dta_conflu %>% distinct(peptide_sequence, uniprot_ac)

regex_aa <- "A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y"
regex_aamod <- "A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y|\\*"

pep_clean <- lapply(pep_n_prot$peptide_sequence, 
                    function(x) paste(unlist(str_extract_all(x, pattern = regex_aamod)), collapse = ""))
pep_unmod <- lapply(pep_clean, function(x) paste(unlist(str_extract_all(x, pattern = regex_aa)), collapse = ""))
nb_mods <- lapply(pep_clean, function(x) length(c(str_match_all(x, "K\\*")[[1]])))





mods_in_pep %>% group_by(uniprot_ac, peptide_seq_clean) %>% 
    mutate(nb_pep = n()) %>% 
    ungroup() %>% 
    filter(nb_pep != 1)

