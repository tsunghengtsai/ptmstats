library(readr)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)

# Prepare data set --------------------------------------------------------
# Load experimental design with data filenames
df_design <- read.table("data/CCCP_confluency/CellDensity_tft_forTH.txt", 
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)

file_name <- paste0("vistaga_", df_design$run_id, ".txt")

sub_cols <- c("Reference", "run_id", "peptide_id", "vista_peak_area_light", 
              "vista_confidence_score", "peptide_m_z", "peptide_sequence", 
              "peptide_trypticity", "peptide_miscleavages", "peptide_validity", 
              "label1", "ms2_rt", "ms2_charge", "is_xq")

# Prepare data set with columns of interest for later use
dta_conflu <- vector("list", length(file_name))
for (i in seq_along(file_name)) {
    file_path <- paste0("data/CCCP_confluency/", file_name[i])
    dta_conflu[[i]] <- read.table(file_path, header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE) %>% 
        tbl_df() %>% select_(.dots = sub_cols)
}
dta_conflu <- bind_rows(dta_conflu)
rm(list = c("i", "file_name", "file_path"))

# save.image("output/dta_conflu.RData")


# Load data set & housekeeping --------------------------------------------
load("output/dta_conflu.RData")

# Confidence threshold
dta_conflu <- dta_conflu %>% filter(vista_confidence_score >= 83)

# Experimental design 
# df_design %>% arrange(id_subject, ann_objective_id, id_injectionset) %>%
#     select(-ann_label)

df_design <- tbl_df(df_design)

dta_conflu <- dta_conflu %>% 
    select(-vista_confidence_score, -peptide_trypticity, -peptide_miscleavages, 
           -peptide_validity)

# OBJ# can be viewed as a batch identity
uniq_subj <- unique(df_design$id_subject)
key_subj <- paste0("S", sprintf("%02d", 1:length(uniq_subj)))
df_design <- df_design %>% 
    mutate(subj = plyr::mapvalues(id_subject, from = uniq_subj, to = key_subj), 
           biorep = ifelse(ann_objective_id == 37122, "B1", "B2"), 
           techrep = paste0("T", id_injectionset), 
           run_bt = paste0(biorep, techrep), 
           run_cbt = paste(id_subject, run_bt, sep = "-"))


# Annotation & data filtering ---------------------------------------------
# [TODO]: need a closer look at `Reference` & `label1`
dta_conflu <- dta_conflu %>% 
    mutate(is_mod = grepl("\\*", peptide_sequence), 
           feature = paste(peptide_sequence, ms2_charge, sep = "_"), 
           log2inty = log2(vista_peak_area_light), 
           cond = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$id_subject), 
           biorep = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$biorep), 
           techrep = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$techrep), 
           run_cbt = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$run_cbt))

# Extract original peptide sequence & other information
df_pep <- dta_conflu %>% 
    distinct(peptide_sequence, is_mod) %>% 
    mutate(unmod_peptide = str_replace_all(peptide_sequence, "\\*", ""), 
           trim_peptide = str_match(unmod_peptide, regex("\\.[ACDEFGHIKLMNPQRSTVWY]+")) %>% 
               str_replace("\\.", ""), 
           nb_aa = str_length(trim_peptide))

# Make sure each peptide was extracted correctly
# df_pep %>% filter(is.na(nb_aa) | nb_aa == 0)

# Retain peptides with length >= 6
df_pep <- df_pep %>% filter(nb_aa >= 6)

# Determine the pairing status of peptide
df_pep <- df_pep %>% group_by(unmod_peptide) %>% 
    mutate(is_paired = ifelse(n_distinct(is_mod) == 2, TRUE, FALSE)) %>% 
    ungroup()

# Features with no charge state or retention time inforamtion
# dta_conflu %>% filter(is.na(ms2_charge) | is.na(ms2_rt)) %>%
#     select(feature, peptide_sequence, ms2_charge, ms2_rt)

# Short matches (peptide length less than 6 AAs)
# dta_conflu %>% anti_join(df_pep) %>% distinct(peptide_sequence, feature)

# Unexpected uniprot ID (non-human proteins in this case)
# regex_uniprot <- ".*([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}).*"
# dta_conflu %>% filter(!grepl(pattern = regex_uniprot, Reference))

# Initial filtering
regex_uniprot <- ".*([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}).*"
dta_conflu <- dta_conflu %>% 
    semi_join(df_pep) %>% 
    filter(!is.na(ms2_charge), !is.na(ms2_rt), grepl(regex_uniprot, Reference))

# Extract uniprot accession number
uniq_ref <- dta_conflu %>% distinct(Reference) %>% .$Reference
uniq_unip <- unlist(lapply(uniq_ref, function(x) str_extract(unlist(str_split(x, '\\|'))[2], regex_uniprot)))
dta_conflu <- dta_conflu %>% 
    mutate(uniprot_ac = plyr::mapvalues(Reference, from = uniq_ref, to = uniq_unip))


# Data exploration --------------------------------------------------------
# Check for duplicate features (4971 in total)
dta_conflu %>% 
    count(run_id, feature) %>% ungroup() %>% 
    filter(n > 1)

# 2151 unique features with duplicates
dta_conflu %>% 
    count(run_id, feature) %>% ungroup() %>% 
    filter(n > 1) %>% 
    distinct(feature)

# 813 with additional criterion for max RT diff > 1.5 min
dta_conflu %>% 
    group_by(run_id, feature) %>% 
    summarise(n = n(), d_rt = max(ms2_rt, na.rm = TRUE) - min(ms2_rt, na.rm = TRUE)) %>% 
    ungroup() %>% 
    filter(n > 1, d_rt > 1.5)

# Protein (Reference)
dta_conflu %>% select(Reference, peptide_sequence) %>% distinct() %>% 
    group_by(peptide_sequence) %>% mutate(nb_prot = n()) %>% ungroup() %>% 
    filter(nb_prot > 1) %>% arrange(peptide_sequence)

# Uniprot ID ([TODO]: address the aliasing pairs)
dta_conflu %>% select(uniprot_ac, peptide_sequence) %>% distinct() %>% 
    group_by(peptide_sequence) %>% mutate(nb_prot = n()) %>% ungroup() %>% 
    filter(nb_prot > 1) %>% arrange(peptide_sequence)

# 592 peptides have modified/unmodified counterparts
df_pep %>% filter(is_paired, !is_mod) %>% distinct(unmod_peptide)


# Data manipulation & transformation --------------------------------------
# Use the max for duplicate features (per run)
# [TODO]: conisder is_xq, ms2_rt
dta_conflu <- dta_conflu %>% 
    group_by(Reference, uniprot_ac, peptide_sequence, is_mod, feature, ms2_charge, cond, biorep, run_cbt, run_id) %>% 
    summarise(log2inty = max(log2inty)) %>% 
    ungroup()

# Remove peptides with ambiguous matches (to >1 uniprot_ac) for now
# pep2uniprot <- dta_conflu %>% 
#     distinct(peptide_sequence, uniprot_ac) %>% 
#     count(peptide_sequence)

# Remove peptides with ambiguous matches (to >1 Reference): more stringent than 
# mapping with uniprot_ac, as multiple Reference may map to the same uniprot_ac
# Example: dta_conflu %>% filter(feature == "K.AFSLK*TSTSAVR.H_3", run_cbt == "CCCP_Light-B1T1")
pep2uniprot <- dta_conflu %>% 
    distinct(peptide_sequence, Reference) %>% 
    count(peptide_sequence)

dta_conflu <- dta_conflu %>% 
    semi_join(filter(pep2uniprot, n == 1))

# Status of peptide pairing & protein normalization
dta_conflu <- dta_conflu %>% 
    left_join(select(df_pep, peptide_sequence, unmod_peptide, is_paired)) %>% 
    group_by(uniprot_ac) %>% 
    mutate(is_ref = any(!is_paired & !is_mod)) %>% ungroup()

# More than 70% of the modified peptides might be eligible for protein normalization
dta_conflu %>% distinct(peptide_sequence, is_mod, is_ref) %>% count(is_mod, is_ref)

# Subset with features observed in both batches
ftr_2b <- dta_conflu %>% 
    distinct(feature, biorep) %>% group_by(feature) %>% 
    filter(all(c("B1", "B2") %in% biorep)) %>% ungroup() %>% 
    distinct(feature)

# Subset with feature consistently absent/present in both batches
ftr_2bc <- dta_conflu %>% 
    distinct(feature, cond, biorep) %>% 
    group_by(feature, cond) %>% 
    mutate(match = all(c("B1", "B2") %in% biorep)) %>% 
    group_by(feature) %>% 
    filter(all(match)) %>% ungroup() %>% 
    distinct(feature)

dta_conflu_2b <- dta_conflu %>% semi_join(ftr_2b)
dta_conflu_2bc <- dta_conflu %>% semi_join(ftr_2bc)

# Normalization (equalizing median) - not appropriate
# med_run <- dta_conflu %>% 
#     group_by(run_id) %>% 
#     summarise(log2inty_med = median(log2inty, na.rm = TRUE)) %>% 
#     mutate(log2inty_ref = median(log2inty_med, na.rm = TRUE))
# 
# dta_conflu <- dta_conflu %>% 
#     left_join(med_run) %>% 
#     mutate(log2inty_eqmed = log2inty - log2inty_med + log2inty_ref) %>% 
#     select(-log2inty_med, -log2inty_ref)


# Visualization -----------------------------------------------------------
# QC plot (boxplot)
dta_conflu %>% 
    mutate(is_mod_fac = factor(ifelse(is_mod, "Modified", "Unmodified"))) %>% 
    ggplot(aes(run_cbt, log2inty)) + 
    geom_boxplot(aes(fill = biorep)) + 
    facet_wrap(~ is_mod_fac) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(x = "Run", y = "Log2-intensity")

dta_conflu_2b %>% 
    mutate(is_mod_fac = factor(ifelse(is_mod, "Modified", "Unmodified"))) %>% 
    ggplot(aes(run_cbt, log2inty)) + 
    geom_boxplot(aes(fill = biorep)) + 
    facet_wrap(~ is_mod_fac) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(x = "Run", y = "Log2-intensity")

dta_conflu_2bc %>% 
    mutate(is_mod_fac = factor(ifelse(is_mod, "Modified", "Unmodified"))) %>% 
    ggplot(aes(run_cbt, log2inty)) + 
    geom_boxplot(aes(fill = biorep)) + 
    facet_wrap(~ is_mod_fac) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(x = "Run", y = "Log2-intensity")

# ggplot(dta_conflu, aes(run_cbt, log2inty_eqmed)) + 
#     geom_boxplot(aes(fill = biorep)) + 
#     facet_wrap(~ is_mod) + 
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Number of modified/unmodified features
# This suggests a discrepancy in the enrichment step across batches
dta_conflu %>% 
    group_by(cond, biorep, run_cbt) %>% 
    summarise(Modified = sum(is_mod), Unmodified = sum(!is_mod)) %>% 
    ungroup() %>% 
    gather(mod, nb_feature, Modified:Unmodified) %>% 
    ggplot(aes(cond, nb_feature)) + 
    geom_jitter(aes(colour = biorep, shape = mod), width = 0.1, size = 4) + 
    labs(x = "Condition", y = "Number of features")

# Number of modified/unmodified peptides
dta_conflu %>% 
    distinct(cond, biorep, run_cbt, peptide_sequence, is_mod) %>% 
    group_by(cond, biorep, run_cbt) %>% 
    summarise(Modified = sum(is_mod), Unmodified = sum(!is_mod)) %>% 
    ungroup() %>% 
    gather(mod, nb_peptide, Modified:Unmodified) %>% 
    ggplot(aes(cond, nb_peptide)) + 
    geom_jitter(aes(colour = biorep, shape = mod), width = 0.1, size = 4) + 
    labs(x = "Condition", y = "Number of peptides")

# Proportion of modified peptide features
# dta_conflu %>% 
#     group_by(cond, biorep, run_cbt) %>% 
#     summarise(p_mod = sum(is_mod) / n()) %>% 
#     ungroup() %>% 
#     ggplot(aes(cond, p_mod)) + 
#     geom_jitter(aes(colour = biorep), width = 0.1, size = 4) + 
#     labs(x = "Condition", y = "Proportion of modified features")

# Number of features
# dta_conflu %>% 
#     group_by(cond, biorep, run_cbt) %>% 
#     summarise(nb_feature = n()) %>% 
#     ungroup() %>% 
#     ggplot(aes(cond, nb_feature)) + 
#     geom_jitter(aes(colour = biorep), width = 0.1, size = 4) + 
#     labs(x = "Condition", y = "Number of features")


# Functions for plotting profiles -----------------------------------------
# Plot feature profiles 
plot_profile <- function(df_allprot, protein, run_level) {
    df_prot <- df_allprot %>% filter(uniprot_ac == protein)
    # Complete possible combinations of peptide features and runs
    mpar <- df_prot %>% distinct(peptide_sequence, feature, is_mod, is_paired)
    df_fill <- df_prot %>% 
        select(feature, run_cbt, log2inty) %>% 
        complete(run_cbt = run_level, feature) %>% 
        left_join(mpar) %>% 
        mutate(is_mod_fac = factor(ifelse(is_mod, "Modified", "Unmodified")), 
               is_par_fac = factor(ifelse(is_paired, "With counterpart", "No counterpart")), 
               run_fac = factor(run_cbt, levels = run_level))
    # Feature profiles categorized by modification and matching status
    df_fill %>% 
        ggplot(aes(run_fac, log2inty, group = feature, colour = peptide_sequence)) + 
        geom_point(size = 4) +
        geom_line() + 
        facet_grid(is_mod_fac ~ is_par_fac) + 
        coord_cartesian(ylim = c(10, 35)) + 
        labs(x = "Run", y = "Log2-intensity", title = protein) + 
        theme(legend.position = "none") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

# Summarization of feature intensities of a protein using Tukey's median polish
sumprot_tmp <- function(df_prot) {
    # Required fields of df_prot: feature, run, log2inty
    inty_wide <- df_prot %>% select(feature, run, log2inty) %>% spread(feature, log2inty)
    inty_mat <- data.matrix(inty_wide[, -1])
    mp_out <- medpolish(inty_mat, na.rm = TRUE, trace.iter = FALSE)
    sum_df <- data_frame(run = inty_wide$run, log2inty = mp_out$overall + mp_out$row)
    
    return(sum_df)
}

# Plot feature profiles with summarization
plot_sumprofile <- function(df_allprot, protein, run_level) {
    df_prot <- df_allprot %>% filter(uniprot_ac == protein)
    # Complete possible combinations of features and runs
    mpar <- df_prot %>% distinct(feature, is_mod, is_paired)
    # mpar <- df_prot %>% distinct(peptide_sequence, feature, is_mod, is_paired)
    df_fill <- df_prot %>% 
        select(feature, run_cbt, log2inty) %>% 
        complete(run_cbt = run_level, feature) %>% 
        left_join(mpar) %>% mutate(quan = "Feature")
    # Median polish summarization
    key_mod <- c(TRUE, TRUE, FALSE, FALSE)
    key_par <- c(TRUE, FALSE, TRUE, FALSE)
    mplist <- vector("list", 4)
    for (i in seq_along(key_mod)) {
        df_sub <- df_prot %>% filter(is_mod == key_mod[i], is_paired == key_par[i])
        if (nrow(df_sub) > 0) {
            mplist[[i]] <- sumprot_tmp(select(df_sub, feature, run = run_cbt, log2inty)) %>% 
                complete(run = run_level) %>% 
                mutate(is_mod = key_mod[i], is_paired = key_par[i])
        }
    }
    # Combine data frames for feature and summarization
    df_sum <- bind_rows(mplist) %>% rename(run_cbt = run) %>% 
        # complete(run_cbt = run_level) %>% 
        mutate(feature = "tmp", quan = "Summarization")
    df_sumfill <- bind_rows(df_fill, df_sum) %>% 
        mutate(is_mod_fac = factor(ifelse(is_mod, "Modified", "Unmodified")), 
               is_par_fac = factor(ifelse(is_paired, "With counterpart", "No counterpart")), 
               run_fac = factor(run_cbt, levels = run_level), quan = factor(quan))
    # Profile plot
    df_sumfill %>% 
        ggplot(aes(run_fac, log2inty, group = feature, colour = quan, size = quan)) + 
        geom_point() + geom_line(size = 0.5) + 
        scale_colour_manual(values = c("lightgray", "darkred")) + 
        scale_size_manual(values = c(1, 3)) + 
        facet_grid(is_mod_fac ~ is_par_fac) + 
        coord_cartesian(ylim = c(10, 35)) + 
        labs(x = "Run", y = "Log2-intensity", title = protein) + 
        theme(legend.title = element_blank(), legend.position = "bottom", legend.box = "horizontal") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


# Profile plots -----------------------------------------------------------
runlvl <- dta_conflu %>% distinct(run_cbt) %>% arrange(run_cbt) %>% .$run_cbt
runlvl_bat <- dta_conflu %>% distinct(biorep, run_cbt) %>% arrange(biorep, run_cbt) %>% .$run_cbt

# USP30 "UBP30_HUMAN"  
pp <- "Q70CQ3"
plot_profile(dta_conflu, pp, runlvl)
plot_profile(dta_conflu, pp, runlvl_bat)
plot_sumprofile(dta_conflu, pp, runlvl_bat)

# MFN1 "MFN1_HUMAN"
pp <- "Q8IWA4"
plot_profile(dta_conflu, pp, runlvl)
plot_profile(dta_conflu, pp, runlvl_bat)
plot_sumprofile(dta_conflu, pp, runlvl_bat)

# Parkin "PRKN2_HUMAN"
pp <- "O60260"
plot_profile(dta_conflu, pp, runlvl)
plot_profile(dta_conflu, pp, runlvl_bat)
plot_sumprofile(dta_conflu, pp, runlvl_bat)

# Tomm20 "TOM20_HUMAN"
pp <- "Q15388"
plot_profile(dta_conflu, pp, runlvl)
plot_profile(dta_conflu, pp, runlvl_bat)
plot_sumprofile(dta_conflu, pp, runlvl_bat)

# VDAC "VDAC3_HUMAN"  
pp <- "Q9Y277"
plot_profile(dta_conflu, pp, runlvl)
plot_profile(dta_conflu, pp, runlvl_bat)
plot_sumprofile(dta_conflu, pp, runlvl_bat)
plot_sumprofile(dta_conflu_2bc, pp, runlvl_bat)

# VDAC "VDAC1_HUMAN"  
pp <- "P21796"
plot_profile(dta_conflu, pp, runlvl)
plot_profile(dta_conflu, pp, runlvl_bat)
plot_sumprofile(dta_conflu, pp, runlvl_bat)
plot_sumprofile(dta_conflu_2bc, pp, runlvl_bat)

# VDAC "VDAC2_HUMAN"  
pp <- "P45880"
plot_profile(dta_conflu, pp, runlvl)
plot_profile(dta_conflu, pp, runlvl_bat)
plot_sumprofile(dta_conflu, pp, runlvl_bat)
plot_sumprofile(dta_conflu_2bc, pp, runlvl_bat)

# GAPDH "G3P_HUMAN"  
pp <- "P04406"
plot_profile(dta_conflu, pp, runlvl)
plot_profile(dta_conflu, pp, runlvl_bat)
plot_sumprofile(dta_conflu, pp, runlvl_bat)
plot_sumprofile(dta_conflu_2bc, pp, runlvl_bat)

# Actin "ACTB_HUMAN"  
pp <- "P60709"
plot_profile(dta_conflu, pp, runlvl)
plot_profile(dta_conflu, pp, runlvl_bat)
plot_sumprofile(dta_conflu, pp, runlvl_bat)
plot_sumprofile(dta_conflu_2bc, pp, runlvl_bat)

# Ubiquitin "RL40_HUMAN" 
pp <- "P62987"
plot_profile(dta_conflu, pp, runlvl)
plot_profile(dta_conflu, pp, runlvl_bat)
plot_sumprofile(dta_conflu, pp, runlvl_bat)
plot_sumprofile(dta_conflu_2bc, pp, runlvl_bat)


## (old) functions for profile plots --------------------------------------
# plot_profiles_pdf <- function(df, protein, runs) {
#     df %>% filter(grepl(pp, Reference)) %>% 
#         select(Reference, feature, run_cbt, log2inty, log2inty_eqmed) %>% 
#         complete(Reference, run_cbt = levels_run, feature) %>% 
#         gather("norm", "log2inty", log2inty:log2inty_eqmed) %>% 
#         ggplot(aes(run_cbt, log2inty, group = feature, colour = feature)) + 
#         geom_point(size = 2) +
#         geom_line() + 
#         facet_wrap(~ norm) + 
#         coord_cartesian(ylim = c(10,35)) + 
#         labs(x = "runs", title = protein) + 
#         theme(legend.position = "none")
# }


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

