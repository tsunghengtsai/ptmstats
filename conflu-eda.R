library(tidyverse)
library(stringr)
library(broom)

# Prepare data set --------------------------------------------------------
# # Load experimental design with data filenames
# df_design <- read.table("data/CCCP_confluency/CellDensity_tft_forTH.txt", 
#                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# 
# file_name <- paste0("vistaga_", df_design$run_id, ".txt")
# 
# sub_cols <- c("Reference", "run_id", "peptide_id", "vista_peak_area_light", 
#               "vista_confidence_score", "peptide_m_z", "peptide_sequence", 
#               "peptide_trypticity", "peptide_miscleavages", "peptide_validity", 
#               "label1", "ms2_rt", "ms2_charge", "is_xq")
# 
# # Prepare data set with columns of interest for later use
# dta_conflu <- vector("list", length(file_name))
# for (i in seq_along(file_name)) {
#     file_path <- paste0("data/CCCP_confluency/", file_name[i])
#     dta_conflu[[i]] <- read.table(file_path, header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE) %>% 
#         tbl_df() %>% select_(.dots = sub_cols)
# }
# dta_conflu <- bind_rows(dta_conflu)
# rm(list = c("i", "file_name", "file_path"))

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
subj_uniq <- unique(df_design$id_subject)
subj_key <- str_c("S", str_pad(seq_along(subj_uniq), width = 2, pad = "0"))
df_design <- df_design %>% 
    mutate(subj = plyr::mapvalues(id_subject, from = subj_uniq, to = subj_key), 
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
           run = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$run_cbt))

# Extract original peptide sequence & other information
df_pep <- dta_conflu %>% 
    distinct(peptide_sequence, is_mod) %>% 
    mutate(unmod_peptide = str_replace_all(peptide_sequence, "\\*", ""), 
           trimmed_pep = str_extract(unmod_peptide, "(?<=\\.)([ACDEFGHIKLMNPQRSTVWY\\*]+)"), 
           trimmed_mod = str_extract(peptide_sequence, "(?<=\\.)([ACDEFGHIKLMNPQRSTVWY\\*]+)"), 
           nb_aa = str_length(trimmed_pep))

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
dta_conflu <- dta_conflu %>% 
    separate(Reference, into = c("p_entry", "uniprot_iso"), sep = "\\|", remove = FALSE)


# Data exploration --------------------------------------------------------
# Check for duplicate features (4971 in all runs)
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

# Mapping peptides to proteins (`Reference` column)
dta_conflu %>% select(Reference, peptide_sequence) %>% distinct() %>% 
    group_by(peptide_sequence) %>% mutate(nb_prot = n()) %>% ungroup() %>% 
    filter(nb_prot > 1) %>% arrange(peptide_sequence)

# Uniprot ID ([TODO]: address the aliasing pairs, search, e.g., P08107)
dta_conflu %>% select(uniprot_iso, peptide_sequence) %>% distinct() %>% 
    group_by(peptide_sequence) %>% mutate(nb_prot = n()) %>% ungroup() %>% 
    filter(nb_prot > 1) %>% arrange(peptide_sequence)

# Among 10696 modified peptides
df_pep %>% filter(is_mod) %>% distinct(peptide_sequence)

# 705 modified peptides have 592 unmodified counterparts
df_pep %>% filter(is_paired, is_mod) %>% distinct(peptide_sequence)
df_pep %>% filter(is_paired) %>% distinct(unmod_peptide)


# Data manipulation & transformation --------------------------------------
# Use the max for duplicate features (per run)
# [TODO]: conisder is_xq, ms2_rt
dta_conflu <- dta_conflu %>% 
    group_by(Reference, uniprot_iso, peptide_sequence, is_mod, feature, cond, biorep, run) %>% 
    summarise(log2inty = max(log2inty)) %>% 
    ungroup()

# Remove peptides with ambiguous matches (to >1 uniprot_iso) for now
# pep2uniprot <- dta_conflu %>% 
#     distinct(peptide_sequence, uniprot_iso) %>% 
#     count(peptide_sequence)

# Remove peptides with ambiguous matches (to >1 Reference): more stringent than 
# mapping with uniprot_iso, as multiple Reference may map to the same uniprot_iso
# Example: dta_conflu %>% filter(feature == "K.AFSLK*TSTSAVR.H_3", run == "CCCP_Light-B1T1")
pep2uniprot <- dta_conflu %>% 
    distinct(peptide_sequence, Reference) %>% 
    count(peptide_sequence)

dta_conflu <- dta_conflu %>% 
    semi_join(filter(pep2uniprot, n == 1))

# Status of peptide pairing & protein normalization
# [TODO]: check the definition of is_ref
dta_conflu <- dta_conflu %>% 
    left_join(select(df_pep, peptide_sequence, unmod_peptide, is_paired)) %>% 
    group_by(uniprot_iso) %>% 
    mutate(is_ref = any(!is_paired & !is_mod)) %>% ungroup()

# More than 70% of the modified peptides might be eligible for protein normalization
# dta_conflu %>% distinct(peptide_sequence, is_mod, is_ref) %>% count(is_mod, is_ref)

# Subset with features observed in both batches
# ftr_2b <- dta_conflu %>% 
#     distinct(feature, biorep) %>% group_by(feature) %>% 
#     filter(all(c("B1", "B2") %in% biorep)) %>% ungroup() %>% 
#     distinct(feature)

# Subset with feature consistently absent/present in both batches
# ftr_2bc <- dta_conflu %>% 
#     distinct(feature, cond, biorep) %>% 
#     group_by(feature, cond) %>% 
#     mutate(match = all(c("B1", "B2") %in% biorep)) %>% 
#     group_by(feature) %>% 
#     filter(all(match)) %>% ungroup() %>% 
#     distinct(feature)

# dta_conflu_2b <- dta_conflu %>% semi_join(ftr_2b)
# dta_conflu_2bc <- dta_conflu %>% semi_join(ftr_2bc)

# Original dataset: 2954 proteins, 21419 features
# dta_conflu %>% distinct(uniprot_iso) %>% nrow()
# dta_conflu %>% distinct(feature) %>% nrow()

# Filtered dataset: 1292 proteins, 4440 features
# dta_conflu_2bc %>% distinct(uniprot_iso) %>% nrow()
# dta_conflu_2bc %>% distinct(feature) %>% nrow()


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
    ggplot(aes(run, log2inty)) + 
    geom_boxplot(aes(fill = biorep)) + 
    facet_wrap(~ is_mod_fac) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(x = "Run", y = "Log2-intensity")

dta_conflu_2b %>% 
    mutate(is_mod_fac = factor(ifelse(is_mod, "Modified", "Unmodified"))) %>% 
    ggplot(aes(run, log2inty)) + 
    geom_boxplot(aes(fill = biorep)) + 
    facet_wrap(~ is_mod_fac) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(x = "Run", y = "Log2-intensity")

dta_conflu_2bc %>% 
    mutate(is_mod_fac = factor(ifelse(is_mod, "Modified", "Unmodified"))) %>% 
    ggplot(aes(run, log2inty)) + 
    geom_boxplot(aes(fill = biorep)) + 
    facet_wrap(~ is_mod_fac) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(x = "Run", y = "Log2-intensity")

# Number of modified/unmodified features
# This suggests a discrepancy in the enrichment step across batches
dta_conflu %>% 
    group_by(cond, biorep, run) %>% 
    summarise(Modified = sum(is_mod), Unmodified = sum(!is_mod)) %>% 
    ungroup() %>% 
    gather(mod, nb_feature, Modified:Unmodified) %>% 
    ggplot(aes(cond, nb_feature)) + 
    geom_jitter(aes(colour = biorep, shape = mod), width = 0.1, size = 4) + 
    labs(x = "Condition", y = "Number of features")

# Number of modified/unmodified peptides
dta_conflu %>% 
    distinct(cond, biorep, run, peptide_sequence, is_mod) %>% 
    group_by(cond, biorep, run) %>% 
    summarise(Modified = sum(is_mod), Unmodified = sum(!is_mod)) %>% 
    ungroup() %>% 
    gather(mod, nb_peptide, Modified:Unmodified) %>% 
    ggplot(aes(cond, nb_peptide)) + 
    geom_jitter(aes(colour = biorep, shape = mod), width = 0.1, size = 4) + 
    labs(x = "Condition", y = "Number of peptides")

# Proportion of modified peptide features
# dta_conflu %>% 
#     group_by(cond, biorep, run) %>% 
#     summarise(p_mod = sum(is_mod) / n()) %>% 
#     ungroup() %>% 
#     ggplot(aes(cond, p_mod)) + 
#     geom_jitter(aes(colour = biorep), width = 0.1, size = 4) + 
#     labs(x = "Condition", y = "Proportion of modified features")

# Number of features
# dta_conflu %>% 
#     group_by(cond, biorep, run) %>% 
#     summarise(nb_feature = n()) %>% 
#     ungroup() %>% 
#     ggplot(aes(cond, nb_feature)) + 
#     geom_jitter(aes(colour = biorep), width = 0.1, size = 4) + 
#     labs(x = "Condition", y = "Number of features")


# Functions for plotting profiles -----------------------------------------
# [TODO]: should use one single function with options for pairing & summarization
# Plot feature profiles 
plot_profile <- function(df_allprot, protein, run_level) {
    df_prot <- df_allprot %>% filter(uniprot_iso == protein)
    # Complete possible combinations of peptide features and runs
    mpar <- df_prot %>% distinct(unmod_peptide, feature, is_mod, is_paired)
    df_fill <- df_prot %>% 
        select(feature, run, log2inty) %>% 
        complete(run = run_level, feature) %>% 
        left_join(mpar) %>% 
        mutate(is_mod_fac = factor(ifelse(is_mod, "Modified", "Unmodified")), 
               is_par_fac = factor(ifelse(is_paired, "With counterpart", "No counterpart")), 
               run_fac = factor(run, levels = run_level))
    # Feature profiles categorized by modification and matching status
    df_fill %>% 
        ggplot(aes(run_fac, log2inty, group = feature, colour = unmod_peptide)) + 
        geom_point(size = 4) +
        geom_line() + 
        facet_grid(is_mod_fac ~ is_par_fac) + 
        coord_cartesian(ylim = c(10, 35)) + 
        labs(x = "Run", y = "Log2-intensity", title = protein) + 
        theme(legend.position = "none") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

# Plot feature profiles with no counterpart
plot_nc_profile <- function(df_allprot, protein, run_level) {
    df_prot <- df_allprot %>% filter(uniprot_iso == protein) %>% 
        filter(!is_paired)
    # Complete possible combinations of peptide features and runs
    mpar <- df_prot %>% 
        # distinct(unmod_peptide, feature, is_mod, is_paired)
        distinct(unmod_peptide, feature, is_mod)
    df_fill <- df_prot %>% 
        select(feature, run, log2inty) %>% 
        complete(run = run_level, feature) %>% 
        left_join(mpar) %>% 
        mutate(is_mod_fac = factor(ifelse(is_mod, "Modified", "Unmodified")), 
               # is_par_fac = factor(ifelse(is_paired, "With counterpart", "No counterpart")), 
               run_fac = factor(run, levels = run_level))
    # Feature profiles categorized by modification and matching status
    df_fill %>% 
        ggplot(aes(run_fac, log2inty, group = feature, colour = unmod_peptide)) + 
        geom_point(size = 4) +
        geom_line() + 
        # facet_grid(is_mod_fac ~ is_par_fac) + 
        facet_grid(is_mod_fac ~ .) + 
        coord_cartesian(ylim = c(10, 35)) + 
        labs(x = "Run", y = "Log2-intensity", title = protein) + 
        theme(legend.position = "none") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


# Plot feature profiles with summarization
plot_sumprofile <- function(df_allprot, protein, run_level) {
    df_prot <- df_allprot %>% filter(uniprot_iso == protein)
    # Complete possible combinations of features and runs
    mpar <- df_prot %>% distinct(feature, is_mod, is_paired)
    df_fill <- df_prot %>% 
        select(feature, run, log2inty) %>% 
        complete(run = run_level, feature) %>% 
        left_join(mpar) %>% mutate(quan = "Feature")
    # Median polish summarization
    key_mod <- c(TRUE, TRUE, FALSE, FALSE)
    key_par <- c(TRUE, FALSE, TRUE, FALSE)
    mplist <- vector("list", 4)
    for (i in seq_along(key_mod)) {
        df_sub <- df_prot %>% filter(is_mod == key_mod[i], is_paired == key_par[i])
        if (nrow(df_sub) > 0) {
            mplist[[i]] <- sumprot_tmp(select(df_sub, feature, run, log2inty)) %>% 
                complete(run = run_level) %>% 
                mutate(is_mod = key_mod[i], is_paired = key_par[i])
        }
    }
    # Combine data frames for feature and summarization
    df_sum <- bind_rows(mplist) %>% mutate(feature = "tmp", quan = "Summarization")
    df_sumfill <- bind_rows(df_fill, df_sum) %>% 
        mutate(is_mod_fac = factor(ifelse(is_mod, "Modified", "Unmodified")), 
               is_par_fac = factor(ifelse(is_paired, "With counterpart", "No counterpart")), 
               run_fac = factor(run, levels = run_level), quan = factor(quan))
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

# Plot feature profiles (no counterpart) with summarization per batch
plot_nc_sumprofile <- function(df_allprot, protein, run_level) {
    df_prot <- df_allprot %>% filter(uniprot_iso == protein) %>% 
        filter(!is_paired)
    # Complete possible combinations of features and runs
    mpar <- df_prot %>% distinct(feature, is_mod)
    df_fill <- df_prot %>% 
        select(feature, run, log2inty) %>% 
        complete(run = run_level, feature) %>% 
        left_join(mpar) %>% mutate(quan = "Feature")
    # Median polish summarization
    key_mod <- c(TRUE, FALSE)
    mplist <- vector("list", 2)
    for (i in seq_along(key_mod)) {
        df_sub <- df_prot %>% 
            filter(is_mod == key_mod[i]) %>% 
            select(feature, run, log2inty, biorep)
        if (nrow(df_sub) > 0) {
            df_sub1 <- df_sub %>% filter(biorep == "B1")
            df_sub2 <- df_sub %>% filter(biorep == "B2")
            if (nrow(df_sub1) > 0) {
                if (nrow(df_sub2) > 0) {
                    df_subtmp <- bind_rows(sumprot_tmp(df_sub1), sumprot_tmp(df_sub2))
                } else {
                    df_subtmp <- sumprot_tmp(df_sub1)
                }
            } else {
                df_subtmp <- sumprot_tmp(df_sub2)
            }
            mplist[[i]] <- df_subtmp %>% 
                complete(run = run_level) %>% 
                mutate(is_mod = key_mod[i])
        }
    }
    # Combine data frames for feature and summarization
    df_sum <- bind_rows(mplist) %>% mutate(feature = "tmp", quan = "Summarization")
    df_sumfill <- bind_rows(df_fill, df_sum) %>% 
        mutate(is_mod_fac = factor(ifelse(is_mod, "Modified", "Unmodified")), 
               run_fac = factor(run, levels = run_level), quan = factor(quan))
    # Profile plot
    df_sumfill %>% 
        ggplot(aes(run_fac, log2inty, group = feature, colour = quan, size = quan)) + 
        geom_point() + geom_line(size = 0.5) + 
        scale_colour_manual(values = c("lightgray", "darkred")) + 
        scale_size_manual(values = c(1, 3)) + 
        facet_grid(is_mod_fac ~ .) + 
        coord_cartesian(ylim = c(10, 35)) + 
        labs(x = "Run", y = "Log2-intensity", title = protein) + 
        theme(legend.title = element_blank(), legend.position = "bottom", legend.box = "horizontal") +
        # theme(legend.position = "none") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


# Plot feature profiles (no counterpart) with 2 ways of summarization
plot_nc_sumprofile2 <- function(df_allprot, protein, run_level) {
    df_prot <- df_allprot %>% filter(uniprot_iso == protein) %>% 
        filter(!is_paired)
    # Complete possible combinations of features and runs
    mpar <- df_prot %>% distinct(feature, is_mod)
    df_fill <- df_prot %>% 
        select(feature, run, log2inty) %>% 
        complete(run = run_level, feature) %>% 
        left_join(mpar) %>% mutate(quan = "Feature")
    # Median polish summarization
    key_mod <- c(TRUE, FALSE)
    mplist <- vector("list", 2)
    mplist2 <- vector("list", 2)
    for (i in seq_along(key_mod)) {
        df_sub <- df_prot %>% 
            filter(is_mod == key_mod[i]) %>% 
            select(feature, run, log2inty, biorep)
        if (nrow(df_sub) > 0) {
            # All batches
            mplist2[[i]] <- sumprot_tmp(df_sub) %>% complete(run = run_level) %>% 
                mutate(is_mod = key_mod[i])
            # Per batch
            df_sub1 <- df_sub %>% filter(biorep == "B1")
            df_sub2 <- df_sub %>% filter(biorep == "B2")
            if (nrow(df_sub1) > 0) {
                if (nrow(df_sub2) > 0) {
                    df_subtmp <- bind_rows(sumprot_tmp(df_sub1), sumprot_tmp(df_sub2))
                } else {
                    df_subtmp <- sumprot_tmp(df_sub1)
                }
            } else {
                df_subtmp <- sumprot_tmp(df_sub2)
            }
            mplist[[i]] <- df_subtmp %>% complete(run = run_level) %>% 
                mutate(is_mod = key_mod[i])
        }
    }
    # Combine data frames for feature and summarization
    df_sum2 <- bind_rows(mplist2) %>% mutate(feature = "tmp2", quan = "Summarization (across batch)")
    df_sum <- bind_rows(mplist) %>% mutate(feature = "tmp", quan = "Summarization (per batch)")
    df_sumfill <- bind_rows(df_fill, df_sum, df_sum2) %>% 
        mutate(is_mod_fac = factor(ifelse(is_mod, "Modified", "Unmodified")), 
               run_fac = factor(run, levels = run_level), quan = factor(quan))
    # Profile plot
    df_sumfill %>% 
        ggplot(aes(run_fac, log2inty, group = feature, colour = quan, size = quan)) + 
        geom_point() + geom_line(size = 0.5) + 
        scale_colour_manual(values = c("lightgray", "deepskyblue4", "darkred")) + 
        scale_size_manual(values = c(1, 2, 3)) + 
        facet_grid(is_mod_fac ~ .) + 
        coord_cartesian(ylim = c(10, 35)) + 
        labs(x = "Run", y = "Log2-intensity", title = protein) + 
        theme(legend.title = element_blank(), legend.position = "bottom", legend.box = "horizontal") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


# TMP summarization functions ---------------------------------------------
# Summarization of feature intensities of a protein using Tukey's median polish
sumprot_tmp <- function(df_prot) {
    # Required fields of df_prot: feature, run, log2inty
    inty_wide <- df_prot %>% select(feature, run, log2inty) %>% spread(feature, log2inty)
    inty_mat <- data.matrix(inty_wide[, -1])
    mp_out <- medpolish(inty_mat, na.rm = TRUE, trace.iter = FALSE)
    sum_df <- data_frame(run = inty_wide$run, log2inty = mp_out$overall + mp_out$row)
    
    return(sum_df)
}


# Run-level summarization (within batch)
sum_feature_bch <- function(df_allprot, proteins) {
    # Housekeeping
    proteins <- unique(proteins)
    df_allprot <- df_allprot %>% filter(uniprot_iso %in% proteins)
    # Nested data frame (protein, mod, paired, batch)
    nest_allprot <- df_allprot %>% 
        group_by(uniprot_iso, is_mod, is_paired, biorep) %>% 
        nest()
    # TMP summarization with function sumprot_tmp
    nest_allprot <- nest_allprot %>% 
        mutate(sumtmp = map(data, sumprot_tmp))
    
    return(unnest(nest_allprot, sumtmp) %>% rename(log2inty_tmp = log2inty))
}


# Run-level summarization (across batches)
sum_feature <- function(df_allprot, proteins) {
    # Housekeeping
    proteins <- unique(proteins)
    df_allprot <- df_allprot %>% filter(uniprot_iso %in% proteins)
    # Nested data frame (protein, mod, paired)
    nest_allprot <- df_allprot %>% 
        group_by(uniprot_iso, is_mod, is_paired) %>% 
        nest()
    # TMP summarization with function sumprot_tmp
    nest_allprot <- nest_allprot %>% 
        mutate(sumtmp = map(data, sumprot_tmp))
    
    return(unnest(nest_allprot, sumtmp) %>% rename(log2inty_tmp = log2inty))
}


# Summarization for all proteins ------------------------------------------
prot2sum <- unique(dta_conflu$uniprot_iso)

df_sum <- sum_feature(dta_conflu, prot2sum)  # across batches
df_sumbch <- sum_feature_bch(dta_conflu, prot2sum)  # per batch


# Comparison between two summarization results ----------------------------
df_sum2 <- df_sum %>% 
    left_join(rename(df_sumbch, log2inty_tmpbch = log2inty_tmp)) %>% 
    mutate(is_mod_fac = factor(ifelse(is_mod, "Modified", "Unmodified")), 
           is_par_fac = factor(ifelse(is_paired, "With counterpart", "No counterpart")), 
           batch = ifelse(grepl("B1", run), "B1", "B2"))

df_sum2 %>% 
    ggplot(aes(log2inty_tmp, log2inty_tmpbch, colour = batch)) + 
    geom_point(alpha = 0.4) + 
    geom_abline(intercept = 0, slope = 1) + 
    labs(x = "Summarization (across batches)", y = "Summarization (per batch)") + 
    facet_grid(is_mod_fac ~ is_par_fac)

df_sum2 %>% 
    mutate(log2inty_diff = log2inty_tmpbch - log2inty_tmp) %>% 
    ggplot(aes(batch, log2inty_diff)) + 
    geom_boxplot() + 
    ylab("Difference of per-batch vs. across-batch summarization") + 
    facet_grid(is_mod_fac ~ is_par_fac)


# Full set for inference --------------------------------------------------
# Number of proteins
df_sum2 %>% distinct(uniprot_iso)

# Number of proteins with fully observed unpaired modifications
df_sum2 %>% filter(!is_paired, is_mod) %>% count(uniprot_iso) %>% filter(n == 16)

# Number of proteins with fully observed unpaired modifications and unmodified peptides
df_sum2 %>% filter(!is_paired) %>% count(uniprot_iso) %>% filter(n == 32)


# Eligibility for protein abundance adjustment ----------------------------
# dta_conflu <- dta_conflu %>% left_join(df_protsum)
# dta_conflu %>% group_by(uniprot_iso) %>% summarise(is_allref = all(is_ref)) %>% count(is_allref)

# df_protsum %>% 
#     mutate(is_ref = !(is_mod | is_paired)) %>% 
#     group_by(uniprot_iso, run) %>% 
#     summarise(run_ref = any(is_ref)) %>% 
#     summarise(prot_ref = all(run_ref)) %>% 
#     count(prot_ref)
# 
# df_protsum %>% filter(grepl("B1", run)) %>% 
#     mutate(is_ref = !(is_mod | is_paired)) %>% 
#     group_by(uniprot_iso, run) %>% 
#     summarise(run_ref = any(is_ref)) %>% 
#     summarise(prot_ref = all(run_ref)) %>% 
#     count(prot_ref)
# 
# df_protsum %>% filter(grepl("B2", run)) %>% 
#     mutate(is_ref = !(is_mod | is_paired)) %>% 
#     group_by(uniprot_iso, run) %>% 
#     summarise(run_ref = any(is_ref)) %>% 
#     summarise(prot_ref = all(run_ref)) %>% 
#     count(prot_ref)


# Modeling with nested data frame -----------------------------------------
# Proteins with fully observed unpaired modifications and unmodified peptides
prot_full <- df_sumbch %>% filter(!is_paired) %>% 
    count(uniprot_iso) %>% filter(n == 32) %>% .$uniprot_iso

# Subset unpaired peptides of fully-observed proteins 
df_full <- df_sumbch %>% 
    filter(!is_paired) %>% 
    filter(uniprot_iso %in% prot_full) %>% 
    separate(run, into = c("group", "biotech"), sep = "-", remove = F) %>% 
    mutate(batch = ifelse(grepl("B1", biotech), "B1", "B2"))


# A batch per model
nested_perbch <- df_full %>% 
    group_by(uniprot_iso, is_mod, batch) %>%
    nest()

nested_perbch <- nested_perbch %>% 
    mutate(fit = map(data, ~ lm(log2inty_tmp ~ 0 + group, data = .))) %>% 
    mutate(param = map(fit, tidy), df_res = map_dbl(fit, df.residual))

param_perbch <- nested_perbch %>% 
    unnest(param) %>% 
    mutate(group = gsub("group", "", term)) %>% 
    select(-term, -statistic, -p.value)


# All bathces in one model
nested_allbch <- df_full %>% 
    group_by(uniprot_iso, is_mod) %>%
    nest()

nested_allbch <- nested_allbch %>% 
    mutate(fit = map(data, ~ lm(log2inty_tmp ~ 0 + group + batch, data = .))) %>% 
    mutate(param = map(fit, tidy), df_res = map_dbl(fit, df.residual))

param_allbch <- nested_allbch %>% 
    unnest(param) %>% 
    mutate(group = gsub("group", "", term)) %>% 
    filter(!grepl("batch", term)) %>% 
    select(-term, -statistic, -p.value)


# Group comparison --------------------------------------------------------
cases = c("CCCP_Satur", "DMSO_Light", "DMSO_Satur")
controls = c("CCCP_Light", "CCCP_Light", "DMSO_Light")

test_null1 <- vector("list", length = length(cases))
test_null2 <- vector("list", length = length(cases))
test_null3 <- vector("list", length = length(cases))
test_null4 <- vector("list", length = length(cases))
for (i in seq_along(cases)) {
    grp_ctrl <- controls[i]
    grp_case <- cases[i]
    # Per-batch model
    diff_perbch <- param_perbch %>% 
        filter(group %in% c(grp_ctrl, grp_case)) %>% 
        mutate(key_grp = ifelse(group == grp_ctrl, "G0", "G1")) %>% 
        complete(uniprot_iso, is_mod, batch, key_grp, df_res) %>% 
        group_by(uniprot_iso, is_mod, batch, df_res) %>% 
        summarise(log2fc = diff(estimate), se2 = sum(std.error ^ 2)) %>% 
        group_by(uniprot_iso, is_mod) %>% 
        # mutate(wt_bch = 1 / n()) %>% 
        ungroup()
    # All-batch model
    diff_allbch <- param_allbch %>% 
        filter(group %in% c(grp_ctrl, grp_case)) %>% 
        mutate(key_grp = ifelse(group == grp_ctrl, "G0", "G1")) %>% 
        complete(uniprot_iso, is_mod, key_grp, df_res) %>% 
        group_by(uniprot_iso, is_mod, df_res) %>% 
        summarise(log2fc = diff(estimate), se2 = sum(std.error ^ 2)) %>% 
        ungroup()
    # H^1
    test_null1[[i]] <- diff_perbch %>% filter(is_mod) %>% 
        group_by(uniprot_iso) %>% 
        mutate(wt_bch = 1 / n(), wse2 = wt_bch ^ 2 * se2) %>% 
        summarise(logFC = mean(log2fc), SE = sqrt(sum(wse2)), 
                  DF = sum(wse2) ^ 2 / sum(wse2 ^ 2 / df_res), 
                  t_stat = logFC / SE, p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)) %>% 
        mutate(ctrx = paste(grp_case, grp_ctrl, sep = " vs "), hyp = "null1")
    # H^2
    test_null2[[i]] <- diff_perbch %>% group_by(uniprot_iso, is_mod) %>% 
        mutate(wt_bch = 1 / n(), wse2 = wt_bch ^ 2 * se2, wlog2fc = wt_bch * log2fc) %>% 
        group_by(uniprot_iso) %>% 
        summarise(logFC = sum(wlog2fc[is_mod]) - sum(wlog2fc[!is_mod]), 
                  SE = sqrt(sum(wse2)), DF = sum(wse2) ^ 2 / sum(wse2 ^ 2 / df_res), 
                  t_stat = logFC / SE, p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)) %>% 
        mutate(ctrx = paste(grp_case, grp_ctrl, sep = " vs "), hyp = "null2")
    # H^3 (grouped summary representation is kept for readability)
    test_null3[[i]] <- diff_allbch %>% filter(is_mod) %>% 
        group_by(uniprot_iso) %>% 
        summarise(logFC = log2fc, SE = sqrt(se2), DF = df_res, 
                  t_stat = logFC / SE, p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)) %>% 
        mutate(ctrx = paste(grp_case, grp_ctrl, sep = " vs "), hyp = "null3")
    # H^4
    test_null4[[i]] <- diff_allbch %>% 
        group_by(uniprot_iso) %>% 
        summarise(logFC = diff(log2fc), SE = sqrt(sum(se2)), DF = sum(se2) ^ 2 / sum(se2 ^ 2 / df_res), 
                  t_stat = logFC / SE, p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)) %>% 
        mutate(ctrx = paste(grp_case, grp_ctrl, sep = " vs "), hyp = "null4")
}

test_all <- bind_rows(test_null1, test_null2, test_null3, test_null4) %>% 
    group_by(hyp) %>% 
    mutate(p_adj = p.adjust(p_val, method = "BH")) %>% 
    ungroup()




# Comparison of the testing results ---------------------------------------
# t-statistic
test_all %>% ggplot(aes(hyp, t_stat)) + 
    geom_boxplot() + geom_point() + 
    facet_wrap(~ ctrx)

test_all %>% ggplot(aes(t_stat, colour = hyp)) + 
    geom_density() + 
    facet_wrap(~ ctrx)

# Degrees of freedom
test_all %>% ggplot(aes(hyp, DF)) + 
    geom_boxplot() + geom_jitter() + 
    facet_wrap(~ ctrx)

# Adjusted p-value
test_all %>% ggplot(aes(hyp, -log10(p_adj))) + 
    geom_boxplot() + geom_jitter() + 
    geom_hline(yintercept = -log10(0.05), linetype = 2) + 
    facet_wrap(~ ctrx)

test_all %>% ggplot(aes(hyp, -log10(p_adj))) + 
    geom_boxplot() + geom_point() + 
    geom_hline(yintercept = -log10(0.05), linetype = 2) + 
    facet_wrap(~ ctrx)



# Significant differences -------------------------------------------------
test_all %>% filter(p_adj < 0.05)

test_sig <- test_all %>% filter(p_adj < 0.05)

test_sig %>% group_by(ctrx, hyp) %>% count()


# Compare H_1 vs H_2 in CCCP_Satur vs CCCP_Light
test_sig %>% 
    filter(hyp %in% c("null1", "null2"), ctrx == "CCCP_Satur vs CCCP_Light") %>% 
    group_by(uniprot_iso) %>% 
    count()

prot_h1h2 <- test_sig %>% 
    filter(hyp %in% c("null1", "null2"), ctrx == "CCCP_Satur vs CCCP_Light") %>% 
    group_by(uniprot_iso) %>% 
    count() %>% 
    filter(n == 2) %>% 
    .$uniprot_iso

prot_h1_nh2 <- test_sig %>% 
    filter(hyp %in% c("null1", "null2"), ctrx == "CCCP_Satur vs CCCP_Light") %>% 
    group_by(uniprot_iso) %>% 
    filter(n() == 1) %>% 
    ungroup() %>% 
    filter(hyp == "null1") %>% 
    .$uniprot_iso

prot_h2_nh1 <- test_sig %>% 
    filter(hyp %in% c("null1", "null2"), ctrx == "CCCP_Satur vs CCCP_Light") %>% 
    group_by(uniprot_iso) %>% 
    filter(n() == 1) %>% 
    ungroup() %>% 
    filter(hyp == "null2") %>% 
    .$uniprot_iso

prot_h1_nh3 <- test_sig %>% 
    filter(hyp %in% c("null1", "null3"), ctrx == "CCCP_Satur vs CCCP_Light") %>% 
    group_by(uniprot_iso) %>% 
    filter(n() == 1) %>% 
    ungroup() %>% 
    filter(hyp == "null1") %>% 
    .$uniprot_iso

# There's no change detected with H_3 but missed with H_1
prot_h3_nh1 <- test_sig %>% 
    filter(hyp %in% c("null1", "null3"), ctrx == "CCCP_Satur vs CCCP_Light") %>% 
    group_by(uniprot_iso) %>% 
    filter(n() == 1) %>% 
    ungroup() %>% 
    filter(hyp == "null3") %>% 
    .$uniprot_iso

# pdf("profile_h1_nh3.pdf")
# for (i in seq_along(prot_h1_nh3)) {
#     print(plot_nc_profile(dta_conflu, prot_h1_nh3[i], runlvl_bat))
#     print(plot_nc_sumprofile(dta_conflu, prot_h1_nh3[i], runlvl_bat))
# }
# dev.off()


# Peptide representation --------------------------------------------------
# Remove singleton (considering to move toward upstream)
dta_conflu <- dta_conflu %>% 
    group_by(uniprot_iso, is_paired, is_mod, feature) %>% 
    filter(n() != 1) %>% ungroup()

# dta_conflu <- dta_conflu %>% 
#     group_by(uniprot_iso, peptide_sequence) %>% 
#     filter(n() != 1) %>% ungroup()

# Reading site data
df_fasmod <- readRDS("output/fasmod.rds")
df_fasmod01 <- df_fasmod %>% 
    filter(nb_mod <= 1) %>% 
    select(uniprot_iso, peptide_sequence, nb_mod, site_str, pps_str)

dta_conflu <- dta_conflu %>% inner_join(df_fasmod01)

pep_full <- dta_conflu %>% filter(!is_paired, is_mod) %>% 
    group_by(uniprot_iso, peptide_sequence) %>% 
    summarise(nb_run = n_distinct(run)) %>% 
    ungroup() %>% 
    filter(nb_run == 16)

unmod_full <- dta_conflu %>% 
    filter(!is_paired, !is_mod) %>% 
    group_by(uniprot_iso) %>% 
    summarise(nb_run = n_distinct(run)) %>% 
    filter(nb_run == 16)

# Keeping only the modified peptides (with unpaired unmodified peptides)
df_fullpep <- dta_conflu %>% semi_join(unmod_full) %>% semi_join(pep_full)

# Adding back the unmodified peptides
df_fullpep <- dta_conflu %>% 
    filter(!is_paired, !is_mod, uniprot_iso %in% df_fullpep$uniprot_iso) %>% 
    bind_rows(df_fullpep)


# prot_fullpep <- unique(df_fullpep$uniprot_iso)
# runlvl_bat <- df_fullpep %>% distinct(biorep, run) %>% arrange(biorep, run) %>% .$run
# pdf("profile_fullpep.pdf")
# for (i in seq_along(prot_fullpep)) {
#     print(plot_nc_profile(df_fullpep, prot_fullpep[i], runlvl_bat))
# }
# dev.off()


# Profile plots for site data ---------------------------------------------
plot_sprofile <- function(df_allprot, protein, run_level) {
    df_prot <- df_allprot %>% filter(uniprot_iso == protein) %>% filter(!is_paired)
    # Complete possible combinations of peptide features and runs
    mpar <- df_prot %>% distinct(site_str, feature, is_mod)
    df_fill <- df_prot %>% 
        select(feature, run, log2inty) %>% 
        complete(run = run_level, feature) %>% 
        left_join(mpar) %>% 
        mutate(is_mod_fac = factor(ifelse(is_mod, "Modified", "Unmodified")), 
               run_fac = factor(run, levels = run_level))
    # Feature profiles categorized by modification and matching status
    filter(df_fill, is_mod) %>% 
        ggplot(aes(run_fac, log2inty, group = feature, colour = site_str)) + 
        geom_point(size = 2, alpha = 0.5) +
        geom_line() + 
        geom_line(data = filter(df_fill, !is_mod), aes(run_fac, log2inty, group = feature), colour = "gray") + 
        geom_point(data = filter(df_fill, !is_mod), aes(run_fac, log2inty), colour = "gray", size = 2) + 
        geom_vline(xintercept = 8.5) + 
        facet_grid(is_mod_fac ~ .) + 
        coord_cartesian(ylim = c(10, 35)) + 
        labs(x = "Run", y = "Log2-intensity", title = protein) + 
        theme_bw() + 
        guides(colour = guide_legend(nrow = 1)) + 
        theme(legend.position = "bottom") + 
        theme(axis.text.x = element_blank())
}


# Profile plots of site data
prot_fullpep <- unique(df_fullpep$uniprot_iso)
runlvl_bat <- df_fullpep %>% distinct(biorep, run) %>% arrange(biorep, run) %>% .$run

pdf("profile_site.pdf", width = 8, height = 6)
for (i in seq_along(prot_fullpep)) {
    print(plot_sprofile(df_fullpep, prot_fullpep[i], runlvl_bat))
}
dev.off()




# Profile plots of the 312 proteins ---------------------------------------
runlvl_bat <- dta_conflu %>% distinct(biorep, run) %>% arrange(biorep, run) %>% .$run

pdf("profile_full.pdf")
for (i in seq_along(prot_full)) {
    print(plot_nc_profile(dta_conflu, prot_full[i], runlvl_bat))
    print(plot_nc_sumprofile(dta_conflu, prot_full[i], runlvl_bat))
}
dev.off()


# Profile plots -----------------------------------------------------------
runlvl <- dta_conflu %>% distinct(run) %>% arrange(run) %>% .$run
runlvl_bat <- dta_conflu %>% distinct(biorep, run) %>% arrange(biorep, run) %>% .$run

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

