# Load libraries and sources ----------------------------------------------
library(tidyverse)
library(stringr)
library(broom)
library(survival)

source("utilities.R")


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

df_work <- dta_conflu
rm(dta_conflu)

# Initial filtering
df_work <- df_work %>% 
    filter(vista_confidence_score >= 83, !is.na(ms2_charge), !is.na(ms2_rt)) %>% 
    select(-vista_confidence_score, -peptide_trypticity, 
           -peptide_miscleavages, -peptide_validity) %>% 
    separate(Reference, into = c("p_entry", "uniprot_iso"), sep = "\\|", remove = FALSE) %>% 
    rename(peptide = peptide_sequence)

# OBJ# can be viewed as a batch identity
subj_uniq <- unique(df_design$id_subject)
subj_key <- str_c("S", str_pad(seq_along(subj_uniq), width = 2, pad = "0"))
df_design <- tbl_df(df_design) %>% 
    mutate(
        group = id_subject, 
        subj = plyr::mapvalues(id_subject, from = subj_uniq, to = subj_key), 
        biorep = if_else(ann_objective_id == 37122, "B1", "B2"), 
        batch = if_else(ann_objective_id == 37122, "BCH1", "BCH2"), 
        techrep = str_c("T", id_injectionset), 
        run_bt = str_c(biorep, techrep), 
        run_cbt = str_c(group, run_bt, sep = "-")
    )


# Annotation & data filtering ---------------------------------------------
df_work <- df_work %>% 
    mutate(
        is_mod = str_detect(peptide, "\\*"), 
        feature = str_c(peptide, ms2_charge, sep = "_"), 
        log2inty = ifelse(vista_peak_area_light <= 1, 0, log2(vista_peak_area_light)), 
        group = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$group), 
        biorep = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$biorep), 
        batch = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$batch), 
        techrep = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$techrep), 
        run = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$run_cbt)
    )

# Use the max for duplicate features (per run)
df_work <- df_work %>% 
    group_by(p_entry, uniprot_iso, peptide, is_mod, feature, group, batch, run) %>% 
    summarise(log2inty = max(log2inty)) %>% 
    ungroup()

# Remove peptides with ambiguous matches (to >1 Reference): more stringent than 
# mapping with uniprot_iso, as multiple Reference may map to the same uniprot_iso
# Example: df_work %>% filter(feature == "K.AFSLK*TSTSAVR.H_3", run == "CCCP_Light-B1T1")
uniq_pep <- df_work %>% 
    distinct(peptide, p_entry, uniprot_iso) %>% 
    count(peptide) %>% 
    filter(n == 1) %>% 
    .$peptide

df_work <- df_work %>% filter(peptide %in% uniq_pep)

# Remove singleton features
df_work <- df_work %>%
    group_by(uniprot_iso, feature, batch) %>%
    filter(n() != 1) %>% ungroup()

# Extract original peptide sequence & other information
df_mod <- df_work %>% 
    distinct(uniprot_iso, peptide) %>% 
    mutate(
        peptide_trimmed = str_extract(peptide, "(?<=\\.)([ACDEFGHIKLMNPQRSTVWY\\*]+)"), 
        peptide_unmod = str_replace_all(peptide, "\\*", ""), 
        peptide_unmod_trimmed = str_extract(peptide_unmod, "(?<=\\.)([ACDEFGHIKLMNPQRSTVWY\\*]+)"), 
        len_peptide = str_length(peptide_unmod_trimmed)
    ) %>% 
    arrange(uniprot_iso)

# Retain peptides with length >= 6
regex_uniprot <- ".*([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}).*"
df_mod <- df_mod %>% 
    # filter(len_peptide >= 6, str_detect(p_entry, "HUMAN")) %>% 
    filter(len_peptide >= 6, str_detect(uniprot_iso, regex_uniprot))


# Site representation for modified peptides -------------------------------
load("output/uniprot_inf.RData")  # Run annotate-seq.R

hs_fasta <- hs_fasta %>% 
    select(uniprot_ac, uniprot_iso, p_entry, mod_idx, mod_res, header, sequence) %>% 
    arrange(uniprot_iso)

df_fasmod <- df_mod %>% 
    inner_join(hs_fasta %>% select(-header)) %>% 
    mutate(pep_idx = str_locate_all(sequence, peptide_unmod_trimmed)) %>% 
    mutate(nb_mch = map_int(pep_idx, ~ nrow(.)))

# Ignoring non-specific matching (peptide mapped to 0 or >1 locations of protein)
df_fasmod <- df_fasmod %>% filter(nb_mch == 1)

# Modification sites
mod_resymb <- str_c(mod_residue, mod_symbol)
df_fasmod <- df_fasmod %>% 
    mutate(
        site_all = map2(peptide_unmod_trimmed, pep_idx, locate_site, mod_residue), 
        site_mod = map2(peptide_trimmed, pep_idx, locate_mod, mod_resymb)
    ) %>% 
    mutate(
        site_null = map2(site_all, site_mod, setdiff),
        nb_site = map_int(site_all, length), 
        nb_mod = map_int(site_mod, length), 
        site_str = map2_chr(site_mod, mod_idx, ant_site, mod_residue), 
        pep_str = map_chr(pep_idx, ~ str_c(., collapse = "-"))
    ) %>% 
    mutate(pps_str = str_c(uniprot_iso, pep_str, site_str, sep = "_"))

# Determine the pairing status of peptide
# [TODO]: determine the pairing status with "site"
df_fasmod <- df_fasmod %>% 
    mutate(is_mod = str_detect(peptide, "\\*")) %>% 
    group_by(peptide_unmod) %>% 
    mutate(is_paired = ifelse(n_distinct(is_mod) == 2, TRUE, FALSE)) %>% 
    ungroup()

# Status of peptide pairing & protein normalization
df_work <- df_work %>% 
    left_join(select(df_fasmod, peptide, peptide_unmod, is_paired))


# Visualization -----------------------------------------------------------
# QC plot (boxplot)
df_work %>% 
    mutate(is_mod_fac = factor(ifelse(is_mod, "Modified", "Unmodified"))) %>% 
    ggplot(aes(run, log2inty)) + 
    geom_boxplot(aes(fill = batch)) + 
    facet_wrap(~ is_mod_fac) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(x = "Run", y = "Log2-intensity")

# Number of modified/unmodified peptides
df_work %>% 
    distinct(group, batch, run, peptide, is_mod) %>% 
    group_by(group, batch, run) %>% 
    summarise(Modified = sum(is_mod), Unmodified = sum(!is_mod)) %>% 
    ungroup() %>% 
    gather(mod, nb_peptide, Modified:Unmodified) %>% 
    ggplot(aes(group, nb_peptide)) + 
    geom_jitter(aes(colour = batch, shape = mod), width = 0.1, size = 4) + 
    labs(x = "Group", y = "Number of peptides")


# Within-batch normalization ----------------------------------------------
# Based on unpaired unmodified peptides
medians <- df_work %>% 
    filter(!is_paired, !is_mod) %>% 
    group_by(batch, run) %>% 
    summarise(log2inty_med = median(log2inty, na.rm = TRUE)) %>% 
    mutate(log2inty_bch = median(log2inty_med)) %>% 
    ungroup()

df_work <- df_work %>% 
    left_join(medians) %>% 
    mutate(log2inty = ifelse(is.na(log2inty_med), log2inty, log2inty - log2inty_med + log2inty_bch)) %>% 
    select(-log2inty_med, -log2inty_bch)


# Potential filtering step ------------------------------------------------
# Ignore multi-site modifications to avoid confounding
df_fasmod01 <- df_fasmod %>% 
    filter(nb_mod <= 1) %>% 
    select(uniprot_iso, peptide, nb_mod, site_str, pps_str)

df_work <- df_work %>% inner_join(df_fasmod01)

# Fully observed sites
site_full <- df_work %>%
    filter(!is_paired, is_mod) %>%
    group_by(uniprot_iso, site_str) %>%
    summarise(nb_run = n_distinct(run)) %>%
    ungroup() %>%
    filter(nb_run == n_distinct(df_design$run_id))

df_site <- df_work %>% 
    filter(!is_mod, !is_paired) %>% 
    semi_join(site_full %>% select(uniprot_iso)) %>% 
    bind_rows(df_work %>% semi_join(site_full %>% select(uniprot_iso, site_str)))

nested_site <- df_site %>%
    select(uniprot_iso, site_str, feature, batch, run, log2inty) %>%
    group_by(uniprot_iso, site_str, batch) %>%
    nest()


# Imputation and summarization --------------------------------------------
nested_site <- nested_site %>% 
    mutate(data = map(data, ~ complete(., feature, run)))

# AFT imputation and MP summarization
nested_site <- nested_site %>%
    mutate(aftdata = map(data, annot_obs)) %>%
    mutate(aftdata = map(aftdata, fill_censored)) %>%
    mutate(mp_aftsum = map(aftdata, sumprot_tmp))


# Fit per-batch and all-batch models --------------------------------------
# Add group information
run2group <- df_work %>% distinct(run, group)
df_sum_aft <- nested_site %>% unnest(mp_aftsum) %>% left_join(run2group)

# One model per batch
nested_perbch_aft <- df_sum_aft %>% 
    group_by(uniprot_iso, site_str, batch) %>% 
    nest() %>% 
    mutate(lm_fit = map(data, lm_perbch)) %>% 
    mutate(
        param = map2(lm_fit, data, tidy_bch), 
        df_res = map_dbl(lm_fit, df.residual)
    )

# One model for all batches (need data from >1 batch)
nested_allbch_aft <- df_sum_aft %>% 
    group_by(uniprot_iso, site_str) %>% 
    nest() %>% 
    mutate(nb_bch = map_int(data, ~ n_distinct(.$batch))) %>% 
    filter(nb_bch > 1) %>% 
    select(-nb_bch)

nested_allbch_aft <- nested_allbch_aft %>% 
    mutate(lm_fit = map(data, lm_allbch)) %>% 
    mutate(
        param = map2(lm_fit, data, tidy_bch), 
        df_res = map_dbl(lm_fit, df.residual)
    )


# Extract estimated parameters --------------------------------------------
site_perbch <- nested_perbch_aft %>% 
    filter(site_str != "UNMOD") %>% 
    unnest(param)

site_allbch <- nested_allbch_aft %>% 
    filter(site_str != "UNMOD") %>% 
    unnest(param)

unmod_perbch <- nested_perbch_aft %>% 
    filter(site_str == "UNMOD") %>% 
    unnest(param) %>% 
    select(-site_str) %>% 
    rename(df_unmod = df_res, est_unmod = estimate, se_unmod = std.error)

unmod_allbch <- nested_allbch_aft %>% 
    filter(site_str == "UNMOD") %>% 
    unnest(param) %>% 
    select(-site_str) %>% 
    rename(df_unmod = df_res, est_unmod = estimate, se_unmod = std.error)

param_perbch <- site_perbch %>% left_join(unmod_perbch) %>% unite(protsite, uniprot_iso, site_str, sep = "--")
param_allbch <- site_allbch %>% left_join(unmod_allbch) %>% unite(protsite, uniprot_iso, site_str, sep = "--")

# Keep parameters with std error for across-batch inference
param_perbch <- param_perbch %>% 
    filter(!is.na(std.error)) %>% 
    group_by(protsite, group) %>% 
    filter(n() == 2) %>% 
    ungroup()
param_allbch <- param_allbch %>% filter(!is.na(std.error))


# Site-level differential analysis ----------------------------------------
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
    test_null1[[i]] <- test_mod_bch(param_perbch, grp_ctrl, grp_case) %>% mutate(hyp = "null1")
    test_null2[[i]] <- test_adjmod_bch(param_perbch, grp_ctrl, grp_case) %>% mutate(hyp = "null2")
    # All-batch model
    test_null3[[i]] <- test_mod(param_allbch, grp_ctrl, grp_case) %>% mutate(hyp = "null3")
    test_null4[[i]] <- test_adjmod(param_allbch, grp_ctrl, grp_case) %>% mutate(hyp = "null4")
}

test_all <- bind_rows(test_null1, test_null2, test_null3, test_null4) %>% 
    group_by(hyp) %>%
    # group_by(hyp, ctrx) %>% 
    mutate(p_adj = p.adjust(p_val, method = "BH")) %>% 
    ungroup()
test_all$p_adj[test_all$logFC %in% c(Inf, -Inf)] <- 0  # Missing in one group

test_all <- test_all %>% 
    mutate(
        model = ifelse(hyp %in% c("null1", "null2"), "per-batch", "all-batch"), 
        protadj = hyp %in% c("null2", "null4")
    )


# Site-level results ------------------------------------------------------
# Volcano plot
test_all %>% 
    filter(abs(logFC) < Inf) %>%
    ggplot(aes(logFC, -log10(p_adj))) + 
    geom_point() + 
    geom_hline(yintercept = -log10(0.05), color = "darkred") + 
    facet_grid(model + protadj ~ ctrx)


# Significant site changes ------------------------------------------------
test_sig <- test_all %>% filter(p_adj < 0.05) %>% 
    separate(protsite, into = c("uniprot_iso", "site_str"), sep = "--", remove = FALSE) %>% 
    select(protsite, uniprot_iso, site_str, ctrx, logFC, SE, DF, t_stat, p_val, p_adj, hyp, model, protadj)

# test_sig %>% count(ctrx, hyp)

# write_csv(test_sig, "output/cccp-full.csv")