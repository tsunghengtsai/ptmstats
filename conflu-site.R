# Load libraries and sources ----------------------------------------------
library(tidyverse)
library(stringr)
library(broom)
library(survival)

source("utilities.R")


# Load data set & housekeeping --------------------------------------------
load("output/dta_conflu.RData")

df_conflu <- dta_conflu
rm(dta_conflu)

# Confidence threshold 83
df_conflu <- df_conflu %>% 
    filter(vista_confidence_score >= 83) %>% 
    select(-vista_confidence_score, -peptide_trypticity, 
           -peptide_miscleavages, -peptide_validity)

# OBJ# can be viewed as a batch identity
subj_uniq <- unique(df_design$id_subject)
subj_key <- str_c("S", str_pad(seq_along(subj_uniq), width = 2, pad = "0"))
df_design <- tbl_df(df_design) %>% 
    mutate(
        subj = plyr::mapvalues(id_subject, from = subj_uniq, to = subj_key), 
        biorep = if_else(ann_objective_id == 37122, "B1", "B2"), 
        batch = biorep, 
        techrep = str_c("T", id_injectionset), 
        run_bt = str_c(biorep, techrep), 
        run_cbt = str_c(id_subject, run_bt, sep = "-")
    )


# Annotation & data filtering ---------------------------------------------
df_conflu <- df_conflu %>% 
    rename(peptide = peptide_sequence) %>% 
    mutate(
        is_mod = str_detect(peptide, "\\*"), 
        feature = str_c(peptide, ms2_charge, sep = "_"), 
        log2inty = log2(vista_peak_area_light), 
        condition = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$id_subject), 
        biorep = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$biorep), 
        batch = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$batch), 
        techrep = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$techrep), 
        run = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$run_cbt)
    )


# Extract original peptide sequence & other information
df_pep <- df_conflu %>% 
    distinct(peptide, is_mod) %>% 
    mutate(
        peptide_unmod = str_replace_all(peptide, "\\*", ""), 
        peptide_unmod_trimmed = str_extract(peptide_unmod, "(?<=\\.)([ACDEFGHIKLMNPQRSTVWY\\*]+)"), 
        peptide_trimmed = str_extract(peptide, "(?<=\\.)([ACDEFGHIKLMNPQRSTVWY\\*]+)"), 
        nb_aa = str_length(peptide_unmod_trimmed)
    )

# Retain peptides with length >= 6
df_pep <- df_pep %>% filter(nb_aa >= 6)

# Determine the pairing status of peptide
df_pep <- df_pep %>% group_by(peptide_unmod) %>% 
    mutate(is_paired = ifelse(n_distinct(is_mod) == 2, TRUE, FALSE)) %>% 
    ungroup()

# Initial filtering
regex_uniprot <- ".*([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}).*"
df_conflu <- df_conflu %>% 
    semi_join(df_pep) %>% 
    filter(!is.na(ms2_charge), !is.na(ms2_rt), grepl(regex_uniprot, Reference))

# Extract uniprot accession number
df_conflu <- df_conflu %>% 
    separate(Reference, into = c("p_entry", "uniprot_iso"), sep = "\\|", remove = FALSE)


# Data manipulation & transformation --------------------------------------
# Use the max for duplicate features (per run)
df_conflu <- df_conflu %>% 
    group_by(Reference, uniprot_iso, peptide, is_mod, feature, condition, biorep, run) %>% 
    summarise(log2inty = max(log2inty)) %>% 
    ungroup()

# Remove peptides with ambiguous matches (to >1 Reference): more stringent than 
# mapping with uniprot_iso, as multiple Reference may map to the same uniprot_iso
# Example: df_conflu %>% filter(feature == "K.AFSLK*TSTSAVR.H_3", run == "CCCP_Light-B1T1")
pep2uniprot <- df_conflu %>% 
    distinct(peptide, Reference) %>% 
    count(peptide)

df_conflu <- df_conflu %>% 
    semi_join(filter(pep2uniprot, n == 1))

# Remove singleton features
df_conflu <- df_conflu %>% 
    group_by(uniprot_iso, feature) %>% 
    filter(n() != 1) %>% ungroup()

# Status of peptide pairing & protein normalization
df_conflu <- df_conflu %>% 
    left_join(select(df_pep, peptide, peptide_unmod, is_paired)) %>% 
    group_by(uniprot_iso) %>% 
    mutate(is_ref = any(!is_paired & !is_mod)) %>% ungroup()


# Visualization -----------------------------------------------------------
# QC plot (boxplot)
df_conflu %>% 
    mutate(is_mod_fac = factor(ifelse(is_mod, "Modified", "Unmodified"))) %>% 
    ggplot(aes(run, log2inty)) + 
    geom_boxplot(aes(fill = biorep)) + 
    facet_wrap(~ is_mod_fac) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(x = "Run", y = "Log2-intensity")

# Number of modified/unmodified peptides
df_conflu %>% 
    distinct(condition, biorep, run, peptide, is_mod) %>% 
    group_by(condition, biorep, run) %>% 
    summarise(Modified = sum(is_mod), Unmodified = sum(!is_mod)) %>% 
    ungroup() %>% 
    gather(mod, nb_peptide, Modified:Unmodified) %>% 
    ggplot(aes(condition, nb_peptide)) + 
    geom_jitter(aes(colour = biorep, shape = mod), width = 0.1, size = 4) + 
    labs(x = "Condition", y = "Number of peptides")


# Peptide representation --------------------------------------------------
# Reading site data
df_fasmod <- readRDS("output/fasmod.rds") %>% 
    rename(
        peptide = peptide_sequence, 
        peptide_unmod = unmod_peptide, 
        peptide_unmod_trimmed = trimmed_pep, 
        peptide_trimmed = trimmed_mod
    )

# Ignore multi-site modifications to avoid confounding
df_fasmod01 <- df_fasmod %>% 
    filter(nb_mod <= 1) %>% 
    select(uniprot_iso, peptide, nb_mod, site_str, pps_str)

df_conflu <- df_conflu %>% inner_join(df_fasmod01)

site_full <- df_conflu %>% 
    filter(!is_paired, is_mod) %>% 
    group_by(uniprot_iso, site_str) %>% 
    summarise(nb_run = n_distinct(run)) %>% 
    ungroup() %>% 
    filter(nb_run == 16)

unmod_full <- df_conflu %>% 
    filter(!is_paired, !is_mod) %>% 
    group_by(uniprot_iso) %>% 
    summarise(nb_run = n_distinct(run)) %>% 
    filter(nb_run == 16)

# Keeping only the modified peptides (with unpaired unmodified peptides)
df_site <- df_conflu %>% semi_join(unmod_full) %>% semi_join(site_full)

# Adding back the unmodified peptides
df_site <- df_conflu %>% 
    filter(!is_paired, !is_mod, uniprot_iso %in% df_site$uniprot_iso) %>% 
    bind_rows(df_site)


# Filtering and imputation ------------------------------------------------
# A feature is considered valid only when it yields > 50% coverage
df_site2 <- df_site %>% 
    group_by(uniprot_iso, feature, biorep) %>% 
    mutate(nb_run = n_distinct(run)) %>% 
    ungroup() %>% 
    filter(nb_run > 4)

nested_site <- df_site2 %>% 
    select(uniprot_iso, site_str, feature, biorep, run, log2inty) %>% 
    group_by(uniprot_iso, site_str, biorep) %>% 
    nest()

nested_site <- nested_site %>% 
    mutate(data = map(data, ~ complete(., feature, run)))

nested_site <- nested_site %>% 
    mutate(
        aftdata = map(data, annot_obs), 
        aftdata = map(aftdata, fill_censored)
    )

df_site_aft <- nested_site %>% 
    unnest(aftdata) %>% 
    mutate(is_mod = site_str != "UNMOD")


# Export profile plots ----------------------------------------------------
# Profile plots of site data
prot_fullpep <- sort(unique(df_site$uniprot_iso))
runlvl_bat <- df_site %>% distinct(biorep, run) %>% arrange(biorep, run) %>% .$run

pdf("profile_site.pdf", width = 8, height = 6)
for (i in seq_along(prot_fullpep)) {
    print(plot_sprofile(df_site, prot_fullpep[i], runlvl_bat))
}
dev.off()


pdf("profile_site_aft.pdf", width = 8, height = 6)
for (i in seq_along(prot_fullpep)) {
    print(plot_sprofile_aft(df_site_aft, prot_fullpep[i], runlvl_bat))
}
dev.off()


# Site-level summarization ------------------------------------------------
prot_fullpep <- unique(df_site$uniprot_iso)
df_sumsite <- sum_siteftr_bch(df_site, prot_fullpep)  # per batch


# Site-level modeling -----------------------------------------------------
# Currently with no adjustment by protein abundance
df_sumsite <- df_sumsite %>% 
    filter(site_str != "UNMOD") %>% 
    # mutate(is_mod = (site_str != "UNMOD")) %>% 
    separate(run, into = c("group", "biotech"), sep = "-", remove = F) %>% 
    mutate(batch = ifelse(grepl("B1", biotech), "B1", "B2")) %>% 
    unite(protsite, uniprot_iso, site_str, sep = "-")

# One model per batch
nested_perbch <- df_sumsite %>% 
    group_by(protsite, batch) %>%
    nest() %>% 
    mutate(fit = map(data, ~ lm(log2inty_tmp ~ 0 + group, data = .))) %>% 
    mutate(param = map(fit, tidy), df_res = map_dbl(fit, df.residual))

param_perbch <- nested_perbch %>% 
    unnest(param) %>% 
    mutate(group = gsub("group", "", term)) %>% 
    select(-term, -statistic, -p.value)

# One model for all batches
nested_allbch <- df_sumsite %>% 
    group_by(protsite) %>%
    nest() %>% 
    mutate(fit = map(data, ~ lm(log2inty_tmp ~ 0 + group + batch, data = .))) %>% 
    mutate(param = map(fit, tidy), df_res = map_dbl(fit, df.residual))

param_allbch <- nested_allbch %>% 
    unnest(param) %>% 
    mutate(group = gsub("group", "", term)) %>% 
    filter(!grepl("batch", term)) %>% 
    select(-term, -statistic, -p.value)


# Site-level differential analysis ----------------------------------------
cases = c("CCCP_Satur", "DMSO_Light", "DMSO_Satur")
controls = c("CCCP_Light", "CCCP_Light", "DMSO_Light")

test_null1 <- vector("list", length = length(cases))
# test_null2 <- vector("list", length = length(cases))
test_null3 <- vector("list", length = length(cases))
# test_null4 <- vector("list", length = length(cases))
for (i in seq_along(cases)) {
    grp_ctrl <- controls[i]
    grp_case <- cases[i]
    # Per-batch model
    diff_perbch <- param_perbch %>% 
        filter(group %in% c(grp_ctrl, grp_case)) %>% 
        mutate(key_grp = ifelse(group == grp_ctrl, "G0", "G1")) %>% 
        complete(protsite, batch, key_grp, df_res) %>% 
        group_by(protsite, batch, df_res) %>% 
        summarise(log2fc = diff(estimate), se2 = sum(std.error ^ 2)) %>% 
        ungroup()
    # All-batch model
    diff_allbch <- param_allbch %>% 
        filter(group %in% c(grp_ctrl, grp_case)) %>% 
        mutate(key_grp = ifelse(group == grp_ctrl, "G0", "G1")) %>% 
        complete(protsite, key_grp, df_res) %>% 
        group_by(protsite, df_res) %>% 
        summarise(log2fc = diff(estimate), se2 = sum(std.error ^ 2)) %>% 
        ungroup()
    # H^1
    test_null1[[i]] <- diff_perbch %>% 
        group_by(protsite) %>% 
        mutate(wt_bch = 1 / n(), wse2 = wt_bch ^ 2 * se2) %>% 
        summarise(logFC = mean(log2fc), SE = sqrt(sum(wse2)), 
                  DF = sum(wse2) ^ 2 / sum(wse2 ^ 2 / df_res), 
                  t_stat = logFC / SE, p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)) %>% 
        mutate(ctrx = paste(grp_case, grp_ctrl, sep = " vs "), hyp = "null1")
    # H^3 (grouped summary representation is kept for readability)
    test_null3[[i]] <- diff_allbch %>% 
        group_by(protsite) %>% 
        summarise(logFC = log2fc, SE = sqrt(se2), DF = df_res, 
                  t_stat = logFC / SE, p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)) %>% 
        mutate(ctrx = paste(grp_case, grp_ctrl, sep = " vs "), hyp = "null3")
}

test_all <- bind_rows(test_null1, test_null3) %>% 
    group_by(hyp) %>% 
    mutate(p_adj = p.adjust(p_val, method = "BH")) %>% 
    ungroup()


# Site-level results ------------------------------------------------------
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
    geom_hline(yintercept = -log10(0.05), colour = "darkred") + 
    facet_wrap(~ ctrx)

test_all %>% ggplot(aes(hyp, -log10(p_adj))) + 
    geom_boxplot() + geom_point() + 
    geom_hline(yintercept = -log10(0.05), colour = "darkred") + 
    facet_wrap(~ ctrx)



# Significant site changes ------------------------------------------------
# test_all %>% filter(p_adj < 0.05)
# test_sig %>% group_by(ctrx, hyp) %>% count()

test_sig <- test_all %>% filter(p_adj < 0.05) %>% 
    separate(protsite, into = c("uniprot_iso", "site_str"), sep = "-", remove = FALSE) 

