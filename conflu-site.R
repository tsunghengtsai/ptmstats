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
        group = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$id_subject), 
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
# [TODO]: postpone this step after filtering
# [TODO]: determine the pairing status with "site"
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
    group_by(Reference, uniprot_iso, peptide, is_mod, feature, group, batch, run) %>% 
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
    group_by(uniprot_iso, feature, batch) %>%
    filter(n() != 1) %>% ungroup()

# Remove non-valid features without sufficient coverage (50%)
# df_conflu <- df_conflu %>% 
#     group_by(uniprot_iso, feature, batch) %>% 
#     filter(n() >= 4) %>% ungroup()

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
    geom_boxplot(aes(fill = batch)) + 
    facet_wrap(~ is_mod_fac) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(x = "Run", y = "Log2-intensity")

# Number of modified/unmodified peptides
df_conflu %>% 
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
med_conflu <- df_conflu %>% 
    filter(!is_paired, !is_mod) %>% 
    group_by(batch, run) %>% 
    summarise(log2inty_med = median(log2inty, na.rm = TRUE)) %>% 
    mutate(log2inty_bch = median(log2inty_med)) %>% 
    ungroup()

df_conflu <- df_conflu %>% 
    left_join(med_conflu) %>% 
    mutate(log2inty = log2inty - log2inty_med + log2inty_bch)

# QC plot (boxplot)
# df_conflu %>% 
#     filter(!is_paired) %>% 
#     mutate(is_mod_fac = factor(ifelse(is_mod, "Modified", "Unmodified"))) %>% 
#     ggplot(aes(run, log2inty)) + 
#     geom_boxplot(aes(fill = batch)) + 
#     facet_wrap(~ is_mod_fac) + 
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
#     labs(x = "Run", y = "Log2-intensity", title = "Unpaired features")


# Site representation -----------------------------------------------------
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

# Focus on conditions with measurements (for a site) in both batches
site_xbch <- df_conflu %>% 
    filter(!is_paired, is_mod) %>% 
    group_by(uniprot_iso, site_str, group) %>% 
    summarise(nb_run = n_distinct(run)) %>% 
    ungroup() %>% 
    filter(nb_run == 4)

unmod_xbch <- df_conflu %>% 
    filter(!is_paired, !is_mod) %>% 
    group_by(uniprot_iso, group) %>% 
    summarise(nb_run = n_distinct(run)) %>% 
    ungroup() %>% 
    filter(nb_run == 4)

# Keeping only the modified peptides (with unpaired unmodified peptides)
df_site <- df_conflu %>% semi_join(unmod_xbch) %>% semi_join(site_xbch)

# Adding back the unmodified peptides
df_site <- df_conflu %>% 
    filter(!is_paired, !is_mod) %>% 
    semi_join(df_site %>% distinct(uniprot_iso, group)) %>% 
    bind_rows(df_site)


# Imputation and summarization --------------------------------------------
nested_site <- df_site %>% 
    select(uniprot_iso, site_str, feature, batch, run, log2inty) %>% 
    group_by(uniprot_iso, site_str, batch) %>% 
    nest()

nested_site <- nested_site %>% 
    mutate(data = map(data, ~ complete(., feature, run)))

# AFT imputation and MP summarization
nested_site <- nested_site %>% 
    mutate(
        mp_sum = map(data, sumprot_tmp), 
        aftdata = map(data, annot_obs)
    ) %>% 
    mutate(aftdata = map(aftdata, fill_censored)) %>% 
    mutate(mp_aftsum = map(aftdata, sumprot_tmp))


# Fit per-batch and all-batch models --------------------------------------
# Add group information
run2group <- df_site %>% distinct(run, group)
df_sum <- nested_site %>% unnest(mp_sum) %>% left_join(run2group)
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


# One model for all batches
nested_allbch_aft <- df_sum_aft %>% 
    group_by(uniprot_iso, site_str) %>% 
    nest() %>% 
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

param_perbch <- site_perbch %>% left_join(unmod_perbch) %>% unite(protsite, uniprot_iso, site_str, sep = "-")
param_allbch <- site_allbch %>% left_join(unmod_allbch) %>% unite(protsite, uniprot_iso, site_str, sep = "-")


# Export profile plots ----------------------------------------------------
# Profile plots of site data
prot_fullpep <- sort(unique(df_site$uniprot_iso))
runlvl_bat <- df_site %>% distinct(batch, run) %>% arrange(batch, run) %>% .$run

pdf("profile_site.pdf", width = 8, height = 6)
for (i in seq_along(prot_fullpep)) {
    print(plot_sprofile(df_site, prot_fullpep[i], runlvl_bat))
}
dev.off()


df_site_aft <- nested_site %>% 
    unnest(aftdata) %>% 
    mutate(is_mod = site_str != "UNMOD")

df_site_aft <- df_site %>% 
    distinct(group, batch, run) %>% 
    right_join(df_site_aft)

df_site_aft <- df_site %>% 
    distinct(uniprot_iso, site_str, feature, peptide, peptide_unmod) %>% 
    right_join(df_site_aft)

pdf("profile_site_aft.pdf", width = 8, height = 6)
for (i in seq_along(prot_fullpep)) {
    print(plot_sprofile_aft(df_site_aft, prot_fullpep[i], runlvl_bat))
}
dev.off()


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
    full_perbch <- param_perbch %>% 
        filter(group %in% c(grp_ctrl, grp_case)) %>% 
        mutate(key_grp = ifelse(group == grp_ctrl, "G0", "G1")) %>% 
        complete(protsite, batch, key_grp)
    diff_perbch <- full_perbch %>% 
        group_by(protsite, batch) %>% 
        filter(!any(is.na(df_res))) %>% 
        group_by(protsite, batch, df_res, df_unmod) %>% 
        summarise(
            log2fc = diff(estimate), se2 = sum(std.error ^ 2), 
            log2fc_unmod = diff(est_unmod), se2_unmod = sum(se_unmod ^ 2)
        ) %>% ungroup()
    part_perbch <- full_perbch %>% anti_join(diff_perbch)
    untest_perbch <- part_perbch %>% 
        filter(!is.na(df_res)) %>% 
        distinct(protsite, key_grp) %>% 
        mutate(
            logFC = ifelse(key_grp == "G1", Inf, -Inf), 
            SE = NA, DF = NA, t_stat = NA, p_val = NA, 
            ctrx = paste(grp_case, grp_ctrl, sep = " vs ")
        ) %>% select(-key_grp)
    # All-batch model
    full_allbch <- param_allbch %>% 
        filter(group %in% c(grp_ctrl, grp_case)) %>% 
        mutate(key_grp = ifelse(group == grp_ctrl, "G0", "G1")) %>% 
        complete(protsite, key_grp)
    diff_allbch <- full_allbch %>% 
        group_by(protsite) %>% 
        filter(!any(is.na(df_res))) %>% 
        group_by(protsite, df_res, df_unmod) %>% 
        summarise(
            log2fc = diff(estimate), se2 = sum(std.error ^ 2), 
            log2fc_unmod = diff(est_unmod), se2_unmod = sum(se_unmod ^ 2)
        ) %>% ungroup()
    part_allbch <- full_allbch %>% anti_join(diff_allbch)
    untest_allbch <- part_allbch %>% 
        filter(!is.na(df_res)) %>% 
        distinct(protsite, key_grp) %>% 
        mutate(
            logFC = ifelse(key_grp == "G1", Inf, -Inf), 
            SE = NA, DF = NA, t_stat = NA, p_val = NA, 
            ctrx = paste(grp_case, grp_ctrl, sep = " vs ")
        ) %>% select(-key_grp)
    # H^1
    test_null1[[i]] <- diff_perbch %>% 
        group_by(protsite) %>% 
        summarise(logFC = mean(log2fc), SE = sqrt(sum(se2)) / n(), 
                  DF = sum(se2) ^ 2 / sum(se2 ^ 2 / df_res), 
                  t_stat = logFC / SE, p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)) %>% 
        mutate(ctrx = paste(grp_case, grp_ctrl, sep = " vs "), hyp = "null1") %>% 
        bind_rows(untest_perbch %>% mutate(hyp = "null1"))
    # H^2
    test_null2[[i]] <- diff_perbch %>% 
        group_by(protsite) %>% 
        summarise(
            logFC = mean(log2fc) - mean(log2fc_unmod), 
            SE = sqrt(sum(se2) + sum(se2_unmod)) / n(), 
            DF = (sum(se2) + sum(se2_unmod)) ^ 2 / sum(se2 ^ 2 / df_res + se2_unmod ^ 2 / df_unmod), 
            t_stat = logFC / SE, p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)
        ) %>% 
        mutate(ctrx = paste(grp_case, grp_ctrl, sep = " vs "), hyp = "null2") %>% 
        bind_rows(untest_perbch %>% mutate(hyp = "null2"))
    # H^3
    test_null3[[i]] <- diff_allbch %>% 
        mutate(logFC = log2fc, SE = sqrt(se2), DF = df_res, 
               t_stat = logFC / SE, p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)) %>% 
        select(protsite, logFC, SE, DF, t_stat, p_val) %>% 
        mutate(ctrx = paste(grp_case, grp_ctrl, sep = " vs "), hyp = "null3") %>% 
        bind_rows(untest_allbch %>% mutate(hyp = "null3"))
    # H^4
    test_null4[[i]] <- diff_allbch %>% 
        mutate(logFC = log2fc - log2fc_unmod, SE = sqrt(se2 + se2_unmod), 
               DF = (se2 + se2_unmod) ^ 2 / (se2 ^ 2 / df_res + se2_unmod ^ 2 / df_unmod), 
               t_stat = logFC / SE, p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)) %>% 
        select(protsite, logFC, SE, DF, t_stat, p_val) %>% 
        mutate(ctrx = paste(grp_case, grp_ctrl, sep = " vs "), hyp = "null4") %>% 
        bind_rows(untest_allbch %>% mutate(hyp = "null4"))
}

test_all <- bind_rows(test_null1, test_null2, test_null3, test_null4) %>% 
    group_by(hyp) %>% 
    mutate(p_adj = p.adjust(p_val, method = "BH")) %>% 
    ungroup()
test_all$p_adj[test_all$logFC %in% c(Inf, -Inf)] <- 0  # Missing in one group

test_all <- test_all %>% 
    mutate(
        model = ifelse(hyp %in% c("null1", "null2"), "per-batch", "all-batch"), 
        protadj = hyp %in% c("null2", "null4")
    )


# Site-level results ------------------------------------------------------
# t-statistic
test_all %>% ggplot(aes(hyp, t_stat)) + 
    geom_boxplot() + geom_point() + 
    facet_wrap(~ ctrx)

test_all %>% ggplot(aes(t_stat, color = hyp)) + 
    geom_density() + 
    facet_wrap(~ ctrx)

test_all %>% ggplot(aes(logFC, color = hyp)) + 
    geom_density() + 
    facet_wrap(~ ctrx)

test_all %>% ggplot(aes(model, t_stat, fill = protadj)) + 
    geom_boxplot() + 
    geom_hline(yintercept = 0, color = "darkred") + 
    facet_wrap(~ ctrx)

test_all %>% ggplot(aes(model, SE, fill = protadj)) + 
    geom_boxplot() + 
    geom_hline(yintercept = 0) + 
    facet_wrap(~ ctrx)


# Degrees of freedom
test_all %>% ggplot(aes(hyp, DF)) + 
    geom_boxplot() + geom_jitter() + 
    facet_wrap(~ ctrx)

test_all %>% ggplot(aes(model, DF, fill = protadj)) + 
    geom_boxplot() + geom_jitter(alpha = 0.25) + 
    facet_wrap(~ ctrx)


# Adjusted p-value
test_all %>% ggplot(aes(hyp, -log10(p_adj))) + 
    geom_boxplot() + geom_jitter() + 
    geom_hline(yintercept = -log10(0.05), color = "darkred") + 
    facet_wrap(~ ctrx)

test_all %>% ggplot(aes(model, -log10(p_adj), fill = protadj)) + 
    geom_boxplot() + geom_jitter(alpha = 0.25) + 
    geom_hline(yintercept = -log10(0.05), color = "darkred") + 
    facet_wrap(~ ctrx)

# Volcano plot
test_all %>% ggplot(aes(logFC, -log10(p_adj))) + 
    geom_point() + 
    geom_hline(yintercept = -log10(0.05), color = "darkred") + 
    facet_grid(protadj ~ model)


# Significant site changes ------------------------------------------------
test_sig <- test_all %>% filter(p_adj < 0.05) %>% 
    separate(protsite, into = c("uniprot_iso", "site_str"), sep = "-", remove = FALSE)

# test_all %>% filter(p_adj < 0.05)
# test_sig %>% group_by(ctrx, hyp) %>% count()
