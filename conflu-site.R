library(tidyverse)
library(stringr)
library(broom)


# Load data set & housekeeping --------------------------------------------
load("output/dta_conflu.RData")

# Confidence threshold
dta_conflu <- dta_conflu %>% filter(vista_confidence_score >= 83)

df_design <- tbl_df(df_design)

dta_conflu <- dta_conflu %>% 
    select(-vista_confidence_score, -peptide_trypticity, -peptide_miscleavages, 
           -peptide_validity)

# OBJ# can be viewed as a batch identity
subj_uniq <- unique(df_design$id_subject)
subj_key <- str_c("S", str_pad(seq_along(subj_uniq), width = 2, pad = "0"))
df_design <- df_design %>% 
    mutate(
        subj = plyr::mapvalues(id_subject, from = subj_uniq, to = subj_key), 
        biorep = if_else(ann_objective_id == 37122, "B1", "B2"), 
        techrep = str_c("T", id_injectionset), 
        run_bt = str_c(biorep, techrep), 
        run_cbt = str_c(id_subject, run_bt, sep = "-")
    )


# Annotation & data filtering ---------------------------------------------
dta_conflu <- dta_conflu %>% 
    mutate(
        is_mod = str_detect(peptide_sequence, "\\*"), 
        feature = str_c(peptide_sequence, ms2_charge, sep = "_"), 
        log2inty = log2(vista_peak_area_light), 
        cond = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$id_subject), 
        biorep = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$biorep), 
        techrep = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$techrep), 
        run = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$run_cbt)
    )

# Extract original peptide sequence & other information
df_pep <- dta_conflu %>% 
    distinct(peptide_sequence, is_mod) %>% 
    mutate(
        unmod_peptide = str_replace_all(peptide_sequence, "\\*", ""), 
        trimmed_pep = str_extract(unmod_peptide, "(?<=\\.)([ACDEFGHIKLMNPQRSTVWY\\*]+)"), 
        trimmed_mod = str_extract(peptide_sequence, "(?<=\\.)([ACDEFGHIKLMNPQRSTVWY\\*]+)"), 
        nb_aa = str_length(trimmed_pep)
    )

# Retain peptides with length >= 6
df_pep <- df_pep %>% filter(nb_aa >= 6)

# Determine the pairing status of peptide
df_pep <- df_pep %>% group_by(unmod_peptide) %>% 
    mutate(is_paired = ifelse(n_distinct(is_mod) == 2, TRUE, FALSE)) %>% 
    ungroup()

# Initial filtering
regex_uniprot <- ".*([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}).*"
dta_conflu <- dta_conflu %>% 
    semi_join(df_pep) %>% 
    filter(!is.na(ms2_charge), !is.na(ms2_rt), grepl(regex_uniprot, Reference))

# Extract uniprot accession number
dta_conflu <- dta_conflu %>% 
    separate(Reference, into = c("p_entry", "uniprot_iso"), sep = "\\|", remove = FALSE)


# Data manipulation & transformation --------------------------------------
# Use the max for duplicate features (per run)
dta_conflu <- dta_conflu %>% 
    group_by(Reference, uniprot_iso, peptide_sequence, is_mod, feature, cond, biorep, run) %>% 
    summarise(log2inty = max(log2inty)) %>% 
    ungroup()

# Remove peptides with ambiguous matches (to >1 Reference): more stringent than 
# mapping with uniprot_iso, as multiple Reference may map to the same uniprot_iso
# Example: dta_conflu %>% filter(feature == "K.AFSLK*TSTSAVR.H_3", run == "CCCP_Light-B1T1")
pep2uniprot <- dta_conflu %>% 
    distinct(peptide_sequence, Reference) %>% 
    count(peptide_sequence)

dta_conflu <- dta_conflu %>% 
    semi_join(filter(pep2uniprot, n == 1))

# Remove singleton features
dta_conflu <- dta_conflu %>% 
    group_by(uniprot_iso, feature) %>% 
    filter(n() != 1) %>% ungroup()

# Status of peptide pairing & protein normalization
dta_conflu <- dta_conflu %>% 
    left_join(select(df_pep, peptide_sequence, unmod_peptide, is_paired)) %>% 
    group_by(uniprot_iso) %>% 
    mutate(is_ref = any(!is_paired & !is_mod)) %>% ungroup()


# Visualization -----------------------------------------------------------
# QC plot (boxplot)
dta_conflu %>% 
    mutate(is_mod_fac = factor(ifelse(is_mod, "Modified", "Unmodified"))) %>% 
    ggplot(aes(run, log2inty)) + 
    geom_boxplot(aes(fill = biorep)) + 
    facet_wrap(~ is_mod_fac) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(x = "Run", y = "Log2-intensity")

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


# Summarization within batch 
#  - require site_str annotation
#  - not conisder unpaired features, multi-site modifications
sum_siteftr_bch <- function(df_allprot, proteins) {
    # Housekeeping
    proteins <- unique(proteins)
    df_allprot <- df_allprot %>% filter(uniprot_iso %in% proteins)
    # Nested data frame (protein, site, batch)
    nest_allprot <- df_allprot %>% 
        group_by(uniprot_iso, site_str, biorep) %>% 
        nest()
    # TMP summarization with function sumprot_tmp
    nest_allprot <- nest_allprot %>% 
        mutate(sumtmp = map(data, sumprot_tmp))
    
    return(unnest(nest_allprot, sumtmp) %>% rename(log2inty_tmp = log2inty))
}


# Peptide representation --------------------------------------------------
# Reading site data
df_fasmod <- readRDS("output/fasmod.rds")

# Ignore multi-site modifications to avoid confounding
df_fasmod01 <- df_fasmod %>% 
    filter(nb_mod <= 1) %>% 
    select(uniprot_iso, peptide_sequence, nb_mod, site_str, pps_str)

dta_conflu <- dta_conflu %>% inner_join(df_fasmod01)

pep_full <- dta_conflu %>% 
    filter(!is_paired, is_mod) %>% 
    group_by(uniprot_iso, site_str) %>% 
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
        geom_line(alpha = 0.75) + 
        geom_line(data = filter(df_fill, !is_mod), aes(run_fac, log2inty, group = feature), colour = "gray") + 
        geom_point(data = filter(df_fill, !is_mod), aes(run_fac, log2inty), colour = "gray", size = 2) + 
        geom_vline(xintercept = 8.5) + 
        facet_grid(is_mod_fac ~ .) + 
        coord_cartesian(ylim = c(10, 35)) + 
        labs(x = "Run", y = "Log2-intensity", title = protein) + 
        theme_bw() + 
        guides(colour = guide_legend(nrow = 1, title = NULL)) + 
        # theme(legend.position = "bottom") + 
        theme(legend.position = c(0.5, 0.065)) + 
        theme(axis.text.x = element_blank())
}



# Export profile plots ----------------------------------------------------
# Profile plots of site data
prot_fullpep <- unique(df_fullpep$uniprot_iso)
runlvl_bat <- df_fullpep %>% distinct(biorep, run) %>% arrange(biorep, run) %>% .$run

pdf("profile_site.pdf", width = 8, height = 6)
for (i in seq_along(prot_fullpep)) {
    print(plot_sprofile(df_fullpep, prot_fullpep[i], runlvl_bat))
}
dev.off()


# Site-level summarization ------------------------------------------------
prot_fullpep <- unique(df_fullpep$uniprot_iso)
df_sumsite <- sum_siteftr_bch(df_fullpep, prot_fullpep)  # per batch


# Site-level modeling -----------------------------------------------------
# Currently with no adjustment by protein abundance
df_sumsite <- df_sumsite %>% 
    filter(site_str != "UNMOD") %>% 
    # mutate(is_mod = (site_str != "UNMOD")) %>% 
    separate(run, into = c("group", "biotech"), sep = "-", remove = F) %>% 
    mutate(batch = ifelse(grepl("B1", biotech), "B1", "B2")) %>% 
    unite(protsite, uniprot_iso, site_str, sep = "-")

# A batch per model
nested_perbch <- df_sumsite %>% 
    group_by(protsite, batch) %>%
    nest() %>% 
    mutate(fit = map(data, ~ lm(log2inty_tmp ~ 0 + group, data = .))) %>% 
    mutate(param = map(fit, tidy), df_res = map_dbl(fit, df.residual))

param_perbch <- nested_perbch %>% 
    unnest(param) %>% 
    mutate(group = gsub("group", "", term)) %>% 
    select(-term, -statistic, -p.value)


# All bathces in one model
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




# Some plots --------------------------------------------------------------

# df_fullpep %>% plot_sprofile(pp, runlvl_bat)

# USP30 "UBP30_HUMAN" - 00
pp <- "Q70CQ3"

# MFN1 "MFN1_HUMAN" - 00
pp <- "Q8IWA4"

# Parkin "PRKN2_HUMAN" - **
pp <- "O60260"
test_sig %>% filter(str_detect(protsite, pp)) %>% arrange(protsite, ctrx, hyp)
df_fullpep %>% plot_sprofile(pp, runlvl_bat)
df_fullpep %>% semi_join(filter(test_sig, uniprot_iso == pp)) %>% plot_sprofile(pp, runlvl_bat)

# Tomm20 "TOM20_HUMAN" - 00
pp <- "Q15388"

# VDAC "VDAC3_HUMAN" - 00
pp <- "Q9Y277"

# VDAC "VDAC1_HUMAN" **
pp <- "P21796"
test_sig %>% filter(str_detect(protsite, pp)) %>% arrange(protsite, ctrx, hyp)
df_fullpep %>% filter(is_mod) %>% plot_sprofile(pp, runlvl_bat)
df_fullpep %>% semi_join(filter(test_sig, uniprot_iso == pp)) %>% plot_sprofile(pp, runlvl_bat)

# VDAC "VDAC2_HUMAN" ** 
pp <- "P45880"
test_sig %>% filter(str_detect(protsite, pp)) %>% arrange(protsite, ctrx, hyp)
df_fullpep %>% filter(is_mod) %>% plot_sprofile(pp, runlvl_bat)
df_fullpep %>% semi_join(filter(test_sig, uniprot_iso == pp)) %>% plot_sprofile(pp, runlvl_bat)

# GAPDH "G3P_HUMAN" **
pp <- "P04406"
test_sig %>% filter(str_detect(protsite, pp)) %>% arrange(protsite, ctrx, hyp)
df_fullpep %>% filter(is_mod) %>% plot_sprofile(pp, runlvl_bat)
df_fullpep %>% semi_join(filter(test_sig, uniprot_iso == pp)) %>% plot_sprofile(pp, runlvl_bat)

# Actin "ACTB_HUMAN" **
pp <- "P60709"
test_sig %>% filter(str_detect(protsite, pp)) %>% arrange(protsite, ctrx, hyp)
df_fullpep %>% filter(is_mod) %>% plot_sprofile(pp, runlvl_bat)
df_fullpep %>% semi_join(filter(test_sig, uniprot_iso == pp)) %>% plot_sprofile(pp, runlvl_bat)

# Ubiquitin "RL40_HUMAN" 
pp <- "P62987"
test_sig %>% filter(str_detect(protsite, pp)) %>% arrange(protsite, ctrx, hyp)
df_fullpep %>% filter(is_mod) %>% plot_sprofile(pp, runlvl_bat)
df_fullpep %>% semi_join(filter(test_sig, uniprot_iso == pp)) %>% plot_sprofile(pp, runlvl_bat)

