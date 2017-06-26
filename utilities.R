# Annotate observation status and add censored values
annot_obs <- function(data) {
    augdata <- data %>% 
        mutate(ind_obs = ifelse(is.na(log2inty), 0, 1)) %>% 
        group_by(feature) %>% 
        mutate(log2inty_aft = ifelse(ind_obs == 0, 0.99 * min(log2inty[ind_obs == 1]), log2inty)) %>% 
        ungroup()
    
    return(augdata)
}

# Fill censored values with AFT
fill_censored <- function(aftdata) {
    if (n_distinct(aftdata$feature) == 1) return(aftdata %>% select(-log2inty_aft))  # only one feature
    if (all(aftdata$ind_obs == 1)) return(aftdata %>% select(-log2inty_aft))  # no missing value
    # AFT model with effects of run and feature
    fit <- survreg(Surv(log2inty_aft, ind_obs, type = "left") ~ run + feature,
                   data = aftdata, dist = "gaussian")
    aftdata <- aftdata %>% 
        mutate(
            log2inty_pred = predict(fit), 
            log2inty = ifelse(ind_obs == 0, pmin(log2inty_pred, log2inty_aft), log2inty_aft)
        ) %>% 
        select(-log2inty_aft, -log2inty_pred)
    
    return(aftdata)
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


# Profile plots for site data 
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
        theme(legend.position = c(0.5, 0.065)) + 
        theme(axis.text.x = element_blank())
}


# Profile plots for site data with annotations for missing values
plot_sprofile_aft <- function(df_allprot, protein, run_level) {
    df_prot <- df_allprot %>% filter(uniprot_iso == protein)
    # Complete possible combinations of peptide features and runs
    mpar <- df_prot %>% distinct(site_str, feature, is_mod)
    df_fill <- df_prot %>% 
        select(feature, run, ind_obs, log2inty) %>% 
        complete(run = run_level, feature) %>% 
        left_join(mpar) %>% 
        mutate(
            is_mod_fac = factor(ifelse(is_mod, "Modified", "Unmodified")), 
            run_fac = factor(run, levels = run_level), 
            ind_obs = factor(ind_obs, levels = c("0", "1"))
        )
    # Feature profiles categorized by modification
    df_fill_mod <- filter(df_fill, is_mod)
    df_fill_unmod <- filter(df_fill, !is_mod)
    ggplot(df_fill_mod, aes(run_fac, log2inty, group = feature, color = site_str, alpha = ind_obs)) + 
        geom_point(size = 2) +
        geom_line() + 
        geom_point(data = df_fill_unmod, aes(run_fac, log2inty, alpha = ind_obs), color = "gray", size = 2) + 
        geom_line(data = df_fill_unmod, aes(run_fac, log2inty, group = feature), color = "gray") + 
        geom_vline(xintercept = 8.5) + 
        scale_alpha_manual(values = c("0" = 0.4, "1" = 0.9), guide = FALSE) + 
        facet_grid(is_mod_fac ~ .) + 
        coord_cartesian(ylim = c(5, 35)) + 
        labs(x = "Run", y = "Log2-intensity", title = protein) + 
        theme_bw() + 
        guides(color = guide_legend(nrow = 1, title = NULL)) + 
        theme(legend.position = c(0.5, 0.065)) + 
        # theme(axis.text.x = element_text(angle = 90, hjust = 1))
        theme(axis.text.x = element_blank())
}


# Linear model with group effect 
lm_perbch <- function(df_perbch) {
    if (n_distinct(df_perbch$group) == 1) {
        fit <- lm(log2inty ~ 1, data = df_perbch)
    } else {
        fit <- lm(log2inty ~ 0 + group, data = df_perbch)
    }
    
    return(fit)
}


# Linear model with group effect and batch effect
lm_allbch <- function(df_allbch) {
    if (n_distinct(df_allbch$batch) == 1) 
        stop("Cannot infer batch effect with a single batch!")
    
    if (n_distinct(df_allbch$group) == 1) {
        fit <- lm(log2inty ~ batch, data = df_allbch)
    } else {
        fit <- lm(log2inty ~ 0 + group + batch, data = df_allbch)
    }
    
    return(fit)
}


# Extract estimate of group effect
tidy_bch <- function(bch_fit, df_bch) {
    params <- tidy(bch_fit) %>% filter(!str_detect(term, "batch"))
    if (n_distinct(df_bch$group) == 1) {
        params <- params %>% 
            mutate(group = unique(df_bch$group)) %>% 
            select(-term, -statistic, -p.value)
    } else {
        params <- params %>% 
            mutate(group = str_replace(term, "group", "")) %>% 
            select(-term, -statistic, -p.value)
    }
    
    return(params)
}


# Testing for differential modification
test_mod <- function(df_mod, grp_ctrl, grp_case) {
    full_mod <- df_mod %>% 
        filter(group %in% c(grp_ctrl, grp_case)) %>% 
        filter(!is.na(df_res)) %>% 
        mutate(key_grp = ifelse(group == grp_ctrl, "G0", "G1")) %>% 
        complete(protsite, key_grp)
    # Model-based inference of log-difference
    diff_mod <- full_mod %>% 
        group_by(protsite) %>% 
        filter(!any(is.na(df_res))) %>% 
        group_by(protsite, df_res) %>% 
        summarise(log2fc = diff(estimate), se2 = sum(std.error ^ 2)) %>% 
        ungroup()
    # Missing in one group
    part_mod <- full_mod %>% anti_join(diff_mod)
    untest_mod <- part_mod %>% 
        filter(!is.na(df_res)) %>% 
        distinct(protsite, key_grp) %>% 
        mutate(
            logFC = ifelse(key_grp == "G1", Inf, -Inf), 
            SE = NA, DF = NA, t_stat = NA, p_val = NA, 
            ctrx = paste(grp_case, grp_ctrl, sep = " vs ")
        ) %>% select(-key_grp)
    res <- diff_mod %>% 
        mutate(logFC = log2fc, SE = sqrt(se2), DF = df_res, 
               t_stat = logFC / SE, p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)) %>% 
        select(protsite, logFC, SE, DF, t_stat, p_val) %>% 
        mutate(ctrx = paste(grp_case, grp_ctrl, sep = " vs ")) %>% 
        bind_rows(untest_mod)
    
    return(res)
}


# Testing for differential modification with adjustment
test_adjmod <- function(df_mod, grp_ctrl, grp_case) {
    full_mod <- df_mod %>% 
        filter(group %in% c(grp_ctrl, grp_case)) %>% 
        filter(!is.na(df_res)) %>% 
        mutate(key_grp = ifelse(group == grp_ctrl, "G0", "G1")) %>% 
        complete(protsite, key_grp)
    # Model-based inference of log-difference
    diff_mod <- full_mod %>% 
        group_by(protsite) %>% 
        filter(all(!is.na(df_res)), all(!is.na(df_unmod))) %>% 
        group_by(protsite, df_res, df_unmod) %>% 
        summarise(
            log2fc = diff(estimate), se2 = sum(std.error ^ 2), 
            log2fc_unmod = diff(est_unmod), se2_unmod = sum(se_unmod ^ 2)
        ) %>% ungroup()
    # Missing in one group
    part_mod <- full_mod %>% anti_join(diff_mod)
    untest_mod <- part_mod %>% 
        filter(!is.na(df_res), !is.na(df_unmod)) %>% 
        distinct(protsite, key_grp) %>% 
        mutate(
            logFC = ifelse(key_grp == "G1", Inf, -Inf), 
            SE = NA, DF = NA, t_stat = NA, p_val = NA, 
            ctrx = paste(grp_case, grp_ctrl, sep = " vs ")
        ) %>% select(-key_grp)
    res <- diff_mod %>% 
        mutate(logFC = log2fc - log2fc_unmod, SE = sqrt(se2 + se2_unmod), 
               DF = (se2 + se2_unmod) ^ 2 / (se2 ^ 2 / df_res + se2_unmod ^ 2 / df_unmod), 
               t_stat = logFC / SE, p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)) %>% 
        select(protsite, logFC, SE, DF, t_stat, p_val) %>% 
        mutate(ctrx = paste(grp_case, grp_ctrl, sep = " vs ")) %>% 
        bind_rows(untest_mod)
    
    return(res)
}


# Batch-aggregated testing for differential modification
test_mod_bch <- function(df_mod, grp_ctrl, grp_case) {
    full_mod <- df_mod %>% 
        filter(group %in% c(grp_ctrl, grp_case)) %>% 
        filter(!is.na(df_res)) %>% 
        mutate(key_grp = ifelse(group == grp_ctrl, "G0", "G1")) %>% 
        complete(protsite, batch, key_grp)
    # Model-based inference of log-difference
    diff_mod <- full_mod %>% 
        group_by(protsite) %>% 
        filter(all(!is.na(df_res))) %>% 
        group_by(protsite, batch, df_res) %>% 
        summarise(log2fc = diff(estimate), se2 = sum(std.error ^ 2)) %>% 
        ungroup()
    # Missing in one group
    part_mod <- full_mod %>% anti_join(diff_mod)
    untest_mod <- part_mod %>% 
        group_by(protsite, key_grp) %>% 
        filter(all(!is.na(df_res))) %>% 
        ungroup() %>% 
        distinct(protsite, key_grp) %>% 
        mutate(
            logFC = ifelse(key_grp == "G1", Inf, -Inf), 
            SE = NA, DF = NA, t_stat = NA, p_val = NA, 
            ctrx = paste(grp_case, grp_ctrl, sep = " vs ")
        ) %>% select(-key_grp)
    res <- diff_mod %>% 
        group_by(protsite) %>% 
        summarise(logFC = mean(log2fc), SE = sqrt(sum(se2)) / n(), 
                  DF = sum(se2) ^ 2 / sum(se2 ^ 2 / df_res), 
                  t_stat = logFC / SE, p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)) %>% 
        mutate(ctrx = paste(grp_case, grp_ctrl, sep = " vs ")) %>% 
        bind_rows(untest_mod)
    
    return(res)
}


# Batch-aggregated testing for differential modification with adjustment
test_adjmod_bch <- function(df_mod, grp_ctrl, grp_case) {
    full_mod <- df_mod %>% 
        filter(group %in% c(grp_ctrl, grp_case)) %>% 
        filter(!is.na(df_res)) %>% 
        mutate(key_grp = ifelse(group == grp_ctrl, "G0", "G1")) %>% 
        complete(protsite, batch, key_grp)
    # Model-based inference of log-difference
    diff_mod <- full_mod %>% 
        group_by(protsite) %>% 
        filter(all(!is.na(df_res)), all(!is.na(df_unmod))) %>% 
        group_by(protsite, batch, df_res, df_unmod) %>% 
        summarise(
            log2fc = diff(estimate), se2 = sum(std.error ^ 2),
            log2fc_unmod = diff(est_unmod), se2_unmod = sum(se_unmod ^ 2)
        ) %>% ungroup()
    # Missing in one group
    part_mod <- full_mod %>% anti_join(diff_mod)
    untest_mod <- part_mod %>% 
        group_by(protsite, key_grp) %>% 
        filter(all(!is.na(df_res)), all(!is.na(df_unmod))) %>% 
        ungroup() %>% 
        distinct(protsite, key_grp) %>% 
        mutate(
            logFC = ifelse(key_grp == "G1", Inf, -Inf), 
            SE = NA, DF = NA, t_stat = NA, p_val = NA, 
            ctrx = paste(grp_case, grp_ctrl, sep = " vs ")
        ) %>% select(-key_grp)
    res <- diff_mod %>% 
        group_by(protsite) %>% 
        summarise(
            logFC = mean(log2fc) - mean(log2fc_unmod), 
            SE = sqrt(sum(se2) + sum(se2_unmod)) / n(), 
            DF = (sum(se2) + sum(se2_unmod)) ^ 2 / sum(se2 ^ 2 / df_res + se2_unmod ^ 2 / df_unmod), 
            t_stat = logFC / SE, p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)
        ) %>% 
        mutate(ctrx = paste(grp_case, grp_ctrl, sep = " vs ")) %>% 
        bind_rows(untest_mod)
    
    return(res)
}

