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
    if (n_distinct(aftdata$feature) == 1) return(aftdata)  # only one feature
    if (all(aftdata$ind_obs == 1)) return(aftdata)  # no missing value
    # AFT model with effects of run and feature
    fit <- survreg(Surv(log2inty_aft, ind_obs, type = "left") ~ run + feature,
                   data = aftdata, dist = "gaussian")
    aftdata <- aftdata %>% 
        mutate(
            log2inty_pred = predict(fit), 
            log2inty_aft = ifelse(ind_obs == 0, log2inty_pred, log2inty_aft)
        ) %>% 
        select(-log2inty_pred)
    
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
        select(feature, run, ind_obs, log2inty_aft) %>% 
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
    ggplot(df_fill_mod, aes(run_fac, log2inty_aft, group = feature, color = site_str, alpha = ind_obs)) + 
        geom_point(size = 2) +
        geom_line() + 
        geom_point(data = df_fill_unmod, aes(run_fac, log2inty_aft, alpha = ind_obs), color = "gray", size = 2) + 
        geom_line(data = df_fill_unmod, aes(run_fac, log2inty_aft, group = feature), color = "gray") + 
        geom_vline(xintercept = 8.5) + 
        scale_alpha_manual(values = c("0" = 0.4, "1" = 0.9), guide = FALSE) + 
        facet_grid(is_mod_fac ~ .) + 
        coord_cartesian(ylim = c(10, 35)) + 
        labs(x = "Run", y = "Log2-intensity", title = protein) + 
        theme_bw() + 
        guides(color = guide_legend(nrow = 1, title = NULL)) + 
        theme(legend.position = c(0.5, 0.065)) + 
        theme(axis.text.x = element_blank())
}


