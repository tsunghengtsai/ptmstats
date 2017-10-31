
# Load required packages & functions --------------------------------------

library(tidyverse)
library(broom)
library(stringr)

source("utilities.R")


# Functions generating data & testing -------------------------------------

# Generate simulation data for site-level modification
gen_sitedata <- function(d_bch, std_bch, d_cnd, nb_sim, nb_rep, nb_grp = 2, log2level = 25) {
    # Error checking
    if (length(d_bch) != length(std_bch)) 
        stop("d_bch & std_bch imply different # of batches")
    if (nb_grp < 2)
        stop("Minimum # of groups for group comparison is 2")
    
    nb_bch <- length(d_bch)
    mu0 <- log2level + d_bch
    mu1 <- mu0 + d_cnd
    
    site_samp <- vector("list", nb_sim)
    for (s in 1:nb_sim) {
        samp_0 <- tibble(
            log2inty = rnorm(nb_bch * nb_rep, mean = rep(mu0, nb_rep), sd = rep(std_bch, nb_rep)), 
            idx_sim = rep(s, nb_bch * nb_rep), 
            group = rep("G0", nb_bch * nb_rep), 
            batch = rep(paste0("B", 1:nb_bch), nb_rep), 
            tech = rep(paste0("T", 1:nb_rep), each = nb_bch), 
            run = paste(group, paste0(batch, tech), sep = "_")
        )
        samp_1 <- tibble(
            log2inty = rnorm(nb_bch * nb_rep, mean = rep(mu1, nb_rep), sd = rep(std_bch, nb_rep)), 
            idx_sim = rep(s, nb_bch * nb_rep), 
            group = rep("G1", nb_bch * nb_rep), 
            batch = rep(paste0("B", 1:nb_bch), nb_rep), 
            tech = rep(paste0("T", 1:nb_rep), each = nb_bch), 
            run = paste(group, paste0(batch, tech), sep = "_")
        )
        if (nb_grp == 2) {
            site_samp[[s]] <- bind_rows(samp_0, samp_1)
        } else {
            # Simulating additional groups of data on the same level as in G0
            add_grp <- nb_grp - 2
            samp_add <- tibble(
                log2inty = rnorm(
                    nb_bch * nb_rep * add_grp, 
                    mean = rep(rep(mu0, nb_rep), add_grp), 
                    sd = rep(rep(std_bch, nb_rep), add_grp)
                ), 
                idx_sim = rep(rep(s, nb_bch * nb_rep), add_grp),
                group = rep(paste0("G", 2:(nb_grp - 1)), each = nb_bch * nb_rep), 
                batch = rep(rep(paste0("B", 1:nb_bch), nb_rep), add_grp), 
                tech = rep(rep(paste0("T", 1:nb_rep), each = nb_bch), add_grp), 
                run = paste(group, paste0(batch, tech), sep = "_")
            )
            site_samp[[s]] <- bind_rows(samp_0, samp_1, samp_add)
        }
    }
    
    return(bind_rows(site_samp))
}


# Generate data for site-level modification & protein abundance
gen_sitenprot <- function(d_bch, std_bch, d_cnd, dprot_cnd, nb_sim, nb_rep, nb_grp = 2, log2level = 25) {
    # Error checking
    if (length(d_bch) != length(std_bch)) 
        stop("d_bch & std_bch imply different # of batches")
    if (nb_grp < 2)
        stop("Minimum # of groups for group comparison is 2")
    
    nb_bch <- length(d_bch)
    mu0 <- log2level + d_bch
    mu1 <- mu0 + d_cnd
    
    mu_prot0 <- rep(log2level, nb_bch)  # assume no batch effect for unmodified peptides
    mu_prot1 <- mu_prot0 + dprot_cnd
    std_prot <- rep(std_bch[1], nb_bch)
    
    site_samp <- prot_samp <- vector("list", nb_sim)
    for (s in 1:nb_sim) {
        # Site-level modification
        samp_0 <- tibble(
            log2inty = rnorm(nb_bch * nb_rep, mean = rep(mu0, nb_rep), sd = rep(std_bch, nb_rep)), 
            idx_sim = rep(s, nb_bch * nb_rep), 
            group = rep("G0", nb_bch * nb_rep), 
            batch = rep(paste0("B", 1:nb_bch), nb_rep), 
            tech = rep(paste0("T", 1:nb_rep), each = nb_bch), 
            run = paste(group, paste0(batch, tech), sep = "_")
        )
        samp_1 <- tibble(
            log2inty = rnorm(nb_bch * nb_rep, mean = rep(mu1, nb_rep), sd = rep(std_bch, nb_rep)), 
            idx_sim = rep(s, nb_bch * nb_rep), 
            group = rep("G1", nb_bch * nb_rep), 
            batch = rep(paste0("B", 1:nb_bch), nb_rep), 
            tech = rep(paste0("T", 1:nb_rep), each = nb_bch), 
            run = paste(group, paste0(batch, tech), sep = "_")
        )
        
        # Protein abundance
        prot_0 <- tibble(
            log2inty = rnorm(nb_bch * nb_rep, mean = rep(mu_prot0, nb_rep), sd = rep(std_prot, nb_rep)), 
            idx_sim = rep(s, nb_bch * nb_rep), 
            group = rep("G0", nb_bch * nb_rep), 
            batch = rep(paste0("B", 1:nb_bch), nb_rep), 
            tech = rep(paste0("T", 1:nb_rep), each = nb_bch), 
            run = paste(group, paste0(batch, tech), sep = "_")
        )
        prot_1 <- tibble(
            log2inty = rnorm(nb_bch * nb_rep, mean = rep(mu_prot1, nb_rep), sd = rep(std_prot, nb_rep)), 
            idx_sim = rep(s, nb_bch * nb_rep), 
            group = rep("G1", nb_bch * nb_rep), 
            batch = rep(paste0("B", 1:nb_bch), nb_rep), 
            tech = rep(paste0("T", 1:nb_rep), each = nb_bch), 
            run = paste(group, paste0(batch, tech), sep = "_")
        )

        if (nb_grp == 2) {
            site_samp[[s]] <- bind_rows(samp_0, samp_1)
            prot_samp[[s]] <- bind_rows(prot_0, prot_1)
        } else {
            # Simulating additional groups of data on the same level as in G0
            add_grp <- nb_grp - 2
            samp_add <- tibble(
                log2inty = rnorm(
                    nb_bch * nb_rep * add_grp, 
                    mean = rep(rep(mu0, nb_rep), add_grp), 
                    sd = rep(rep(std_bch, nb_rep), add_grp)
                ), 
                idx_sim = rep(rep(s, nb_bch * nb_rep), add_grp),
                group = rep(paste0("G", 2:(nb_grp - 1)), each = nb_bch * nb_rep), 
                batch = rep(rep(paste0("B", 1:nb_bch), nb_rep), add_grp), 
                tech = rep(rep(paste0("T", 1:nb_rep), each = nb_bch), add_grp), 
                run = paste(group, paste0(batch, tech), sep = "_")
            )
            site_samp[[s]] <- bind_rows(samp_0, samp_1, samp_add)
            
            prot_add <- tibble(
                log2inty = rnorm(
                    nb_bch * nb_rep * add_grp, 
                    mean = rep(rep(mu_prot0, nb_rep), add_grp), 
                    sd = rep(rep(std_prot, nb_rep), add_grp)
                ), 
                idx_sim = rep(rep(s, nb_bch * nb_rep), add_grp),
                group = rep(paste0("G", 2:(nb_grp - 1)), each = nb_bch * nb_rep), 
                batch = rep(rep(paste0("B", 1:nb_bch), nb_rep), add_grp), 
                tech = rep(rep(paste0("T", 1:nb_rep), each = nb_bch), add_grp), 
                run = paste(group, paste0(batch, tech), sep = "_")
            )
            prot_samp[[s]] <- bind_rows(prot_0, prot_1, prot_add)
        }
    }
    
    # Combine site-level & protein data
    full_samp <- bind_rows(
        bind_rows(site_samp) %>% mutate(site = "MOD"), 
        bind_rows(prot_samp) %>% mutate(site = "UNMOD")
    )
    
    return(full_samp)
}



# Generate two-group data with multiple sites per protein
gen_msites <- function(d_bch, std_bch, d_cnd, dprot_cnd, nb_sim, nb_rep, nb_site = 1, log2level = 25) {
    # Error checking
    if (length(d_bch) != length(std_bch)) 
        stop("d_bch & std_bch imply different # of batches")
    if (nb_site < 1)
        stop("Number of sites should be greater or equal to 1")

    nb_grp <- 2
    
    nb_bch <- length(d_bch)
    mu0 <- log2level + d_bch
    mu1 <- mu0 + d_cnd
    
    mu_prot0 <- rep(log2level, nb_bch)  # assume no batch effect for unmodified peptides
    mu_prot1 <- mu_prot0 + dprot_cnd
    std_prot <- rep(std_bch[1], nb_bch)
    
    site_samp <- prot_samp <- vector("list", nb_sim)
    for (s in 1:nb_sim) {
        # Site-level modification
        samp_0 <- samp_1 <- vector("list", nb_site)
        for (ss in 1:nb_site) {
            samp_0[[ss]] <- tibble(
                log2inty = rnorm(nb_bch * nb_rep, mean = rep(mu0, nb_rep), sd = rep(std_bch, nb_rep)), 
                idx_sim = rep(s, nb_bch * nb_rep), 
                group = rep("G0", nb_bch * nb_rep), 
                batch = rep(paste0("B", 1:nb_bch), nb_rep), 
                tech = rep(paste0("T", 1:nb_rep), each = nb_bch), 
                run = paste(group, paste0(batch, tech), sep = "_"), 
                site = rep(paste0("S", ss), nb_bch * nb_rep)
            )
            samp_1[[ss]] <- tibble(
                log2inty = rnorm(nb_bch * nb_rep, mean = rep(mu1, nb_rep), sd = rep(std_bch, nb_rep)), 
                idx_sim = rep(s, nb_bch * nb_rep), 
                group = rep("G1", nb_bch * nb_rep), 
                batch = rep(paste0("B", 1:nb_bch), nb_rep), 
                tech = rep(paste0("T", 1:nb_rep), each = nb_bch), 
                run = paste(group, paste0(batch, tech), sep = "_"), 
                site = rep(paste0("S", ss), nb_bch * nb_rep)
            )
        }

        # Protein abundance
        prot_0 <- tibble(
            log2inty = rnorm(nb_bch * nb_rep, mean = rep(mu_prot0, nb_rep), sd = rep(std_prot, nb_rep)), 
            idx_sim = rep(s, nb_bch * nb_rep), 
            group = rep("G0", nb_bch * nb_rep), 
            batch = rep(paste0("B", 1:nb_bch), nb_rep), 
            tech = rep(paste0("T", 1:nb_rep), each = nb_bch), 
            run = paste(group, paste0(batch, tech), sep = "_"), 
            site = rep("UNMOD", nb_bch * nb_rep)
        )
        prot_1 <- tibble(
            log2inty = rnorm(nb_bch * nb_rep, mean = rep(mu_prot1, nb_rep), sd = rep(std_prot, nb_rep)), 
            idx_sim = rep(s, nb_bch * nb_rep), 
            group = rep("G1", nb_bch * nb_rep), 
            batch = rep(paste0("B", 1:nb_bch), nb_rep), 
            tech = rep(paste0("T", 1:nb_rep), each = nb_bch), 
            run = paste(group, paste0(batch, tech), sep = "_"), 
            site = rep("UNMOD", nb_bch * nb_rep)
        )
        
        site_samp[[s]] <- bind_rows(bind_rows(samp_0), bind_rows(samp_1))
        prot_samp[[s]] <- bind_rows(prot_0, prot_1)
    }
    
    # Combine site-level & protein data
    full_samp <- bind_rows(bind_rows(site_samp), bind_rows(prot_samp))
    
    return(full_samp)
}


test_perbch_sim <- function(simdata) {
    # Fit per-batch model
    nested_perbch <- simdata %>% 
        group_by(idx_sim, site, batch) %>% 
        # group_by(idx_sim, batch) %>% 
        nest() %>% 
        mutate(lm_fit = map(data, lm_perbch)) %>% 
        mutate(
            param = map2(lm_fit, data, tidy_bch), 
            df_res = map_dbl(lm_fit, df.residual)
        )
    # Extract parameters
    param_perbch <- nested_perbch %>% 
        select(idx_sim, site, batch, param, df_res) %>% 
        # select(idx_sim, batch, param, df_res) %>% 
        unnest(param)
    # Model-based inference of log-difference
    diff_mod <- param_perbch %>% 
        filter(group %in% c("G0", "G1")) %>% 
        group_by(idx_sim, site, batch, df_res) %>% 
        # group_by(idx_sim, batch, df_res) %>% 
        summarise(log2fc = diff(estimate), se2 = sum(std.error ^ 2)) %>% 
        ungroup()
    res <- diff_mod %>% 
        group_by(idx_sim, site) %>% 
        # group_by(idx_sim) %>% 
        summarise(
            logFC = mean(log2fc), SE = sqrt(sum(se2)) / n(), 
            DF = sum(se2) ^ 2 / sum(se2 ^ 2 / df_res), t_stat = logFC / SE, 
            p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)
        ) %>% 
        ungroup() %>% select(-SE)
        # select(-SE)

    return(res)
}


test_allbch_sim <- function(simdata) {
    # Fit all-batch model
    nested_allbch <- simdata %>% 
        group_by(idx_sim, site) %>% 
        # group_by(idx_sim) %>% 
        nest() %>% 
        mutate(lm_fit = map(data, lm_allbch)) %>% 
        mutate(
            param = map2(lm_fit, data, tidy_bch), 
            df_res = map_dbl(lm_fit, df.residual)
        )
    # Extract parameters
    param_allbch <- nested_allbch %>% 
        select(idx_sim, site, param, df_res) %>% 
        # select(idx_sim, param, df_res) %>% 
        unnest(param)
    # Model-based inference of log-difference
    diff_mod <- param_allbch %>% 
        filter(group %in% c("G0", "G1")) %>% 
        group_by(idx_sim, site, df_res) %>% 
        # group_by(idx_sim, df_res) %>% 
        summarise(log2fc = diff(estimate), se2 = sum(std.error ^ 2)) %>% 
        ungroup()
    res <- diff_mod %>% 
        group_by(idx_sim, site) %>% 
        mutate(
            logFC = log2fc, SE = sqrt(se2), DF = df_res, t_stat = logFC / SE, 
            p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)
        ) %>% 
        ungroup() %>% select(idx_sim, site, logFC, DF, t_stat, p_val)
        # select(idx_sim, logFC, DF, t_stat, p_val)
    
    return(res)
}


test_adjperbch_sim <- function(simdata) {
    # Fit per-batch model
    nested_perbch <- simdata %>% 
        group_by(idx_sim, site, batch) %>% 
        nest() %>% 
        mutate(lm_fit = map(data, lm_perbch)) %>% 
        mutate(
            param = map2(lm_fit, data, tidy_bch), 
            df_res = map_dbl(lm_fit, df.residual)
        )
    # Extract parameters
    nested_param <- nested_perbch %>% 
        select(idx_sim, site, batch, param, df_res)
    param_mod <- nested_param %>% 
        filter(site != "UNMOD") %>% 
        unnest(param)
    param_unmod <- nested_param %>% 
        filter(site == "UNMOD") %>% 
        select(-site) %>% 
        unnest(param) %>% 
        rename(df_unmod = df_res, est_unmod = estimate, se_unmod = std.error)
    param_perbch <- param_mod %>% left_join(param_unmod)
    # Model-based inference of log-difference
    diff_mod <- param_perbch %>% 
        filter(group %in% c("G0", "G1")) %>% 
        group_by(idx_sim, site, batch, df_res, df_unmod) %>% 
        # group_by(idx_sim, batch, df_res, df_unmod) %>% 
        summarise(
            log2fc = diff(estimate), se2 = sum(std.error ^ 2),
            log2fc_unmod = diff(est_unmod), se2_unmod = sum(se_unmod ^ 2)
        ) %>% 
        ungroup()
    res <- diff_mod %>% 
        group_by(idx_sim, site) %>% 
        # group_by(idx_sim) %>% 
        summarise(
            logFC = mean(log2fc) - mean(log2fc_unmod), 
            SE = sqrt(sum(se2) + sum(se2_unmod)) / n(), 
            DF = (sum(se2) + sum(se2_unmod)) ^ 2 / sum(se2 ^ 2 / df_res + se2_unmod ^ 2 / df_unmod), 
            t_stat = logFC / SE, p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)
        ) %>% 
        ungroup() %>% select(-SE)
        # select(-SE)

    return(res)
}

test_adjallbch_sim <- function(simdata) {
    # Fit all-batch model
    nested_allbch <- simdata %>% 
        group_by(idx_sim, site) %>% 
        nest() %>% 
        mutate(lm_fit = map(data, lm_allbch)) %>% 
        mutate(
            param = map2(lm_fit, data, tidy_bch), 
            df_res = map_dbl(lm_fit, df.residual)
        )
    # Extract parameters
    nested_param <- nested_allbch %>% 
        select(idx_sim, site, param, df_res)
    param_mod <- nested_param %>% 
        filter(site != "UNMOD") %>% 
        unnest(param)
    param_unmod <- nested_param %>% 
        filter(site == "UNMOD") %>% 
        select(-site) %>% 
        unnest(param) %>% 
        rename(df_unmod = df_res, est_unmod = estimate, se_unmod = std.error)
    param_allbch <- param_mod %>% left_join(param_unmod)
    # Model-based inference of log-difference
    diff_mod <- param_allbch %>% 
        filter(group %in% c("G0", "G1")) %>% 
        group_by(idx_sim, site, df_res, df_unmod) %>% 
        # group_by(idx_sim, df_res, df_unmod) %>% 
        summarise(
            log2fc = diff(estimate), se2 = sum(std.error ^ 2), 
            log2fc_unmod = diff(est_unmod), se2_unmod = sum(se_unmod ^ 2)
        ) %>% 
        ungroup()
    res <- diff_mod %>% 
        group_by(idx_sim, site) %>% 
        mutate(logFC = log2fc - log2fc_unmod, SE = sqrt(se2 + se2_unmod), 
               DF = (se2 + se2_unmod) ^ 2 / (se2 ^ 2 / df_res + se2_unmod ^ 2 / df_unmod), 
               t_stat = logFC / SE, p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)) %>% 
        select(idx_sim, site, logFC, DF, t_stat, p_val) %>% 
        ungroup()
        # select(idx_sim, logFC, DF, t_stat, p_val)
    
    return(res)
}



# t-test across batches (ignore batch effects)
ttest_allbch_sim <- function(simdata) {
    # test across batches
    test_allbch <- simdata %>% 
        filter(group %in% c("G0", "G1")) %>% 
        group_by(idx_sim, site) %>% 
        # group_by(idx_sim) %>% 
        nest() %>% 
        mutate(ttest = map(data, ~t.test(log2inty ~ group, data = .))) %>% 
        mutate(param = map(ttest, tidy))
    
    res <- test_allbch %>% 
        select(idx_sim, site, param) %>% 
        # select(idx_sim, param) %>% 
        unnest(param) %>% 
        mutate(logFC = -estimate) %>% 
        select(idx_sim, site, logFC, t_stat = statistic, DF = parameter, p_val = p.value)
        # select(idx_sim, logFC, t_stat = statistic, DF = parameter, p_val = p.value)
    
    return(res)
}


# Smallest p-values among batches from t-test (significant in ANY batches)
ttest_minpval_sim <- function(simdata) {
    # test within batches
    test_perbch <- simdata %>% 
        filter(group %in% c("G0", "G1")) %>% 
        group_by(idx_sim, site, batch) %>% 
        # group_by(idx_sim, batch) %>% 
        nest() %>% 
        mutate(ttest = map(data, ~t.test(log2inty ~ group, data = .))) %>% 
        mutate(param = map(ttest, tidy))
    # Take the most significant batch
    res <- test_perbch %>% 
        select(idx_sim, site, batch, param) %>% 
        # select(idx_sim, batch, param) %>% 
        unnest(param) %>% 
        group_by(idx_sim, site) %>% 
        # group_by(idx_sim) %>% 
        filter(p.value == min(p.value)) %>% 
        ungroup() %>% 
        mutate(logFC = -estimate) %>% 
        select(idx_sim, site, logFC, t_stat = statistic, DF = parameter, p_val = p.value)
        # select(idx_sim, logFC, t_stat = statistic, DF = parameter, p_val = p.value)
    
    return(res)
}


# Largest p-values among batches from t-test (significant in ALL batches)
ttest_maxpval_sim <- function(simdata) {
    # test within batches
    test_perbch <- simdata %>% 
        filter(group %in% c("G0", "G1")) %>% 
        group_by(idx_sim, site, batch) %>% 
        # group_by(idx_sim, batch) %>% 
        nest() %>% 
        mutate(ttest = map(data, ~t.test(log2inty ~ group, data = .))) %>% 
        mutate(param = map(ttest, tidy))
    # logFC based on average over batches; the rest takes the least significant batch
    res <- test_perbch %>% 
        select(idx_sim, site, batch, param) %>% 
        # select(idx_sim, batch, param) %>% 
        unnest(param) %>% 
        group_by(idx_sim, site) %>% 
        # group_by(idx_sim) %>% 
        mutate(logFC = -mean(estimate)) %>% 
        filter(p.value == max(p.value)) %>% 
        ungroup() %>% 
        select(idx_sim, site, logFC, t_stat = statistic, DF = parameter, p_val = p.value)
        # select(idx_sim, logFC, t_stat = statistic, DF = parameter, p_val = p.value)
    
    return(res)
}


alltest_sim <- function(simdata) {
    simdata <- simdata %>% mutate(site = "MOD")
    res <- bind_rows(
        simdata %>% ttest_allbch_sim() %>% mutate(method = "t-test (no batch)"),
        simdata %>% ttest_maxpval_sim() %>% mutate(method = "t-test (all batch)"),
        simdata %>% ttest_minpval_sim() %>% mutate(method = "t-test (any batch)"),
        simdata %>% test_perbch_sim() %>% mutate(method = "Proposed (per-batch)"),
        simdata %>% test_allbch_sim() %>% mutate(method = "Proposed (all-batch)")
    )
    
    return(res)
}


alladjtest_sim <- function(simdata) {
    res <- bind_rows(
        simdata %>% test_adjperbch_sim() %>% mutate(method = "Proposed (adj. per-batch)"),
        simdata %>% test_adjallbch_sim() %>% mutate(method = "Proposed (adj. all-batch)"),
        simdata %>% filter(site != "UNMOD") %>% test_perbch_sim() %>% mutate(method = "Proposed (per-batch)"),
        simdata %>% filter(site != "UNMOD") %>% test_allbch_sim() %>% mutate(method = "Proposed (all-batch)")
    )
    
    return(res)
}


# Simulation 1: batch effects ---------------------------------------------

nb_sim <- 500
nb_bch <- 2
s_bch1 <- 0.2

# Parameters & forms of batch effects
params_rep <- c(2, 3, 5)  # number of replicates
params_ngrp <- c(2, 3, 4)  # number of conditions
params_dgrp <- c(0, 0.5, 0.75, 1)  # change between conditions
params_dbch <- c(0, 1, 2)  # difference in signal intensities
params_s_bch2 <- c(0.2, 0.3)  # difference in variability

nb_paramset <- length(params_rep) * length(params_dbch) * length(params_dgrp) * 
    length(params_ngrp) * length(params_s_bch2)
syn_null_res <- syn_pos_res <- syn_neg_res <- vector("list", nb_paramset)

idx_res <- 1
set.seed(2017)
for (b in seq_along(params_dbch)) {
    for (d in seq_along(params_dgrp)) {
        for (s in seq_along(params_s_bch2)) {
            for (n in seq_along(params_ngrp)) {
                for (r in seq_along(params_rep)) {
                    cat("param set", idx_res, "out of", nb_paramset, "\n")
                    
                    # Difference across batches
                    bch_effect <- c(0, params_dbch[b])
                    s_bch <- c(s_bch1, params_s_bch2[s])
                    
                    # Interaction between batch and condition
                    cnd_bch_null <- rep(params_dgrp[d], nb_bch)
                    cnd_bch_poss <- params_dgrp[d] * c(1, 1 + 0.25)
                    cnd_bch_negs <- params_dgrp[d] * c(1, 1 - 0.25)
                    
                    # Generate simulation data & evaluate methods
                    null_change <- gen_sitedata(bch_effect, s_bch, cnd_bch_null, nb_sim, params_rep[r], params_ngrp[n])
                    poss_change <- gen_sitedata(bch_effect, s_bch, cnd_bch_poss, nb_sim, params_rep[r], params_ngrp[n])
                    negs_change <- gen_sitedata(bch_effect, s_bch, cnd_bch_negs, nb_sim, params_rep[r], params_ngrp[n])
                    
                    syn_null_res[[idx_res]] <- null_change %>% 
                        alltest_sim() %>% 
                        mutate(
                            bchLFC = params_dbch[b], trueLFC = params_dgrp[d], 
                            nb_rep = params_rep[r], nb_grp = params_ngrp[n], 
                            sd_bch2 = params_s_bch2[s]
                        )
                    syn_pos_res[[idx_res]] <- poss_change %>% 
                        alltest_sim() %>% 
                        mutate(
                            bchLFC = params_dbch[b], trueLFC = params_dgrp[d], 
                            nb_rep = params_rep[r], nb_grp = params_ngrp[n], 
                            sd_bch2 = params_s_bch2[s]
                        )
                    syn_neg_res[[idx_res]] <- negs_change %>% 
                        alltest_sim() %>% 
                        mutate(
                            bchLFC = params_dbch[b], trueLFC = params_dgrp[d], 
                            nb_rep = params_rep[r], nb_grp = params_ngrp[n], 
                            sd_bch2 = params_s_bch2[s]
                        )
                    
                    idx_res <- idx_res + 1
                }
            }
        }
    }
}
syn_null_res <- bind_rows(syn_null_res)
syn_pos_res <- bind_rows(syn_pos_res)
syn_neg_res <- bind_rows(syn_neg_res)

saveRDS(syn_null_res, "sim-synnull-171011.RDS")
saveRDS(syn_pos_res, "sim-synpos-171011.RDS")
saveRDS(syn_neg_res, "sim-synneg-171011.RDS")


# Simulation 2a: protein-level adjustment ---------------------------------

nb_sim <- 250
nb_bch <- 2
s_bch1 <- 0.2

# Parameters
params_rep <- c(2, 3, 5)  # number of replicates
params_ngrp <- c(2, 3, 4)  # number of conditions
params_dgrp <- c(0, 0.5)  # change between conditions
params_dprot <- c(0, 0.5)  # change in protein abundance
params_s_bch2 <- c(0.2, 0.3)  # difference in variability

nb_paramset <- length(params_rep) * length(params_dprot) * length(params_dgrp) * 
    length(params_ngrp) * length(params_s_bch2)
prot_res <- vector("list", nb_paramset)

idx_res <- 1
set.seed(2017)
for (b in seq_along(params_dprot)) {
    for (d in seq_along(params_dgrp)) {
        for (s in seq_along(params_s_bch2)) {
            for (n in seq_along(params_ngrp)) {
                for (r in seq_along(params_rep)) {
                    cat("param set", idx_res, "out of", nb_paramset, "\n")
                    
                    # Parameters (no shift across batches, no interaction)
                    bch_effect <- c(0, nb_bch)  # [TODO]: this should be c(0, 0)
                    s_bch <- c(s_bch1, params_s_bch2[s])
                    cnd_bch <- rep(params_dgrp[d], nb_bch)
                    prot_bch <- rep(params_dprot[b], nb_bch)
                    
                    # Generate simulation data & evaluate methods
                    protdata <- gen_sitenprot(bch_effect, s_bch, cnd_bch, prot_bch, nb_sim, params_rep[r], params_ngrp[n])
                    prot_res[[idx_res]] <- protdata %>% 
                        alladjtest_sim() %>% 
                        mutate(
                            protLFC = params_dprot[b], trueLFC = params_dgrp[d], 
                            nb_rep = params_rep[r], nb_grp = params_ngrp[n], 
                            sd_bch2 = params_s_bch2[s]
                        )

                    idx_res <- idx_res + 1
                }
            }
        }
    }
}
prot_res <- bind_rows(prot_res)

saveRDS(prot_res, "sim-prot-171016.RDS")



# Simulation 2b: FDR with protein-level adjustment ------------------------

nb_sim <- 1000
nb_bch <- 2
nb_grp <- 2
s_bch1 <- 0.2

# Parameters
params_rep <- c(2, 3, 5)  # number of replicates
params_nst <- c(1, 3, 7)  # number of sites
params_dgrp <- c(0.5, 0.75, 1)  # change between conditions
params_pnl <- c(0.5, 0.8)

nb_paramset <- length(params_rep) * length(params_nst) * length(params_dgrp) * length(params_pnl)
mix_res <- vector("list", nb_paramset)

idx_res <- 1
set.seed(2017)
for (p in seq_along(params_pnl)) {
    p_null <- params_pnl[p]
    nb_null <- as.integer(nb_sim * p_null)
    nb_chg <- nb_sim - nb_null
    
    for (d in seq_along(params_dgrp)) {
        for (s in seq_along(params_nst)) {
            for (r in seq_along(params_rep)) {
                cat("param set", idx_res, "out of", nb_paramset, "\n")
                
                # Parameters (no shift across batches, no interaction)
                bch_effect <- rep(0, nb_bch)
                s_bch <- rep(s_bch1, nb_bch)
                cnd_bch <- rep(params_dgrp[d], nb_bch)
                cnd_null <- rep(0, nb_bch)
                prot_bch <- rep(0, nb_bch)
                
                # Generate simulation data & evaluate methods
                protdata_chg <- gen_msites(bch_effect, s_bch, cnd_bch, prot_bch, nb_chg, params_rep[r], params_nst[s])
                protdata_null <- gen_msites(bch_effect, s_bch, cnd_null, prot_bch, nb_null, params_rep[r], params_nst[s])
                protdata <- bind_rows(protdata_chg, protdata_null %>% mutate(idx_sim = idx_sim + nb_chg))
                
                mix_res[[idx_res]] <- protdata %>% 
                    alladjtest_sim() %>% 
                    mutate(
                        p_null = params_pnl[p], trueLFC = params_dgrp[d], 
                        nb_rep = params_rep[r], nb_site = params_nst[s]
                    )

                idx_res <- idx_res + 1
            }
        }
    }
}
mix_res <- bind_rows(mix_res)

saveRDS(mix_res, "sim-mix-171031.RDS")


# Summarize & visualize results for Simulation 2b -------------------------

mix_res <- mix_res %>% 
    mutate(is_de = idx_sim <= nb_sim * (1 - p_null)) %>% 
    group_by(p_null, trueLFC, nb_rep, nb_site, method) %>% 
    mutate(p_adj = p.adjust(p_val, method = "BH")) %>%
    ungroup()

summary_mix <- mix_res %>% 
    group_by(p_null, trueLFC, nb_rep, nb_site, method) %>% 
    summarise(
        n_pos = sum(is_de), n_neg = sum(!is_de), n_ppos = sum(p_adj <= 0.05), n_tp = sum(p_adj <= 0.05 & is_de), 
        n_tn = sum(p_adj > 0.05 & !is_de), n_fp = sum(p_adj <= 0.05 & !is_de)
    ) %>% 
    mutate(
        sens = n_tp / n_pos, spec = n_tn / n_neg, 
        ppv = n_tp / (n_tp + n_fp), fdr = n_fp / (n_tp + n_fp)
    ) %>% 
    ungroup()


summary_mix %>% 
    filter(str_detect(method, "adj")) %>% 
    mutate(nb_site = ifelse(nb_site == 1, paste(nb_site, "site"), paste(nb_site, "sites"))) %>% 
    mutate(trueLFC = paste("Log2-FC of", trueLFC)) %>% 
    ggplot(aes(fdr, sens, color = method, shape = factor(p_null))) + 
    geom_point(aes(size = factor(nb_rep)), alpha = 0.8) + 
    geom_vline(xintercept = 0.05) + 
    facet_grid(trueLFC ~ nb_site) + 
    labs(title = "Proposed methods with protein-level adjustment were able to controll FDR across various cases", x = "FDR", y = "Power") +
    theme(legend.position = "bottom", legend.box = "horizontal")
ggsave("prot_fdr.png", width = 9, height = 6)


mix_res %>% 
    filter(p_null == 0.5) %>% 
    filter(nb_site == 3) %>% 
    filter(method == "Proposed (adj. per-batch)") %>% 
    ggplot(aes(p_val, fill = (is_de), color = (is_de))) + 
    geom_histogram(binwidth = 0.025) + 
    geom_vline(xintercept = 0.05) + 
    facet_grid(trueLFC ~ nb_rep)


# Summarize & visualize results for Simulation 2a -------------------------

prot_res %>% 
    filter(protLFC == 0, trueLFC == 0)

null_res <- prot_res %>% 
    filter(protLFC == trueLFC) %>% 
    group_by(protLFC, trueLFC, sd_bch2, nb_rep, nb_grp, method) %>% 
    summarise(fpr = sum(p_val < 0.05) / n()) %>% 
    ungroup()

alt_res <- prot_res %>% 
    filter(protLFC != trueLFC) %>% 
    group_by(protLFC, trueLFC, sd_bch2, nb_rep, nb_grp, method) %>% 
    summarise(power = sum(p_val < 0.05) / n()) %>% 
    ungroup()

null_res <- null_res %>% 
    mutate(
        prot_lvl = paste0("Protein level change: ", protLFC), 
        bch_sd = paste0("Batch sd (0.2, ", sd_bch2, ")"), 
        str_grp = paste0(nb_grp, " conditions"), 
        n_rep = factor(nb_rep)
    )

null_res %>% 
    ggplot(aes(method, fpr, size = n_rep, color = method)) + 
    geom_point(alpha = 0.7) + 
    facet_grid(str_grp ~ prot_lvl + bch_sd) + 
    geom_hline(yintercept = 0.05, linetype = "dashed") + 
    labs(x = "Method", y = "FPR", 
         title = "Analysis without accounting for changes in protein abundance led to wrong conclusions") + 
    theme(axis.text.x = element_blank(), legend.position = "bottom", legend.box = "horizontal")
ggsave("prot_fpr.png", width = 9, height = 6)


alt_res <- alt_res %>% 
    mutate(
        prot_lvl = paste0("Protein level change: ", protLFC), 
        bch_sd = paste0("Batch sd (0.2, ", sd_bch2, ")"), 
        str_grp = paste0(nb_grp, " conditions"), 
        n_rep = factor(nb_rep)
    )

alt_res %>% 
    ggplot(aes(trueLFC - protLFC, power, color = method, shape = n_rep)) + 
    # geom_line(aes(linetype = n_rep)) +
    geom_point() + 
    facet_grid(str_grp ~ prot_lvl + sd_bch2) + 
    scale_x_continuous(breaks = c(-0.5, 0.5), limits = c(-0.45, 0.55)) + 
    labs(title = "Proposed models improved power", x = "Adjusted log (base 2) fold change", y = "Power") +
    theme(legend.position = "bottom", legend.box = "horizontal")

alt_res %>% 
    filter(nb_grp == 2) %>% 
    ggplot(aes(trueLFC - protLFC, power, color = method)) + 
    geom_point() + 
    facet_grid(n_rep ~ prot_lvl + sd_bch2) + 
    scale_x_continuous(breaks = c(-0.5, 0.5), limits = c(-0.45, 0.55)) + 
    labs(title = "Proposed models improved power", x = "Adjusted log (base 2) fold change", y = "Power") +
    theme(legend.position = "bottom", legend.box = "horizontal")


alt_res %>% 
    filter(nb_grp == 2) %>% 
    ggplot(aes(method, power, fill = method)) + 
    geom_bar(stat = "identity") + 
    facet_grid(n_rep ~ prot_lvl + bch_sd) + 
    labs(x = "Method", y = "Power") +
    theme(axis.text.x = element_blank(), legend.position = "bottom", legend.box = "horizontal")
ggsave("prot_pwr.png", width = 9, height = 6)


# Summarize & visualize results for Simulation 1 --------------------------

sim_res <- readRDS("sim-synnull-171011.RDS")
sim_res <- readRDS("sim-synpos-171011.RDS")
sim_res <- readRDS("sim-synneg-171011.RDS")

null_res <- sim_res %>% 
    filter(trueLFC == 0) %>% 
    group_by(bchLFC, sd_bch2, nb_rep, nb_grp, method) %>% 
    summarise(fpr = sum(p_val < 0.05) / n()) %>% 
    ungroup()

alt_res <- sim_res %>% 
    filter(trueLFC != 0) %>% 
    group_by(trueLFC, bchLFC, sd_bch2, nb_rep, nb_grp, method) %>% 
    summarise(power = sum(p_val < 0.05) / n()) %>% 
    ungroup()

# Estimation
sim_res <- sim_res %>% 
    mutate(
        bch_lvl = paste0("Batch level (0, ", bchLFC, ")"), 
        bch_sd = paste0("Batch sd (0.2, ", sd_bch2, ")"), 
        str_lfc = paste0("Log2-FC of ", trueLFC)
    )

sim_res %>% 
    filter(nb_grp == 2) %>% 
    filter(bchLFC != 2, trueLFC != 1) %>% 
    ggplot(aes(x = factor(nb_rep), y = logFC - trueLFC, color = method)) + 
    geom_boxplot(position = position_dodge(0.75)) + 
    geom_hline(yintercept = 0) + 
    facet_grid(str_lfc ~ bch_lvl + bch_sd) + 
    labs(x = "Method", y = "Estimation error", subtitle = "No synergy, 2 conditions", 
         title = "Estimation based on most significant batch with t-test was biased") + 
    theme(legend.position = "bottom", legend.box = "horizontal")
ggsave("synnull_est_2.png", width = 9, height = 6)

null_res <- null_res %>% 
    mutate(
        bch_lvl = paste0("Batch level (0, ", bchLFC, ")"), 
        bch_sd = paste0("Batch sd (0.2, ", sd_bch2, ")"), 
        str_grp = paste0(nb_grp, " conditions"), 
        n_rep = factor(nb_rep)
    )

null_res %>% 
    ggplot(aes(method, fpr, size = n_rep, color = method)) + 
    geom_point(alpha = 0.7) + 
    facet_grid(str_grp ~ bch_lvl + bch_sd) + 
    geom_hline(yintercept = 0.05, linetype = "dashed") + 
    labs(x = "Method", y = "FPR", subtitle = "No synergy", 
         title = "Proposed per-batch model better calibrated Type I error rate") + 
    theme(axis.text.x = element_blank(), legend.position = "bottom", legend.box = "horizontal")
ggsave("synnull_fpr.png", width = 9, height = 6)
# ggsave("synpos_fpr.png", width = 9, height = 6)
# ggsave("synneg_fpr.png", width = 9, height = 6)


alt_res <- alt_res %>% 
    mutate(
        bch_lvl = paste0("Batch level (0, ", bchLFC, ")"), 
        bch_sd = paste0("Batch sd (0.2, ", sd_bch2, ")"), 
        str_grp = paste0(nb_grp, " conditions")
    )

nrep <- 2
# sttl <- paste0("No synergy, ", nrep, " replicates")
# sttl <- paste0("Positive synergy, ", nrep, " replicates")
sttl <- paste0("Negative synergy, ", nrep, " replicates")
alt_res %>% 
    filter(nb_rep == nrep) %>% 
    ggplot(aes(trueLFC, power, group = method, color = method)) + 
    geom_line() + 
    geom_point() + 
    facet_grid(str_grp ~ bch_lvl + bch_sd) + 
    scale_x_continuous(breaks = c(0.5, 0.75, 1), limits = c(0.45, 1.05)) + 
    labs(title = "Proposed models improved power", x = "Log (base 2) fold change", y = "Power", subtitle = sttl) +
    # labs(title = "Ignoring batch effects lost power", x = "Log (base 2) fold change", y = "Power", subtitle = sttl) +
    theme(legend.position = "bottom", legend.box = "horizontal")
# ggsave(paste0("synnull_pwr_", nrep,".png"), width = 9, height = 6)
# ggsave(paste0("synpos_pwr_", nrep,".png"), width = 9, height = 6)
ggsave(paste0("synneg_pwr_", nrep,".png"), width = 9, height = 6)


