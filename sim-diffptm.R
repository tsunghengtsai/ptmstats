
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
        # site_samp[[s]] <- bind_rows(samp0, samp1)
    }
    
    return(bind_rows(site_samp))
}


test_perbch_sim <- function(simdata) {
    # Fit per-batch model
    nested_perbch <- simdata %>% 
        group_by(idx_sim, batch) %>% 
        nest() %>% 
        mutate(lm_fit = map(data, lm_perbch)) %>% 
        mutate(
            param = map2(lm_fit, data, tidy_bch), 
            df_res = map_dbl(lm_fit, df.residual)
        )
    # Extract parameters
    param_perbch <- nested_perbch %>% 
        select(idx_sim, batch, param, df_res) %>% 
        unnest(param)
    # Model-based inference of log-difference
    diff_mod <- param_perbch %>% 
        filter(group %in% c("G0", "G1")) %>% 
        group_by(idx_sim, batch, df_res) %>% 
        summarise(log2fc = diff(estimate), se2 = sum(std.error ^ 2)) %>% 
        ungroup()
    res <- diff_mod %>% 
        group_by(idx_sim) %>% 
        summarise(
            logFC = mean(log2fc), SE = sqrt(sum(se2)) / n(), 
            DF = sum(se2) ^ 2 / sum(se2 ^ 2 / df_res), t_stat = logFC / SE, 
            p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)
        ) %>% 
        select(-SE)

    return(res)
}


test_allbch_sim <- function(simdata) {
    # Fit all-batch model
    nested_allbch <- simdata %>% 
        group_by(idx_sim) %>% 
        nest() %>% 
        mutate(lm_fit = map(data, lm_allbch)) %>% 
        mutate(
            param = map2(lm_fit, data, tidy_bch), 
            df_res = map_dbl(lm_fit, df.residual)
        )
    # Extract parameters
    param_allbch <- nested_allbch %>% 
        select(idx_sim, param, df_res) %>% 
        unnest(param)
    # Model-based inference of log-difference
    diff_mod <- param_allbch %>% 
        filter(group %in% c("G0", "G1")) %>% 
        group_by(idx_sim, df_res) %>% 
        summarise(log2fc = diff(estimate), se2 = sum(std.error ^ 2)) %>% 
        ungroup()
    res <- diff_mod %>% 
        mutate(
            logFC = log2fc, SE = sqrt(se2), DF = df_res, t_stat = logFC / SE, 
            p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)
        ) %>% 
        select(idx_sim, logFC, DF, t_stat, p_val)
    
    return(res)
}


# t-test across batches (ignore batch effects)
ttest_allbch_sim <- function(simdata) {
    # test across batches
    test_allbch <- simdata %>% 
        filter(group %in% c("G0", "G1")) %>% 
        group_by(idx_sim) %>% 
        nest() %>% 
        mutate(ttest = map(data, ~t.test(log2inty ~ group, data = .))) %>% 
        mutate(param = map(ttest, tidy))
    
    res <- test_allbch %>% 
        select(idx_sim, param) %>% 
        unnest(param) %>% 
        mutate(logFC = -estimate) %>% 
        select(idx_sim, logFC, t_stat = statistic, DF = parameter, p_val = p.value)
    
    return(res)
}


# Smallest p-values among batches from t-test (significant in ANY batches)
ttest_minpval_sim <- function(simdata) {
    # test within batches
    test_perbch <- simdata %>% 
        filter(group %in% c("G0", "G1")) %>% 
        group_by(idx_sim, batch) %>% 
        nest() %>% 
        mutate(ttest = map(data, ~t.test(log2inty ~ group, data = .))) %>% 
        mutate(param = map(ttest, tidy))
    # Take the most significant batch
    res <- test_perbch %>% 
        select(idx_sim, batch, param) %>% 
        unnest(param) %>% 
        group_by(idx_sim) %>% 
        filter(p.value == min(p.value)) %>% 
        ungroup() %>% 
        mutate(logFC = -estimate) %>% 
        select(idx_sim, logFC, t_stat = statistic, DF = parameter, p_val = p.value)
    
    return(res)
}


# Largest p-values among batches from t-test (significant in ALL batches)
ttest_maxpval_sim <- function(simdata) {
    # test within batches
    test_perbch <- simdata %>% 
        filter(group %in% c("G0", "G1")) %>% 
        group_by(idx_sim, batch) %>% 
        nest() %>% 
        mutate(ttest = map(data, ~t.test(log2inty ~ group, data = .))) %>% 
        mutate(param = map(ttest, tidy))
    # logFC based on average over batches; the rest takes the least significant batch
    res <- test_perbch %>% 
        select(idx_sim, batch, param) %>% 
        unnest(param) %>% 
        group_by(idx_sim) %>% 
        mutate(logFC = -mean(estimate)) %>% 
        filter(p.value == max(p.value)) %>% 
        ungroup() %>% 
        select(idx_sim, logFC, t_stat = statistic, DF = parameter, p_val = p.value)
    
    return(res)
}


alltest_sim <- function(simdata) {
    res <- bind_rows(
        simdata %>% ttest_allbch_sim() %>% mutate(method = "t-test (no batch)"),
        simdata %>% ttest_maxpval_sim() %>% mutate(method = "t-test (all batch)"),
        simdata %>% ttest_minpval_sim() %>% mutate(method = "t-test (any batch)"),
        simdata %>% test_perbch_sim() %>% mutate(method = "Proposed (per-batch)"),
        simdata %>% test_allbch_sim() %>% mutate(method = "Proposed (all-batch)")
    )
    
    return(res)
}



# First investigation: batch effects --------------------------------------

nb_sim <- 500
nb_bch <- 2
s_bch1 <- 0.2

# Parameters & forms of batch effects
params_rep <- c(2, 3, 5)  # sample size
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


# Summarize & visualize results for first investigation -------------------

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


