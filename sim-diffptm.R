
# Load required packages & functions --------------------------------------

library(tidyverse)
library(broom)
library(stringr)

source("utilities.R")


# Functions generating data & testing -------------------------------------

# generate simulation data for site-level modification
gen_sitedata <- function(d_bch, std_bch, d_cnd, nb_sim, nb_rep, log2level = 25) {
    # Error checking
    if (length(d_bch) != length(std_bch)) 
        stop("d_bch & std_bch imply different # of batches")
    
    nb_bch <- length(d_bch)
    mu0 <- log2level + d_bch
    mu1 <- mu0 + d_cnd
    
    site_samp <- vector("list", nb_sim)
    for (s in 1:nb_sim) {
        samp1 <- tibble(
            log2inty = rnorm(nb_bch * nb_rep, mean = rep(mu0, nb_rep), sd = rep(std_bch, nb_rep)), 
            idx_sim = rep(s, nb_bch * nb_rep), 
            group = rep("G0", nb_bch * nb_rep), 
            batch = rep(paste0("B", 1:nb_bch), nb_rep), 
            tech = rep(paste0("T", 1:nb_rep), each = nb_bch), 
            run = paste(group, paste0(batch, tech), sep = "_")
        )
        samp2 <- tibble(
            log2inty = rnorm(nb_bch * nb_rep, mean = rep(mu1, nb_rep), sd = rep(std_bch, nb_rep)), 
            idx_sim = rep(s, nb_bch * nb_rep), 
            group = rep("G1", nb_bch * nb_rep), 
            batch = rep(paste0("B", 1:nb_bch), nb_rep), 
            tech = rep(paste0("T", 1:nb_rep), each = nb_bch), 
            run = paste(group, paste0(batch, tech), sep = "_")
        )
        site_samp[[s]] <- bind_rows(samp1, samp2)
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


# No (or same) changes across conditions & batches ------------------------

nb_sim <- 1000
nb_bch <- 2
s_bch1 <- 0.2

params_rep <- c(2, 3, 5)
params_dbch <- c(0, 1)
params_dcnd <- c(0, 0.5, 1, 2)
params_s_bch2 <- c(0.2, 0.5)

same_change_res <- vector("list", length(params_rep) * length(params_dbch) * 
                              length(params_dcnd) * length(params_s_bch2))
idx_res <- 1
set.seed(2017)
for (b in seq_along(params_dbch)) {
    for (d in seq_along(params_dcnd)) {
        for (s in seq_along(params_s_bch2)) {
            for (r in seq_along(params_rep)) {
                cat("param set", idx_res, "\n")
                
                bch_effect <- c(0, params_dbch[b])
                cnd_bch <- rep(params_dcnd[d], nb_bch)
                s_bch <- c(s_bch1, params_s_bch2[s])
                same_change <- gen_sitedata(bch_effect, s_bch, cnd_bch, nb_sim, params_rep[r])
                same_change_res[[idx_res]] <- same_change %>% 
                    alltest_sim() %>% 
                    mutate(
                        bchLFC = params_dbch[b], trueLFC = params_dcnd[d], 
                        nb_rep = params_rep[r], sd_bch2 = params_s_bch2[s]
                    )
                
                idx_res <- idx_res + 1
            }
        }
    }
}

df_same_change <- bind_rows(same_change_res)

df_same_change %>% 
    ggplot(aes(x = factor(nb_rep), y = logFC - trueLFC, color = method)) + 
    geom_boxplot(position = position_dodge(0.75)) + 
    geom_hline(yintercept = 0) + 
    facet_grid(trueLFC ~ bchLFC + sd_bch2)

same_null <- df_same_change %>% 
    filter(trueLFC == 0) %>% 
    group_by(bchLFC, sd_bch2, nb_rep, method) %>% 
    summarise(fpr = sum(p_val < 0.05) / n()) %>% 
    ungroup()

same_alt <- df_same_change %>% 
    filter(trueLFC != 0) %>% 
    group_by(trueLFC, bchLFC, sd_bch2, nb_rep, method) %>% 
    summarise(power = sum(p_val < 0.05) / n()) %>% 
    ungroup()

same_null %>% ggplot(aes(method, fpr, color = factor(nb_rep))) + 
    geom_point(position = position_dodge(0.75)) + 
    geom_hline(yintercept = 0.05) + 
    facet_grid(. ~ bchLFC + sd_bch2)

same_alt %>% ggplot(aes(method, power, color = factor(nb_rep))) + 
    geom_point(position = position_dodge(0.75)) + 
    facet_grid(trueLFC ~ bchLFC + sd_bch2)

