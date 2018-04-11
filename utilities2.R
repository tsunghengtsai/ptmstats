# Sequence annotation -----------------------------------------------------

# Annotate modification site
annotate_site <- function(index, index_full, residue) {
    if (is_empty(index)) {
        ant_str <- "None"
    } else {
        ant_len <- max(str_length(index_full))
        ant_idx <- str_pad(index, width = ant_len, pad = "0")
        ant_str <- str_c(residue, ant_idx, collapse = "-")
    }
    
    return(ant_str)
}


# Locate potential modification sites
locate_site <- function(pep_seq, pep_idx, mod_res) {
    pep_start <- unname(pep_idx[, "start"])
    site_relidx <- str_locate_all(pep_seq, mod_res)[[1]]
    site_relstart <- unname(site_relidx[, "start"])
    site_idx <- site_relstart + pep_start - 1
    
    return(site_idx)
}


# Locate modified sites
locate_mod <- function(pep_seq, pep_idx, mod_res_symbol) {
    pep_start <- unname(pep_idx[, "start"])
    site_relidx <- str_locate_all(pep_seq, mod_res_symbol)[[1]]
    site_relstart <- unname(site_relidx[, "start"])
    mod_idx <- site_relstart - seq_along(site_relstart) + pep_start
    
    return(mod_idx)
}


# Data processing & summarization -----------------------------------------

# Normalization with unpaired unmodified peptides
normalize_ptm <- function(data) {
    cols <- names(data)
    # Based on unmodified peptides
    if ("batch" %in% cols) {
        medians <- data %>% 
            filter(!is_mod) %>% 
            group_by(batch, run) %>% 
            summarise(log2inty_med = median(log2inty, na.rm = TRUE)) %>% 
            mutate(log2inty_bch = median(log2inty_med)) %>% 
            ungroup()
    } else {
        medians <- data %>% 
            filter(!is_mod) %>% 
            group_by(run) %>% 
            summarise(log2inty_med = median(log2inty, na.rm = TRUE)) %>% 
            mutate(log2inty_bch = median(log2inty_med)) %>% 
            ungroup()
    }
    normdata <- data %>%
        left_join(medians) %>%
        mutate(log2inty = ifelse(is.na(log2inty_med), log2inty, log2inty - log2inty_med + log2inty_bch)) %>%
        select(-log2inty_med, -log2inty_bch)
    
    return(normdata)
}


# Annotate censored values and fill with AFT
fill_censored_aft <- function(data) {
    if (n_distinct(data$feature) == 1) return(data)  # only one feature
    if (all(!is.na(data$log2inty))) return(data)  # no missing value
    # Annotate observation status and add censored values
    aftdata <- data %>% 
        mutate(ind_obs = ifelse(is.na(log2inty), 0L, 1L)) %>% 
        group_by(feature) %>% 
        mutate(log2inty_aft = ifelse(ind_obs == 0, 0.99 * min(log2inty[ind_obs == 1]), log2inty)) %>% 
        ungroup()
    # AFT model with effects of run and feature
    fit <- tryCatch(
        survreg(Surv(log2inty_aft, ind_obs, type = "left") ~ run + feature, data = aftdata, dist = "gaussian"), 
        error = function(e) e, warning = function(w) w
    )
    if (is(fit, "warning")) return(data)  # not converged
    aftdata <- aftdata %>% 
        mutate(
            log2inty_pred = predict(fit), 
            log2inty = ifelse(ind_obs == 0, pmin(log2inty_pred, log2inty_aft), log2inty_aft)
        ) %>% 
        select(-log2inty_aft, -log2inty_pred, -ind_obs)
    
    return(aftdata)
}


# Summarization of feature intensities at run level using Tukey's median polish
summarize_feature <- function(df_prot, method = "tmp") {
    if (method == "tmp") {
        inty_wide <- df_prot %>% 
            select(feature, run, log2inty) %>% 
            spread(feature, log2inty)
        inty_mat <- data.matrix(inty_wide[, -1])
        mp_out <- medpolish(inty_mat, na.rm = TRUE, trace.iter = FALSE)
        df_sum <- data_frame(run = inty_wide$run, log2inty = mp_out$overall + mp_out$row)
    } else {
        df_sum <- df_prot %>% 
            group_by(run) %>% 
            summarise(log2inty = log2(sum(2 ^ log2inty, na.rm = TRUE)))
    }
    
    return(df_sum)
}


# Whole-plot modeling -----------------------------------------------------

# Linear model with group effect (and potentially batch effect)
lm_group <- function(df_onesite, w_batch = FALSE) {
    if (w_batch) {
        # Model with batch effect
        if (n_distinct(df_onesite$batch) == 1) 
            stop("Cannot estimate batch effect with a single batch!")
        
        if (n_distinct(df_onesite$group) == 1) {
            fit <- lm(log2inty ~ batch, data = df_onesite)
        } else {
            fit <- lm(log2inty ~ 0 + group + batch, data = df_onesite)
        }
    } else {
        if (n_distinct(df_onesite$group) == 1) {
            fit <- lm(log2inty ~ 1, data = df_onesite)
        } else {
            fit <- lm(log2inty ~ 0 + group, data = df_onesite)
        }
    }
    
    return(fit)
}

nest_site <- function(df_sum, w_batch = FALSE) {
    if (w_batch) {
        if (!("batch" %in% names(df_sum)))
            stop("There is no information about batch!")
        
        # One model for all batches (need data with >1 batches)
        nested_sum <- df_sum %>% 
            group_by(uniprot_iso, site_str) %>% 
            nest() %>% 
            mutate(nb_bch = map_int(data, ~ n_distinct(.$batch))) %>% 
            filter(nb_bch > 1) %>% 
            select(-nb_bch)
        # Linear model and associated parameter estimates
        nested_sum <- nested_sum %>% 
            mutate(lm_fit = map(data, lm_group, w_batch)) %>% 
            mutate(
                param = map2(lm_fit, data, tidy_bch), 
                df_res = map_dbl(lm_fit, df.residual)
            )
        # Remove cases not eligible for hypothesis testing (SE is NA)
        # [TODO]: consider another option to report FC alone 
        nested_sum <- nested_sum %>% 
            filter(df_res > 0, !map_lgl(param, ~any(is.na(.$std.error))))
    } else {
        # One model per site (and potentially batch)
        if ("batch" %in% names(df_sum)) {
            nested_sum <- df_sum %>%
                group_by(uniprot_iso, site_str, batch) %>%
                nest()
        } else {
            nested_sum <- df_sum %>%
                group_by(uniprot_iso, site_str) %>%
                nest()
        }
        # Linear model and associated parameter estimates
        nested_sum <- nested_sum %>% 
            mutate(lm_fit = map(data, lm_group, w_batch)) %>%
            mutate(
                param = map2(lm_fit, data, tidy_bch),
                df_res = map_dbl(lm_fit, df.residual)
            )
        # Remove cases not eligible for hypothesis testing (SE is NA)
        # [TODO]: consider another option to report FC alone 
        nested_sum <- nested_sum %>% 
            filter(df_res > 0, !map_lgl(param, ~any(is.na(.$std.error))))
    }
    
    return(nested_sum)
}



# # Linear model with group effect 
# lm_perbch <- function(df_perbch) {
#     if (n_distinct(df_perbch$group) == 1) {
#         fit <- lm(log2inty ~ 1, data = df_perbch)
#     } else {
#         fit <- lm(log2inty ~ 0 + group, data = df_perbch)
#     }
#     
#     return(fit)
# }
# 
# 
# # Linear model with group effect and batch effect
# lm_allbch <- function(df_allbch) {
#     if (n_distinct(df_allbch$batch) == 1) 
#         stop("Cannot infer batch effect with a single batch!")
#     
#     if (n_distinct(df_allbch$group) == 1) {
#         fit <- lm(log2inty ~ batch, data = df_allbch)
#     } else {
#         fit <- lm(log2inty ~ 0 + group + batch, data = df_allbch)
#     }
#     
#     return(fit)
# }


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


# Extract estimated parameters
extract_param <- function(nested_data) {
    if ("batch" %in% names(nested_data)) {
        # per-batch model
        nested_param <- nested_data %>% 
            select(uniprot_iso, site_str, batch, param, df_res)
    } else {
        # all-batch model
        nested_param <- nested_data %>% 
            select(uniprot_iso, site_str, param, df_res)
    }
    param_mod <- nested_param %>% 
        filter(site_str != "None") %>% 
        unnest(param)
    param_unmod <- nested_param %>% 
        filter(site_str == "None") %>% 
        select(-site_str) %>% 
        unnest(param) %>% 
        rename(df_unmod = df_res, est_unmod = estimate, se_unmod = std.error)
    
    return(left_join(param_mod, param_unmod) %>% unite(protsite, uniprot_iso, site_str, sep = "--"))
}


# Differential analysis ---------------------------------------------------

compare_mod <- function(df_mod, grp_ctrl, grp_case, protadj = TRUE) {
    if ("batch" %in% names(df_mod)) {
        # Batch-aggregated testing (with multiple difference estimates)
        if (protadj) {
            # With protein-level adjustment
            full_mod <- df_mod %>% 
                filter(group %in% c(grp_ctrl, grp_case)) %>% 
                filter(!is.na(df_res) | !is.na(df_unmod)) %>% 
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
            res <- diff_mod %>% 
                group_by(protsite) %>% 
                summarise(
                    log2FC = mean(log2fc) - mean(log2fc_unmod), 
                    std_error = sqrt(sum(se2) + sum(se2_unmod)) / n(), 
                    DF = (sum(se2) + sum(se2_unmod)) ^ 2 / sum(se2 ^ 2 / df_res + se2_unmod ^ 2 / df_unmod), 
                    statistic = log2FC / std_error, p_value = 2 * pt(abs(statistic), df = DF, lower.tail = FALSE)
                ) %>% 
                mutate(contrast = paste(grp_case, grp_ctrl, sep = " vs "))
            # Missing in one group
            part_mod <- full_mod %>% 
                anti_join(diff_mod %>% select(protsite)) %>% 
                group_by(protsite) %>% 
                filter(all(is.na(df_res) == is.na(df_unmod))) %>% 
                filter(
                    n_distinct(is.na(df_res[key_grp == "G0"])) == 1, 
                    n_distinct(is.na(df_res[key_grp == "G1"])) == 1
                ) %>% 
                ungroup() %>% 
                filter(!is.na(df_res)) %>% 
                distinct(protsite, key_grp)
            if (nrow(part_mod) > 0) {
                untest_mod <- part_mod %>% 
                    mutate(
                        log2FC = ifelse(key_grp == "G1", Inf, -Inf), 
                        std_error = NA, DF = NA, statistic = NA, p_value = NA, 
                        contrast = paste(grp_case, grp_ctrl, sep = " vs ")
                    ) %>% select(-key_grp)
                res <- res %>% bind_rows(untest_mod)
            }
        } else {
            # Without protein-level adjustment
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
            res <- diff_mod %>% 
                group_by(protsite) %>% 
                summarise(log2FC = mean(log2fc), std_error = sqrt(sum(se2)) / n(), 
                          DF = sum(se2) ^ 2 / sum(se2 ^ 2 / df_res), 
                          statistic = log2FC / std_error, p_value = 2 * pt(abs(statistic), df = DF, lower.tail = FALSE)) %>% 
                mutate(contrast = paste(grp_case, grp_ctrl, sep = " vs "))
            # Missing in one group
            part_mod <- full_mod %>% 
                anti_join(diff_mod %>% select(protsite)) %>% 
                group_by(protsite) %>% 
                filter(
                    n_distinct(is.na(df_res[key_grp == "G0"])) == 1, 
                    n_distinct(is.na(df_res[key_grp == "G1"])) == 1
                ) %>% 
                ungroup() %>% 
                filter(!is.na(df_res)) %>% 
                distinct(protsite, key_grp)
            if (nrow(part_mod) > 0) {
                untest_mod <- part_mod %>% 
                    mutate(
                        log2FC = ifelse(key_grp == "G1", Inf, -Inf), 
                        std_error = NA, DF = NA, statistic = NA, p_value = NA, 
                        contrast = paste(grp_case, grp_ctrl, sep = " vs ")
                    ) %>% select(-key_grp)
                res <- res %>% bind_rows(untest_mod)
            }
        }
    } else {
        # Single difference estimate between groups
        if (protadj) {
            # With protein-level adjustment
            full_mod <- df_mod %>% 
                filter(group %in% c(grp_ctrl, grp_case)) %>% 
                filter(!is.na(df_res) | !is.na(df_unmod)) %>% 
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
            res <- diff_mod %>% 
                mutate(log2FC = log2fc - log2fc_unmod, std_error = sqrt(se2 + se2_unmod), 
                       DF = (se2 + se2_unmod) ^ 2 / (se2 ^ 2 / df_res + se2_unmod ^ 2 / df_unmod), 
                       statistic = log2FC / std_error, p_value = 2 * pt(abs(statistic), df = DF, lower.tail = FALSE)) %>% 
                select(protsite, log2FC, std_error, DF, statistic, p_value) %>% 
                mutate(contrast = paste(grp_case, grp_ctrl, sep = " vs "))
            # Missing in one group
            part_mod <- full_mod %>% 
                anti_join(diff_mod %>% select(protsite)) %>% 
                group_by(protsite) %>% 
                filter(all(is.na(df_res) == is.na(df_unmod))) %>% 
                ungroup() %>% 
                filter(!is.na(df_res)) %>% 
                distinct(protsite, key_grp)
            if (nrow(part_mod) > 0) {
                untest_mod <- part_mod %>% 
                    mutate(
                        log2FC = ifelse(key_grp == "G1", Inf, -Inf), 
                        std_error = NA, DF = NA, statistic = NA, p_value = NA, 
                        contrast = paste(grp_case, grp_ctrl, sep = " vs ")
                    ) %>% select(-key_grp)
                res <- res %>% bind_rows(untest_mod)
            }
        } else {
            # Without protein-level adjustment
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
            res <- diff_mod %>% 
                mutate(log2FC = log2fc, std_error = sqrt(se2), DF = df_res, 
                       statistic = log2FC / std_error, p_value = 2 * pt(abs(statistic), df = DF, lower.tail = FALSE)) %>% 
                select(protsite, log2FC, std_error, DF, statistic, p_value) %>% 
                mutate(contrast = paste(grp_case, grp_ctrl, sep = " vs "))
            # Missing in one group
            part_mod <- full_mod %>% 
                anti_join(diff_mod %>% select(protsite)) %>% 
                filter(!is.na(df_res)) %>% 
                distinct(protsite, key_grp)
            if (nrow(part_mod) > 0) {
                untest_mod <- part_mod %>% 
                    mutate(
                        log2FC = ifelse(key_grp == "G1", Inf, -Inf), 
                        std_error = NA, DF = NA, statistic = NA, p_value = NA, 
                        contrast = paste(grp_case, grp_ctrl, sep = " vs ")
                    ) %>% select(-key_grp)
                res <- res %>% bind_rows(untest_mod)
            }
        }
    }
    
    return(res)
}


# # Testing for differential modification
# test_mod <- function(df_mod, grp_ctrl, grp_case) {
#     full_mod <- df_mod %>% 
#         filter(group %in% c(grp_ctrl, grp_case)) %>% 
#         filter(!is.na(df_res)) %>% 
#         mutate(key_grp = ifelse(group == grp_ctrl, "G0", "G1")) %>% 
#         complete(protsite, key_grp)
#     # Model-based inference of log-difference
#     diff_mod <- full_mod %>% 
#         group_by(protsite) %>% 
#         filter(!any(is.na(df_res))) %>% 
#         group_by(protsite, df_res) %>% 
#         summarise(log2fc = diff(estimate), se2 = sum(std.error ^ 2)) %>% 
#         ungroup()
#     res <- diff_mod %>% 
#         mutate(logFC = log2fc, SE = sqrt(se2), DF = df_res, 
#                t_stat = logFC / SE, p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)) %>% 
#         select(protsite, logFC, SE, DF, t_stat, p_val) %>% 
#         mutate(ctrx = paste(grp_case, grp_ctrl, sep = " vs "))
#     # Missing in one group
#     part_mod <- full_mod %>% 
#         anti_join(diff_mod %>% select(protsite)) %>% 
#         filter(!is.na(df_res)) %>% 
#         distinct(protsite, key_grp)
#     if (nrow(part_mod) > 0) {
#         untest_mod <- part_mod %>% 
#             mutate(
#                 logFC = ifelse(key_grp == "G1", Inf, -Inf), 
#                 SE = NA, DF = NA, t_stat = NA, p_val = NA, 
#                 ctrx = paste(grp_case, grp_ctrl, sep = " vs ")
#             ) %>% select(-key_grp)
#         res <- res %>% bind_rows(untest_mod)
#     }
#     
#     return(res)
# }
# 
# 
# # Testing for differential modification with adjustment
# test_adjmod <- function(df_mod, grp_ctrl, grp_case) {
#     full_mod <- df_mod %>% 
#         filter(group %in% c(grp_ctrl, grp_case)) %>% 
#         filter(!is.na(df_res) | !is.na(df_unmod)) %>% 
#         mutate(key_grp = ifelse(group == grp_ctrl, "G0", "G1")) %>% 
#         complete(protsite, key_grp)
#     # Model-based inference of log-difference
#     diff_mod <- full_mod %>% 
#         group_by(protsite) %>% 
#         filter(all(!is.na(df_res)), all(!is.na(df_unmod))) %>% 
#         group_by(protsite, df_res, df_unmod) %>% 
#         summarise(
#             log2fc = diff(estimate), se2 = sum(std.error ^ 2), 
#             log2fc_unmod = diff(est_unmod), se2_unmod = sum(se_unmod ^ 2)
#         ) %>% ungroup()
#     res <- diff_mod %>% 
#         mutate(logFC = log2fc - log2fc_unmod, SE = sqrt(se2 + se2_unmod), 
#                DF = (se2 + se2_unmod) ^ 2 / (se2 ^ 2 / df_res + se2_unmod ^ 2 / df_unmod), 
#                t_stat = logFC / SE, p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)) %>% 
#         select(protsite, logFC, SE, DF, t_stat, p_val) %>% 
#         mutate(ctrx = paste(grp_case, grp_ctrl, sep = " vs "))
#     # Missing in one group
#     part_mod <- full_mod %>% 
#         anti_join(diff_mod %>% select(protsite)) %>% 
#         group_by(protsite) %>% 
#         filter(all(is.na(df_res) == is.na(df_unmod))) %>% 
#         ungroup() %>% 
#         filter(!is.na(df_res)) %>% 
#         distinct(protsite, key_grp)
#     if (nrow(part_mod) > 0) {
#         untest_mod <- part_mod %>% 
#             mutate(
#                 logFC = ifelse(key_grp == "G1", Inf, -Inf), 
#                 SE = NA, DF = NA, t_stat = NA, p_val = NA, 
#                 ctrx = paste(grp_case, grp_ctrl, sep = " vs ")
#             ) %>% select(-key_grp)
#         res <- res %>% bind_rows(untest_mod)
#     }
#     
#     return(res)
# }
# 
# 
# # Batch-aggregated testing for differential modification
# test_mod_bch <- function(df_mod, grp_ctrl, grp_case) {
#     full_mod <- df_mod %>% 
#         filter(group %in% c(grp_ctrl, grp_case)) %>% 
#         filter(!is.na(df_res)) %>% 
#         mutate(key_grp = ifelse(group == grp_ctrl, "G0", "G1")) %>% 
#         complete(protsite, batch, key_grp)
#     # Model-based inference of log-difference
#     diff_mod <- full_mod %>% 
#         group_by(protsite) %>% 
#         filter(all(!is.na(df_res))) %>% 
#         group_by(protsite, batch, df_res) %>% 
#         summarise(log2fc = diff(estimate), se2 = sum(std.error ^ 2)) %>% 
#         ungroup()
#     res <- diff_mod %>% 
#         group_by(protsite) %>% 
#         summarise(logFC = mean(log2fc), SE = sqrt(sum(se2)) / n(), 
#                   DF = sum(se2) ^ 2 / sum(se2 ^ 2 / df_res), 
#                   t_stat = logFC / SE, p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)) %>% 
#         mutate(ctrx = paste(grp_case, grp_ctrl, sep = " vs "))
#     # Missing in one group
#     part_mod <- full_mod %>% 
#         anti_join(diff_mod %>% select(protsite)) %>% 
#         group_by(protsite) %>% 
#         filter(
#             n_distinct(is.na(df_res[key_grp == "G0"])) == 1, 
#             n_distinct(is.na(df_res[key_grp == "G1"])) == 1
#         ) %>% 
#         ungroup() %>% 
#         filter(!is.na(df_res)) %>% 
#         distinct(protsite, key_grp)
#     if (nrow(part_mod) > 0) {
#         untest_mod <- part_mod %>% 
#             mutate(
#                 logFC = ifelse(key_grp == "G1", Inf, -Inf), 
#                 SE = NA, DF = NA, t_stat = NA, p_val = NA, 
#                 ctrx = paste(grp_case, grp_ctrl, sep = " vs ")
#             ) %>% select(-key_grp)
#         res <- res %>% bind_rows(untest_mod)
#     }
#     
#     return(res)
# }
# 
# 
# # Batch-aggregated testing for differential modification with adjustment
# test_adjmod_bch <- function(df_mod, grp_ctrl, grp_case) {
#     full_mod <- df_mod %>% 
#         filter(group %in% c(grp_ctrl, grp_case)) %>% 
#         filter(!is.na(df_res) | !is.na(df_unmod)) %>% 
#         mutate(key_grp = ifelse(group == grp_ctrl, "G0", "G1")) %>% 
#         complete(protsite, batch, key_grp)
#     # Model-based inference of log-difference
#     diff_mod <- full_mod %>% 
#         group_by(protsite) %>% 
#         filter(all(!is.na(df_res)), all(!is.na(df_unmod))) %>% 
#         group_by(protsite, batch, df_res, df_unmod) %>% 
#         summarise(
#             log2fc = diff(estimate), se2 = sum(std.error ^ 2),
#             log2fc_unmod = diff(est_unmod), se2_unmod = sum(se_unmod ^ 2)
#         ) %>% ungroup()
#     res <- diff_mod %>% 
#         group_by(protsite) %>% 
#         summarise(
#             logFC = mean(log2fc) - mean(log2fc_unmod), 
#             SE = sqrt(sum(se2) + sum(se2_unmod)) / n(), 
#             DF = (sum(se2) + sum(se2_unmod)) ^ 2 / sum(se2 ^ 2 / df_res + se2_unmod ^ 2 / df_unmod), 
#             t_stat = logFC / SE, p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)
#         ) %>% 
#         mutate(ctrx = paste(grp_case, grp_ctrl, sep = " vs "))
#     # Missing in one group
#     part_mod <- full_mod %>% 
#         anti_join(diff_mod %>% select(protsite)) %>% 
#         group_by(protsite) %>% 
#         filter(all(is.na(df_res) == is.na(df_unmod))) %>% 
#         filter(
#             n_distinct(is.na(df_res[key_grp == "G0"])) == 1, 
#             n_distinct(is.na(df_res[key_grp == "G1"])) == 1
#         ) %>% 
#         ungroup() %>% 
#         filter(!is.na(df_res)) %>% 
#         distinct(protsite, key_grp)
#     
#     if (nrow(part_mod) > 0) {
#         untest_mod <- part_mod %>% 
#             mutate(
#                 logFC = ifelse(key_grp == "G1", Inf, -Inf), 
#                 SE = NA, DF = NA, t_stat = NA, p_val = NA, 
#                 ctrx = paste(grp_case, grp_ctrl, sep = " vs ")
#             ) %>% select(-key_grp)
#         res <- res %>% bind_rows(untest_mod)
#     }
#     
#     return(res)
# }


