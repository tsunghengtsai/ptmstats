# Load libraries and sources ----------------------------------------------

library(tidyverse)
library(stringr)
library(broom)
library(survival)

library(devtools)
load_all(pkg = "~/Projects/MSstatsPTM/")


# Prepare data set --------------------------------------------------------

# Load experimental design with data filenames
df_design <- read.table(
    "data/RIP2/OBJ002798_tft_20180307113542.tsv",
    header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

file_name <- paste0("vistaga_", df_design$run_id, ".tsv")

sub_cols <- c("Reference", "run_id", "peptide_id", "vista_peak_area_light",
              "vista_confidence_score", "peptide_m_z", "peptide_sequence",
              "peptide_trypticity", "peptide_miscleavages", "peptide_validity",
              "label1", "ms2_rt", "ms2_charge", "is_xq")

# Prepare data set with columns of interest for later use
dta_rip2 <- vector("list", length(file_name))
for (i in seq_along(file_name)) {
    file_path <- paste0("data/RIP2/", file_name[i])
    dta_rip2[[i]] <- read.table(
        file_path, header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE
    ) %>% 
        tbl_df() %>% 
        select_(.dots = sub_cols)
}
dta_rip2 <- bind_rows(dta_rip2)
rm(list = ls()[!(ls() %in% c("df_design", "dta_rip2"))])

# save.image("output/dta_rip2.RData")


# Load data set & housekeeping --------------------------------------------

# Parameters
site_spec <- TRUE  # Site-level analysis; default TRUE
only_protadj <- FALSE  # Ignore unadjustable PTM; default FALSE
min_len_peptide <- 6  # Minimum acceptable length of peptide

# Modifications of interest
mod_residue <- "K"
mod_symbol <- "\\*"

# Load prepared dataset
load("output/dta_rip2.RData")

df_work <- dta_rip2
# rm(dta_rip2)

# Initial filtering
df_work <- df_work %>% 
    filter(vista_confidence_score >= 83, !is.na(ms2_charge), !is.na(ms2_rt)) %>%
    select(-vista_confidence_score, -peptide_trypticity, 
           -peptide_miscleavages, -peptide_validity) %>% 
    rename(peptide = peptide_sequence)

# Group information is organized in slightly different ways across GNE datasets
df_design <- df_design %>% 
    mutate(group = str_replace_all(
        cond_treatment, 
        c("min|_1hr|_only" = "", "MDP" = "mdp", "none" = "ctrl")
    ))

# OBJ# can be viewed as a batch identity
subj_uniq <- unique(df_design$id_subject)
subj_key <- str_c("S", str_pad(subj_uniq, width = 2, pad = "0"))
bio_uniq <- unique(df_design$ann_bioreplicate)
bio_key <- str_c("B", bio_uniq)
df_design <- tbl_df(df_design) %>% 
    mutate(
        subj = plyr::mapvalues(id_subject, from = subj_uniq, to = subj_key), 
        biorep = plyr::mapvalues(ann_bioreplicate, from = bio_uniq, to = bio_key), 
        batch = ifelse(biorep == "B1", "BCH1", "BCH2"), 
        techrep = str_c("T", id_injectionset), 
        run_bt = str_c(biorep, techrep), 
        run_cbt = str_c(group, run_bt, sep = "-")
    )


# Annotation & data filtering ---------------------------------------------

# Read fasta annotation
hs_fasta <- tidy_fasta("data/Sequence/homo_sapiens_all_20160725.fasta")

df_work <- df_work %>% 
    rename(entry_name = Reference) %>% 
    left_join(hs_fasta %>% distinct(uniprot_ac, entry_name)) %>% 
    rename(uniprot_iso = uniprot_ac)

df_work <- df_work %>% 
    mutate(
        is_mod = str_detect(peptide, mod_symbol),
        feature = str_c(peptide, ms2_charge, sep = "_"), 
        log2inty = ifelse(vista_peak_area_light <= 1, 0, log2(vista_peak_area_light)), 
        group = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$group), 
        batch = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$batch), 
        run = plyr::mapvalues(run_id, from = df_design$run_id, to = df_design$run_cbt)
    )

# Use the max for duplicate features (per run)
df_work <- df_work %>% 
    group_by(uniprot_iso, peptide, is_mod, feature, group, batch, run) %>% 
    summarise(log2inty = max(log2inty)) %>% 
    ungroup()

# Remove peptides with ambiguous matches (to >1 uniprot isoforms)
# uniq_pep <- df_work %>% 
#     distinct(peptide, uniprot_iso) %>% 
#     count(peptide) %>% 
#     filter(n == 1) %>% 
#     .$peptide
# df_work <- df_work %>% filter(peptide %in% uniq_pep)

# Remove features exclusively observed in one batch
# df_work <- df_work %>%
#     group_by(feature, batch) %>%
#     filter(n() != 1) %>% 
#     ungroup()


# Site representation for modified peptides -------------------------------

# Annotate potential modification sites on protein sequences
hs_fasta <- hs_fasta %>% 
    mutate(
        rng_site_all = str_locate_all(sequence, mod_residue), 
        idx_site_all = map(rng_site_all, ~.[, "start"]), 
        aa_site_all = str_match_all(sequence, mod_residue), 
        aa_site_all = map(aa_site_all, ~.[, 1])
    ) %>% 
    select(uniprot_ac, uniprot_iso, entry_name, idx_site_all, aa_site_all, header, sequence) %>% 
    arrange(uniprot_iso)

# hs_fasta <- hs_fasta %>% 
#     mutate(
#         rng_site_all = str_locate_all(sequence, mod_residue), 
#         aa_site_all = str_match_all(sequence, mod_residue)
#     ) %>% 
#     mutate(
#         idx_site_all = map(rng_site_all, ~.[, "start"]), 
#         aa_site_all = map(aa_site_all, ~.[, 1])
#     )
# 
# hs_fasta <- hs_fasta %>% 
#     select(uniprot_ac, uniprot_iso, entry_name, idx_site_all, aa_site_all, header, sequence) %>% 
#     arrange(uniprot_iso)

# Observed peptide sequence & other information
df_mod <- df_work %>% 
    distinct(uniprot_iso, peptide, is_mod) %>%
    mutate(
        peptide_trimmed = str_extract(peptide, "(?<=\\.)([ACDEFGHIKLMNPQRSTVWY\\*]+)"), 
        peptide_unmod = str_replace_all(peptide, mod_symbol, ""), 
        peptide_unmod_trimmed = str_extract(peptide_unmod, "(?<=\\.)([ACDEFGHIKLMNPQRSTVWY\\*]+)"), 
        len_peptide = str_length(peptide_unmod_trimmed)
    ) %>% 
    arrange(uniprot_iso)

# Retain peptides with length >= 6
regex_uniprot_iso <- regex("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})([-]\\d{1,}){0,1}")
# regex_uniprot <- ".*([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}).*"
df_mod <- df_mod %>% 
    filter(len_peptide >= min_len_peptide, str_detect(uniprot_iso, regex_uniprot_iso))

# Integrate fasta information with observed peptides
# [TODO]: use extended AAs for more specific matching
df_fasmod <- df_mod %>% 
    inner_join(hs_fasta %>% select(uniprot_iso, sequence, idx_site_all)) %>% 
    mutate(rng_peptide = str_locate_all(sequence, peptide_unmod_trimmed)) %>% 
    mutate(nb_mch = map_int(rng_peptide, ~ nrow(.)))

# Ignoring non-specific matching (peptide mapped to 0 or >1 locations of protein)
df_fasmod <- df_fasmod %>% 
    filter(nb_mch == 1) %>% 
    mutate(aa_start = map_int(rng_peptide, ~.[, "start"]))

# Observed modification sites
mod_resymb <- str_c(mod_residue, mod_symbol)
df_fasmod <- df_fasmod %>% 
    mutate(
        idx_site = pmap(list(peptide_unmod_trimmed, aa_start, mod_residue), locate_site), 
        idx_mod = pmap(list(peptide_trimmed, aa_start, mod_resymb), locate_mod)
    ) %>% 
    mutate(
        nb_site = map_int(idx_site, length), 
        nb_mod = map_int(idx_mod, length), 
        len_site = map(idx_site_all, ~str_length(.[length(.)])), 
        site = pmap_chr(list(idx_mod, mod_residue, len_site), annot_site), 
        # site = map2_chr(idx_mod, idx_site_all, annotate_site, mod_residue), 
        peptide_str = map_chr(rng_peptide, ~ str_c(., collapse = "-")), 
        full_site_str = str_c(uniprot_iso, peptide_str, site, sep = "_")
    ) %>% 
    select(uniprot_iso, peptide, peptide_unmod, nb_site, nb_mod, is_mod, 
           idx_site, idx_mod, site, peptide_str, full_site_str)

# Check if unmodified peptides have a site modified elsewhere
# Firstly, find protein with at least one unmodified peptide/site (to be checked)
nested_w_unmod <- df_fasmod %>% 
    group_by(uniprot_iso) %>% 
    filter(any(!is_mod)) %>% 
    nest()
# Secondly, compare the unmodified peptides with modified sites
nested_w_unmod <- nested_w_unmod %>% 
    mutate(unmod_data = map(data, ~filter(., !is_mod))) %>% 
    mutate(site_mod_prot = map(data, ~unlist(.$idx_mod))) %>% 
    mutate(unmod_data = map2(unmod_data, site_mod_prot, ~mutate(.x, w_mod = idx_site %in% .y)))

# Restrict on peptides modified on 1 site and unmodified peptides with no site modified elsewhere
# (Another option) Restrict on peptides modified on 1 site and unmodified peptides with no site
df_unconfound <- bind_rows(
    nested_w_unmod %>% unnest(unmod_data) %>% filter(!w_mod) %>% 
        select(uniprot_iso, peptide, is_mod, site), 
    df_fasmod %>% filter(nb_mod == 1) %>% 
        select(uniprot_iso, peptide, is_mod, site)
)

# Proteins eligible for protein-level adjustment
# if (only_protadj) {
#     df_unconfound <- df_unconfound %>%
#         group_by(uniprot_iso) %>%
#         filter(n_distinct(is_mod) == 2) %>%
#         ungroup()
# }

# df_work <- df_work %>% inner_join(df_unconfound)
df_work <- df_work %>% 
    inner_join(df_unconfound) %>% 
    rename(protein = uniprot_iso)

# Visualization -----------------------------------------------------------

# Key for ordering
run_fac <- str_c(
    rep(c("ctrl", "mdp_30", "mdp_60", "cmpd89", "cmpd89_mdp_30", "cmpd89_mdp_60"), each = 4), 
    rep(c("B1T1", "B1T2", "B2T1", "B2T2"), 6), 
    sep = "-"
)
grp_fac <- c("ctrl", "mdp_30", "mdp_60", "cmpd89", "cmpd89_mdp_30", "cmpd89_mdp_60")

# QC plot (boxplot)
df_work %>% 
    mutate(is_mod_fac = factor(ifelse(is_mod, "Modified", "Unmodified"))) %>% 
    mutate(run = factor(run, levels = run_fac, labels = run_fac)) %>% 
    ggplot(aes(run, log2inty)) + 
    geom_boxplot(aes(fill = batch)) + 
    geom_vline(xintercept = c(4, 8, 12, 16, 20) + 0.5, linetype = "dashed") + 
    facet_wrap(~ is_mod_fac) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(x = "Run", y = "Log2-intensity", 
         title = "Batch effect in terms of intensity level")
# ggsave("rip2_box.png", width = 9, height = 6)
# ggsave("rip2_box_norm.png", width = 9, height = 6)

# Number of modified/unmodified peptides
df_work %>% 
    distinct(group, batch, run, peptide, is_mod) %>% 
    group_by(group, batch, run) %>% 
    summarise(Modified = sum(is_mod), Unmodified = sum(!is_mod)) %>% 
    ungroup() %>% 
    mutate(group = factor(group, levels = grp_fac, labels = grp_fac)) %>% 
    gather(mod, nb_peptide, Modified:Unmodified) %>% 
    ggplot(aes(group, nb_peptide)) + 
    geom_jitter(aes(colour = batch, shape = mod), width = 0.1, size = 4) + 
    labs(x = "Group", y = "Number of peptides", 
         title = "Batch effect in terms of # identified peptides")
# ggsave("rip2_nbpep.png", width = 6, height = 4)


# Normalization & filtering -----------------------------------------------

# Based on unpaired unmodified peptides
df_work <- normalize_ptm(df_work)

# Annotation for either site-level analysis or protein-level analysis
if (!site_spec) {
    df_work <- df_work %>% mutate(site = ifelse(site == "None", "None", "MOD"))
}

# Fully observed sites (in at least one conditions)
site_full <- df_work %>%
    group_by(protein, site, group) %>%
    summarise(nb_run = n_distinct(run)) %>%
    ungroup() %>%
    left_join(df_design %>% group_by(group) %>% summarise(nb_run_design = n())) %>% 
    filter(nb_run == nb_run_design) %>% 
    distinct(protein, site)

# if (only_protadj) {
#     site_full <- site_full %>% 
#         group_by(protein) %>% 
#         filter(n_distinct(site == "None") == 2) %>% 
#         ungroup()
# }

df_site <- df_work %>% 
    semi_join(site_full)


# Imputation and summarization --------------------------------------------

if ("batch" %in% names(df_site)) {
    nested_site <- df_site %>%
        group_by(protein, site, batch) %>%
        nest()
} else {
    nested_site <- df_site %>%
        group_by(protein, site) %>%
        nest()
}

nested_site <- nested_site %>% 
    mutate(data = map(data, ~ complete(., feature, nesting(run, group))))

nested_site <- nested_site %>% 
    mutate(aftdata = map(data, fill_censored_aft))

nested_site <- nested_site %>% 
    mutate(sumdata = map(aftdata, summarize_feature, "tmp"))


# Fit per-batch and all-batch models --------------------------------------
# Add group information
run2group <- df_site %>% 
    distinct(run, group)
df_sum <- nested_site %>% 
    unnest(sumdata) %>% 
    left_join(run2group)

# # Whole-plot modeling
# nested_perbch <- nest_site(df_sum, w_batch = FALSE)
# nested_allbch <- nest_site(df_sum, w_batch = TRUE)
# 
# # Extract estimated model parameters
# param_perbch <- extract_param(nested_perbch)
# param_allbch <- extract_param(nested_allbch)

# Whole-plot modeling
param_perbch <- model_ptm(df_sum, w_batch = FALSE)
param_allbch <- model_ptm(df_sum, w_batch = TRUE)


# Site-level differential analysis ----------------------------------------

cases <- c("cmpd89_mdp_30")
controls <- c("mdp_30")

# cases <- c("mdp_30")
# controls <- c("ctrl")
# cases <- c("cmpd89_mdp_30", "mdp_30")
# controls <- c("mdp_30", "ctrl")

test_1 <- test_2 <- test_3 <- test_4 <- vector("list", length = length(cases))
for (i in seq_along(cases)) {
    grp_ctrl <- controls[i]
    grp_case <- cases[i]
    # Per-batch model
    test_1[[i]] <- compare_mod(param_perbch, grp_ctrl, grp_case, protadj = FALSE) %>% 
        mutate(model = "per-batch", protadj = "no adjustment")
    test_2[[i]] <- compare_mod(param_perbch, grp_ctrl, grp_case, protadj = TRUE) %>% 
        mutate(model = "per-batch", protadj = "protein adjustment")
    # All-batch model
    test_3[[i]] <- compare_mod(param_allbch, grp_ctrl, grp_case, protadj = FALSE) %>% 
        mutate(model = "all-batch", protadj = "no adjustment")
    test_4[[i]] <- compare_mod(param_allbch, grp_ctrl, grp_case, protadj = TRUE) %>% 
        mutate(model = "all-batch", protadj = "protein adjustment")
    
    # Append unadjustable results (EXPERIMENTAL)
    test_2[[i]] <- bind_rows(
        test_2[[i]], 
        test_1[[i]] %>% 
            anti_join(test_2[[i]] %>% select(protein, site)) %>% 
            mutate(model = "per-batch", protadj = "protein adjustment")
    )
    test_4[[i]] <- bind_rows(
        test_4[[i]], 
        test_3[[i]] %>% 
            anti_join(test_4[[i]] %>% select(protein, site)) %>% 
            mutate(model = "all-batch", protadj = "protein adjustment")
    )

} 

test_res <- bind_rows(test_1, test_2, test_3, test_4) %>% 
    group_by(model, protadj) %>%
    mutate(p_adjusted = p.adjust(p_value, method = "BH")) %>% 
    ungroup()
test_res$p_adjusted[test_res$log2FC %in% c(Inf, -Inf)] <- 0  # Missing in one group


# Site-level results ------------------------------------------------------

library(ggrepel)

# Volcano plot
test_res %>% 
    group_by(protein, site, contrast) %>% filter(n() == 4) %>% ungroup() %>% 
    filter(abs(log2FC) < Inf) %>%
    ggplot(aes(log2FC, -log10(p_adjusted))) + 
    geom_point(color = "darkgray") + 
    geom_hline(yintercept = 0) + 
    geom_vline(xintercept = 0) + 
    geom_hline(yintercept = -log10(0.05), color = "darkred") + 
    geom_vline(xintercept = c(2, -2), color = "darkred", linetype = "dashed") + 
    facet_grid(model ~ protadj, scales = "free_y")
# ggsave("rip2_volcanos.png", width = 16, height = 8)

test_res_prop <- test_res %>% 
    # filter(model == "per-batch", protadj == "protein adjustment") %>%
    filter(model == "all-batch", protadj == "protein adjustment") %>%
    mutate(
        col_sig = ifelse(p_adjusted > 0.05 | abs(log2FC) < 2, "nonsig", ifelse(log2FC > 0, "upreg", "dnreg"))
    )

test_res_prop <- test_res_prop %>% 
    # separate(protsite, c("protein", "site"), sep = "--", remove = FALSE) %>% 
    left_join(hs_fasta %>% distinct(uniprot_iso, entry_name) %>% rename(protein = uniprot_iso)) %>% 
    mutate(protsite2 = str_c(entry_name, " (", site, ")"))

test_res_prop %>% 
    ggplot(aes(log2FC, -log10(p_adjusted), color = col_sig)) + 
    geom_point(aes(color = col_sig), alpha = 0.75, size = 3) +
    geom_hline(yintercept = -log10(0.05), color = "gray20") + 
    scale_color_manual(
        values = c("gray65", "blue", "red"), 
        limits = c("nonsig", "dnreg", "upreg"), 
        labels = c("No significant change", "Down regulation", "Up regulation")
    ) + 
    # geom_text_repel(
    #     data = filter(test_res_prop, !is.na(p_value), p_adjusted < 0.05, abs(log2FC) > 2),
    #     aes(label = protsite2)
    # ) +
    geom_text_repel(
        data = filter(test_res_prop, !is.na(p_value), str_detect(entry_name, "RIPK2|IKBE|NFKB1")),
        aes(label = protsite2), color = "black"
    ) +
    labs(x = "log (base 2) FC", y = "-log (base 10) adjusted p-value", color = "", 
         title = str_c(test_res_prop$contrast[1], " (", test_res_prop$model[1], ")")) + 
    theme(legend.position = "bottom", legend.box = "horizontal")
# ggsave("rip2_volcano_all_key.png", width = 12, height = 8)
# ggsave("rip2_volcano_all_sig.png", width = 12, height = 8)
# ggsave("rip2_volcano_per_key.png", width = 12, height = 8)
# ggsave("rip2_volcano_per_sig.png", width = 12, height = 8)

# Histogram of p-values
test_res %>% 
    group_by(protein, site, contrast) %>% filter(n() == 4) %>% ungroup() %>% 
    ggplot(aes((p_value))) + 
    geom_histogram(binwidth = 0.01) + 
    facet_grid(model + protadj ~ ., scales = "free_y")


# Evaluate modeling strategies --------------------------------------------

test_allbch <- test_res %>% 
    filter(model == "all-batch")

test_allbch_pair <- test_allbch %>% 
    select(protsite, protadj, p_adjusted) %>% 
    spread(protadj, p_adjusted)

sig_null <- test_allbch_pair %>% 
    filter(`no adjustment` < 0.05, `protein adjustment` > 0.05) %>% 
    separate(protsite, into = c("protein", "site"), sep = "--")
sig_adj <- test_allbch_pair %>% 
    filter(`no adjustment` > 0.05, `protein adjustment` < 0.05) %>% 
    separate(protsite, into = c("protein", "site"), sep = "--")
sig_both <- test_allbch_pair %>% 
    filter(`no adjustment` < 0.05, `protein adjustment` < 0.05) %>% 
    separate(protsite, into = c("protein", "site"), sep = "--")

runlvl_bat <- str_c(
    rep(rep(c("ctrl", "mdp_30", "mdp_60", "cmpd89", "cmpd89_mdp_30", "cmpd89_mdp_60"), each = 2), 2), 
    c(rep(c("B1T1", "B1T2"), 6), rep(c("B2T1", "B2T2"), 6)), 
    sep = "-"
)

unnested_site <- nested_site %>% unnest(aftdata) 
# unnested_site <- nested_site %>% unnest(aftdata) %>% filter(is_mod)
subset_site <- function(df_full, test_ex) {
    df_sub <- df_full %>%
        semi_join(test_ex %>% select(protein, site)) %>%
        bind_rows(df_full %>% filter(site == "None", protein %in% test_ex$protein)) %>%
        mutate(is_mod = str_detect(feature, "\\*"))
    
    return(df_sub)
}
sub_perbch <- subset_site(
    unnested_site, 
    test_allbch %>% separate(protsite, into = c("protein", "site"))
)

plot_sprofile <- function(df_allprot, protein, run_level) {
    df_prot <- df_allprot %>% filter(protein == protein) 
    # Complete possible combinations of peptide features and runs
    mpar <- df_prot %>% distinct(site, feature, is_mod)
    df_fill <- df_prot %>% 
        select(feature, run, log2inty) %>% 
        complete(run = run_level, feature) %>% 
        left_join(mpar) %>% 
        mutate(is_mod_fac = factor(ifelse(is_mod, "Modified", "Unmodified")), 
               run_fac = factor(run, levels = run_level))
    # Feature profiles categorized by modification and matching status
    filter(df_fill, is_mod) %>% 
        ggplot(aes(run_fac, log2inty, group = feature, colour = site)) + 
        geom_point(size = 2, alpha = 0.5) +
        geom_line(alpha = 0.75) + 
        geom_line(data = filter(df_fill, !is_mod), aes(run_fac, log2inty, group = feature), colour = "gray") + 
        geom_point(data = filter(df_fill, !is_mod), aes(run_fac, log2inty), colour = "gray", size = 2) + 
        geom_vline(xintercept = 12.5) + 
        facet_grid(is_mod_fac ~ .) + 
        coord_cartesian(ylim = c(10, 35)) + 
        labs(x = "Run", y = "Log2-intensity", title = protein) + 
        theme_bw() + 
        guides(colour = guide_legend(nrow = 1, title = NULL)) + 
        theme(legend.position = c(0.5, 0.065)) + 
        theme(axis.text.x = element_blank())
}

# subprot <- sort(unique(sig_adj$protein))
# pdf("profile_rip2_adj.pdf")
# subprot <- sort(unique(sig_null$protein))
# pdf("profile_rip2_null.pdf")
subprot <- sort(unique(sig_both$protein))
pdf("profile_rip2_both.pdf")
for (i in seq_along(subprot)) {
    print(plot_sprofile(sub_perbch, subprot[i], runlvl_bat))
}
dev.off()


# Significant site changes ------------------------------------------------

test_res <- test_res %>% 
    mutate(protsite = str_replace(protsite, "--", " (") %>% str_c(")")) %>%
    select(protsite, contrast, log2FC, std_error, DF, statistic, p_value, p_adjusted, model, protadj)

if (site_spec) {
    write_csv(test_res, "output/rip2-tmp-4mdl.csv")
} else {
    write_csv(test_res, "output/rip2-tmp-4mdl-prot.csv")
}


test_res_sig <- test_res %>% 
    filter(p_adjusted < 0.05, !is.na(p_value)) %>% 
    filter(abs(log2FC) >= 2)

test_res_sig %>% count(contrast, model, protadj)

# test_resfair_sig <- test_res %>% 
#     group_by(protsite, contrast) %>% filter(n() == 4) %>% ungroup() %>% 
#     filter(p_adjusted < 0.05) %>% 
#     separate(protsite, into = c("protein", "site"), sep = "--", remove = FALSE) %>% 
#     select(protsite, protein, site, contrast, log2FC, std_error, DF, statistic, p_value, p_adjusted, model, protadj)
# 
# test_resfair_sig %>% count(contrast, model, protadj)

# write_csv(test_res_sig, "output/usp30-full.csv")

