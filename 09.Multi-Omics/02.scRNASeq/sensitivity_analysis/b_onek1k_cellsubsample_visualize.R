##### 01 - library #####
library(dplyr)
library(data.table)
library(ggplot2)
library(ggh4x)

gc()
rm(list=ls())


setwd("~/linke/chip_gwas_rev/")


##### 02 - files #####
# =============================================================================
# 1. Load main result (full OneK1K, full gene set)
# =============================================================================

load_downstream <- function(path, trait, rep = NA, subsample_type = NA) {
    f <- file.path(path, paste0(trait, ".scdrs_group.predicted.celltype.l2"))
    
    ## readin ing
    df <- fread(f, data.table = FALSE)
    df$trait          <- trait
    df$rep            <- rep
    df$subsample_type <- subsample_type
    # rename index column to cell_type if needed
    if (!"group" %in% names(df)) {
        df <- df %>% rename(group = 1)
    }
    df
}

traits      <- c("CHIP", "DNMT3A", "TET2")
main_dir    <- "Results/onek1k_sc_filtered/downstream/"
cell_sub_dir <- "Results/sensitivity_analysis/onek1k_50kcells_sub/"
n_reps      <- 20

# Main result — one row per cell type per trait
main_results <- do.call(rbind, lapply(traits, function(trait) {
    load_downstream(main_dir, trait,
                    rep = NA, subsample_type = "main")
}))


# =============================================================================
# 2. Load cell subsampling results (20 reps)
# =============================================================================

cell_sub_results <- do.call(rbind, lapply(0:(n_reps - 1), function(i_rep) {
    rep_dir <- file.path(cell_sub_dir, paste0("rep", i_rep), "downstream")
    do.call(rbind, lapply(traits, function(trait) {
        load_downstream(rep_dir, trait,
                        rep = i_rep, subsample_type = "cell_subsample")
    }))
}))


# =============================================================================
# 3. FDR correction — within each atlas × trait, matching your main analysis
# =============================================================================

# Main result FDR
main_results <- main_results %>%
    group_by(trait) %>%
    mutate(assoc_fdr = p.adjust(assoc_mcp, method = "BH")) %>%
    ungroup() %>%
    mutate(
        sig_label = case_when(
            assoc_fdr < 0.05 ~ "FDR < 0.05",
            assoc_fdr < 0.20 ~ "FDR < 0.20",
            assoc_mcp < 0.01 ~ "P < 0.01",
            TRUE             ~ "NS"
        )
    )

# Cell subsample FDR — per rep per trait
cell_sub_results <- cell_sub_results %>%
    group_by(trait, rep) %>%
    mutate(assoc_fdr = p.adjust(assoc_mcp, method = "BH")) %>%
    ungroup()


write.table(main_results, "Results/sensitivity_analysis/onek1k_cellsubsample_main_results.tsv", 
            sep="\t", row.names = F, quote = F)
write.table(cell_sub_results, "Results/sensitivity_analysis/onek1k_cellsubsample_n20sensi_results.tsv", 
            sep="\t", row.names = F, quote = F)

# =============================================================================
# 4. Order cell types by main result Z-score (within each trait)
# =============================================================================

# For each trait, order cell types by main result assoc_mcz descending
get_celltype_order <- function(trait_name) {
    main_results %>%
        filter(trait == trait_name) %>%
        arrange(assoc_fdr) %>%
        pull(group)
}


# =============================================================================
# 5. Figure 5b equivalent — boxplot of -log10(p) across 20 reps
#    with main result overlaid as red dot
# =============================================================================

plot_sensitivity <- function(trait_name, subsample_label = "Cell subsampling") {
    
    # Cell types ordered by main result for this trait
    ct_order <- get_celltype_order(trait_name)
    
    # Subsample data for this trait
    sub_df <- cell_sub_results %>%
        filter(trait == trait_name) %>%
        mutate(
            neg_log10_p = -log10(assoc_fdr + 1e-10),
            group = factor(group, levels = rev(ct_order))
        )
    
    # Main result for this trait
    main_df <- main_results %>%
        filter(trait == trait_name) %>%
        mutate(
            neg_log10_p = -log10(assoc_fdr + 1e-10),
            group = factor(group, levels = rev(ct_order))
        )
    
    # FDR 0.2 threshold line — same as paper
    # Martin used FDR < 0.2 across all cell types
    fdr_line_02 <- -log10(0.2)
    fdr_line_005 <- -log10(0.05)
    
    ggplot() +
        # Boxplot from 20 subsampling reps
        geom_boxplot(
            data    = sub_df,
            aes(x = group, y = neg_log10_p),
            fill    = "lightblue",
            color   = "grey40",
            outlier.size  = 0.5,
            outlier.alpha = 0.5,
            width   = 0.6
        ) +
        # Main result as red dot
        geom_point(
            data  = main_df,
            aes(x = group, y = neg_log10_p,
                shape = sig_label),
            color = "red",
            size  = 2.5,
            alpha = 0.9
        ) +
        # FDR threshold line
        geom_hline(
            yintercept = fdr_line_02,
            linetype   = "dashed",
            color      = "darkred",
            linewidth  = 0.5
        ) +
        # FDR threshold line
        geom_hline(
            yintercept = fdr_line_005,
            linetype   = "dashed",
            color      = "darkred",
            linewidth  = 0.5
        ) +
        annotate("text",
                 x     = 1,
                 y     = fdr_line_02-0.02,
                 label = "FDR = 0.2",
                 color = "darkred",
                 size  = 3,
                 hjust = 0) +
        annotate("text",
                 x     = 1,
                 y     = fdr_line_005-0.02,
                 label = "FDR = 0.05",
                 color = "darkred",
                 size  = 3,
                 hjust = 0) +
        scale_shape_manual(
            values = c("FDR < 0.05" = 16,
                       "FDR < 0.20" = 17,
                       "P < 0.01"   = 15,
                       "NS"         = 4),
            name   = "Main result"
        ) +
        coord_flip() +
        labs(
            title    = paste0("OneK1K — ", trait_name),
            subtitle = paste0(subsample_label, " (n=", n_reps, " reps × 50K cells)"),
            x        = "Cell type",
            y        = expression(-log[10](p))
        ) +
        theme_classic() +
        theme(
            axis.text.y     = element_text(size = 8),
            axis.text.x     = element_text(size = 8),
            plot.title      = element_text(face = "bold", size = 11),
            plot.subtitle   = element_text(size = 8, color = "grey40"),
            legend.position = "bottom",
            legend.text     = element_text(size = 8),
            panel.grid.minor = element_blank()
        )
}

# Generate one plot per trait (just need to run once)
# dir.create("Results/figures/onek1k_sensitivity", recursive = TRUE, showWarnings = FALSE)

for (trait in traits) {
    p <- plot_sensitivity(trait)
    ggsave(
        paste0("Results/figures/onek1k_sensitivity/onek1k_cell_subsample_", trait, ".png"),
        p, width = 7, height = 9, dpi = 300
    )
    print(p)
    cat("Saved:", trait, "\n")
}


# =============================================================================
# 6. Combined three-trait panel (matches Martin's Figure 5b layout)
# =============================================================================

# Order by mean Z across all traits
global_order <- main_results %>%
    group_by(group) %>%
    summarise(mean_z = mean(assoc_mcz, na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_z) %>%
    pull(group)

combined_sub <- cell_sub_results %>%
    mutate(
        neg_log10_p = -log10(assoc_fdr + 1e-10),
        group = factor(group, levels = global_order),
        trait = factor(trait, levels = traits)
    )

combined_main <- main_results %>%
    mutate(
        neg_log10_p = -log10(assoc_fdr + 1e-10),
        group = factor(group, levels = global_order),
        trait = factor(trait, levels = traits)
    )

fdr_line_02 <- -log10(0.2)
fdr_line_005 <- -log10(0.05)

p_combined <- ggplot() +
    geom_boxplot(
        data  = combined_sub,
        aes(x = group, y = neg_log10_p),
        fill  = "lightblue",
        color = "grey40",
        outlier.size  = 0.3,
        outlier.alpha = 0.4,
        width = 0.6
    ) +
    geom_point(
        data  = combined_main,
        aes(x = group, y = neg_log10_p,
            shape = sig_label),
        color = "red",
        size  = 1.8
    ) +
    geom_hline(
        yintercept = fdr_line_02,
        linetype   = "dashed",
        color      = "darkred",
        linewidth  = 0.4
    ) +
    geom_hline(
        yintercept = fdr_line_005,
        linetype   = "dashed",
        color      = "darkred",
        linewidth  = 0.4
    ) +
    scale_shape_manual(
        values = c("FDR < 0.05" = 16,
                   "FDR < 0.20" = 17,
                   "P < 0.01"   = 15,
                   "NS"         = 4),
        name   = "Main result"
    ) +
    coord_flip() +
    ggh4x::facet_grid2(.~ trait, scales = "free_x") +
    labs(
        # title    = "OneK1K sensitivity analysis — cell subsampling",
        # subtitle = paste0(n_reps, " repetitions × 50K cells | Red = main result (full atlas)"),
        x        = "Cell type",
        y        = expression(-log[10](p))
    ) +
    theme_classic() +
    theme(
        axis.text.y      = element_text(size = 7),
        axis.text.x      = element_text(size = 7),
        strip.text       = element_text(face = "bold", size = 9),
        legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        panel.spacing    = unit(0.5, "lines")
    )

print(p_combined)
ggsave("Results/figures/onek1k_sensitivity/onek1k_cell_subsample_combined.png",
       p_combined, width = 14, height = 8, dpi = 300)



