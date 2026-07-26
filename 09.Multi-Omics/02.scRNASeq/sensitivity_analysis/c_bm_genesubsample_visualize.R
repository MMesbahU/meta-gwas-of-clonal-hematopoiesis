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
# MAIN RESULTS — load results
# =============================================================================
traits <- c("CHIP", "DNMT3A", "TET2")

load_downstream <- function(path, trait, rep = NA, subsample_type = NA) {
    f <- file.path(path, paste0(trait, ".scdrs_group.anno"))
    
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

main_dir    <- "Results/BM_standard_design/downstream/"
# Main result — one row per cell type per trait
main_results <- do.call(rbind, lapply(traits, function(trait) {
    load_downstream(main_dir, trait,
                    rep = NA, subsample_type = "main")
}))

main_results <- main_results %>%
    group_by(trait, rep) %>%
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


# =============================================================================
# GENE SUBSAMPLING — load results
# =============================================================================
gene_sub_dir <- "Results/sensitivity_analysis/bm_500genes_sub/"
n_batches    <- 12
batch_size   <- 5

# Trait naming in gene subsampling: CHIP_subsample500_rep0, etc.
# Need to map back to original trait names
get_base_trait <- function(trait_name) {
    # "CHIP_subsample500_rep3" -> "CHIP"
    gsub("_subsample500genes_rep\\d+$", "", trait_name)
}

get_rep_num <- function(trait_name) {
    # "CHIP_subsample500_rep3" -> 3
    as.integer(sub(".*_subsample500genes_rep(\\d+)$", "\\1", trait_name))
}

# Load from batch folders — each batch has multiple traits
gene_sub_results <- do.call(rbind, lapply(0:(n_batches - 1), function(i_batch) {
    batch_dir <- file.path(gene_sub_dir,
                           paste0("batch", i_batch),
                           "downstream")
    
    # Find all score files in this batch's downstream folder
    # Trait names look like: CHIP_subsample500_rep0.scdrs_group.predicted.celltype.l2
    score_files <- list.files(
        batch_dir,
        pattern = "\\.scdrs_group\\.anno$",
        full.names = FALSE
    )
    
    if (length(score_files) == 0) {
        warning("No score files in: ", batch_dir)
        return(NULL)
    }
    
    do.call(rbind, lapply(score_files, function(fname) {
        # Extract trait name from filename
        trait_full <- sub("\\.scdrs_group\\.anno$", "", fname)
        
        df <- fread(file.path(batch_dir, fname), data.table = FALSE)
        if (!"group" %in% names(df)) df <- df %>% rename(group = 1)
        
        df$trait_full    <- trait_full
        df$trait         <- get_base_trait(trait_full)
        df$rep           <- get_rep_num(trait_full)
        df
    }))
}))


cat("Gene subsample rows:", nrow(gene_sub_results), "\n")
cat("Unique traits:", unique(gene_sub_results$trait), "\n")
cat("Reps per trait:", table(gene_sub_results$trait[!duplicated(
    paste(gene_sub_results$trait, gene_sub_results$rep))]), "\n")


# =============================================================================
# FDR correction for gene subsampling — per trait per rep
# =============================================================================
gene_sub_results <- gene_sub_results %>%
    group_by(trait, rep) %>%
    mutate(assoc_fdr = p.adjust(assoc_mcp, method = "BH")) %>%
    ungroup()


# =============================================================================
# function section for plotting
# =============================================================================
get_celltype_order <- function(trait_name) {
    main_results %>%
        filter(trait == trait_name) %>%
        arrange(assoc_fdr) %>%
        pull(group)
}


# =============================================================================
# Plot — gene subsampling version (same style as cell subsampling)
# =============================================================================
plot_sensitivity_gene <- function(trait_name) {
    
    ct_order <- get_celltype_order(trait_name)
    
    sub_df <- gene_sub_results %>%
        filter(trait == trait_name) %>%
        mutate(
            neg_log10_p = -log10(assoc_fdr + 1e-10),
            group = factor(group, levels = rev(ct_order))
        )
    
    main_df <- main_results %>%
        filter(trait == trait_name) %>%
        mutate(
            neg_log10_p = -log10(assoc_fdr + 1e-10),
            group = factor(group, levels = rev(ct_order))
        )
    
    fdr_line_02 <- -log10(0.2)
    fdr_line_005 <- -log10(0.05)
    
    ggplot() +
        geom_boxplot(
            data    = sub_df,
            aes(x = group, y = neg_log10_p),
            fill    = "lightyellow",    # different color from cell subsample
            color   = "grey40",
            outlier.size  = 0.5,
            outlier.alpha = 0.5,
            width   = 0.6
        ) +
        geom_point(
            data  = main_df,
            aes(x = group, y = neg_log10_p, shape = sig_label),
            color = "red",
            size  = 2.5,
            alpha = 0.9
        ) +
        geom_hline(yintercept = -log10(0.2),
                   linetype = "dashed", color = "darkred", linewidth = 0.5) +
        geom_hline(yintercept = -log10(0.05),
                   linetype = "dashed", color = "darkred", linewidth = 0.5) +
        annotate("text", x = 1, y = -log10(0.2) - 0.02,
                 label = "FDR = 0.2",  color = "darkred", size = 3, hjust = 0) +
        annotate("text", x = 1, y = -log10(0.05) - 0.02,
                 label = "FDR = 0.05", color = "darkred", size = 3, hjust = 0) +
        scale_shape_manual(
            values = c("FDR < 0.05" = 16, "FDR < 0.20" = 17,
                       "P < 0.01" = 15, "NS" = 4),
            name   = "Main result"
        ) +
        coord_flip() +
        labs(
            title    = paste0("BoneMarrow — ", trait_name),
            subtitle = paste0("Gene subsampling (n=", n_reps,
                              " reps × 500 genes from top 1000)"),
            x = "Cell type",
            y = expression(-log[10](p))
        ) +
        theme_classic() +
        theme(
            axis.text.y     = element_text(size = 8),
            axis.text.x     = element_text(size = 8),
            plot.title      = element_text(face = "bold", size = 11),
            plot.subtitle   = element_text(size = 8, color = "grey40"),
            legend.position = "bottom"
        )
}

## variable section
n_reps <- 20

## starts to plot
for (trait in traits) {
    p <- plot_sensitivity_gene(trait)
    ggsave(
        paste0("Results/figures/bm_sensitivity/bm_gene_subsample_",
               trait, ".png"),
        p, width = 7, height = 9, dpi = 300
    )
    print(p)
    cat("Saved gene subsample:", trait, "\n")
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

combined_sub <- gene_sub_results %>%
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
        fill  = "lightyellow",
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
        # title    = "BoneMarrow sensitivity analysis — cell subsampling",
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
ggsave("Results/figures/bm_sensitivity/bm_gene_subsample_combined.png",
       p_combined, width = 14, height = 8, dpi = 300)



