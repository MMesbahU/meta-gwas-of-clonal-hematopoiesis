##### 01 - library #####
library(data.table)
library(dplyr)
library(ggplot2)
library(ggh4x)

gc()
rm(list = ls())


setwd("~/linke/chip_gwas_rev/")

##### 02 - files #####
onk1k_main     <- fread("Results/sensitivity_analysis/onek1k_cellsubsample_main_results.tsv",  data.table = F)
onk1k_cellsensi <- fread("Results/sensitivity_analysis/onek1k_cellsubsample_n20sensi_results.tsv", data.table = F)
bm_main        <- fread("Results/sensitivity_analysis/bm_cellsubsample_main_results.tsv",  data.table = F)
bm_cellsensi   <- fread("Results/sensitivity_analysis/bm_cellsubsample_n20sensi_results.tsv", data.table = F)

##### 03 - restrict to 3 traits #####
traits_keep <- c("CHIP", "DNMT3A", "TET2")

##### 04 - z-score equivalents of FDR thresholds #####
# Two-sided: a p-value of 0.2 / 0.05 corresponds to +/- these z magnitudes
z_thresh_02  <- qnorm(1 - 0.2 / 2)
z_thresh_005 <- qnorm(1 - 0.05 / 2)


##### 05 - reusable combined-panel plotting function #####
plot_combined_beta <- function(main_df, sub_df, atlas_label, out_file,
                               width = 14, height = 8) {
    
    main_df <- main_df %>% filter(trait %in% traits_keep)
    sub_df  <- sub_df  %>% filter(trait %in% traits_keep)
    
    # Order cell types by mean z across the 3 traits (based on main result)
    global_order <- main_df %>%
        filter(trait == "CHIP") %>%
        group_by(group) %>%
        summarise(mean_z = mean(assoc_mcz, na.rm = TRUE), .groups = "drop") %>%
        arrange(mean_z) %>%
        pull(group)
    
    combined_sub <- sub_df %>%
        mutate(
            group = factor(group, levels = global_order),
            trait = factor(trait, levels = traits_keep)
        )
    
    combined_main <- main_df %>%
        mutate(
            group = factor(group, levels = global_order),
            trait = factor(trait, levels = traits_keep)
        )
    
    p_combined <- ggplot() +
        geom_boxplot(
            data  = combined_sub,
            aes(x = group, y = assoc_mcz),
            fill  = "lightblue",
            color = "grey40",
            outlier.size  = 0.3,
            outlier.alpha = 0.4,
            width = 0.6
        ) +
        geom_point(
            data  = combined_main,
            aes(x = group, y = assoc_mcz, shape = sig_label),
            color = "red",
            size  = 1.8
        ) +
        geom_hline(yintercept = 0,  linetype = "solid", color = "darkred", linewidth = 0.4) +
        # # FDR = 0.2 lines (both directions, since z can be negative)
        # geom_hline(yintercept = z_thresh_02,  linetype = "dashed", color = "darkred", linewidth = 0.4) +
        # geom_hline(yintercept = -z_thresh_02, linetype = "dashed", color = "darkred", linewidth = 0.4) +
        # # FDR = 0.05 lines
        # geom_hline(yintercept = z_thresh_005,  linetype = "dashed", color = "darkred", linewidth = 0.4) +
        # geom_hline(yintercept = -z_thresh_005, linetype = "dashed", color = "darkred", linewidth = 0.4) +
        scale_shape_manual(
            values = c("FDR < 0.05" = 16,
                       "FDR < 0.20" = 17,
                       "P < 0.01"   = 15,
                       "NS"         = 4),
            name   = "Main result"
        ) +
        coord_flip() +
        ggh4x::facet_grid2(. ~ trait, scales = "free_x") +
        labs(
            x = "Cell type",
            y = "scDRS association z-score (beta)"
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
    ggsave(out_file, p_combined, width = width, height = height, dpi = 300)

}

##### 06 - run for BM and OneK1K #####
plot_combined_beta(
    main_df  = bm_main,
    sub_df   = bm_cellsensi,
    out_file = "Results/figures/bm_sensitivity/bm_cell_subsample_combined_beta.png"
)

plot_combined_beta(
    main_df  = onk1k_main,
    sub_df   = onk1k_cellsensi,
    out_file = "Results/figures/onek1k_sensitivity/onek1k_cell_subsample_combined_beta.png"
)