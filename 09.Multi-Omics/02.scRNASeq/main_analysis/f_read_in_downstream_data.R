##### 01 - libraries #####
library(data.table)
library(dplyr)
library(glue)
library(ggplot2)

gc()
rm(list=ls())

setwd("~/linke/chip_gwas_rev/")



##### 02 - files #####
results_full <- data.frame()

manifest <- fread("Data/sc_annot_files/file_list.txt", data.table = F, header = F)
files <- gsub("/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/Data/sc_annot_files/", "", manifest$V1)
atlas_name <- gsub(".h5ad","",files)[c(2,4,8)]



for (i in seq_along(atlas_name)){
    folder_name <- glue("Results/{atlas_name[i]}/downstream")
    
    if (i != 3){
        chip <- fread(glue("{folder_name}/CHIP.scdrs_group.anno"), data.table = F) %>%
            mutate(chip = "CHIP", .before=1)
        dnmt3a <- fread(glue("{folder_name}/DNMT3A.scdrs_group.anno"), data.table = F) %>%
            mutate(chip = "DNMT3A", .before=1)
        tet2 <- fread(glue("{folder_name}/TET2.scdrs_group.anno"), data.table = F) %>%
            mutate(chip = "TET2", .before=1)
        
        result <- rbind(chip, dnmt3a, tet2) %>%
            mutate(atlas = atlas_name[i], .before=1)
        results_full <- rbind(results_full, result)
        
    } else if (i == 3){
        chip <- fread(glue("{folder_name}/CHIP.scdrs_group.predicted.celltype.l2"), data.table = F) %>%
            mutate(chip = "CHIP", .before=1)
        dnmt3a <- fread(glue("{folder_name}/DNMT3A.scdrs_group.predicted.celltype.l2"), data.table = F) %>%
            mutate(chip = "DNMT3A", .before=1)
        tet2 <- fread(glue("{folder_name}/TET2.scdrs_group.predicted.celltype.l2"), data.table = F) %>%
            mutate(chip = "TET2", .before=1)
        
        result <- rbind(chip, dnmt3a, tet2) %>%
            mutate(atlas = atlas_name[i], .before=1)
        results_full <- rbind(results_full, result)
    }
    
}

# =============================================================================
# 1. Prep data
# =============================================================================

# FDR correction within each trait
results_plot <- results_full %>%
    # group_by(atlas,chip) %>%
    mutate(
        assoc_fdr   = p.adjust(assoc_mcp,  method = "BH"),
        hetero_fdr  = p.adjust(hetero_mcp, method = "BH")
    ) %>%
    ungroup() %>%
    mutate(
        neg_log10_fdr = -log10(assoc_fdr + 1e-10),   # add small value to avoid -Inf
        sig_label     = case_when(
            assoc_fdr < 0.05 ~ "FDR < 0.05",
            assoc_fdr < 0.20 ~ "FDR < 0.20",
            TRUE             ~ "NS"
        ),
        pct_sig_cells = n_fdr_0.05 / n_cell * 100
    )


# =============================================================================
# 2. Dot plot: Z-score × cell type × trait
# =============================================================================

# Order cell types by mean Z-score across traits
celltype_order <- results_plot %>%
    group_by(group) %>%
    summarise(mean_z = mean(assoc_mcz, na.rm = TRUE)) %>%
    arrange(mean_z) %>%
    pull(group)

results_plot <- results_plot %>%
    mutate(group = factor(group, levels = celltype_order))

p_dot <- ggplot(results_plot,
                aes(x = chip, y = group,
                    color = assoc_mcz,
                    size  = pct_sig_cells,
                    shape = sig_label)) +
    geom_point(alpha = 0.85) +
    scale_color_gradient2(
        low      = "blue",
        mid      = "white",
        high     = "red",
        midpoint = 0,
        name     = "Association\nZ-score"
    ) +
    scale_size_continuous(
        range = c(1, 8),
        name  = "% sig cells\n(FDR < 0.05)"
    ) +
    scale_shape_manual(
        values = c("FDR < 0.05" = 19,   # filled circle
                   "FDR < 0.20" = 21,   # open circle
                   "NS"         = 4),   # X
        name   = "Significance"
    ) +
    facet_wrap(~ atlas, scales = "free_y") +
    labs(
        title = "scDRS enrichment: CHIP GWAS gene sets",
        x     = "Trait",
        y     = "Cell type"
    ) +
    theme_bw() +
    theme(
        axis.text.x  = element_text(angle = 45, hjust = 1),
        axis.text.y  = element_text(size  = 8),
        strip.text   = element_text(face  = "bold"),
        legend.position = "right"
    )


# ggsave("Results/figures/main_analysis/scDRS_dotplot.png", p_dot, width = 14, height = 10, dpi = 300)
print(p_dot)

