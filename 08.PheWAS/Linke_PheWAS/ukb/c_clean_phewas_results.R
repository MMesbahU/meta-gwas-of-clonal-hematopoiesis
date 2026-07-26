###### 01 - library ######
library(dplyr)
library(data.table)
library(glue)
library(ggplot2)
library(ggrepel)
library(ggh4x)

gc()
rm(list=ls())

setwd("/medpop/esp2/lli/chip_gwas_rev/")


##### 02 - files #####
phecodes_list <- fread("Data/ukb_phewas/phecode_list.tsv", data.table = F, header=F)$V1

# i=1
results_full <- data.frame()
for (i in seq_along(phecodes_list)){
    phecode <- phecodes_list[i]
    file_name <- glue("Results/ukb_phewas/HR_ukbb_mesbahCHIP_Phewas_CoxPH.PheCode{phecode}.tsv")
    
    # only EUR
    dat <- fread(file_name, data.table = F) %>%
        mutate(phecode_cat = phecode, .before=1)
    
    results_full <- rbind(results_full, dat)
}



write.table(results_full, "Results/ukb_phewas/ukbb_phewas_results_full.tsv", sep="\t", row.names = F, quote = F)







##### 03 - plotting #####
gc()
rm(list=ls())

results_full <- fread("Results/ukb_phewas/ukbb_phewas_results_full.tsv", data.table = FALSE)

# ---- 1. Restrict to CHIP exposures ----
plot_df <- results_full %>%
    mutate(
        pval_adj     = p.adjust(pval, method = "BH"),
        neg_log10_p  = -log10(pval_adj),
        direction    = ifelse(hr > 1, "Risk", "Protective")
    )

# ---- 2. Order phenotypes within category for a clean x-axis ----
plot_df <- plot_df %>%
    arrange(phecode_cat, outcomes) %>%
    mutate(
        outcomes = gsub("_", " ", outcomes),
        outcomes  = factor(outcomes, levels = unique(outcomes)),
        knn       = factor(knn, levels = c("EUR", "AFR", "AMR", "EAS", "SAS")),
        exposures = factor(exposures, levels = c("CHIP", "expandedCHIP", "DNMT3A", "expandedDNMT3A",
                                                 "TET2", "expandedTET2", "ASXL1", "expandedASXL1",
                                                 "JAK2", "expandedJAK2", "DDR", "expandedDDR",
                                                 "Splice", "expandedSplice"),
                                      labels = c("CHIP", "Large CHIP", "DNMT3A", "Large DNMT3A",
                                                 "TET2", "Large TET2", "ASXL1", "Large ASXL1",
                                                 "JAK2", "Large JAK2", "DDR", "Large DDR",
                                                 "Splice", "Large Splice")
        )
    )

# ---- 3. Significance threshold ----
sig_threshold <- -log10(0.05)

# ---- 4. Build one plot per ancestry, faceted by CHIP subtype ----
make_phewas_plot <- function(df, ancestry_label) {
    
    sub_df   <- df %>% filter(knn == ancestry_label)
    label_df <- sub_df %>% 
        filter(neg_log10_p > sig_threshold) %>%
        group_by(exposures) %>%
        slice_max(neg_log10_p, n = 5) %>%
        ungroup()
    
    ggplot(sub_df, aes(x = outcomes, y = neg_log10_p, fill = phecode_cat)) +
        geom_point(aes(shape = direction), color = "black", stroke = 0.3, alpha = 0.9) +
        geom_hline(yintercept = sig_threshold, linetype = "dashed", color = "grey40") +
        geom_text_repel(
            data = label_df,
            aes(label = outcomes),
            size = 2.8, max.overlaps = 15, show.legend = FALSE, inherit.aes = TRUE, fontface = "bold",
            min.segment.length = 0
        ) +
        facet_nested_wrap(.~exposures, scales = "free_y", ncol = 1) +
        scale_shape_manual(values = c("Risk" = 24, "Protective" = 25), breaks = c("Risk", "Protective")) + 
        labs(
            x = NULL, y = expression(-log[10](p[BH])),
            fill = "Phecode Category", shape = "Direction",
            title = glue("{ancestry_label} PheWAS")
        ) +
        theme_classic() +
        guides(
            fill  = guide_legend(nrow = 2, 
                                 override.aes = list(shape = 24, size = 3)),
            shape = guide_legend(nrow = 1)
        ) +
        theme(
            axis.text.x  = element_blank(),
            axis.ticks.x = element_blank(),
            strip.text.y = element_text(angle = 0),
            legend.position = "bottom"
        )
}

ancestries <- levels(plot_df$knn)

plots <- lapply(ancestries, function(anc) make_phewas_plot(plot_df, anc))
names(plots) <- ancestries

# ---- 5. Save each ancestry's plot ----
for (anc in ancestries) {
    ggsave(
        filename = glue("Results/ukb_phewas/figures/phewas_plot_{anc}.png"),
        plot = plots[[anc]],
        width = 12, height = 17, dpi = 300
    )
}

