##### 01 - library #####
library(dplyr)
library(data.table)

gc()
rm(list=ls())


setwd("~/linke/chip_gwas_rev/")



##### 02 - files #####
gs_full <- fread("Data/scDRS/munged/chip_dnmt3a_tet2_munged_n1000_ensg.gs",
                 data.table = FALSE, sep = "\t")

print(names(gs_full))   # should be TRAIT, GENESET
print(nrow(gs_full))    # should be 3 (CHIP, DNMT3A, TET2)

n_sub  <- 500
n_reps <- 20

set.seed(613)

gs_subsample <- do.call(rbind, lapply(1:nrow(gs_full), function(i) {
    trait    <- gs_full$TRAIT[i]
    # Gene set is comma-separated "GENE:SCORE" pairs
    gene_vec <- strsplit(gs_full$GENESET[i], ",")[[1]]
    cat(trait, "— total genes:", length(gene_vec), "\n")
    
    do.call(rbind, lapply(0:(n_reps - 1), function(i_rep) {
        set.seed(i_rep)   # matches Martin's np.random.seed(i_rep)
        genes_sub  <- sample(gene_vec, size = n_sub, replace = FALSE)
        trait_name <- paste0(trait, "_subsample500genes_rep", i_rep)
        data.frame(
            TRAIT   = trait_name,
            GENESET = paste(genes_sub, collapse = ","),
            stringsAsFactors = FALSE
        )
    }))
}))


cat("Subsampled gs rows:", nrow(gs_subsample), "\n")
# Should be 60: 3 traits × 20 reps
print(head(gs_subsample$TRAIT, 10))

# Write full subsampled gs (all 60 rows)
dir.create("Data/scDRS/munged/sensitivity_analysis/", recursive = TRUE, showWarnings = FALSE)
write.table(gs_subsample,
            "Data/scDRS/munged/sensitivity_analysis/chip_dnmt3a_tet2_subsample500genes.gs",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Split into batches of 5 for parallel jobs
# 60 rows / 5 per batch = 12 batches
batch_size <- 5
n_batches  <- ceiling(nrow(gs_subsample) / batch_size)
batch_dir  <- "Data/scDRS/munged/sensitivity_analysis/gene_subsample_batches"
dir.create(batch_dir, recursive = TRUE, showWarnings = FALSE)

for (i_batch in 0:(n_batches - 1)) {
    idx_start <- i_batch * batch_size + 1
    idx_end   <- min((i_batch + 1) * batch_size, nrow(gs_subsample))
    batch_df  <- gs_subsample[idx_start:idx_end, ]
    out_file  <- file.path(batch_dir,
                           sprintf("batch%02d.gs", i_batch))
    write.table(batch_df, out_file,
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Batch", i_batch, ":", nrow(batch_df), "traits →", out_file, "\n")
}
cat("Done writing", n_batches, "batch files\n")


















####### BM subsampling
gs_full <- fread("Data/scDRS/munged/chip_dnmt3a_tet2_munged_n1000.gs",
                 data.table = FALSE, sep = "\t")

print(names(gs_full))   # should be TRAIT, GENESET
print(nrow(gs_full))    # should be 3 (CHIP, DNMT3A, TET2)

n_sub  <- 500
n_reps <- 20

set.seed(613)

gs_subsample <- do.call(rbind, lapply(1:nrow(gs_full), function(i) {
    trait    <- gs_full$TRAIT[i]
    # Gene set is comma-separated "GENE:SCORE" pairs
    gene_vec <- strsplit(gs_full$GENESET[i], ",")[[1]]
    cat(trait, "— total genes:", length(gene_vec), "\n")
    
    do.call(rbind, lapply(0:(n_reps - 1), function(i_rep) {
        set.seed(i_rep)   # matches Martin's np.random.seed(i_rep)
        genes_sub  <- sample(gene_vec, size = n_sub, replace = FALSE)
        trait_name <- paste0(trait, "_subsample500genes_rep", i_rep)
        data.frame(
            TRAIT   = trait_name,
            GENESET = paste(genes_sub, collapse = ","),
            stringsAsFactors = FALSE
        )
    }))
}))


cat("Subsampled gs rows:", nrow(gs_subsample), "\n")
# Should be 60: 3 traits × 20 reps
print(head(gs_subsample$TRAIT, 10))

# Write full subsampled gs (all 60 rows)
write.table(gs_subsample,
            "Data/scDRS/munged/sensitivity_analysis/chip_dnmt3a_tet2_subsample500genes_genesymbol.gs",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Split into batches of 5 for parallel jobs
# 60 rows / 5 per batch = 12 batches
batch_size <- 5
n_batches  <- ceiling(nrow(gs_subsample) / batch_size)
batch_dir  <- "Data/scDRS/munged/sensitivity_analysis/gene_subsample_batches_genesymbol"
dir.create(batch_dir, recursive = TRUE, showWarnings = FALSE)

for (i_batch in 0:(n_batches - 1)) {
    idx_start <- i_batch * batch_size + 1
    idx_end   <- min((i_batch + 1) * batch_size, nrow(gs_subsample))
    batch_df  <- gs_subsample[idx_start:idx_end, ]
    out_file  <- file.path(batch_dir,
                           sprintf("batch%02d.gs", i_batch))
    write.table(batch_df, out_file,
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Batch", i_batch, ":", nrow(batch_df), "traits →", out_file, "\n")
}
cat("Done writing", n_batches, "batch files\n")