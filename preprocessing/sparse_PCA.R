# Gene module discovery from single-cell transcriptomes using sparse PCA.


# Import presets, libraries, and data --------------
source("presets.R")
library(PMA)
library(tibble)

# Immune gene lists downloaded from ImmPort (https://www.immport.org/shared/genelists)
immune_genes <- read_csv("data/immune_genes.csv")
seurat <- readRDS("data/Processed_PNOIT2_seurat.rds")


# Filter down the genes and cells ----------

# Immune + variable genes
genelist <- unique(c(immune_genes$Symbol, seurat@var.genes))
write.csv(genelist, "results/Supplementary_Table_3.csv")

genes <- t(as.matrix(seurat@data[which(rownames(seurat@data) %in% genelist),]))
metadata <- seurat@meta.data
genes_sub <- genes[which(metadata$condition %in% c("CD154","CD137")),]
metadata_sub <- metadata[which(metadata$condition %in% c("CD154","CD137")),]
rm(seurat)

# Remove the lowest-variance genes (lowest 1%)
# This is to remove immune genes that are basically not detected in our data
variances <- apply(genes_sub, 2, var)
genes_sub <- genes_sub[,(variances > quantile(variances, 0.01))]

# Downsample the top 85 (out of 86) samples to have the same number of cells.
# Cuts down on compute time and prevents the PCA from being dominated by a few samples
# (the 86th sample has very few cells so it is excluded from the downsampling)
sample_counts <- metadata_sub %>% count(orig.ident, patient) %>% arrange(n)

set_aside_indices <- which((metadata_sub$orig.ident == sample_counts[1,]$orig.ident) & (metadata_sub$patient == sample_counts[1,]$patient))

downsample <- metadata_sub %>%
  filter(!((orig.ident == sample_counts[1,]$orig.ident) & (patient == sample_counts[1,]$patient))) %>%
  rownames_to_column('rowname') %>%
  group_by(orig.ident, patient) %>%
  sample_n(size = sample_counts$n[2]) %>%
  column_to_rownames('rowname')
  
downsample_indices <- match(rownames(downsample), rownames(metadata_sub))
indices <- c(set_aside_indices, downsample_indices)
genes_sub <- genes_sub[indices,]
metadata_sub <- metadata_sub[indices,]


# Run sparse PCA --------------
# Using methods described in Witten et al, Biostatistics 2009

# Arg choices:
# sumabsv: Penalty to achieve desired level of sparsity. 1.8 selected by tuning.
# orth: set to TRUE to enforce a pseudo-orthogonal constraint on the components.
# niter: # iterations run per component. Default 20, chose 50 to ensure convergence.
out.orth <- SPC(x = scale(genes_sub), sumabsv=1.8, K=100, orth=TRUE, niter=50)
out.orth$genenames <- colnames(genes_sub)

# Invert any majority-negative components to be majority-positive
signs <- apply(out.orth$v, 2, function(x) sum(x<0)/sum(x!=0))
out.orth$v.adj <- out.orth$v
out.orth$v.adj[,(signs>0.5)] <- -out.orth$v.adj[,(signs>0.5)]


# Export processed data ------------

scores <- scale(genes[,match(out.orth$genenames, colnames(genes))]) %*% (out.orth$v.adj)
saveRDS(scores, "data/Allcells_sparsePCA_scores.rds")
saveRDS(out.orth, "data/CD137andCD154downsampled_sparsePCA.rds")

# TODO: export seurat object with gene module scores integrated into metadata

