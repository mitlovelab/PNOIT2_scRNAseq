# Figure 6: Treg phenotypes are not significantly modulated by OIT.


# Pre-subsetting down the right cells ----

# library(Seurat)
# library(dplyr)
# seurat = readRDS('data/Processed_PNOIT2_seurat.RDS')
# data_with_modules = readRDS('data/tcr_scores_genes_gates.RDS')
# cells = as.data.frame(as.matrix(t(seurat@data)))
# cell_metadata = seurat@meta.data
# gene_names = rownames(data_with_modules@data)
# gene_names = c(gene_names, "IL10", "TGFB1", "TGFB2", "TGFB3")
# cells2 = cells[,which(colnames(cells) %in% gene_names)]
# overall_cells = cbind(cells2, cell_metadata)
# final_cells = left_join(overall_cells, data_with_modules@meta.data, by=c('Barcode','patient','timept','condition','group','outcome'))
# final_cells = final_cells %>% rename(TRB_CDR3 = TRB_CDR3.y) %>% select(-TRB_CDR3.x)
# write.csv(final_cells, "data/cells_100k_with_modules_and_IL10.csv")
# final_cells_treg = final_cells[which(final_cells$MODULE_1 > 1),]
# write.csv(final_cells_treg, "../data/treg_cells_100k_with_modules_and_IL10.csv")
# # For 6A and D, keep all clonotypes with at least one Treg cell
# final_cells = final_cells %>% filter(condition != 'DblNeg', outcome != 'placebo')
# df_treg = final_cells %>% filter(MODULE_1 > 2)
# treg_clonotypes <- levels(factor(df_treg$TRB_CDR3))
# treg_clones = final_cells %>% filter(TRB_CDR3 %in% treg_clonotypes)
# treg_clones = treg_clones %>% filter(!is.na(MODULE_1))
# write.csv(treg_clones, 'data/treatment_treg_clones_100k.csv')



# Load libraries and presets ----

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(Seurat)
library(pheatmap)
library(dendsort)

source("utils.R")
source("presets.R")

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 3)/3, symbols = c("****", "***", "**", "*", "ns"))


# Figure 6A: Make box+swarmplots for the Treg figure ----

df = read.csv('data/treatment_treg_clones_100k.csv')

df_by_patient = df %>% 
  group_by(outcome, patient, timept, condition) %>%
  summarize(meantreg = mean(MODULE_1),
            meanFOXP3 = mean(FOXP3),
            meanIL10 = mean(IL10)) %>%
  mutate(timept = factor(timept, levels=c('BL','BU','MN','AV')))

df_CD154 = df_by_patient %>% filter(condition == 'CD154')
df_CD137 = df_by_patient %>% filter(condition == 'CD137')

df_overall = df %>% 
  filter(patient != 'P105') %>% # remove the patient we didn't have CD137 for
  group_by(outcome, patient, timept) %>%
  summarize(meantreg = mean(MODULE_1),
            meanFOXP3 = mean(FOXP3),
            meanIL10 = mean(IL10)) %>%
  mutate(timept = factor(timept, levels=c('BL','BU','MN','AV')))

timept = 'timept'

for (gene in c('meantreg','meanFOXP3','meanIL10')){
  ggplot(data = df_overall, aes_string(x = timept, y = gene)) + 
    geom_boxplot(lwd = 0.1, outlier.shape = NA) + 
    geom_beeswarm(aes(color = outcome), stroke=0, pch=16, size=1.8, cex=4.5) +
    stat_compare_means(comparisons = list(c(1,2), c(1,3), c(1,4)),  
                       method = 'wilcox', paired = FALSE, size = 2.1, symnum.args = symnum.args) +
    labs(x = "Timepoint", y = gene) + 
    scale_color_manual(values = patientpalette[c(2,5,8)]) +
    ggpreset + 
    theme(text = element_text(size = 8), axis.text = element_text(size = 8))
  ggsave(paste(gene,"_6A_by_patient.pdf", sep=""), width = 2, height = 2.2, unit = "in", family = "ArialMT", useDingbats=FALSE)
}


# Cluster the Tregs - only need to run this once, can skip to next section ----

df_treg = read.csv('data/treg_cells_100k_with_modules_and_IL10.csv')

treg = CreateSeuratObject(setNames(data.frame(t(df_treg[,2:556])), df_treg[,1]))
treg = FindVariableFeatures(treg)
treg@meta.data$patient = df_treg$patient
treg@meta.data$TRB_CDR3 = df_treg$TRB_CDR3
treg@meta.data$nUMI = df_treg$nUMI
treg@meta.data$MODULE_1 = df_treg$MODULE_1
treg@meta.data$timept = df_treg$timept
treg@meta.data$outcome = df_treg$outcome
treg@meta.data$Barcode = df_treg$Barcode
treg = ScaleData(treg, vars.to.regress = 'nUMI', model.use = 'poisson')
treg = RunPCA(treg, pcs.compute = 50, do.print = FALSE)
treg = RunUMAP(treg, dims=1:10)
treg@meta.data = cbind(treg@meta.data, treg@reductions$umap@cell.embeddings)

ggplot(treg@meta.data, aes(x=UMAP_1, y=UMAP_2, col=seurat_clusters)) + 
  scale_color_manual(values=c('#e6ab02','#a6761d','#386cb0','#666666')) + 
  geom_point(size=0.6, stroke=0) + ggpreset
ggsave("results/Figure6_Treg_seuratclusters_UMAP.pdf", width=3, height=3, unit = "in", family = "ArialMT", useDingbats =FALSE)

treg <- FindNeighbors(treg, dims = 1:20)
treg <- FindClusters(treg, resolution = 0.2)

a = FindMarkers(treg, 2, test.use = 't', logfc.threshold = 0.5)
write.csv(a, 'results/treg_cluster0_DEgenes.csv')

treg@meta.data$cluster = "Conventional Treg"
treg@meta.data$cluster[treg@meta.data$seurat_clusters == 1] = "Tfh-like Treg"
treg@meta.data$cluster[treg@meta.data$seurat_clusters == 2] = "IL7R+ Treg"

treg@meta.data$TGFB1 = df_treg$TGFB1
treg@meta.data$IL10 = df_treg$IL10
treg@meta.data$FOXP3 = df_treg$FOXP3
# only need to save this once:
#saveRDS(treg, 'data/treg_seurat_new.RDS')

# Figure 6B: UMAP of Treg clusters ----
treg = readRDS('data/treg_seurat_new.RDS')

ggplot(treg@meta.data, aes(x=UMAP_1, y=UMAP_2, col=cluster)) + 
  scale_color_manual(values=c('#e6ab02','#386cb0','#666666')) + 
  geom_point(size=0.5, stroke=0) + ggpreset
ggsave("results/Figure6_Treg_clusters_UMAP.pdf", width=2.5, height=2.5, unit = "in", family = "ArialMT", useDingbats =FALSE)


# Figure 6C: Heatmap of DE genes across Treg clusters ----

# Manually curated from DE genes list + Treg genes, in desired plotting order
genes = c(
  # Conventional Tregs
  "ANTXR2","CCR5","HLA.DRB1","FOXP3","IL2RA","TGFB1",
  # Thf-like Tregs
  "TNFRSF4","TFRC","IL10","CTLA4","ICOS","NR4A3","REL","TNFRSF9",
  # IL7R+ Tregs
  "CCR7","SESN3","IL7R","GIMAP5","TCF7","ITGA4"
)


df = treg@assays$RNA[genes,] %>% as.matrix() %>% t() %>% as.data.frame()
df$cluster = factor(treg@meta.data$cluster, levels=c("Conventional Treg","Tfh-like Treg","IL7R+ Treg"))
df$patient = treg@meta.data$patient
df = df %>% group_by(cluster, patient) %>% summarize_all(~mean(.)) %>% as.data.frame() 
rownames(df) = paste(df$cluster, df$patient)
meta = data.frame(Cluster = df$cluster, row.names = rownames(df))
df = df %>% select(-c("cluster","patient"))
df = scale(df)
# cap outlier z-scores:
df[df > 3] = 3
df[df < -2] = -2
ann_colors = list(Cluster = c("Conventional Treg" = "#e6ab02",
                              "Tfh-like Treg" = "#666666",
                              "IL7R+ Treg" = "#386cb0"))

pheatmap(t(df), color = inferno(100), 
         annotation_col = meta, annotation_colors = ann_colors, 
         cluster_cols= FALSE, cluster_rows= FALSE, show_colnames = FALSE, 
         fontsize = 8, treeheight_row = 6, 
         width = 3.3, height = 2.9, 
         border_color = NA, filename = "results/Treg_DE_genes.pdf")


# Figure 6D: Trends within Treg clusters over time ----

# Use the revised definition of Tregs that includes clones and modules:
df = read.csv('data/treatment_treg_clones_100k.csv')
df$cluster = 'none'

# Assign each clonotype to a Treg cluster based on the plurality of its module-expressing cells
treg_meta = readRDS('data/treg_seurat_new.RDS')@meta.data

treg_clonotypes = levels(factor(df$TRB_CDR3))
df$TRB_CDR3 = as.character(df$TRB_CDR3)
treg_meta$TRB_CDR3 = as.character(treg_meta$TRB_CDR3)

for (clone in treg_clonotypes){
  treg_meta_clone = treg_meta %>% filter(TRB_CDR3 == clone)
  if (nrow(treg_meta_clone) > 0){
    max_cluster = names(sort(-table(treg_meta_clone$cluster)))[1]
    df$cluster[which(df$TRB_CDR3 == clone)] = max_cluster
  }
}

# Now make the box plots:

grouped_treg = df %>% 
  filter(patient != 'P105') %>% # remove the patient we didn't have CD137 for
  group_by(patient, outcome, timept, cluster) %>% 
  summarize(MeanTreg=mean(MODULE_1),
            MeanIL10 = mean(IL10),
            MeanTGFB1 = mean(TGFB1),
            MeanFOXP3 = mean(FOXP3))

re_order = c('BL','BU','MN','AV')

ggplot(grouped_treg, aes(x=factor(timept, level=re_order), y=MeanIL10)) + 
  geom_boxplot(lwd = 0.1, outlier.size=0) + 
  stat_compare_means(comparisons = list(c(1,2), c(1,3), c(1,4)),  
                     method = 'wilcox', paired = FALSE, size = 2.1, symnum.args = symnum.args)+ 
  geom_beeswarm(aes(color=outcome), size=0.7, cex=3) +
  scale_color_manual(values=patientpalette[c(2,5,8)]) +
  facet_grid(cols = vars(cluster)) + ggpreset
ggsave("results/Figure6_Treg_clusters_IL10.pdf", width=3.6, height=1.8, unit = "in", family = "ArialMT", useDingbats =FALSE)

ggplot(grouped_treg, aes(x=factor(timept, level=re_order), y=MeanFOXP3)) + 
  geom_boxplot(lwd = 0.1, outlier.size=0) + 
  stat_compare_means(comparisons = list(c(1,2), c(1,3), c(1,4)),  
                     method = 'wilcox', paired = FALSE, size = 2.1, symnum.args = symnum.args)+ 
  geom_beeswarm(aes(color=outcome), size=0.7, cex=3) +
  scale_color_manual(values=patientpalette[c(2,5,8)]) +
  facet_grid(cols = vars(cluster)) + ggpreset
ggsave("results/Figure6_Treg_clusters_FOXP3.pdf", width=3.6, height=1.8, unit = "in", family = "ArialMT", useDingbats =FALSE)

ggplot(grouped_treg, aes(x=factor(timept, level=re_order), y=MeanTreg)) + 
  geom_boxplot(lwd = 0.1, outlier.size=0) + 
  stat_compare_means(comparisons = list(c(1,2), c(1,3), c(1,4)),  
                     method = 'wilcox', paired = FALSE, size = 2.1, symnum.args = symnum.args)+ 
  geom_beeswarm(aes(color=outcome), size=0.7, cex=3) +
  scale_color_manual(values=patientpalette[c(2,5,8)]) +
  facet_grid(cols = vars(cluster)) + ggpreset
ggsave("results/Figure6_Treg_clusters_Module1.pdf", width=3.6, height=1.8, unit = "in", family = "ArialMT", useDingbats =FALSE)



