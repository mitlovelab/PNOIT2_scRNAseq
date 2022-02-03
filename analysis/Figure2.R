# Figure 2: Peanut-reactive T cells have diverse and distinct transcriptional signatures.

# Import presets, libraries, functions, and processed data --------------

source("analysis/presets.R")
Total_meta_data <- readRDS("data/Total_meta_data.rds")
seurat <-  readRDS("data/tcr_scores_genes_gates.RDS")

# Figure 2A -------------

# get patients in the order that matches color palette
# TODO: move this to preprocessing section
Total_meta_data <- Total_meta_data %>% 
  mutate(patient = factor(patient, levels = c("P105","P106","P111","P33","P90","P93","P69","P95","P97","P84","P96","P107")))

plot_data = Total_meta_data %>% filter(downsampled) %>% sample_frac()

ggplot(plot_data, aes(x = UMAP1_ds, y = UMAP2_ds, color = patient)) +
  geom_point_rast(size=0.25, stroke=0, raster.dpi = 1000) + 
  xlab(label="UMAP 1") + ylab(label="UMAP 2") + 
  scale_color_manual(values = patientpalette) + 
  ggpreset + UMAPpreset 

ggsave("results/Fig2A_UMAP_patient.pdf", width=2.5, height=2.5, useDingbats=FALSE)


ggplot(plot_data, aes(x = UMAP1_ds,y = UMAP2_ds, color = orig.ident)) +
  geom_point_rast(size=0.25, stroke=0, raster.dpi = 1000) + 
  xlab(label="UMAP 1") + ylab(label="UMAP 2") +
  scale_color_manual(values = subsetpalette) + 
  ggpreset + UMAPpreset

ggsave("results/Fig2A_UMAP_subset.pdf", width=2.5, height=2.5, useDingbats=FALSE)


# Figure 2B ------------

# Plot differentially expressed genes between CD154+, CD137+, and double-negative subsets

# First examine the top genes using ROC test in seurat
# Randomly subset the data to speed things up
# TODO: replace outdated seurat functions
#seurat <- SetAllIdent(object=seurat, id="condition")
#seurat2 <- SubsetData(object=seurat, max.cells.per.ident = 2000)
#FindMarkers(seurat2, ident.1 = "CD137", ident.2 = c("CD154","DblNeg"), test.use = "roc", max.cells.per.ident=2000)
#FindMarkers(seurat2, ident.1 = "CD154", ident.2 = c("CD137","DblNeg"), test.use = "roc", max.cells.per.ident=2000)
#FindMarkers(seurat2, ident.1 = "DblNeg", ident.2 = c("CD154","CD137"), test.use = "roc", max.cells.per.ident=2000)

# Final set of genes for plotting:
# TODO: not all genes are present in the truncated data, fix
genes <- c(
  "TNFRSF9","HLA-A","IKZF2","FOXP3","TIGIT","IL2RA","TNFRSF1B","BATF", # enriched in CD137
  "CD40LG","BHLHE40","LCP1","VIM","SOS1","ACTG1","ACTB","PKM", # enriched in CD154
  "TPT1","IL7R","RPL31","RPL30","TCF7","PABPC1","TMEM66" # enriched in DblNeg
)

# Get mean expression of each gene in each sorted subset/patient
genes_summary <- as.data.frame(t(as.matrix(seurat@data[genes,])))
genes_summary$patient = seurat@meta.data$patient
genes_summary$condition = seurat@meta.data$condition

genes_summary <- genes_summary %>%
  group_by(condition, patient) %>%
  summarise_all(mean)

# Sort to keep the genes in the desired order, and scale the columns
genes_plot <- genes_summary %>%
  ungroup() %>%
  select(genes) %>%
  mutate_each(funs(scale))

pheatmap(t(as.matrix(genes_plot)), color = viridis(100, option='inferno'), 
         cluster_rows = FALSE, cluster_cols= FALSE, 
         fontsize = 6, width = 4, height = 4.5, border_color = NA, 
         filename = "results/Fig2B_DEgenes.pdf")


# Figure 2C ------------

# 6 selected modules:
mods <- c(7, 1, 38, 9, 2, 4)
mods_list <- paste("MODULE",mods,sep = "_")
# And corresponding descriptions:
Interpret <- c("TH2", "Treg", "TH1", "TH17", "Costimulatory markers", "Interferon response")

# Set bounds on outlier values, for plotting purposes
Set_maxmin <- function(meta_data, module){
  upperbound <- mean(meta_data[[module]]) + 4*sd(meta_data[[module]])
  lowerbound <- mean(meta_data[[module]]) - 4*sd(meta_data[[module]])
  module <- sym(module)
  meta_data <- meta_data %>% mutate((!!module) := ifelse((!!module)>upperbound,upperbound,(!!module))) %>%
    mutate((!!module) := ifelse((!!module)<lowerbound,lowerbound,(!!module)))
  return(meta_data)
}


Module_umap <- map(mods_list, function(mods_list) ggplot(data=Set_maxmin(Total_meta_data %>% 
                                                    filter(downsampled) %>% sample_frac(), mods_list),
                                aes(x = UMAP1_ds,y = UMAP2_ds, color = !!sym(mods_list))) +
                                geom_point_rast(size=0.8, stroke=0) + 
                                scale_color_gradientn(colors=modulepalette) + 
                                scale_x_continuous(expand=c(0,0)) + 
                                scale_y_continuous(expand=c(0,0)) + 
                                ggpreset + UMAPpreset + 
                                theme(axis.title = element_blank(), panel.border = element_blank(),
                                plot.margin = unit(c(0.1,0.4,0.1,0.1), "cm")))

names(Module_umap) <- mods_list

plot_grid(plotlist = Module_umap, nrow=4, rel_widths = rep(c(2.6),2))

ggsave("results/Fig2C_Selected_modules.pdf", width=2.8, height=4.9) 


#TODO: Pull in sparse PCA bar plot code from supplement files

