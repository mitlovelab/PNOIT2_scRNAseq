# Figure 3: T-helper gene modules are associated with clonal expansion and expression in activated cells.


# Source functions and libraries---------------------------------------------------------------

source("analysis/presets.R")


# import dataset ----------------------------------------------------------

# 300k cells, more inclusive
Total_meta_data <- read_rds("data/Total_meta_data.rds")

# 100k cells, higher-quality subset
HQ_meta_data <- read_rds("data/Processed_PNOIT2_cellmetadata.rds")

# Functions ---------------------------------------------------------------

Add_TRAB_count <- function(meta.data, group){
  require(tidyverse)
  require(viridis)
  group <- enquo(group)
  meta.data %>% group_by(!!group) %>% dplyr::count(TRB_CDR3)%>%
    filter(!is.na(TRB_CDR3)) %>% arrange(desc(n)) -> seurat_CDR3B_count
  colnames(seurat_CDR3B_count) <- c(quo_name(group),"TRB_CDR3","clones")
  meta.data %>% group_by(!!group) %>% dplyr::count(TRA_CDR3) %>% 
    filter(!is.na(TRA_CDR3)) %>% arrange(desc(n)) -> seurat_CDR3A_count
  colnames(seurat_CDR3A_count) <- c(quo_name(group),"TRA_CDR3","clones")
  
  left_join(meta.data, seurat_CDR3B_count, by = c("TRB_CDR3",quo_name(group))) %>% .$clones -> meta.data$TRB_count_group
  
  left_join(meta.data, seurat_CDR3A_count, by = c("TRA_CDR3",quo_name(group))) %>% .$clones -> meta.data$TRA_count_group
  
  return(meta.data)
}


shift_fn <- function(x){
  if(min(x)>=0){
    return(x)
  } else if(min(x)<0){
    return(x+abs(min(x)))
  }
}

# Figure 3A ---------------------------------------------------------------

# Join the 300k and 100k datasets
#TODO: clean this up to be more concise
HQ_meta_data <- HQ_meta_data %>% rename(umap1_100k = umap1, umap2_100k = umap2)

Total_meta_data <- Total_meta_data %>% 
  left_join(HQ_meta_data %>% 
              select(Barcode, patient, condition, umap1_100k, umap2_100k), 
              by = c("Barcode","patient","condition"))

# Calculate the TCR counts:

Total_CDR3B_count <- Total_meta_data %>% 
  count(TRB_CDR3) %>% 
  filter(!is.na(TRB_CDR3)) %>% 
  arrange(desc(n))
colnames(Total_CDR3B_count) <- c("TRB_CDR3","clones")

Total_CDR3A_count <- Total_meta_data %>% 
  count(TRA_CDR3) %>% 
  filter(!is.na(TRA_CDR3)) %>% 
  arrange(desc(n))
colnames(Total_CDR3A_count) <- c("TRA_CDR3","clones")

left_join(Total_meta_data, Total_CDR3B_count, by = "TRB_CDR3") %>% .$clones -> Total_meta_data$Total_TRB_count
left_join(Total_meta_data, Total_CDR3A_count, by = "TRA_CDR3") %>% .$clones -> Total_meta_data$Total_TRA_count

# Plot total TRA and TRB count. Shuffle the order of data points first.
Total_meta_data <- Total_meta_data %>% sample_frac()

Total_meta_temp <- Total_meta_data %>% 
  select(umap1_100k,umap2_100k,Total_TRB_count,Total_TRA_count) %>% 
  filter(!is.na(umap1_100k))

TRB_breaks <- log2(2^(c(0,seq(7))))
TRB_breaks_label <- 2^(c(0,seq(7)))

# Filter this by NAs, then bind the rows back
Total_meta_b <- Total_meta_temp %>% filter(is.na(Total_TRB_count))
Total_meta_b_1 <- Total_meta_temp %>% filter(!is.na(Total_TRB_count))
Total_meta_b <- bind_rows(Total_meta_b,Total_meta_b_1)

p <- ggplot(Total_meta_b, aes(umap1_100k, umap2_100k, color = log2(Total_TRB_count))) +
  geom_point_rast(size=0.25, stroke=0, raster.dpi = 600) + 
  labs(x="UMAP 1", y ="UMAP 2") + 
  scale_color_gradientn(colors=rev(viridis(20)), na.value="gray90",name = "Clonal size \n(TCRB)",
                        breaks = TRB_breaks,labels = TRB_breaks_label, limits = c(0,log2(128))) + 
  TSNE_theme + 
  theme(legend.position = c(0.2,0.2), legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size=6), legend.background = element_blank())

ggsave("results/Fig3A_UMAP_TRBclonal_size.pdf", p, width=2.3, height=2.3, family = 'ArialMT')

# Now do the same for TRA

TRA_breaks <- log2(2^(c(0,seq(8))))
TRA_breaks_label <- 2^(c(0,seq(8)))

# Filter this by NAs, then bind the rows back

Total_meta_a <- Total_meta_temp %>% filter(is.na(Total_TRA_count))
Total_meta_a_1 <- Total_meta_temp %>% filter(!is.na(Total_TRA_count))
Total_meta_a <- bind_rows(Total_meta_a,Total_meta_a_1)

p_a <- ggplot(Total_meta_a, aes(umap1_100k, umap2_100k, color = log2(Total_TRA_count))) +
  geom_point_rast(size=0.25, stroke=0, raster.dpi = 600) + 
  labs(x="UMAP 1", y ="UMAP 2") + 
  scale_color_gradientn(colors=rev(viridis(20, option="plasma")), na.value="gray90",name = "Clonal size \n(TCRA)",
                        breaks = TRA_breaks,labels = TRA_breaks_label, limits = c(0,log2(256))) + 
  TSNE_theme + 
  theme(legend.position = c(0.2,0.2), legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size=6), legend.background = element_blank())

ggsave("results/Fig3A_UMAP_TRAclonal_size.pdf", p_a, width=2.3, height=2.3, family = 'ArialMT')


# Figure 3B ---------------------------------------------------------------

df_norm_tp <- Total_meta_data %>% group_by(outcome, patient, timept,condition, TRB_CDR3) %>%
  filter(!is.na(TRB_CDR3)) %>%
  summarize(n = n()) %>%
  summarize(div = diversity(n, index = 'shannon', base = n()))

#double check on the recovery
recovery_trb_summary <- Total_meta_data  %>% 
  group_by(outcome, patient, timept, condition) %>% 
  summarise(cell_n = n(), trb = sum(!is.na(TRB_CDR3)), recovery = trb/cell_n)
#there is a double negative sample that is an outlier. P111 DblNeg, at baseline.


ggplot(df_norm_tp, aes(x = condition, y = div, color = condition)) + geom_beeswarm(size = 0.5,cex=2, priority = 'random')  + 
  stat_compare_means(comparisons = list(c(1,3),c(2,3)), size = 2.5, method = "wilcox.test") + 
  scale_color_manual(values = subsetpalette[c(2,6,10)]) + 
  labs(x = '', y = 'Normalized Shannon index \nof TCRb repertoire') +
  No_legend + Axis_themes + ggpreset#+ coord_cartesian(ylim = c(0.95, 1.03))

ggsave('results/Fig3B_TRB_entropy_norm_all_tp_v4.pdf',width = 2.2, height = 2.8, family = 'ArialMT', useDingbats =FALSE)

# Figure 3C ---------------------------------------------------------------

# now make the plot of clonal sizes. downsample first so it's not skewed by condition:
cell_minimum <- Total_meta_data %>% filter(!is.na(TRB_CDR3)) %>% count(condition) %>% .$n %>% min()
Total_meta_data_ds <- Total_meta_data %>% filter(!is.na(TRB_CDR3)) %>% group_by(condition) %>% sample_n(cell_minimum)


TRB_count_by_stim <- Total_meta_data_ds %>% filter(!is.na(TRB_CDR3)) %>% group_by(condition) %>% count(TRB_CDR3)
TRB_count_by_stim <- TRB_count_by_stim %>% rename(count = n)
TRB_count_by_stim <- TRB_count_by_stim %>% mutate(count = as.numeric(count)) %>% ungroup()

TRB_breaks <- log2(2^(c(0,seq(6))))
TRB_breaks_label <- 2^(c(0,seq(6)))

ggplot(TRB_count_by_stim, aes(x = condition, y = log2(count))) + 
  geom_quasirandom(varwidth = FALSE, size = 0.25, method ="tukey") +
  stat_compare_means(comparisons = list(c(1,3),c(2,3)),size = 2.5, method = "wilcox.test") + 
  Axis_themes + labs(x = '', y = 'TRb clonal size') + 
  scale_y_continuous(breaks = TRB_breaks, labels = TRB_breaks_label) 

ggsave("results/Fig3C_clonal_sizes_entropy.pdf",width = 2.8, height = 3, units = "in", family = 'ArialMT')



# Figure 3D ---------------------------------------------------------------

# Level of sharing of TCRb clonotypes between conditions
# "condition": a sorted subset (CD137+, CD154+, CD137-CD154-) at a timepoint

# Downsample to the patients for which we have TCR data from all 3 sorted subsets
TCR_data <- Total_meta_data %>%
  # TODO: update with latest list of patients with TCR
  filter(patient %in% c("P106","P111","P33","P95")) %>%
  filter(!is.na(TRB_CDR3))

conditions <- levels(droplevels(TCR_data$orig.ident))

# Get % shared clonotypes between each pair of conditions
pairwiseTCR <- matrix(0, nrow=length(conditions), ncol=length(conditions))
rownames(pairwiseTCR) <- conditions
colnames(pairwiseTCR) <- conditions
for (i in 1:nrow(pairwiseTCR)){
  clones1 <- TCR_data[which(TCR_data$orig.ident == rownames(pairwiseTCR)[i]),]
  clones1 <- as.character(unique(clones1$TRB_CDR3))
  for (j in 1:nrow(pairwiseTCR)){
    clones2 <- TCR_data[which(TCR_data$orig.ident == rownames(pairwiseTCR)[j]),]
    clones2 <- as.character(unique(clones2$TRB_CDR3))
    # Normalize by geometric mean of number of clonotypes
    geom_mean <- sqrt(length(clones1)*length(clones2))
    pairwiseTCR[i,j] <- 100*sum(table(c(clones1, clones2)) > 1) / geom_mean
  }
  pairwiseTCR[i,i] <- NA
}

pheatmap(pairwiseTCR, color = viridis(100), 
         cluster_rows = FALSE, cluster_cols= FALSE, 
         fontsize = 6, width = 4.2, height = 4, border_color = NA, 
         filename = "results/Fig3D_TCRsharing.pdf")



# Figure 3E ---------------------------------------------------------------


# Names of the module-gate yes/no columns
Module_gates <- Total_meta_data %>% select(ends_with("GATED")) %>% names()
Module_names <- Total_meta_data %>% select(starts_with("MODULE")) %>% select(-ends_with("GATED")) %>% names()

Total_meta_data <- Add_TRAB_count(Total_meta_data,condition) # this adds count by condition
#previously count was added by total (i.e. count in whole dataset, not in each sorted condition)
Total_meta_data <- Total_meta_data %>% mutate_at(c(Module_names),~shift_fn(.))

Total_filtered <- Total_meta_data %>% filter(condition != "DblNeg") %>%
  filter(!is.na(TRB_CDR3)) %>% filter(HighQ == TRUE) # remove lower quality cells

Total_filtered_all_cond <- Total_meta_data %>% filter(!is.na(TRB_CDR3)) %>% filter(HighQ == TRUE) #use high quality cells

# OK now, we need to get the average clonal size for all module-scoring cells
Module_gates <- Module_gates[1:50]

expansion <- data.frame(matrix(0, nrow=50, ncol=1))
colnames(expansion) <- c("Module")
rownames(expansion) <- Module_gates

# Calculate clonal size of module-expressing cells, and enrichment in module expression above double-negative cells, for every gene module
expansion$Expansion_154 <- 0
for (j in seq_along(Module_gates)){
  Cellsmod <- Total_filtered %>% filter(!!sym(Module_gates[j]) == TRUE) %>% filter(condition == "CD154")
  expansion$Expansion_154[j] <- mean(Cellsmod$TRB_count_group)
}

expansion$Expansion_137 <- 0
for (j in seq_along(Module_gates)){
  Cellsmod <- Total_filtered %>% filter(!!sym(Module_gates[j]) == TRUE) %>% filter(condition == "CD137")
  expansion$Expansion_137[j] <- mean(Cellsmod$TRB_count_group)
}


expansion$Enrichment_154 <- 0
for (j in seq_along(Module_gates)){
  Cellsmod <- Total_filtered_all_cond  %>% 
    select(!!sym(Module_names[j]), !!sym(Module_gates[j]), condition)
  pos <- Cellsmod %>% filter(condition %in% c("CD154")) %>% filter(!!sym(Module_gates[j])) %>% 
    pull(!!sym(Module_names[j]))
  neg <- Cellsmod %>% filter(condition %in% c("DblNeg")) %>% #filter(!!sym(Module_gates[j])) %>% 
    pull(!!sym(Module_names[j]))
  expansion$Enrichment_154[j] <- (mean(pos) / mean(neg))
}
expansion$Enrichment_137 <- 0
for (j in seq_along(Module_gates)){
  Cellsmod <- Total_filtered_all_cond  %>% 
    select(!!sym(Module_names[j]), !!sym(Module_gates[j]), condition)
  pos <- Cellsmod %>% filter(condition %in% c("CD137")) %>% filter(!!sym(Module_gates[j])) %>% 
    pull(!!sym(Module_names[j]))
  neg <- Cellsmod %>% filter(condition %in% c("DblNeg")) %>% #filter(!!sym(Module_gates[j])) %>% 
    pull(!!sym(Module_names[j]))
  expansion$Enrichment_137[j] <- (mean(pos) / mean(neg))
}

ggplot(expansion, aes(x=Expansion_154, y=Enrichment_154)) + 
  geom_point(size=1.5, stroke = 0, pch = 16) + 
  geom_text_repel(data = expansion %>% top_n(n = 6, wt = Enrichment_154), aes(label = Module), size = 2,force = 0.1) + 
  geom_text_repel(data = expansion %>% top_n(n = 6, wt = Expansion_154), aes(label = Module), size = 2,force = 0.1) +
  labs(y = "Fold Enrichment (CD154)", x = "Average clonal size (CD154)") + ggpreset
ggsave("results/Fig3E_expansion_vs_enrichment_CD154_labeled_V3.pdf", width=2.6, height=2.6)

ggplot(expansion, aes(x=Expansion_137, y=Enrichment_137)) + 
  geom_point(size=1.5, stroke = 0, pch = 16) + 
  geom_text_repel(data = expansion %>% top_n(n = 6, wt = Enrichment_137), aes(label = Module), size = 2,force = 0.1) + 
  geom_text_repel(data = expansion %>% top_n(n = 6, wt = Expansion_137), aes(label = Module), size = 2,force = 0.1) +
  labs(y = "Fold Enrichment (CD137)", x = "Average clonal size (CD137)") + ggpreset
ggsave("results/Fig3_expansion_vs_enrichment_CD137_labeled_V3.pdf", width=2.6, height=2.6)

ggplot(expansion, aes(x=log10(Enrichment_154), y=log10(Enrichment_137))) + 
  geom_point(size=1.5, stroke = 0, pch = 16) + 
  geom_text_repel(aes(label = Module), size = 2,force = 2, segment.size = 0.2) +
  ggpreset + labs(x = "Fold Enrichment (CD154)", y = "Fold Enrichment (CD137)")
ggsave("results/Fig3_enrichment_vs_enrichment_CD137v154_labeled_v3.pdf", width=2.6, height=2.6, useDingbats=FALSE)
