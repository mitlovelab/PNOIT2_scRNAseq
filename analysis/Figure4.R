# Figure 4: T-helper subtypes are clonally distinct and exhibit TCR convergence.

# libraries ---------------------------------------------------------------

source('analysis/presets.R')


# Functions ---------------------------------------------------------------
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
sort_hclust_avg <- function(...) as.hclust(dendsort(as.dendrogram(...),type = "average"))


assignPheno = function(trb, seurat) {
  matches = seurat@meta.data[seurat@meta.data$TRB_CDR3 == trb,'pheno']
  matches = matches[!is.na(matches)]
  phenos = table(matches)
  a = names(phenos)[which.max(phenos)]
  a
}
# give the clone heatmap a shot
assignPatient = function(seurat, tcr) {
  pheno = table(seurat@meta.data[seurat@meta.data$TRB_CDR3 == tcr, 'patient'])
  pheno = pheno[order(pheno, decreasing = TRUE)]
  a= names(pheno)[1]
  a
}
# Setting up data objects -------------------------------------------------

th2 <- readRDS('th2mostcurrent.RDS')
th1_sub <-  readRDS('th1_sub.RDS')
th1 <-  readRDS('th1.RDS')
th17 <-  readRDS('th17.RDS')
all <-  readRDS('tcr_scores_genes_gates.RDS')
andy <-  readRDS('tcr_scores_genes_v2.RDS')


# produce TH1 and Th17 umap objects ------------------------------------------------

th1_names <-  th1_sub@cell.names[!th1_sub@cell.names %in% th2@cell.names]
cols <- intersect(colnames(th1@meta.data), colnames(th2@meta.data))
cols <-  cols[-3]
th1_add <-  SubsetData(th1_sub, th1_names)
cols <-  cols[cols %in% colnames(th1_add@meta.data)]
# make sure the meta data are the same before merging
th1_add@meta.data = th1_add@meta.data[,cols]
th2@meta.data = th2@meta.data[,cols]
th2@meta.data$group <- as.character(th2@meta.data$group)
seurat <-  MergeSeurat(th2, th1_add)
# same with th17 as well
th17_names <-  th17@cell.names[!th17@cell.names %in% seurat@cell.names]
th17_add <-  SubsetData(th17, th17_names)
seurat <-  MergeSeurat(seurat, th17_add) 

# Next, call the cluster of each cell directly

seurat@meta.data$cluster = NA
modmat = seurat@meta.data[,c("MODULE_7", "MODULE_9", 'MODULE_38')]
seurat@meta.data$cluster = apply(modmat, 1, which.max)
table(seurat@meta.data$cluster)

# Produce the umaps for each Th function

th1_sub = SetAllIdent(th1_sub, 'cluster')
th1_sub@meta.data$UMAP1 = th1_sub@dr$umap@cell.embeddings[,1]
th1_sub@meta.data$UMAP2 = th1_sub@dr$umap@cell.embeddings[,2]
th1_umap <- ggplot(th1_sub@meta.data, aes(x = UMAP1, y = UMAP2, color = cluster)) + 
  geom_point(stroke = 0, pch = 16, size = 0.5) +
  ggpreset + scale_color_brewer(palette = 'Set2') + No_axis_labels#+ guides(color = FALSE)
ggsave('results/Fig4_Th1_umap.pdf',th1_umap + TSNE_theme, height = 2, width = 2, units = 'in')


th17@meta.data$UMAP1 = th17@dr$umap@cell.embeddings[,1]
th17@meta.data$UMAP2 = th17@dr$umap@cell.embeddings[,2]
th17@meta.data$cluster = 'Th17'
df = th17@meta.data[,c('UMAP1', 'UMAP2', 'cluster')]
th17_umap <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = factor(cluster))) + 
  geom_point(stroke = 0, pch = 16, size = 0.5) +
  ggpreset + scale_color_manual(values = brewer.pal(3, 'Set2')[3]) + No_axis_labels
rm(df)
ggsave('results/Fig4_th17_umap.pdf', th17_umap + TSNE_theme, height = 2, width = 2, units = 'in')


th2@meta.data$UMAP1 = th2@dr$umap@cell.embeddings[,1]
th2@meta.data$UMAP2 = th2@dr$umap@cell.embeddings[,2]
th2_umap <- ggplot(th2@meta.data, aes(x = UMAP1, y = UMAP2, color = cluster)) + 
  geom_point(stroke = 0, pch = 16,size = 0.5) + TSNE_theme + 
  scale_color_manual(values = brewer.pal(6, 'Set2')[c(5,4,6)])
ggsave('results/Fig4_th2_umap.pdf', th2_umap, width = 2, height = 2, units = 'in', family = 'ArialMT')

ggplot(th2@meta.data, aes(x = UMAP1, y = UMAP2, color = cluster)) + 
  geom_point(stroke = 0, pch = 16,size = 1) + 
  scale_color_manual(values = brewer.pal(6, 'Set2')[c(5,4,6)])
ggsave('results/Fig4_th_legend.pdf', width = 4, height = 4, units = 'in', family = 'ArialMT')


# Fig 4B: IgE titers vs Th2 subset gene expression --------------------------------

# TODO: move this here from supplement code


# Fig 4C: TRB sharing between Th subtypes -----------------------------------------

seurat@meta.data$MODULE_7_GATED = all@meta.data[seurat@cell.names, 'MODULE_7_GATED']
seurat@meta.data$MODULE_9_GATED = all@meta.data[seurat@cell.names, 'MODULE_9_GATED']
seurat@meta.data$MODULE_38_GATED = all@meta.data[seurat@cell.names, 'MODULE_38_GATED']

# If cells score positive for multiple Th modules (Th1/2/17), assign each cell 
# exclusively to one by the max of its rank in the Th1, 2, and 17 module scores.

seurat@meta.data$cluster = NA
all@meta.data$th1_rank = rank(all@meta.data$MODULE_38)
all@meta.data$th2_rank = rank(all@meta.data$MODULE_7)
all@meta.data$th17_rank = rank(all@meta.data$MODULE_9)

modmat = all@meta.data[seurat@cell.names,c("th1_rank", 'th2_rank', 'th17_rank')]

seurat@meta.data$cluster = apply(modmat, 1, which.max)
seurat@meta.data$cluster = as.character(seurat@meta.data$cluster)
seurat@meta.data$cluster[seurat@meta.data$cluster == '1'] = 'Th1'
seurat@meta.data$cluster[seurat@meta.data$cluster == '2'] = 'Th2'
seurat@meta.data$cluster[seurat@meta.data$cluster == '3'] = 'Th17'

tab <- table(seurat@meta.data$TRB_CDR3, seurat@meta.data$cluster)
rsums = apply(tab, 1, sum)
tab = tab/rsums
tab = tab[rsums > 5,]# filtering for at least 4 cells.

#plots of the TRB sharing between th1/2/17 types.
hmp_cluster_row <- hclust(dist((tab)), method = "ward.D2")
hmp_cluster_row <- sort_hclust(hmp_cluster_row)

pheatmap(tab, color = viridis(50),cluster_rows = hmp_cluster_row, cluster_cols= TRUE,
         fontsize = 6, treeheight_row = 6, treeheight_col = 6,
         width = 2.5, height = 6, border_color = NA, filename = "results/Fig4C_TRB_Th_sharing_hmp_test.pdf") 
dev.off()

# move the cluster information to the pheno column.
seurat@meta.data$pheno = ''
seurat@meta.data$pheno[seurat@meta.data$cluster == 'Th2'] = th2@meta.data[seurat@cell.names[seurat@meta.data$cluster == 'Th2'], 'cluster']
# the issue is that the th2 meta data doesn't have cluster info, I think.
seurat@meta.data$pheno[seurat@meta.data$cluster == 'Th1'] = th1_sub@meta.data[seurat@cell.names[seurat@meta.data$cluster == 'Th1'], 'cluster']
seurat@meta.data$pheno[seurat@meta.data$cluster == 'Th17'] = 'Th17'


seurat@meta.data$th1 = as.numeric(seurat@cell.names %in% th1_sub@cell.names[th1_sub@meta.data$cluster == 'Th1'])
seurat@meta.data$th1fh = as.numeric(seurat@cell.names %in% th1_sub@cell.names[th1_sub@meta.data$cluster == 'Tfh1-like'])
seurat@meta.data$th2 = as.numeric(seurat@cell.names %in% th2@cell.names[th2@meta.data$cluster == 'GATA3hi'])
seurat@meta.data$th2fh = as.numeric(seurat@cell.names %in% th2@cell.names[th2@meta.data$cluster == 'Tfh-like'])
seurat@meta.data$th2reg = as.numeric(seurat@cell.names %in% th2@cell.names[th2@meta.data$cluster == 'Treg-like'])
seurat@meta.data$th17 = as.numeric(seurat@cell.names %in% th17@cell.names)
seurat@meta.data$outcome = all@meta.data[seurat@cell.names, 'outcome']
sub = seurat@meta.data[seurat@meta.data$outcome != 'placebo',c('TRB_CDR3','th1', 'th1fh', 'th2', 'th2fh', 'th2reg', 'th17')]
sub = sub %>% group_by(TRB_CDR3) %>% filter(!is.na(TRB_CDR3)) %>% summarize_all(.funs =sum) %>% as.data.frame()
trb_count = table(seurat@meta.data$TRB_CDR3)
head(sub)
tab = sub
rownames(tab) = tab[,1]
rsums = trb_count[tab$TRB_CDR3]
tab = tab[,-1]
tab = tab/rsums
tab = tab[rsums > 3,]
pheatmap(tab, color = viridis(100))


module_trees <- cutree(hmp_cluster_row, k = 5)
module_dgram_table <- data.frame(module = module_trees)
module_dgram_table <- module_dgram_table %>% 
  rownames_to_column(var = "TRB_CDR3") %>% 
  mutate(module = factor(module)) %>% 
  group_by(module) %>%
  sample_n(ifelse(n()>=20,20, n()))
sampled_TRB <- module_dgram_table$TRB_CDR3

tab_test <- tab[rownames(tab) %in% sampled_TRB,]

hmp_cluster_row <- hclust(dist((tab_test)), method = "complete")
hmp_cluster_row <- sort_hclust(hmp_cluster_row)

hmp_cluster_col <-hclust(dist(t(tab_test)), method = "complete")
hmp_cluster_col <- sort_hclust(hmp_cluster_col)

# organize final TRB_heatmap
pheatmap(tab_test, color = viridis(50),cluster_rows = hmp_cluster_row, cluster_cols= TRUE,
         fontsize = 5, treeheight_row = 6, treeheight_col = 6,
         width = 3, height = 5, border_color = NA)

# order the columns to match with earlier legend:
tab_test <- tab_test[,c(6,1,3,4,2,5)]
Annotation_col <- data.frame(Subset = c("Th17","Th1","Th2A-like","Tfh2-like","Tfh1-like","Treg-like"))
rownames(Annotation_col) <- colnames(tab_test)

module_dgram_temp <- module_dgram_table %>% filter(TRB_CDR3 %in% rownames(tab_test)) 
Annotation_row <- data.frame(module = module_dgram_temp$module)
rownames(Annotation_row) <- module_dgram_temp$TRB_CDR3

clusterpal <- c(brewer.pal(6, 'Set2'))
names(clusterpal) <-  c("Tfh1-like", 'Th1', 'Th17', 'Tfh2-like', 'Th2A-like', 'Treg-like')
my_colour <-  list(Subset = clusterpal)

pheatmap(tab_test, color = viridis(50),cluster_rows = as.hclust(rev(as.dendrogram(hmp_cluster_row))), 
         cluster_cols= FALSE,annotation_colors = my_colour,annotation_col = Annotation_col,
         fontsize = 5, treeheight_row = 8, show_colnames = FALSE, annotation_names_row = FALSE,
         width = 3, height = 5, border_color = NA, annotation_row = Annotation_row,
         filename = "results/TRB_Th_subsets_sharing_hmp_v3.pdf") #this is closer,

# Also label the rows, so we can see the order of the leaves in case we want to flip any:
plot(hmp_cluster_row)
pheatmap(tab_test, color = viridis(50),cluster_rows = hmp_cluster_row, cluster_cols= FALSE,
         fontsize = 5, treeheight_row = 6, treeheight_col = 6,
         width = 3, height = 5, border_color = NA, annotation_row = Annotation_label)

# Fig 4D: Gene expression heatmap -------------------------------------------------

# Genes pre-selected by ROC test in Seurat, and TH2/1/17/Tfh genes of interest:
genes = c('ICOS', 'CXCR5', "NFKB1", 'DUSP4', 'REL', 'NR4A3', 'TNF', 'IL2', 'GZMB', 'IFNG', "DPP4", 
          "IL22", 'CYTIP', 'CXCR6','IL17A', "IL17F", 'PTPN13', 'TOX', "ZEB2", "SOS1",
          "CHDH", 'PTGDR2', 'KLRB1', 'ITGA4', 'ANXA1', 'SOS1', 'TNFSF10', 'IL16', 'IL17RB', 'GATA3', 
          'TNFRSF9', 'CTLA4', 'IL5', 'IL9', 'IL13', 'NFKB1', 'ICOS', 'IL4', 'NR4A3', 'PDCD1', 'CXCR5', 
          'CD27', 'TIGIT', "FOXP3", "IKZF2", "IL1R1", 'IL1R2')
genes = unique(genes)
df = seurat@data[genes,] %>% as.matrix() %>% t() %>% as.data.frame()
df$cluster = seurat@meta.data$pheno
df$patient = seurat@meta.data$patient
df$outcome = all@meta.data[rownames(df), 'outcome']
df$cluster = factor(as.character(df$cluster), levels = c('Th17', 'Th1', 'GATA3hi', 'Tfh-like', 'Tfh1-like', 'Treg-like'))
df = df %>% filter(!is.na(cluster)) %>% group_by(cluster, outcome, patient) %>% summarize_all(~mean(.)) %>% as.data.frame() 
rownames(df) = paste(df$cluster, df$patient)
meta = data.frame(Cluster = df$cluster, row.names = rownames(df))
df = df[,-c(1:3)]
df = scale(df)
df[df > 2.5] = 2.5
df[df < -1] = -1
clusterpal = brewer.pal(6, 'Set2')
names(clusterpal) = c("Tfh1-like", 'Th1', 'Th17', 'Tfh-like', 'GATA3hi', 'Treg-like')
outcomepal = patientpalette[c(1,4,7,10)]
names(patientpalette) = c("P105", "P106", "P111", "P33", 'P90', 'P93', 
                          "P69", "P95", 'P97', 'P84', 'P96', 'P107')
#names(outcomepal) = c('Full Tolerance', 'Partial Tolerance', "Treatment Failure", "Placebo")

hmp_cluster_row <- hclust(dist(t(df)), method = "ward.D2")
hmp_cluster_row <- sort_hclust_avg(hmp_cluster_row)

pheatmap(t(df), color = inferno(100), cluster_rows =hmp_cluster_row, annotation_col = meta, 
         annotation_colors= list(Cluster = clusterpal, Patient = patientpalette), cluster_cols= FALSE, 
         show_colnames = FALSE, fontsize = 6, treeheight_row = 6,
         width = 3, height = 5, border_color = NA, filename = "results/Fig4D_Th_subsets_gene_hmp_v3.pdf")


# Fig 4E: TCRdist -----------------------------------------------------------------

# Output of TCRdist python script:
distmat = read.csv('th17atrbout.txt', stringsAsFactors = FALSE, row.names = 1) %>% as.data.frame()

# TODO: comment this plotting section better
df = seurat@meta.data %>% filter(!is.na(pheno)) %>% filter(!is.na(TRB_CDR3))
distmat = distmat[rownames(distmat) %in% df$TRB_CDR3,rownames(distmat) %in% df$TRB_CDR3]
a = sapply(rownames(distmat),function(x) assignPheno(x, seurat))
distmat$pheno1 = a
distmat = distmat %>% as.data.frame()
hash = data.frame(row.names = rownames(distmat), Pheno = distmat$pheno1)

distmat$TCR = rownames(distmat)
df = melt(distmat, value.name = 'dist', id = c('pheno1', 'TCR'))
df$variablepheno = hash[df$variable, 'Pheno']
df$bin = cut(df$dist, c(0,4,8,12,1000))
df$same = as.numeric(df$pheno1 == df$variablepheno) 
hash2 = data.frame(row.names = rownames(distmat), patient = sapply(rownames(distmat), function(x) assignPatient(seurat, x)))
df$pat1 = hash2[df$TCR, 'patient']
df$pat2 = hash2[df$variable, 'patient']

result = df %>% filter(pat1 == pat2) %>% group_by(pheno1, bin, same) %>%
  summarize(count = n()) %>% mutate(total = sum(count), p = count/total) %>%
  filter(same == 1) %>% mutate(sd = sqrt(p*(1-p)/total))
result$ratio = NA
for (a in unique(result$pheno1)) {
  temp = result[result$pheno1 == a & is.na(result$bin),]
  result$ratio[result$pheno1 == a] =temp$count / temp$total
}
result = result %>% mutate(prat = p/ratio, sd = sd/ratio)
lv = levels(result$bin)

result2 = result
result = result[!is.na(result$bin),]

blank = result[1,]
blank = rbind(result, result, result, result, result)
blank$pheno1 = c('Tfh-like', 'Tfh-like','Treg-like', 'Treg-like', 'Treg-like')
blank$prat = 0
#TODO: annotate where the 1.44 comes from
result2 = result %>% group_by(pheno1) %>% mutate(sd_rel = sd/prat*1.44, p_back = tail(p, n =1),prat = p/tail(p, n= 1), alternative = 'greater')
result2$p_val = NA
for (i in 1:dim(result2)[1]) {
  result2$p_val[i] = prop.test(result2$count[i], result2$total[i], p = result2$p_back[i], alternative = 'greater')$p.value
}
result2[result2$p_val < 0.01,]
ggplot(result2, aes(x = bin, y = log2(prat), fill = pheno1)) + geom_col(position = 'dodge') + scale_fill_brewer(palette = 'Set2') + 
  geom_errorbar(aes(ymin = log2(prat)-sd_rel, ymax =log2(prat) + sd_rel), position = 'dodge')


blank = result[1,]
blank = rbind(blank, blank, blank, blank)
blank$pheno1 = c('Tfh1-like', 'Tfh-like', 'Treg-like', 'Treg-like')
blank$prat = 1
blank$sd = 0
blank$bin = levels(result$bin)[c(result2$bin[1], result2$bin[3], result2$bin[1:2])]
blank
blank$pheno1 = as.character(blank$pheno1)
result2$pheno1 = as.character(result2$pheno1)
temp = rbind(result2,blank)


ggplot(temp, aes(x = bin, y = log2(prat), fill = pheno1)) + geom_col(position = 'dodge') + scale_fill_brewer(palette = 'Set2') + 
  geom_errorbar(aes(ymin = log2(prat)-sd_rel, ymax =log2(prat) + sd_rel), position = 'dodge', size = .5) + scale_x_discrete(limits = levels(result$bin)[c(1:2, 3:4,5,6)])

temp <- temp %>% mutate(pval_symbol = case_when(p_val < 5e-2 & p_val > 5e-3 ~ "*",
                                        p_val < 5e-3 & p_val > 5e-4 ~ "**",
                                        p_val < 5e-4 & p_val > 5e-5 ~ "***",
                                        p_val < 5e-5 ~ "****"))

temp <- temp %>% 
  ungroup() %>%
  mutate(bin = factor(bin, levels = c('(0,4]','(4,8]','(8,12]','(12,1e+03]'), labels = c('1-4','5-8','9-12','>13'))) %>%
  mutate(pheno1 = factor(pheno1)) %>% 
  mutate(pheno1 = fct_recode(pheno1, Th2Alike = 'GATA3hi', Tfh2like = 'Tfh-like'))

th_subpalette <- brewer.pal(6, 'Set2')
names(th_subpalette) <- c('Tfh1-like','Th1-conv','Th17','Tfh2-like','Th2A-like','Treg-like')

temp <- temp %>% 
  mutate(pheno1 = fct_recode(pheno1, `Th2A-like` = "Th2Alike",
                                     `Th1-conv` = "Th1",
                                     `Tfh2-like` = "Tfh2like")) %>% 
  mutate(pheno1 = fct_relevel(pheno1, 'Tfh1-like','Tfh2-like','Th1-conv','Th17','Th2A-like','Treg-like')) %>% 
  filter(!is.na(sd_rel))

ggplot(temp, aes(x = bin, y = log2(prat), color = pheno1)) + 
  geom_point(size = 1, shape = 15) +
  scale_color_manual(values = th_subpalette) +
  geom_errorbar(aes(ymin = log2(prat)-sd_rel, ymax =log2(prat) + sd_rel), color = 'black',
                                                      position = 'dodge', size = .25) + 
  geom_text(aes(x = bin, y = log2(prat) + sd_rel, label = pval_symbol), 
            position = position_dodge(width = 1),vjust = -0.05, size = 2.1, color = 'black') + 
  geom_text(aes(x = bin, y = log2(prat) + sd_rel, label = total), 
            position = position_dodge(width = 1),vjust = -1.5, size = 2.1, color = 'black') +
  labs(x = "Distance betwen TCR pair (TCRdist)", y = "log2(likelihood of TCR pair in same subset)") +
  ggpreset + theme(legend.position = 'right',legend.key.size = unit(0.5,"line")) + facet_wrap(~ pheno1, nrow = 1) +
  theme(strip.background = element_rect(color=NA, fill='grey87', size=0.5),strip.text.x = element_text(size=6)) +
  No_legend

ggsave("results/Fig4E_tcrdist.pdf",width = 5.5, height = 1.7, unit = "in")

save.image("Figure4_v3.Rdata")
