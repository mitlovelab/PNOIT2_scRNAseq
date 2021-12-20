# Figure 5: Th1 and Th2 effector, but not Tfh-like, subsets are suppressed by OIT.

# libraries ---------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(factoextra)
library(viridis)
library(pheatmap)
library(RColorBrewer)
library(ggbeeswarm)

# Figure configuration ----------------------------------------------------

source("presets.R")

# data setup --------------------------------------------------------------

th1_all = readRDS('th1_all.RDS')
th17_all = readRDS('th17_all.RDS')
th2_all = readRDS('th2_all.RDS')
all_clones <-  readRDS('all_clones_orig.RDS')

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 3)/3, symbols = c("****", "***", "**", "*", "ns"))


# Figure 5A ---------------------------------------------------------------

df = th2_all %>% 
  filter(outcome != 'placebo') %>% 
  group_by(outcome, patient, timept) %>%
  summarize(mean = mean(MODULE_7), count = n(), sd = sd(MODULE_7)/sqrt(sum(count))) %>% as.data.frame()
df$th = 'Th2'
master = df

df = th1_all %>% 
  filter(outcome != 'placebo') %>% 
  group_by(outcome, patient, timept) %>%
  summarize(mean = mean(MODULE_38), count = n(), sd = sd(MODULE_38)/sqrt(sum(count))) %>% as.data.frame()
df$th = 'Th1'
master = rbind(master, df)

df = th17_all %>% 
  filter(outcome != 'placebo') %>% 
  group_by(outcome, patient, timept) %>%
  summarize(mean = mean(MODULE_9), count = n(), sd = sd(MODULE_9)/sqrt(sum(count))) %>% as.data.frame()
df$th = 'Th17'
master = rbind(master, df)
master$th = factor(master$th, levels=c('Th2','Th1','Th17'))

ggplot(data = master, aes(x = timept, y = mean)) + 
  geom_boxplot(lwd = 0.1, outlier.shape = NA) + 
  geom_beeswarm(aes(color = outcome), stroke=0, pch=16, size=1.7, cex=6) +
  labs(x = "Timepoint", y = 'Mean Module Expression') + 
  scale_color_manual(values = patientpalette[c(2,5,8)]) +
  stat_compare_means(comparisons = list(c(1,3), c(1,4), c(3,4)),  
                     method = 'wilcox', paired = TRUE, size = 2.1, symnum.args = symnum.args)+ 
  facet_wrap('th') + ggpreset + 
  theme(text = element_text(size = 8), axis.text = element_text(size = 8)) + 
  theme(strip.background = element_rect(color=NA, fill='grey87', size=0.5),
        strip.text.x = element_text(size=8))

ggsave("results/Figure5A_by_patient_v3.pdf",width = 5.3, height = 2, unit = "in", family = "ArialMT")



# Figure 5B - filling out specific subtypes overtime ----------------------------------

all_clones$MODULE = 0
all_clones$MODULE[all_clones$con %in% c("GATA3hi", "Tfh-like", 'Treg-like')] = all_clones$MODULE_7_GATED[all_clones$con %in% c("GATA3hi", "Tfh-like", 'Treg-like')]
all_clones$MODULE[all_clones$con %in% c("Tfh1-like", 'Th1')] = all_clones$MODULE_38_GATED[all_clones$con %in%  c("Tfh1-like", 'Th1')]
all_clones$MODULE[all_clones$con %in% c("Th17")] = all_clones$MODULE_9_GATED[all_clones$con %in%  c("Th17")]
df = all_clones %>% group_by(outcome, patient, timept,con, TRB_CDR3) %>% filter(outcome != 'placebo') %>% summarize(n = n(), on = sum(MODULE)) %>% mutate(rat = on/n) %>%
  filter(n > 1)

df <- df %>% ungroup() %>% mutate(con = factor(con, levels = c("GATA3hi","Tfh-like","Treg-like","Th1","Tfh1-like","Th17")))


ggplot(df, aes(x = timept, y = rat))  + facet_wrap('con') + 
  stat_compare_means(comparisons = list(c(1,3), c(1,4)), size = 2.1,symnum.args = symnum.args) + 
  geom_jitter(size=0.65, stroke=0, col='gray70', width=0.35, height=0.018) + 
  geom_boxplot(aes(fill = con), outlier.shape=NA, lwd=0.4, width=0.18, fatten=3) + 
  scale_fill_manual(values = brewer.pal(6,'Set2')[c(5,4,6,2,1,3)]) + No_legend + ggpreset + 
  labs(x = "Timepoint", y = "Fractional clonal expression") +
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), limits=c(-0.02,1.2)) +
  theme(strip.background = element_rect(color=NA, fill='grey87', size=0.5),strip.text.x = element_text(size=6))

ggsave("results/Figure5B_th_subset_ratios_v4.pdf", width=5.9, height=3.2, unit = "in", family = "ArialMT", useDingbats =FALSE)


# Figure 5C: Degree of Th2 suppression by clinical group -------------------------

# TODO: pull this plot from supplement code


# Figure 5D and E: PCA vs clinical outcome in CD154+ cells -----------------------

Total_meta_data <- readRDS("data/Total_meta_data.rds")

# Remove 7 of the 50 gene modules that are related to B cells, sex-linked genes, and other confounders.
PCA_meta_data <- Total_meta_data %>% 
  select(patient, condition, timept, outcome, starts_with('MODULE'), -ends_with('GATED')) %>% 
  select(-MODULE_19, -MODULE_15, -MODULE_28, -MODULE_5, -MODULE_11, -MODULE_39, -MODULE_34) %>% 
  filter(condition == "CD154") %>% 
  group_by(patient, condition, timept, outcome) %>% 
  summarise_all(~mean(.)) %>% 
  ungroup()

PCA_input <- PCA_meta_data %>%
  select(-c(patient, condition, timept, outcome)) %>% 
  select(c(1:43))

PCA_input_BL <- PCA_meta_data %>%
  filter(timept == "BL") %>%
  select(-c(patient, condition, timept, outcome)) %>% 
  select(c(1:43))

# Generate the PCA with just the baseline data
PCA_result_BL <-  prcomp(PCA_input_BL, scale. = TRUE)

# Then predict on all timepoints
PCA_prediction <-  as.matrix(PCA_input) %*% as.matrix(PCA_result_BL$rotation)
PCA_prediction <- as_tibble(PCA_prediction)
PCA_prediction$patient <- PCA_meta_data$patient
PCA_prediction$outcome <- PCA_meta_data$outcome
PCA_prediction$timept <- PCA_meta_data$timept

ggplot(PCA_prediction, aes(x = outcome, y = PC1, fill = outcome)) + 
  geom_jitter(data = PCA_prediction %>% filter(timept != 'BL'), color = 'white',pch = 21,size = 1, width = 0.2) + 
  geom_jitter(data = PCA_prediction %>% filter(timept == 'BL'), color = 'black',pch = 21,size = 1, width = 0.2) +
  Axis_themes + theme(legend.position = 'none') + 
  labs(x = '') + theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  stat_compare_means(comparisons = list(c(1,2), c(1,3), c(2,3)), size = 2.1)
ggsave("Figure5D_PCA_BL_v2.pdf", width = 3, height = 3, units = 'in', family = 'ArialMT')


# Get the loadings in PC1
PC1 <- PCA_result_BL$rotation[,"PC1"]
PC1 <- PC1[order(PC1)]
PC1_tb <- enframe(PC1) %>% rename(Loading = value)
PC1_tb <- PC1_tb %>% mutate(Module = str_extract(name, "[0-9]+"))
PC1_tb <- PC1_tb %>% mutate(Module = fct_reorder(factor(Module),Loading))

ggplot(PC1_tb, aes(x = Module, y = Loading)) + 
  geom_bar(stat = 'identity', color = 'black',fill = 'grey87', width = 0.8, size = 0.15) + 
  theme(axis.title.x = element_blank()) + labs(x = "PC1 Modules (38.3%)") +
  coord_flip() + Axis_themes

ggsave("Baseline_PC1_loadings.pdf",height = 3, width = 3, family = "ArialMT")


# Pick out the highest few loadings for plotting

rank_to_pick <- 5
PC1_tb_filtered <- PC1_tb %>% 
  top_n(rank_to_pick, wt = Loading) %>% 
  mutate(Module = fct_reorder(factor(Module),Loading))

ggplot(PC1_tb_filtered, aes(x = Module, y = Loading)) + 
  geom_bar(stat = 'identity', color = 'black',fill = 'grey87', width = 0.8, size = 0.15) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  coord_flip() + Axis_themes

ggsave("Fig5E_PC1_loadings_filtered.pdf",height = 2, width = 1.25, family = "ArialMT")


