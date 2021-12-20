# Figure 2: Peanut-reactive T cells have diverse and distinct transcriptional signatures.

# Import presets, libraries, functions, and processed data --------------
source("presets.R")
Total_meta_data <- readRDS("data/Total_meta_data.rds")


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
#TODO: CD154 vs CD137 genes heatmap, need to pull from extended figure 1 code


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

