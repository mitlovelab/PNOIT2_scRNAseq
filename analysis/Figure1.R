# Figure 1: Peanut-reactive T cells decrease in frequency during OIT.


source("analysis/presets.R")

# FIG 1C: Flow data -----------------
# Flow data summary:
Flowdata <- read_csv("data/metadata/PNOIT2_flow_sort_summary.csv") %>%
  mutate(Patient = factor(paste("P", Patient, sep=""),
                          levels=c("P105","P106","P111","P33","P90","P93","P69","P95","P97","P84","P96","P107")),
         Timept = factor(Timept, levels=c("BL","BU","MN","AV")),
         group = ifelse(Patient %in% c("P84","P96","P107"), "Placebo", "Treatment"),
         outcome = case_when(Patient %in% c("P105","P106","P111") ~ "Tolerance",
                             Patient %in% c("P33","P90","P93") ~ "Partial tolerance",
                             Patient %in% c("P69","P95","P97") ~ "Treatment failure",
                             Patient %in% c("P84","P96","P107") ~ "Placebo")) %>%
  mutate(group = factor(group,levels=c("Treatment","Placebo")), grouptimept = paste(group,Timept))

#Flowdata$group <- factor(Flowdata$group, levels=c("Treatment","Placebo"))


# CD137
toplot <- "CD137+"
# Plot treatment and placebo in separate charts but on the same axes:
lims <- c(0, 3.7)

Flowdata$toplot <- Flowdata[,which(names(Flowdata) == toplot)][[1]]

Treat <- Flowdata[which(Flowdata$group == "Treatment"),]
ggplot(Treat, aes(x=Timept, y=toplot, fill=Timept)) +
  geom_violin(trim=FALSE, fill="gray95") + geom_boxplot(width = 0.2, fill = 'white') +
  geom_quasirandom(dodge.width = 0.5, size=1.9, aes(col=Patient)) + 
  scale_color_manual(values=patientpalette) + 
  ylab(paste("Percent ", toplot, " in CD4 memory T cells", sep="")) + 
  scale_x_discrete(labels=c("BL","BU","MN","AV")) + 
  scale_y_continuous(limits=lims) + ggpreset + 
ggsave(paste("results/", substr(toplot,1,5), "_Treatment_flow_summary.pdf",sep="_"), 
       width=2.8, height=2.4, useDingbats=FALSE, family = "ArialMT")

Placebo <- Flowdata[which(Flowdata$group == "Placebo"),]
ggplot(Placebo, aes(x=Timept, y=toplot, fill=Timept)) +
  geom_violin(trim=FALSE, fill="gray95") + geom_boxplot(width = 0.2, fill = 'white') +
  geom_quasirandom(dodge.width = 0.5, size=1.9, aes(col=Patient)) + 
  scale_color_manual(values=patientpalette[10:12]) + 
  ylab(paste("Percent ", toplot, " in CD4 memory T cells", sep="")) + 
  scale_x_discrete(labels=c("BL","BU","MN","AV")) + 
  scale_y_continuous(limits=lims) + ggpreset
ggsave(paste("results/", substr(toplot,1,5), "_Placebo_flow_summary.pdf",sep="_"), 
       width=2.2, height=2.4, useDingbats=FALSE, family = "ArialMT")

#CD154 equivalent
toplot <- "CD154+"
lims <- c(0, 7.5)

Flowdata$toplot <- Flowdata[,which(names(Flowdata) == toplot)][[1]]

Treat <- Flowdata[which(Flowdata$group == "Treatment"),]
ggplot(Treat, aes(x=Timept, y=toplot, fill=Timept)) +
  geom_violin(trim=FALSE, fill="gray95") + geom_boxplot(width = 0.2, fill = 'white') +
  geom_quasirandom(dodge.width = 0.5, size=1.9, aes(col=Patient)) + 
  scale_color_manual(values=patientpalette) + 
  ylab(paste("Percent ", toplot, " in CD4 memory T cells", sep="")) + 
  scale_x_discrete(labels=c("BL","BU","MN","AV")) + 
  scale_y_continuous(limits=lims) + ggpreset
ggsave(paste("results/", substr(toplot,1,5), "_Treatment_flow_summary.pdf",sep="_"), 
       width=2.8, height=2.4, useDingbats=FALSE, family = "ArialMT")

Placebo <- Flowdata[which(Flowdata$group == "Placebo"),]
ggplot(Placebo, aes(x=Timept, y=toplot, fill=Timept)) +
  geom_violin(trim=FALSE, fill="gray95") + geom_boxplot(width = 0.2, fill = 'white') +
  geom_quasirandom(dodge.width = 0.5, size=1.9, aes(col=Patient)) + 
  scale_color_manual(values=patientpalette[10:12]) + 
  ylab(paste("Percent ", toplot, " in CD4 memory T cells", sep="")) + 
  scale_x_discrete(labels=c("BL","BU","MN","AV")) + 
  scale_y_continuous(limits=lims) + ggpreset
ggsave(paste("results/", substr(toplot,1,5), "_Placebo_flow_summary.pdf",sep="_"), 
       width=2.2, height=2.4, useDingbats=FALSE, family = "ArialMT")
