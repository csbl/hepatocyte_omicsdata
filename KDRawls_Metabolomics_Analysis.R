# Rawls Metabolomics Data Analysis for Liver Data
# Uses Fold Changes and Mann-Whitney U Test to evaluate significance
# Analyze Primary Metabolites, Lipids, and Biogenic Amines
# Imports Raw Metabolomics Data
# Outputs part of Supplementary Data 3 and Figure 3. 
# 04/17/19

# Load Pacakges 
library("openxlsx")
library("tidyverse")
library("matrixStats")
library("gplots")
library("reshape2")
library("ggdendro")
library("grid")
library("gridExtra")
library("ggrepel")
library("ggpubr")

## -------- Create Normalized Dataframes for Metabolomics Data -----------------------

# Read in Raw Metabolomics Datasheets
hep.raw.metabolomics <- read.xlsx("Rawls_Supplementary_data3.xlsx", sheet = 2)
annote1 <- data.frame(hep.raw.metabolomics[,c(1:2)]) %>% dplyr::rename(met = "met_name") %>% dplyr::rename(Annotation = "westcoast_metname")
hep.raw.metabolomics <- hep.raw.metabolomics[,c(2:26)] %>% dplyr::rename(met = "met_name")
mets <- hep.raw.metabolomics$met

hep.raw.lipids <- read.xlsx("Rawls_Supplementary_data3.xlsx", sheet = 3)
hep.raw.lipids <- na.omit(hep.raw.lipids) # Remove lipids that don't have peaks for each sample
annote2 <- data.frame(hep.raw.lipids[,c(1:3)])
hep.raw.lipids <- hep.raw.lipids[,c(3,17:40)]
mets.lp <- hep.raw.lipids$met

hep.raw.amines <- read.xlsx("Rawls_Supplementary_data3.xlsx", sheet = 4)
hep.raw.amines <- na.omit(hep.raw.amines) # Remove biogenic amines that don't have peaks for each sample
annote3 <- data.frame(hep.raw.amines[,c(1:2)])
hep.raw.amines <- hep.raw.amines[,c(2,9:32)]
mets.bg <- hep.raw.amines$met

# Combine Data Frame 
hep.combined.metabolomics <- rbind(hep.raw.metabolomics,hep.raw.lipids,hep.raw.amines)

# Find Row Standard Deviations
hep.combined.metabolomics$stdevs <- rowSds(data.matrix(hep.combined.metabolomics[,c(2:25)]))

# Calculate t stats for each metabolite 
n <- length(hep.combined.metabolomics$stdevs)
t_apap1 <- matrix(0,n,2)
t_apap2 <- matrix(0,n,2)
t_ccl41 <- matrix(0,n,2)
t_ccl42 <- matrix(0,n,2)
t_tcdd1 <- matrix(0,n,2)
t_tcdd2 <- matrix(0,n,2)
t_tce1 <-  matrix(0,n,2)
t_tce2 <-  matrix(0,n,2)


for (i in 1:n) {
  # Compare Conditions with Control Data
  t_apap1[i,1] <-  log2(mean(data.matrix(hep.combined.metabolomics[i,c(6:9)]))/mean(data.matrix(hep.combined.metabolomics[i,c(2:5)])))
  t_apap1[i,2] <-  wilcox.test(data.matrix(hep.combined.metabolomics[i,c(6:9)]),data.matrix(hep.combined.metabolomics[i,c(2:5)]))$p.value  
  
  t_ccl41[i,1] <-  log2(mean(data.matrix(hep.combined.metabolomics[i,c(10:13)]))/mean(data.matrix(hep.combined.metabolomics[i,c(2:5)])))
  t_ccl41[i,2] <-  wilcox.test(data.matrix(hep.combined.metabolomics[i,c(10:13)]),data.matrix(hep.combined.metabolomics[i,c(2:5)]))$p.value 
  
  t_tcdd1[i,1] <-  log2(mean(data.matrix(hep.combined.metabolomics[i,c(14:17)]))/mean(data.matrix(hep.combined.metabolomics[i,c(2:5)])))
  t_tcdd1[i,2] <-  wilcox.test(data.matrix(hep.combined.metabolomics[i,c(14:17)]),data.matrix(hep.combined.metabolomics[i,c(2:5)]))$p.value 
  
  t_tce1[i,1] <-  log2(mean(data.matrix(hep.combined.metabolomics[i,c(18:21)]))/mean(data.matrix(hep.combined.metabolomics[i,c(2:5)])))
  t_tce1[i,2] <-  wilcox.test(data.matrix(hep.combined.metabolomics[i,c(18:21)]),data.matrix(hep.combined.metabolomics[i,c(2:5)]))$p.value  
}

# Find Metabolites Produced or Consumed in the control condition (Control vs Blank Media)
t_ctl1 <- matrix(0,n,2)

for (i in 1:n) {
  t_ctl1[i,1] <-  log2(mean(data.matrix(hep.combined.metabolomics[i,c(2:5)]))/(mean(data.matrix(hep.combined.metabolomics[i,c(22:25)]))))
  t_ctl1[i,2] <-  wilcox.test(data.matrix(hep.combined.metabolomics[i,c(2:5)]),data.matrix(hep.combined.metabolomics[i,c(22:25)]))$p.value   
}

# Filter for Significance 
z_apap1 <- t_apap1
z_ccl41 <- t_ccl41 
z_tcdd1 <- t_tcdd1 
z_tce1 <- t_tce1
z_ctl1 <- t_ctl1 

for (k in 1:n) {
  if (z_apap1[k,2] > 0.05) {
    z_apap1[k,1] = 0
  }
}

for (k in 1:n) {
  if (z_ccl41[k,2] > 0.05) {
    z_ccl41[k,1] = 0
  }
}

for (k in 1:n) {
  if (z_tcdd1[k,2] > 0.05) {
    z_tcdd1[k,1] = 0
  }
}

for (k in 1:n) {
  if (z_tce1[k,2] > 0.05) {
    z_tce1[k,1] = 0
  }
}

for (k in 1:n) {
  if (z_ctl1[k,2] > 0.05) {
    z_ctl1[k,1] = 0
  }
}

# Reassign to new data frame
normal.hep.metabolomics <- data.frame(t_apap1[,1],t_ccl41[,1],t_tcdd1[,1],t_tce1[,1])
colnames(normal.hep.metabolomics) <- c("t_apap1","t_ccl41","t_tcdd1","t_tce1")
normal.hep.metabolomics$met <- c(mets,mets.lp,mets.bg)

filtered.hep.metabolomics <- data.frame(z_apap1[,1],z_ccl41[,1],z_tcdd1[,1],z_tce1[,1])
colnames(filtered.hep.metabolomics) <- c("t_apap1","t_ccl41","t_tcdd1","t_tce1")
filtered.hep.metabolomics$met <- c(mets,mets.lp,mets.bg)

control.hep.metabolomics <- data.frame(t_ctl1[,1])
colnames(control.hep.metabolomics) <- c("t_ctl1")
control.hep.metabolomics$met <- c(mets,mets.lp,mets.bg)
combined.annotation <- rbind(annote1,annote2[,c(2:3)],annote3)

## -------- Create Heatmap of Metabolomics Data -----------------------

## Create Metabolomics Heatmap 
data_hmap <- normal.hep.metabolomics[,c(1:4)]

rownames(data_hmap) <- normal.hep.metabolomics$met
hr <- hclust(dist(data_hmap, method = "euclidean"), method = "complete")
hc <- hclust(as.dist(1-cor(data_hmap, method = "spearman")), method = "complete")

# ggplot heatmap 
melted_metabolites <- melt(normal.hep.metabolomics[,c(1:5)], id.vars = "met", 
    variable.name = "condition",value.name = "fold_change") %>% rename(metabolite = "met")
metabolite_order <- as.list(normal.hep.metabolomics[hr$order,5]) %>% unlist()
condition_unorder <- colnames(normal.hep.metabolomics[,c(1:4)])
condition_order <- condition_unorder[hc$order]
melted_metabolites$metabolite <- factor(melted_metabolites$metabolite,levels = metabolite_order, ordered = TRUE)
melted_metabolites$condition <- factor(melted_metabolites$condition,levels = condition_order, ordered = TRUE)

melted_metabolites$fold_change[melted_metabolites$fold_change < -2] <- -2
melted_metabolites$fold_change[melted_metabolites$fold_change > 2] <- 2

f3.1 <- ggplot(melted_metabolites, aes(x = condition,y = metabolite)) + geom_tile(aes(fill = as.numeric(melted_metabolites$fold_change)),size = 0.01) +
  scale_fill_gradient2(low = "blue",mid = "white",high = "red",limits = c(-2,2), labels = c("Increased levels","Decreased levels"), breaks = c(2,-2)) + 
  scale_x_discrete(position = "bottom", labels = c(t_apap1 = expression(paste("APAP 6hr")),t_ccl41 = expression(paste("CCl"[4]*" 6hr")), 
      t_tcdd1 = expression(paste("TCDD 6hr")), t_tce1 = expression(paste("TCE 6hr")))) + 
  theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank()) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(family = "Arial",face = "bold",size = 30, 
        angle = 30, vjust = 0.8, hjust = 1)) +
  theme(axis.ticks.y = element_blank()) +
  theme(plot.margin = unit(c(0.5,-0.5,0.5,0.5),"cm")) +
  theme(legend.position = "top") + theme(legend.title=element_text(size=26), legend.text=element_text(size=24)) +
  guides(fill=guide_legend(title = expression(paste("log"[2]*" metabolite changes")))) + 
  labs(x = "", y = "")


# Create GG Plot for dendrogram of genes
f3.2 <- ggdendrogram(hr, rotate = TRUE, theme_dendro = TRUE) + 
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) + 
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank())

f3.3 <- ggdendrogram(hc, rotate = FALSE, theme_dendro = TRUE) +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

## -------- Create Scatter Plots of Metabolomics Data -----------------------


# Calculate absolute fold changes
n <- length(hep.combined.metabolomics$stdevs)
m_apap1 <- matrix(0,n,2)
m_apap2 <- matrix(0,n,2)
m_ccl41 <- matrix(0,n,2)
m_ccl42 <- matrix(0,n,2)
m_tcdd1 <- matrix(0,n,2)
m_tcdd2 <- matrix(0,n,2)
m_tce1 <-  matrix(0,n,2)
m_tce2 <-  matrix(0,n,2)


for (i in 1:n) {
  m_apap1[i,1] <-  (mean(data.matrix(hep.combined.metabolomics[i,c(6:9)]))/mean(data.matrix(hep.combined.metabolomics[i,c(2:5)])))
  m_apap1[i,2] <-  wilcox.test(as.numeric(data.matrix(hep.combined.metabolomics[i,c(6:9)])),as.numeric(data.matrix(hep.combined.metabolomics[i,c(2:5)])))$p.value  

  m_ccl41[i,1] <-  (mean(data.matrix(hep.combined.metabolomics[i,c(10:13)]))/ mean(data.matrix(hep.combined.metabolomics[i,c(2:5)])))
  m_ccl41[i,2] <-  wilcox.test(data.matrix(hep.combined.metabolomics[i,c(10:13)]),data.matrix(hep.combined.metabolomics[i,c(2:5)]))$p.value 
  
  m_tcdd1[i,1] <- (mean(data.matrix(hep.combined.metabolomics[i,c(14:17)]))/mean(data.matrix(hep.combined.metabolomics[i,c(2:5)])))
  m_tcdd1[i,2] <-  wilcox.test(data.matrix(hep.combined.metabolomics[i,c(14:17)]),data.matrix(hep.combined.metabolomics[i,c(2:5)]))$p.value 
  
  m_tce1[i,1] <- (mean(data.matrix(hep.combined.metabolomics[i,c(18:21)]))/mean(data.matrix(hep.combined.metabolomics[i,c(2:5)])))
  m_tce1[i,2] <-  wilcox.test(data.matrix(hep.combined.metabolomics[i,c(18:21)]),data.matrix(hep.combined.metabolomics[i,c(2:5)]))$p.value 
  
  m_apap2[i,1] <-  (mean(data.matrix(hep.combined.metabolomics[i,c(6:9)]))/mean(data.matrix(hep.combined.metabolomics[i,c(22:25)])))
  m_apap2[i,2] <-  wilcox.test(as.numeric(data.matrix(hep.combined.metabolomics[i,c(6:9)])),as.numeric(data.matrix(hep.combined.metabolomics[i,c(22:25)])))$p.value  
  
  m_ccl42[i,1] <-  (mean(data.matrix(hep.combined.metabolomics[i,c(10:13)]))/ mean(data.matrix(hep.combined.metabolomics[i,c(22:25)])))
  m_ccl42[i,2] <-  wilcox.test(data.matrix(hep.combined.metabolomics[i,c(10:13)]),data.matrix(hep.combined.metabolomics[i,c(22:25)]))$p.value 
  
  m_tcdd2[i,1] <- (mean(data.matrix(hep.combined.metabolomics[i,c(14:17)]))/ mean(data.matrix(hep.combined.metabolomics[i,c(22:25)])))
  m_tcdd2[i,2] <-  wilcox.test(data.matrix(hep.combined.metabolomics[i,c(14:17)]),data.matrix(hep.combined.metabolomics[i,c(22:25)]))$p.value 
  
  m_tce2[i,1] <- (mean(data.matrix(hep.combined.metabolomics[i,c(18:21)]))/ mean(data.matrix(hep.combined.metabolomics[i,c(22:25)])))
  m_tce2[i,2] <-  wilcox.test(data.matrix(hep.combined.metabolomics[i,c(18:21)]),data.matrix(hep.combined.metabolomics[i,c(22:25)]))$p.value 
}



# Construct the data frame for subsequent analyses
apap1 <- data.frame(m_apap1,m_apap2)
ccl41 <- data.frame(m_ccl41,m_ccl42)
tcdd1 <- data.frame(m_tcdd1,m_tcdd2)
tce1 <- data.frame(m_tce1,m_tce2)

colnames(apap1) <- c("control_FC","control_pval","blank_FC","blank_pval")
colnames(ccl41) <- c("control_FC","control_pval","blank_FC","blank_pval")
colnames(tcdd1) <- c("control_FC","control_pval","blank_FC","blank_pval")
colnames(tce1) <- c("control_FC","control_pval","blank_FC","blank_pval")

apap1$mets <- c(mets,mets.lp,mets.bg)
ccl41$mets <- c(mets,mets.lp,mets.bg)
tcdd1$mets <- c(mets,mets.lp,mets.bg)
tce1$mets <- c(mets,mets.lp,mets.bg)

# Create a separate variable for significance
apap1 <- apap1 %>% mutate(identified = grepl("\\D",mets)) %>% mutate(control_sig = as.integer(control_pval < 0.05)) %>% mutate(blank_sig = as.integer(blank_pval < 0.05))
ccl41 <- ccl41 %>% mutate(identified = grepl("\\D",mets)) %>% mutate(control_sig = as.integer(control_pval < 0.05)) %>% mutate(blank_sig = as.integer(blank_pval < 0.05))
tcdd1 <- tcdd1 %>% mutate(identified = grepl("\\D",mets)) %>% mutate(control_sig = as.integer(control_pval < 0.05)) %>% mutate(blank_sig = as.integer(blank_pval < 0.05))
tce1 <- tce1 %>% mutate(identified = grepl("\\D",mets)) %>% mutate(control_sig = as.integer(control_pval < 0.05)) %>% mutate(blank_sig = as.integer(blank_pval < 0.05))

apap1 <- apap1 %>% mutate(control_compare = as.integer(control_FC > 1)) %>% mutate(blank_compare = as.integer(blank_FC > 1)) %>% mutate(color = "black")
ccl41 <- ccl41 %>% mutate(control_compare = as.integer(control_FC > 1)) %>% mutate(blank_compare = as.integer(blank_FC > 1)) %>% mutate(color = "black")
tcdd1 <- tcdd1 %>% mutate(control_compare = as.integer(control_FC > 1)) %>% mutate(blank_compare = as.integer(blank_FC > 1)) %>% mutate(color = "black")
tce1 <- tce1 %>% mutate(control_compare = as.integer(control_FC > 1)) %>% mutate(blank_compare = as.integer(blank_FC > 1)) %>% mutate(color = "black")

# Label metabolite production/consumption behavior
apap1[apap1$control_sig == 0 & apap1$blank_sig == 0,"color"] <- "black"
apap1[(apap1$control_sig == 1 | apap1$blank_sig == 1) & apap1$control_FC > 1 & apap1$blank_FC > 1,"color"] <- "Increased Produciton"
apap1[(apap1$control_sig == 1 | apap1$blank_sig == 1) & apap1$control_FC > 1 & apap1$blank_FC < 1,"color"] <- "Decreased Consumption"
apap1[(apap1$control_sig == 1 | apap1$blank_sig == 1) & apap1$control_FC < 1 & apap1$blank_FC > 1 ,"color"] <- "Decreased Production"
apap1[(apap1$control_sig == 1 | apap1$blank_sig == 1) & apap1$control_FC < 1 & apap1$blank_FC < 1 ,"color"] <- "Increased Consumption"

ccl41[ccl41$control_sig == 0 & ccl41$blank_sig == 0,"color"] <- "black"
ccl41[(ccl41$control_sig == 1 | ccl41$blank_sig == 1) & ccl41$control_FC > 1 & ccl41$blank_FC > 1,"color"] <- "Increased Production"
ccl41[(ccl41$control_sig == 1 | ccl41$blank_sig == 1) & ccl41$control_FC > 1 & ccl41$blank_FC < 1,"color"] <- "Decreased Consumption"
ccl41[(ccl41$control_sig == 1 | ccl41$blank_sig == 1) & ccl41$control_FC < 1 & ccl41$blank_FC > 1 ,"color"] <- "Decreased Production"
ccl41[(ccl41$control_sig == 1 | ccl41$blank_sig == 1) & ccl41$control_FC < 1 & ccl41$blank_FC < 1 ,"color"] <- "Increased Consumption"

tcdd1[tcdd1$control_sig == 0 & tcdd1$blank_sig == 0,"color"] <- "black"
tcdd1[(tcdd1$control_sig == 1 | tcdd1$blank_sig == 1) & tcdd1$control_FC > 1 & tcdd1$blank_FC > 1,"color"] <- "Increased Production"
tcdd1[(tcdd1$control_sig == 1 | tcdd1$blank_sig == 1) & tcdd1$control_FC > 1 & tcdd1$blank_FC < 1,"color"] <- "Decreased Consumption"
tcdd1[(tcdd1$control_sig == 1 | tcdd1$blank_sig == 1) & tcdd1$control_FC < 1 & tcdd1$blank_FC > 1 ,"color"] <- "Decreased Production"
tcdd1[(tcdd1$control_sig == 1 | tcdd1$blank_sig == 1) & tcdd1$control_FC < 1 & tcdd1$blank_FC < 1 ,"color"] <- "Increased Consumption"

tce1[tce1$control_sig == 0 & tce1$blank_sig == 0,"color"] <- "black"
tce1[(tce1$control_sig == 1 | tce1$blank_sig == 1) & tce1$control_FC > 1 & tce1$blank_FC > 1,"color"] <- "Increased Production"
tce1[(tce1$control_sig == 1 | tce1$blank_sig == 1) & tce1$control_FC > 1 & tce1$blank_FC < 1,"color"] <- "Decreased Consumption"
tce1[(tce1$control_sig == 1 | tce1$blank_sig == 1) & tce1$control_FC < 1 & tce1$blank_FC > 1 ,"color"] <- "Decreased Production"
tce1[(tce1$control_sig == 1 | tce1$blank_sig == 1) & tce1$control_FC < 1 & tce1$blank_FC < 1 ,"color"] <- "Increased Consumption"

f3.4 <- apap1 %>% filter(control_sig == 1 | blank_sig == 1) %>% filter(color != "black") %>%  
  ggplot(aes(x=log2(blank_FC), y=log2(control_FC))) + 
  geom_point(alpha=0.6, size=4, aes(color = factor(color))) + 
  theme_bw() + 
  scale_color_manual(name = "", values = c("gold","lightblue","orange","purple")) +
  annotate("text",label = "Increased Production (175)", x = 2.8, y = 2.8, size = 9, color = "black", family = "sans") +
  annotate("text",label = "Decreased Production (145)", x = 2.8, y = -2.8, size = 9, color = "black", family = "sans") +  
  annotate("text",label = "Increased Consumption (38)", x = -2.65, y = -2.8, size = 9, color = "black", family = "sans") +
  annotate("text",label = "Decreased Consumption (13)", x = -2.65, y = 2.8, size = 9, color = "black", family = "sans") +
  ggtitle("APAP") + 
  xlab("Fold change (treated/blank)") + ylab("Fold change (treated/control)") +
  theme(axis.text = element_text(family = "Arial",face = "bold",size = 20)) + 
  theme(axis.title = element_text(family = "Arial",face = "bold",size = 24)) + 
  theme(plot.title = element_text(family = "Arial",face = "bold",size = 26)) + 
  geom_vline(xintercept = 0, linetype="dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme(legend.position="none") +
  xlim(-4,4) + ylim(-3,3)

f3.5 <- ccl41 %>% filter(control_sig == 1 | blank_sig == 1) %>% filter(color != "black") %>% 
  ggplot(aes(x=log2(blank_FC), y=log2(control_FC))) + 
  geom_point(alpha=0.6, size=4, aes(color = factor(color))) + 
  theme_bw() + 
  scale_color_manual(name = "", values = c("gold","lightblue","orange","purple")) +
  annotate("text",label = "Increased Production (105)", x = 2.8, y = 2.8, size = 9, color = "black", family = "sans") +
  annotate("text",label = "Decreased Production (208)", x = 2.8, y = -2.8, size = 9, color = "black", family = "sans") +  
  annotate("text",label = "Increased Consumption (32)", x = -2.65, y = -2.8, size = 9, color = "black", family = "sans") +
  annotate("text",label = "Decreased Consumption (19)", x = -2.65, y = 2.8, size = 9, color = "black", family = "sans") +
  ggtitle(expression(bold(paste("CCl"[4]*"")))) + 
  xlab("Fold change (treated/blank)") + ylab("Fold change (treated/control)") +
  theme(axis.text = element_text(family = "Arial",face = "bold",size = 20)) + 
  theme(axis.title = element_text(family = "Arial",face = "bold",size = 24)) + 
  theme(plot.title = element_text(family = "Arial",face = "bold",size = 26)) + 
  geom_vline(xintercept = 0, linetype="dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme(legend.position="none") +
  xlim(-4,4) + ylim(-3,3)

f3.6 <- tcdd1 %>% filter(control_sig == 1 | blank_sig == 1) %>% filter(color != "black") %>% 
  ggplot(aes(x=log2(blank_FC), y=log2(control_FC))) + 
  geom_point(alpha=0.6, size=4, aes(color = factor(color))) + 
  theme_bw() + 
  scale_color_manual(name = "", values = c("gold","lightblue","orange","purple")) +
  annotate("text",label = "Increased Production (127)", x = 2.8, y = 2.8, size = 9, color = "black") +
  annotate("text",label = "Decreased Production (209)", x = 2.8, y = -2.8, size = 9, color = "black") +  
  annotate("text",label = "Increased Consumption (53)", x = -2.65, y = -2.8, size = 9, color = "black") +
  annotate("text",label = "Decreased Consumption (15)", x = -2.65, y = 2.8, size = 9, color = "black") +
  ggtitle("TCDD") + 
  xlab("Fold change (treated/blank)") + ylab("Fold change (treated/control)") +
  theme(axis.text = element_text(family = "Arial",face = "bold",size = 20)) + 
  theme(axis.title = element_text(family = "Arial",face = "bold",size = 24)) + 
  theme(plot.title = element_text(family = "Arial",face = "bold",size = 26)) + 
  geom_vline(xintercept = 0, linetype="dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme(legend.position="none") +
  xlim(-4,4) + ylim(-3,3)

f3.7 <- tce1 %>% filter(control_sig == 1 | blank_sig == 1) %>% filter(color != "black") %>% 
  ggplot(aes(x=log2(blank_FC), y=log2(control_FC))) + 
  theme_bw() + 
  geom_point(alpha=0.6, size=4, aes(color = factor(color))) + 
  scale_color_manual(name = "", values = c("gold","lightblue","orange","purple")) +
  annotate("text",label = "Increased Production (111)", x = 2.8, y = 2.8, size = 9, color = "black") +
  annotate("text",label = "Decreased Production (190)", x = 2.8, y = -2.8, size = 9, color = "black") +  
  annotate("text",label = "Increased Consumption (38)", x = -2.65, y = -2.8, size = 9, color = "black") +
  annotate("text",label = "Decreased Consumption (18)", x = -2.65, y = 2.8, size = 9, color = "black") +
  ggtitle("TCE") + 
  xlab("Fold change (treated/blank)") + ylab("Fold change (treated/control)") +
  theme(axis.text = element_text(family = "Arial",face = "bold",size = 20)) + 
  theme(axis.title = element_text(family = "Arial",face = "bold",size = 24)) + 
  theme(plot.title = element_text(family = "Arial",face = "bold",size = 26)) + 
  geom_vline(xintercept = 0, linetype="dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme(legend.position="none", legend.text=element_text(size=16)) +
  xlim(-4,4) + ylim(-3,3)

scatter_figs <- ggarrange(f3.4,f3.5,f3.6,f3.7,nrow = 4,ncol = 1,labels = c("A","B","C","D"), font.label = list(size = 30), common.legend = FALSE)
ann_scatter <- annotate_figure(scatter_figs, top = text_grob("Hepatocyte Metabolomics by Drug",face= "bold", size = 24, family = "Arial"))
figure3 <- ggarrange(ann_scatter,f3.1,nrow = 1,ncol = 2, labels = c("","E"), font.label = list(size = 30))
figure3

ggsave(figure3,filename = "Figure3_Metabolomicsfigure.png",width = 26,height = 25)

