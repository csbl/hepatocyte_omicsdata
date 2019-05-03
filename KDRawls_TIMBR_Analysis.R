# Rawls TIMBR Score Analysis for liver Data
# Imports TIMBR Production Scores
# Analyze similarities and differences between compounds
# Outputs Figures 4 and 5
# This script is adapted from blais_timbr_analysis.R
# Script comes from "Reconciled rat and human metabolic networks for comparative toxicogenomics analyses and biomarker predictions" by Blais et al.
# Article and data files are published at https://www.nature.com/articles/ncomms14250 and github.com/csbl/ratcon1
# 04/17/19

# Load Pacakges 
library(tidyverse)
library(openxlsx)
library(gplots)
library(reshape2)
library(gridExtra)
library(matrixStats)
library(ggpubr)
library(VennDiagram)
library(ggplotify)

# Source Files
# Read in helper file from "Reconciled rat and human metabolic networks for comparative toxicogenomics analyses and biomarker predictions" by Blais et al.
# Article and data files are published at https://www.nature.com/articles/ncomms14250
source("ncomm_helper.R")

# Read in workbook from "Reconciled rat and human metabolic networks for comparative toxicogenomics analyses and biomarker predictions" by Blais et al.
# Article and data files are published at https://www.nature.com/articles/ncomms14250
rxn.info.load = readWorkbook("ncomm_blais_supplementary_data3.xlsx",startRow = 2) %>% as.tbl


# Load in TIMBR Data from Supplementary Data 4
timbr.predictions.load = read.xlsx("Rawls_Supplementary_data4.xlsx",sheet = 3) %>% as.tbl
expression.summary.load = read.xlsx("Rawls_Supplementary_data4.xlsx",sheet = 2) %>% as.tbl

rxn.info = rxn.info.load %>% 
  mutate(enabled = ifelse(!is.na(enabled), enabled == "true", F)) %>% 
  mutate(enabled_hsa = ifelse(!is.na(enabled_hsa), enabled_hsa == "true", F)) %>% 
  mutate(enabled_rno = ifelse(!is.na(enabled_rno), enabled_rno == "true", F))%>% 
  mutate(n_hsa = as.numeric(as.character(n_gene_hsa)),
         n_rno = as.numeric(as.character(n_gene_rno)))  %>%
  mutate(met = gsub(" exchange| demand","",trimws(rxn_name))) %>%
  mutate(met = ifelse(rxn_id == "RCR30122","vitamin C",met))

timbr.predictions = timbr.predictions.load %>%
  mutate(production_score = as.numeric(as.character(production_score))) %>% 
  left_join(rxn.info %>% select(rxn_id, rxn_name, met)) %>% 
  mutate(drug_id = paste0(drug_id,"_",time_id))

expression.summary = expression.summary.load %>% 
  mutate(drug_selected = ifelse(!is.na(drug_selected), grepl("true", tolower(drug_selected)), NA)) %>%
  mutate(dose_selected = ifelse(!is.na(dose_selected), grepl("true", tolower(dose_selected)), NA)) %>% 
  mutate(drug_name = drug_id) %>% 
  mutate(drug_id = paste0(drug_id,"_",time_id))


# Determine rank by TIMBR Score
timbr.compare = timbr.predictions %>%  
  semi_join(expression.summary %>% filter(drug_selected, dose_selected) %>% 
              select(organism_id, drug_id, time_id, dose_id) %>% distinct) %>%
  left_join(rxn.info %>% select(rxn_id, rxn_name, met)) %>% 
  dcast(drug_id + rxn_id + rxn_name + met ~ organism_id, 
        value.var = "production_score", fill = NA) %>% as.tbl %>% 
  arrange(rno) %>% group_by(drug_id) %>%  mutate(rno_rank = 1:n(), rno_pctile = 100 * rno_rank / n()) %>% ungroup %>%
  mutate(rno_qrtile = ifelse(rno_pctile > 75, 1, ifelse(rno_pctile < 25, -1, 0))) %>% 
  mutate(rno_border = ifelse(rno_pctile > 75 | rno_pctile < 25, "#000000", ifelse(rno_pctile < 25, "#000000", NA))) %>% 
  left_join(expression.summary %>% select(drug_id, drug_name) %>% mutate(drug_abbrev = drug_id) %>% distinct)

timbr.drug.summary = timbr.compare %>%
  group_by(drug_id, drug_name, drug_abbrev) %>% 
  summarize(drug_sd_rno = sd(rno),
            drug_ave_rno = mean(rno)) %>% ungroup %>% 
  mutate(drug_ave = drug_ave_rno) %>% 
  mutate(drug_index = 1:n())

timbr.rxn.summary = timbr.compare %>%
  group_by(rxn_id, rxn_name, met) %>% 
  summarize(rxn_ave_rno = mean(rno)) %>% ungroup %>% 
  mutate(rxn_ave = rxn_ave_rno) %>% 
  mutate(rxn_index = 1:n())

# Plot TIMBR Summary Figure (Figure 4A)
f4.1 <- timbr.compare %>% filter(drug_id == "acetaminophen_t1" | drug_id == "carbontetrachloride_t1" | drug_id == "TCDD_t1" | drug_id == "trichloroethylene_t1") %>%  
  ggplot(aes(x = drug_id, y = rno)) + geom_jitter(width = 0.2,aes(color = drug_name, alpha = 0.4, size = 3)) + 
  theme_bw() +
  scale_x_discrete(position = "bottom", labels = c(acetaminophen_t1 = expression(paste("APAP")),carbontetrachloride_t1 = expression(paste("CCl"[4]*"")), 
                                                   TCDD_t1 = expression(paste("TCDD")), trichloroethylene_t1 = expression(paste("TCE")))) + 
  scale_color_manual(name = "",values = c("purple2","grey40","springgreen","orange")) +
  ggtitle("Hepatocyte 6hr TIMBR Production Scores") + 
  ylab("Production Score") + xlab("Experimental Condition") +
  geom_hline(yintercept = 1, linetype = "dotted", color = "red", size = 3) + 
  geom_hline(yintercept = -1, linetype = "dotted", color = "red", size = 3) + 
  theme(plot.title = element_text(family = "Arial",face = "bold",size = 45)) +
  theme(axis.text = element_text(family = "Arial",size = 40)) + 
  theme(plot.margin = unit(c(2,2,2,2), "lines")) +
  theme(axis.title = element_text(family = "Arial",face = "bold",size = 40)) +
  theme(legend.position = "none") + 
  guides(fill=FALSE, color=FALSE)

## -------------------- Create TIMBR Venn Diagrams -------------------------------
apap_up <- timbr.compare %>% filter(drug_id == "acetaminophen_t1") %>% filter(rno > 0) %>% select(met) %>% unlist()
apap_dwn <- timbr.compare %>% filter(drug_id == "acetaminophen_t1") %>% filter(rno < 0) %>% select(met) %>% unlist()

ccl4_up <- timbr.compare %>% filter(drug_id == "carbontetrachloride_t1") %>% filter(rno > 0) %>% select(met) %>% unlist()
ccl4_dwn <- timbr.compare %>% filter(drug_id == "carbontetrachloride_t1") %>% filter(rno < 0) %>% select(met) %>% unlist()

tcdd_up <- timbr.compare %>% filter(drug_id == "TCDD_t1") %>% filter(rno > 0) %>% select(met) %>% unlist()
tcdd_dwn <- timbr.compare %>% filter(drug_id == "TCDD_t1") %>% filter(rno < 0) %>% select(met) %>% unlist()

tce_up <- timbr.compare %>% filter(drug_id == "trichloroethylene_t1") %>% filter(rno > 0) %>% select(met) %>% unlist()
tce_dwn <- timbr.compare %>% filter(drug_id == "trichloroethylene_t1") %>% filter(rno < 0) %>% select(met) %>% unlist()

# Draw Venn Diagram 
x <- calculate.overlap(list(apap_up,ccl4_up,tcdd_up,tce_up))
y <- calculate.overlap(list(apap_dwn,ccl4_dwn,tcdd_dwn,tce_dwn))

up.area1 <- length(apap_up)
up.area2 <- length(ccl4_up)
up.area3 <- length(tcdd_up)
up.area4 <- length(tce_up)
up.n12 <- length(intersect(apap_up,ccl4_up))
up.n13 <- length(intersect(apap_up,tcdd_up))
up.n14 <- length(intersect(apap_up,tce_up))
up.n23 <- length(intersect(ccl4_up,tcdd_up))
up.n24 <- length(intersect(ccl4_up,tce_up))
up.n34 <- length(intersect(tcdd_up,tce_up))
up.n123 <- length(x$a12) + length(x$a6)
up.n124 <- length(x$a11) + length(x$a6)
up.n134 <- length(x$a5) + length(x$a6)
up.n234 <- length(x$a7) + length(x$a6)
up.n1234 <- length(x$a6)

# Create Figure 4B, TIMBR Up Scores
f4.11 <- as_ggplot(draw.quad.venn(area1 = up.area1, area2 = up.area2, area3 = up.area3, area4 = up.area4, 
        n12 = up.n12, n13 = up.n13, n14 = up.n14, n23 = up.n23, n24 = up.n24, n34 = up.n34, n123 = up.n123, 
        n124 = up.n124, n134 = up.n134, n234 = up.n234, n1234 = up.n1234,
        category = c("APAP",expression(paste("CCl"[4]*"")),"TCDD","TCE"),
        fill = c("Purple","grey","Green","orange"), alpha = 0.4, fontface = rep("plain",15), fontfamily = rep("sans",15),  
        cex = rep(4,15), cat.cex = c(3,3,3,3), cat.fontface = rep("plain",4), cat.fontfamily = rep("sans",4))) 

f4.12 <- as.ggplot(f4.11) + ggtitle("Metabolites predicted to increase")+ 
  theme(plot.margin = unit(c(3,3,3,3), "lines")) +
  theme(plot.title = element_text(hjust = 0.5, size = 36, family = "Arial", face = "bold"))
f4.12


dwn.area1 <- length(apap_dwn)
dwn.area2 <- length(ccl4_dwn)
dwn.area3 <- length(tcdd_dwn)
dwn.area4 <- length(tce_dwn)
dwn.n12 <- length(intersect(apap_dwn,ccl4_dwn))
dwn.n13 <- length(intersect(apap_dwn,tcdd_dwn))
dwn.n14 <- length(intersect(apap_dwn,tce_dwn))
dwn.n23 <- length(intersect(ccl4_dwn,tcdd_dwn))
dwn.n24 <- length(intersect(ccl4_dwn,tce_dwn))
dwn.n34 <- length(intersect(tcdd_dwn,tce_dwn))
dwn.n123 <- length(y$a12) + length(y$a6)
dwn.n124 <- length(y$a11) + length(y$a6)
dwn.n134 <- length(y$a5) + length(y$a6)
dwn.n234 <- length(y$a7) + length(y$a6)
dwn.n1234 <- length(y$a6) 

# Create Figure 4C, TIMBR Down Scores
f4.21 <- as_ggplot(draw.quad.venn(area1 = dwn.area1, area2 = dwn.area2, area3 = dwn.area3, area4 = dwn.area4, 
        n12 = dwn.n12, n13 = dwn.n13, n14 = dwn.n14, n23 = dwn.n23, n24 = dwn.n24, n34 = dwn.n34, 
        n123 = dwn.n123, n124 = dwn.n124, n134 = dwn.n134, n234 = dwn.n234, n1234 = dwn.n1234,
        category = c("APAP",expression(paste("CCl"[4]*"")),"TCDD","TCE"),
        fill = c("Purple","grey","Green","orange"), alpha = 0.4, fontface = rep("plain",15), fontfamily = rep("sans",15),  
        cex = rep(4,15), cat.cex = c(3,3,3,3), cat.fontface = rep("plain",4), cat.fontfamily = rep("sans",4))) 

f4.22 <- as.ggplot(f4.21) + ggtitle("Metabolites predicted to decrease")+ 
  theme(plot.margin = unit(c(3,3,3,3), "lines")) +
  theme(plot.title = element_text(hjust = 0.5, size = 36, family = "Arial", face = "bold"))
f4.22


## -------------------------- Create TIMBR Bar charts -----------------------

# Read in annotation file
data <- read.xlsx("Rawls_Supplementary_data5.xlsx", sheet = 2) %>% 
  select(-KEGG_ID) %>% na.omit() %>% arrange(Drug,Class,Subclass)

# Filter names to make plots legible
data$Subclass[data$Subclass == "Carbohydrates and carbohydrate conjugates"] <- "Carbohydrates and conjugates"
data$Subclass[data$Subclass == "Carbohydrates and carbohydrate derivatives"] <- "Carbohydrates and derivatives"

# Filter out names to make plotting easier:
data <-  data %>% 
  filter(Subclass != "Pyrimidine ribonucleotides" & Subclass != "Phenols" & 
         Subclass != "Other non-metal oxides" & Subclass != "Estrane steroids")


# Combined
f4.4 <- data %>%  
  filter(Drug == "Combined") %>% 
  ggplot(aes(x = reorder(Subclass, Subclass, FUN = length), fill = Production_comparison)) + 
  geom_bar(stat = "count", color = "black") +
  geom_text(stat='count', aes(label = ..count..),position = position_stack(),
            hjust = 0.9,color = "Black",size = 12) +
  guides(fill=guide_legend(reverse = TRUE, title = "")) +
  scale_fill_manual(values = c("goldenrod1","yellow")) +
  theme_classic() + coord_flip() + ylim(0,17) +
  ggtitle("Overlap TIMBR Scores") + 
  ylab("") + xlab("Metabolite Classification") + 
  theme(legend.text = element_text(family = "Arial", size = 35)) +
  theme(legend.position = c(0.4,0.2)) +
  theme(plot.title = element_text(family = "Arial",face = "bold",size = 45)) +
  theme(axis.text = element_text(family = "Arial",size = 40)) + 
  theme(axis.title = element_text(family = "Arial",face = "bold",size = 40))

# APAP 
f4.5 <- data %>%  
  filter(Drug == "APAP") %>% 
  ggplot(aes(x = reorder(Subclass, Subclass, FUN = length), fill = Production_comparison)) + 
  geom_bar(stat = "count", color = "black") +
  geom_text(stat='count', aes(label = ..count..),position = position_stack(),
            hjust = 0.9,color = "White",size = 12) +
  guides(fill=guide_legend(reverse = TRUE, title = "")) +
  scale_fill_manual(values = c("darkviolet","orchid1")) +
  theme_classic() + coord_flip() + ylim(0,17) +
  ggtitle("APAP TIMBR Scores") + 
  ylab("") + xlab("Metabolite Classification") + 
  theme(legend.text = element_text(family = "Arial", size = 35)) +
  theme(legend.position = c(0.4,0.2)) +
  theme(plot.title = element_text(family = "Arial",face = "bold",size = 45)) +
  theme(axis.text = element_text(family = "Arial",size = 40)) + 
  theme(axis.title = element_text(family = "Arial",face = "bold",size = 40))

# CCl4 
f4.6 <- data %>%  
  filter(Drug == "CCl4") %>% 
  ggplot(aes(x = reorder(Subclass, Subclass, FUN = length), fill = Production_comparison)) + 
  geom_bar(stat = "count", color = "black") +
  geom_text(stat='count', aes(label = ..count..),position = position_stack(),
            hjust = 0.9,color = "Black",size = 12) +
  guides(fill=guide_legend(reverse = TRUE, title = "")) +
  scale_fill_manual(values = c("grey60","grey90")) +
  theme_classic() + coord_flip() + ylim(0,17) +
  ggtitle(expression(bold(paste("CCl"[4]*" TIMBR Scores")))) + 
  ylab("") + xlab("Metabolite Classification") + 
  theme(legend.text = element_text(family = "Arial", size = 35)) +
  theme(legend.position = c(0.4,0.2)) +
  theme(plot.title = element_text(family = "Arial",face = "bold",size = 45)) +
  theme(axis.text = element_text(family = "Arial",size = 40)) + 
  theme(axis.title = element_text(family = "Arial",face = "bold",size = 40))


# TCDD Increased
f4.7 <- data %>%  
  filter(Drug == "TCDD") %>% 
  ggplot(aes(x = reorder(Subclass, Subclass, FUN = length), fill = Production_comparison)) + 
  geom_bar(stat = "count", color = "black") +
  geom_text(stat='count', aes(label = ..count..),position = position_stack(),
            hjust = 0.9,color = "White",size = 12) +
  guides(fill=guide_legend(reverse = TRUE, title = "")) +
  scale_fill_manual(values = c("darkgreen","springgreen4")) +
  theme_classic() + coord_flip() + ylim(0,17) +
  ggtitle("TCDD TIMBR Scores") + 
  ylab("") + xlab("Metabolite Classification") + 
  theme(legend.text = element_text(family = "Arial", size = 35)) +
  theme(legend.position = c(0.4,0.2)) +
  theme(plot.title = element_text(family = "Arial",face = "bold",size = 45)) +
  theme(axis.text = element_text(family = "Arial",size = 40)) + 
  theme(axis.title = element_text(family = "Arial",face = "bold",size = 40))

# TCE Increased
f4.8 <- data %>%  
  filter(Drug == "TCE") %>% 
  ggplot(aes(x = reorder(Subclass, Subclass, FUN = length), fill = Production_comparison)) + 
  geom_bar(stat = "count", color = "black") +
  geom_text(stat='count', aes(label = ..count..),position = position_stack(),
            hjust = 0.9,color = "White",size = 12) +
  guides(fill=guide_legend(reverse = TRUE, title = "")) +
  scale_fill_manual(values = c("darkorange4","orange")) +
  theme_classic() + coord_flip() + ylim(0,17) +
  ggtitle("TCE TIMBR Scores") + 
  ylab("") + xlab("Metabolite Classification") + 
  theme(legend.text = element_text(family = "Arial", size = 35)) +
  theme(legend.position = c(0.4,0.2)) +
  theme(plot.title = element_text(family = "Arial",face = "bold",size = 45)) +
  theme(axis.text = element_text(family = "Arial",size = 40)) + 
  theme(axis.title = element_text(family = "Arial",face = "bold",size = 40))

left_figs <- ggarrange(f4.1,f4.12,f4.22,f4.4, nrow = 4, ncol = 1, labels = c("A","B","C","D"), font.label = list(size = 40))
right_figs <- ggarrange(f4.5,f4.6,f4.7,f4.8, nrow = 4, ncol = 1, labels = c("E","F","G","H"), font.label = list(size = 40))

figure4 <- ggarrange(left_figs,right_figs)

ggsave(figure4, filename = "Figure4_timbrsummaryfigure.png", width = 42, height = 45)


## -------------------------- Create Validation Heatmap -----------------------

timbr.compare <- timbr.compare %>% mutate(prod_score = rno)
timbr.compare[timbr.compare$rno_pctile < 75 & timbr.compare$rno_pctile > 25,"prod_score"] <- 0

# Read in analzyed metabolomics data
filtered_hep_metabolomics <- read.xlsx("Rawls_Supplementary_data3.xlsx", sheet = 6) %>% select(-X1)
control_hep_metabolomics <- read.xlsx("Rawls_Supplementary_data3.xlsx", sheet = 7) %>% select(-X1)

# Filter to only contain metabolites that are produced in the control case. 
avail_mets <- control_hep_metabolomics %>% filter(t_ctl1 > 0) %>% select(met)

new_filtered_mets <- filtered_hep_metabolomics[filtered_hep_metabolomics$met %in% avail_mets$met,]

metab.timbr.melt <- melt(new_filtered_mets, id.vars = "met", variable.name = "drug_id")

timbr.validation.1 <- left_join(metab.timbr.melt,timbr.compare, by = "met") %>% 
  filter(drug_id.x == "t_apap1") %>% filter(drug_id.y == "acetaminophen_t1")

timbr.validation.2 <- left_join(metab.timbr.melt,timbr.compare, by = "met") %>% 
  filter(drug_id.x == "t_ccl41") %>% filter(drug_id.y == "carbontetrachloride_t1")

timbr.validation.3 <- left_join(metab.timbr.melt,timbr.compare, by = "met") %>% 
  filter(drug_id.x == "t_tcdd1") %>% filter(drug_id.y == "TCDD_t1")

timbr.validation.4 <- left_join(metab.timbr.melt,timbr.compare, by = "met") %>% 
  filter(drug_id.x == "t_tce1") %>% filter(drug_id.y == "trichloroethylene_t1")

full_validation <- rbind(timbr.validation.1,timbr.validation.2,timbr.validation.3,timbr.validation.4) %>% 
  select(met,drug_id.y,value,prod_score) %>% 
  mutate(agreement = sign(value)*sign(prod_score)) %>% 
  mutate(stats = "increase.increase")

full_validation[full_validation$value == 0 & full_validation$prod_score == 0,"agreement"] <- 1

full_validation[full_validation$agreement == 1 & full_validation$prod_score < 0,"stats"] <- "decrease.decrease"
full_validation[full_validation$agreement == 1 & full_validation$prod_score == 0,"stats"] <- "no_change.no_change"
full_validation[full_validation$agreement == -1 & full_validation$prod_score > 0,"stats"] <- "decrease.increase"
full_validation[full_validation$agreement == -1 & full_validation$prod_score < 0,"stats"] <- "increase.decrease"
full_validation[full_validation$agreement == 0 & full_validation$prod_score > 0,"stats"] <- "no_change.increase"
full_validation[full_validation$agreement == 0 & full_validation$prod_score < 0,"stats"] <- "no_change.decrease"
full_validation[full_validation$agreement == 0 & full_validation$value > 0,"stats"] <- "increase.no_change"
full_validation[full_validation$agreement == 0 & full_validation$value < 0,"stats"] <- "decrease.no_change"


full_validation$stats <- factor(full_validation$stats, levels = c("no_change.no_change",
                                                                  "increase.increase","decrease.decrease","decrease.increase","increase.decrease",
                                                                  "increase.no_change","decrease.no_change","no_change.decrease","no_change.increase"))
full_validation <- full_validation %>% mutate(correct = agreement)
full_validation[full_validation$correct != 1,"correct"] <- 0

x.triangle.size = 0.45
drug_hindex <-  1:length(unique(full_validation$drug_id.y))
met_hindex <- 1:length(unique(full_validation$met))

position_xvector <- data.frame(p1 = -1, p2 = 1, p3 = 1, p4 = -1, p5 = -1, p6 = 1)  
position_yvector <- data.frame(p1 = -1, p2 = -1, p3 = 1, p4 = 1, p5 = -1, p6 = 1)

drug_pos <- data.frame(matrix(0,ncol = 6,nrow = length(unique(full_validation$drug_id.y))))
for (i in 1:length(drug_hindex)) {
  if (i == 1) {
    drug_pos[i,] <-  position_xvector
  } else {
    drug_pos[i,] <-  drug_hindex[i-1]  + (2*(i-1)) + position_xvector
  }
}
colnames(drug_pos) <- c("xmet.1","xmet.2","xmet.3","xtim.1","xtim.2","xtim.3")
met_pos <- data.frame(matrix(0,ncol = 6,nrow = length(unique(full_validation$met)))) # Change 16 for unique metabolites 
for (i in 1:length(met_hindex)) {
  if (i == 1) {
    met_pos[i,] <-  position_yvector
  } else {
    met_pos[i,] <-  met_hindex[i-1]  + (2*(i-1)) + position_yvector
  }
}
colnames(met_pos) <- c("ymet.1","ymet.2","ymet.3","ytim.1","ytim.2","ytim.3")
rownames(met_pos) <- unique(sort(full_validation$met))

full_combo.1 <- cbind(met_pos,drug_pos[1,])
full_combo.2 <- cbind(met_pos,drug_pos[2,])
full_combo.3 <- cbind(met_pos,drug_pos[3,])
full_combo.4 <- cbind(met_pos,drug_pos[4,])


full_combo <- rbind(full_combo.1,full_combo.2,full_combo.3,full_combo.4)
condition_vector <- c("acetaminophen_t1","carbontetrachloride_t1","TCDD_t1","trichloroethylene_t1") 
newdrugs <- data.frame(drug_id.x = rep(condition_vector,each = length(unique(full_validation$met)))) # change to number of unique metabolites
full_combo <- cbind(full_combo,newdrugs)
newmets <- unique(sort(full_validation$met, decreasing = TRUE))
met_vectors <- data.frame(met = rep(newmets,length(unique(full_validation$drug_id.y))))
full_combo <- cbind(full_combo,met_vectors) %>% rename(drug_id.y = "drug_id.x")

full_combination <- inner_join(full_combo,full_validation, by = "drug_id.y") %>% filter(met.x == met.y) %>% 
  rename(exper_val = "value") 

x1 <- cbind(full_combination$xmet.1,full_combination$xmet.2,full_combination$xmet.3)
x2 <- cbind(full_combination$xtim.1,full_combination$xtim.2,full_combination$xtim.3)
y1 <- cbind(full_combination$ymet.1,full_combination$ymet.2,full_combination$ymet.3)
y2 <- cbind(full_combination$ytim.1,full_combination$ytim.2,full_combination$ytim.3)

x1 <- data.frame(x1)
y1 <- data.frame(y1)
x2 <- data.frame(x2)
y2 <- data.frame(y2)


end_val <- length(full_combination$prod_score)
ids <- factor(c(1:end_val))
ids2 <- factor(c((end_val+1):(end_val+end_val)))

full_combination <- full_combination %>% rename(drug_id = "drug_id.y", met = "met.x")

values <- data.frame(full_combination,id = ids) %>% select(drug_id,met,prod_score,id) %>% rename(filler = "prod_score")
values2 <- data.frame(full_combination,id = ids2) %>% select(drug_id,met,exper_val,id) %>% rename(filler = "exper_val")

value_join <- rbind(values,values2)

positions <- data.frame(id = rep(ids,each = 3),x = c(t(x1)), y = c(t(y1)))
positions2 <- data.frame(id = rep(ids2,each = 3), x = c(t(x2)), y = c(t(y2)))

mega_values <- rbind(values,values2)
mega_position <- rbind(positions,positions2)

metaframe <- merge(value_join,mega_position, by = c("id"))


metaframe$filler[metaframe$filler > 0] <- 1
metaframe$filler[metaframe$filler < 0] <- -1
metaframe$filler <- as.factor(metaframe$filler)

# Heatmap figure
f5.1 <- ggplot(metaframe,aes(x = x, y = y)) +
  geom_point(alpha = 0) + 
  scale_x_continuous(breaks = seq(0,9,3), labels = c(acetaminophen_t1 = expression(paste("APAP 6hr")),
                                                     carbontetrachloride_t1 = expression(paste("CCl"[4]*" 6hr")),TCDD_t1 = expression(paste("TCDD 6hr")),
                                                     trichloroethylene_t1 = expression(paste("TCE 6hr"))), expand = c(0,0)) + 
  scale_y_continuous(limits = c(-1,59),breaks = seq(0,59,3), labels = sort(unique(full_validation$met), decreasing = TRUE),expand = c(0,0)) + 
  geom_polygon(color = "gray1",aes(fill = metaframe$filler, group = metaframe$id), alpha = 1) +
  theme_bw() + theme(panel.border = element_blank()) + theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 30))+ theme(axis.text = element_text(face = "bold", size = 30, family = "Arial")) +
  theme(legend.text = element_text(size = 30, family = "Arial")) + labs(x = "", y = "") +
  scale_fill_manual(labels = c("-1" = "Reduced Production","0" = "No change","1" = "Elevated Production"), values = c("-1" = "blue","0" = "white","1" = "red"))  

## -------------------------- Create Statistical Table (Figure 5B) -----------------------

f5.2 <- full_validation %>% 
  ggplot(aes(x = stats)) + geom_bar(width = 0.8, aes(fill = factor(correct, levels = c(1,0), ordered = TRUE))) + 
  geom_text(stat='count', aes(label = ..count..), hjust = -0.5, size = 10) + ylim(0,35) +
  scale_x_discrete(limits = c("decrease.increase","increase.decrease","decrease.no_change", 
    "increase.no_change","no_change.increase","no_change.decrease","increase.increase", 
    "decrease.decrease","no_change.no_change"), labels = c(no_change.no_change = expression(paste("E: NC\nP: NC")),
    decrease.decrease = expression(paste("E:  ↓  \nP: ↓ ")),increase.increase = expression(paste("E:  ↑  \nP: ↑ ")),
    no_change.decrease = expression(paste("E: NC\nP: ↓ ")),no_change.increase = expression(paste("E: NC\nP: ↑ ")),
    increase.no_change = expression(paste("E: ↑  \nP: NC")),decrease.no_change = expression(paste("E: ↓  \nP: NC")),
    increase.decrease = expression(paste("E:  ↑  \nP: ↓ ")),decrease.increase = expression(paste("E:  ↓  \nP: ↑ ")))) + 
  coord_flip(clip = "on") + theme_bw() + 
  scale_fill_manual(labels = c("1" = "Agreement","0" = "Disagreement"), 
                    values = c("1" = "grey50","0" = "grey90")) + 
  labs(x = "", y = "") + ggtitle("Summary of TIMBR Predictions") + 
  theme(legend.title = element_blank()) + 
  theme(legend.text = element_text(face = "bold",family = "Arial", size = 30)) + 
  theme(legend.position = c(0.8,0.1)) + 
  theme(plot.title = element_text(family = "Arial",face = "bold",size = 35)) +
  theme(axis.text = element_text(family = "Arial",face = "bold",size = 27)) 

timbr_figs <- ggarrange(f5.1,f5.2, ncol = 1, nrow = 2, heights = c(2,1), widths = c(1,2),labels = c("A","B"), font.label = list(size = 30))
timbr_figs
ggsave(timbr_figs,filename = "Figure5_timbrheatmapfigure.png",width = 25,height = 30)
