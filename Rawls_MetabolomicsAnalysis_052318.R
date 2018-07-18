# Source helper files for analysis 
source("ncomm_helper.R")
source("https://bioconductor.org/biocLite.R")

# Load necessary libraries for analysis and plotting 
library("reshape2")
library("dplyr")
library("ggplot2")
library("xlsx")
library("openxlsx")
library("pheatmap")
library("ggthemes")
library("matrixStats")
library("tibble")
library("gplots")
library("ggpubr")


# Load reaction info, TIMBR Scores, and Expression Data
rxn.info.load = readWorkbook("ncomm_blais_supplementary_data3.xlsx",startRow = 2) %>% as.tbl
timbr.predictions.load = readWorkbook("Rawls_supplement_productionscores.xlsx",startRow = 1) %>% as.tbl
expression.summary.load = readWorkbook("Rawls_supplement_expressionsummary.xlsx",startRow = 1) %>% as.tbl

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

# Rank TIMBR scores by their values for each drug  
timbr.compare = timbr.predictions %>% filter(drug_name != "gentamicin") %>% 
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


# Change production scores to 0, if they are between 1 and -1 
kris.timbr.compare.1 <- timbr.compare %>%
  filter(rno < 1 & rno > -1) %>%
  mutate(prod_score = 0)

kris.timbr.compare.2 <- timbr.compare %>%
  filter(rno > 1 | rno < -1) %>%
  mutate(prod_score = rno)

kris.timbr.predictions <- rbind(kris.timbr.compare.1,kris.timbr.compare.2) %>% arrange(drug_id,rxn_id)

# Read in metabolomics 
hep.raw.metabolomics <- readWorkbook("Rawls_supplement_metabolomics.xlsx")
annote1 <- data.frame(hep.raw.metabolomics[,c(1:2)]) %>% dplyr::rename(met = "met_name") %>% dplyr::rename(Annotation = "westcoast_metname")
hep.raw.metabolomics <- hep.raw.metabolomics[,c(2:65)] %>% dplyr::rename(met = "met_name")
mets <- hep.raw.metabolomics$met

hep.raw.lipids <- readWorkbook("Lipid_Data_Filtered_Hep_080317.xlsx")
hep.raw.lipids <- na.omit(hep.raw.lipids)
annote2 <- data.frame(hep.raw.lipids[,c(1:3)])
hep.raw.lipids <- hep.raw.lipids[,c(3,17:79)]
mets.lp <- hep.raw.lipids$met

hep.raw.amines <- readWorkbook("Biogenic_Amine_Data_Filtered_Hep_080317.xlsx")
hep.raw.amines <- na.omit(hep.raw.amines)
annote3 <- data.frame(hep.raw.amines[,c(1:2)])
hep.raw.amines <- hep.raw.amines[,c(2,9:71)]
mets.bg <- hep.raw.amines$met

## Filter The Data
# Remove excess apap and ccl4 data 
hep.raw.metabolomics <- hep.raw.metabolomics[,c(c(2:25),c(34:64))]
hep.raw.lipids <- hep.raw.lipids[,c(c(2:25),c(34:64))]
hep.raw.amines <- hep.raw.amines[,c(c(2:25),c(34:64))]

# Combine Data Frame 
hep.combined.metabolomics <- rbind(hep.raw.metabolomics,hep.raw.lipids,hep.raw.amines)

# Find Row Standard Deviations
hep.combined.metabolomics$stdevs <- rowSds(data.matrix(hep.combined.metabolomics[,c(2:55)]))

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
  t_apap1[i,1] <-  t.test(data.matrix(hep.combined.metabolomics[i,c(5:8)]),data.matrix(hep.combined.metabolomics[i,c(1:4)]), alternative = "two.sided", var.equal = FALSE)$statistic %>% unlist() %>% unname()  
  t_apap1[i,2] <-  t.test(data.matrix(hep.combined.metabolomics[i,c(5:8)]),data.matrix(hep.combined.metabolomics[i,c(1:4)]), alternative = "two.sided", var.equal = FALSE)$p.value %>% unlist() %>% unname()  
  
  t_apap2[i,1] <-  t.test(data.matrix(hep.combined.metabolomics[i,c(46:49)]),data.matrix(hep.combined.metabolomics[i,c(37:40)]), alternative = "two.sided", var.equal = FALSE)$statistic %>% unlist() %>% unname()
  t_apap2[i,2] <-  t.test(data.matrix(hep.combined.metabolomics[i,c(46:49)]),data.matrix(hep.combined.metabolomics[i,c(37:40)]), alternative = "two.sided", var.equal = FALSE)$p.value %>% unlist() %>% unname()
  
  t_ccl41[i,1] <-  t.test(data.matrix(hep.combined.metabolomics[i,c(9:12)]),data.matrix(hep.combined.metabolomics[i,c(1:4)]), alternative = "two.sided", var.equal = FALSE)$statistic %>% unlist() %>% unname() 
  t_ccl41[i,2] <-  t.test(data.matrix(hep.combined.metabolomics[i,c(9:12)]),data.matrix(hep.combined.metabolomics[i,c(1:4)]), alternative = "two.sided", var.equal = FALSE)$p.value %>% unlist() %>% unname() 
  
  t_ccl42[i,1] <-  t.test(data.matrix(hep.combined.metabolomics[i,c(50:55)]),data.matrix(hep.combined.metabolomics[i,c(41:45)]), alternative = "two.sided", var.equal = FALSE)$statistic %>% unlist() %>% unname()
  t_ccl42[i,2] <-  t.test(data.matrix(hep.combined.metabolomics[i,c(50:55)]),data.matrix(hep.combined.metabolomics[i,c(41:45)]), alternative = "two.sided", var.equal = FALSE)$p.value %>% unlist() %>% unname()
  
  t_tcdd1[i,1] <-  t.test(data.matrix(hep.combined.metabolomics[i,c(13:16)]),data.matrix(hep.combined.metabolomics[i,c(1:4)]), alternative = "two.sided", var.equal = FALSE)$statistic %>% unlist() %>% unname()
  t_tcdd1[i,2] <-  t.test(data.matrix(hep.combined.metabolomics[i,c(13:16)]),data.matrix(hep.combined.metabolomics[i,c(1:4)]), alternative = "two.sided", var.equal = FALSE)$p.value %>% unlist() %>% unname()
  
  t_tcdd2[i,1] <-  t.test(data.matrix(hep.combined.metabolomics[i,c(25:28)]),data.matrix(hep.combined.metabolomics[i,c(21:24)]), alternative = "two.sided", var.equal = FALSE)$statistic %>% unlist() %>% unname()
  t_tcdd2[i,2] <-  t.test(data.matrix(hep.combined.metabolomics[i,c(25:28)]),data.matrix(hep.combined.metabolomics[i,c(21:24)]), alternative = "two.sided", var.equal = FALSE)$p.value %>% unlist() %>% unname()
  
  t_tce1[i,1] <-  t.test(data.matrix(hep.combined.metabolomics[i,c(17:20)]),data.matrix(hep.combined.metabolomics[i,c(1:4)]), alternative = "two.sided", var.equal = FALSE)$statistic %>% unlist() %>% unname()  
  t_tce1[i,2] <-  t.test(data.matrix(hep.combined.metabolomics[i,c(17:20)]),data.matrix(hep.combined.metabolomics[i,c(1:4)]), alternative = "two.sided", var.equal = FALSE)$p.value %>% unlist() %>% unname()  
  
  t_tce2[i,1] <-  t.test(data.matrix(hep.combined.metabolomics[i,c(29:32)]),data.matrix(hep.combined.metabolomics[i,c(21:24)]), alternative = "two.sided", var.equal = FALSE)$statistic %>% unlist() %>% unname()  
  t_tce2[i,2] <-  t.test(data.matrix(hep.combined.metabolomics[i,c(29:32)]),data.matrix(hep.combined.metabolomics[i,c(21:24)]), alternative = "two.sided", var.equal = FALSE)$p.value %>% unlist() %>% unname()  
}


# Find Metabolites Produced or Consumed in the control condition 
t_ctl1 <- matrix(0,n,2)
t_ctl2 <- matrix(0,n,2)
t_actl2 <-  matrix(0,n,2)
t_cctl2 <-  matrix(0,n,2)

for (i in 1:n) {
  t_ctl1[i,1] <-  t.test(data.matrix(hep.combined.metabolomics[i,c(1:4)]),data.matrix(hep.combined.metabolomics[i,c(33:36)]), alternative = "two.sided", var.equal = FALSE)$statistic %>% unlist() %>% unname()  
  t_ctl1[i,2] <-  t.test(data.matrix(hep.combined.metabolomics[i,c(1:4)]),data.matrix(hep.combined.metabolomics[i,c(33:36)]), alternative = "two.sided", var.equal = FALSE)$p.value %>% unlist() %>% unname()  
  
  t_ctl2[i,1] <-  t.test(data.matrix(hep.combined.metabolomics[i,c(21:24)]),data.matrix(hep.combined.metabolomics[i,c(33:36)]), alternative = "two.sided", var.equal = FALSE)$statistic %>% unlist() %>% unname()  
  t_ctl2[i,2] <-  t.test(data.matrix(hep.combined.metabolomics[i,c(21:24)]),data.matrix(hep.combined.metabolomics[i,c(33:36)]), alternative = "two.sided", var.equal = FALSE)$p.value %>% unlist() %>% unname()  
  
  t_actl2[i,1] <-  t.test(data.matrix(hep.combined.metabolomics[i,c(37:40)]),data.matrix(hep.combined.metabolomics[i,c(33:36)]), alternative = "two.sided", var.equal = FALSE)$statistic %>% unlist() %>% unname()  
  t_actl2[i,2] <-  t.test(data.matrix(hep.combined.metabolomics[i,c(37:40)]),data.matrix(hep.combined.metabolomics[i,c(33:36)]), alternative = "two.sided", var.equal = FALSE)$p.value %>% unlist() %>% unname()  
  
  t_cctl2[i,1] <-  t.test(data.matrix(hep.combined.metabolomics[i,c(41:45)]),data.matrix(hep.combined.metabolomics[i,c(33:36)]), alternative = "two.sided", var.equal = FALSE)$statistic %>% unlist() %>% unname() 
  t_cctl2[i,2] <-  t.test(data.matrix(hep.combined.metabolomics[i,c(41:45)]),data.matrix(hep.combined.metabolomics[i,c(33:36)]), alternative = "two.sided", var.equal = FALSE)$p.value %>% unlist() %>% unname()  

}

# Filtering for significance.  
for (j in 1:n) {
  if (t_apap1[j,2] > 0.05) {
    t_apap1[j,1] = 0
  }
  
  if (t_apap2[j,2] > 0.05) {
    t_apap2[j,1] = 0
  }
  
  if (t_ccl41[j,2] > 0.05) {
    t_ccl41[j,1] = 0
  }
  
  if (t_ccl42[j,2] > 0.05) {
    t_ccl42[j,1] = 0
  }
  
  if (t_tcdd1[j,2] > 0.05) {
    t_tcdd1[j,1] = 0
  }
  
  if (t_tcdd2[j,2] > 0.05) {
    t_tcdd2[j,1] = 0
  }
  
  if (t_tce1[j,2] > 0.05) {
    t_tce1[j,1] = 0
  }
  
  if (t_tce2[j,2] > 0.05) {
    t_tce2[j,1] = 0
  }
  
  if (t_ctl1[j,2] > 0.05) {
    t_ctl1[j,1] = 0
  }
  
  if (t_ctl2[j,2] > 0.05) {
    t_ctl2[j,1] = 0
  }
  
  if (t_actl2[j,2] > 0.05) {
    t_actl2[j,1] = 0
  }
  
  if (t_cctl2[j,2] > 0.05) {
    t_cctl2[j,1] = 0
  }
}

# Reassign to new data frame
normal.hep.metabolomics <- data.frame(t_apap1[,1],t_apap2[,1],t_ccl41[,1],t_ccl42[,1],t_tcdd1[,1],t_tcdd2[,1],t_tce1[,1],t_tce2[,1])
colnames(normal.hep.metabolomics) <- c("t_apap1","t_apap2","t_ccl41","t_ccl42","t_tcdd1","t_tcdd2","t_tce1","t_tce2")
normal.hep.metabolomics$met <- c(mets,mets.lp,mets.bg)

compared_blank.hep.metabolomics <- data.frame(t_ctl1[,1],t_ctl2[,1],t_actl2[,1],t_cctl2[,1])
colnames(compared_blank.hep.metabolomics) <- c("t_ctl1","t_ctl2","t_actl2","t_cctl2")
compared_blank.hep.metabolomics$met <- c(mets,mets.lp,mets.bg)

combined.annotation <- rbind(annote1,annote2[,c(2:3)],annote3)

normal.hep.metabolomics %>% filter(t_apap1 != 0) %>% dplyr::select(met) %>% inner_join(combined.annotation, by = "met") %>% write.csv("apap1_mets.csv")
normal.hep.metabolomics %>% filter(t_apap2 != 0) %>% dplyr::select(met) %>% inner_join(combined.annotation, by = "met") %>% write.csv("apap2_mets.csv")
normal.hep.metabolomics %>% filter(t_ccl41 != 0) %>% dplyr::select(met) %>% inner_join(combined.annotation, by = "met") %>% write.csv("ccl41_mets.csv")
normal.hep.metabolomics %>% filter(t_ccl42 != 0) %>% dplyr::select(met) %>% inner_join(combined.annotation, by = "met") %>% write.csv("ccl42_mets.csv")
normal.hep.metabolomics %>% filter(t_tcdd1 != 0) %>% dplyr::select(met) %>% inner_join(combined.annotation, by = "met") %>% write.csv("tcdd1_mets.csv")
normal.hep.metabolomics %>% filter(t_tcdd2 != 0) %>% dplyr::select(met) %>% inner_join(combined.annotation, by = "met") %>% write.csv("tcdd2_mets.csv")
normal.hep.metabolomics %>% filter(t_tce1 != 0) %>% dplyr::select(met) %>% inner_join(combined.annotation, by = "met") %>% write.csv("tce1_mets.csv")
normal.hep.metabolomics %>% filter(t_tce2 != 0) %>% dplyr::select(met) %>% inner_join(combined.annotation, by = "met") %>% write.csv("tce2_mets.csv")


# Write out normalized metabolomics (Added 02/06/18)
write.xlsx2(normal.hep.metabolomics, "C:/Users/kr2up/Desktop/020618_normalizedmetabolomics.xlsx",row.names = FALSE)

# Melt the data frame 
normalized.melt.metabolomics <- melt(normal.hep.metabolomics, id.vars = "met")

normalized.melt.metabolomics <- normalized.melt.metabolomics %>% 
  mutate(time_id = gsub("\\D","",gsub("[4]","",normalized.melt.metabolomics$variable))) %>% 
  mutate(drug_id = gsub("\\d","",gsub("t_","",normalized.melt.metabolomics$variable))) %>% 
  select(met,value,time_id,drug_id)

normalized.melt.metabolomics.2 <- normalized.melt.metabolomics %>% 
  filter(drug_id == "ccl") %>% mutate(drug_id = paste0(drug_id,"4"))

normalized.melt.metabolomics.1 <- normalized.melt.metabolomics %>% filter(drug_id == "apap" | drug_id == "tcdd" | drug_id == "tce")

normal.melt.metabolomics <- rbind(normalized.melt.metabolomics.1,normalized.melt.metabolomics.2) %>% arrange(drug_id)

filt.melt.data <- normal.melt.metabolomics %>% 
  filter(drug_id == "apap" | drug_id == "ccl4" | drug_id == "tcdd" | 
           drug_id == "tce") %>% filter(time_id == "1")

filt.melt.data$drug_id <- as.factor(filt.melt.data$drug_id)

# Find Standard Deviation 
normal.bin.metabolomics <- normal.hep.metabolomics
normal.bin.metabolomics[normal.bin.metabolomics > 5] <- 5
normal.bin.metabolomics[normal.bin.metabolomics < -5] <- -5
normal.bin.metabolomics$met <- normal.hep.metabolomics$met

trans_bin.metabolomics <- t(normal.bin.metabolomics)
conds <- c("apap1","apap2","ccl41","ccl42",
           "tcdd1", "tcdd2","tce1", "tce2")

trans_metabolomics <- trans_bin.metabolomics[c(1:8),]
rownames(trans_metabolomics) <- conds
class(trans_metabolomics) <- "numeric"
normalized.met.hmap <- data.matrix(trans_metabolomics)

# Compare timbr predictions with metabolomics data
metabolomics.timbr.compare <-  normal.bin.metabolomics 

metab.timbr.melt <- melt(metabolomics.timbr.compare, id.vars = "met") 
metab.timbr.melt$drug_id <- metab.timbr.melt$variable

# Match metabolite from TIMBR scores with metabolites from metabolomics data
timbr.validation.1 <- left_join(metab.timbr.melt,kris.timbr.predictions, by = "met") %>% 
  filter(drug_id.x == "t_apap1") %>% filter(drug_id.y == "acetaminophen_t1")

timbr.validation.2 <- left_join(metab.timbr.melt,kris.timbr.predictions, by = "met") %>% 
  filter(drug_id.x == "t_apap2") %>% filter(drug_id.y == "acetaminophen_t2")

timbr.validation.3 <- left_join(metab.timbr.melt,kris.timbr.predictions, by = "met") %>%
  filter(drug_id.x == "t_ccl41") %>% filter(drug_id.y == "carbontetrachloride_t1")

timbr.validation.4 <- left_join(metab.timbr.melt,kris.timbr.predictions, by = "met") %>% 
  filter(drug_id.x == "t_ccl42") %>% filter(drug_id.y == "carbontetrachloride_t2")

timbr.validation.5 <- left_join(metab.timbr.melt,kris.timbr.predictions, by = "met") %>% 
  filter(drug_id.x == "t_tcdd1") %>% filter(drug_id.y == "TCDD_t1")

timbr.validation.6 <- left_join(metab.timbr.melt,kris.timbr.predictions, by = "met") %>%
  filter(drug_id.x == "t_tcdd2") %>% filter(drug_id.y == "TCDD_t2")

timbr.validation.7 <- left_join(metab.timbr.melt,kris.timbr.predictions, by = "met") %>% 
  filter(drug_id.x == "t_tce1") %>% filter(drug_id.y == "trichloroethylene_t1")

timbr.validation.8 <- left_join(metab.timbr.melt,kris.timbr.predictions, by = "met") %>% 
  filter(drug_id.x == "t_tce2") %>% filter(drug_id.y == "trichloroethylene_t2")


full.validation <- rbind(timbr.validation.1,timbr.validation.2,timbr.validation.3,
                         timbr.validation.4,timbr.validation.5,timbr.validation.6,
                         timbr.validation.7,timbr.validation.8) 

full.validation <- full.validation %>% select(met,drug_id.y,value,prod_score) %>% 
  mutate(agreement = value*prod_score) 

# Separate each prediction match in data frames
full.validation_1 <- full.validation %>% 
  filter(sign(value) == sign(prod_score)) %>% mutate(agreement = 1)

full.validation_1.1 <- full.validation_1 %>% 
  filter(sign(value) > 0) %>% mutate(stats = "increase.increase")

full.validation_1.2 <- full.validation_1 %>% 
  filter(sign(value) == 0) %>% mutate(stats = "no_change.no_change")

full.validation_1.3 <- full.validation_1 %>% 
  filter(sign(value) < 0) %>% mutate(stats = "decrease.decrease")

full.validation_1 <- rbind(full.validation_1.1,full.validation_1.2,full.validation_1.3)

full.validation_2 <- full.validation %>% 
  filter((value == 0 & prod_score != 0) |(value != 0 & prod_score == 0) ) %>% mutate(agreement = 0)

full.validation_2.1 <- full.validation_2 %>% 
  filter(sign(value) > 0 & prod_score == 0) %>% mutate(stats = "increase.no_change")

full.validation_2.2 <- full.validation_2 %>% 
  filter(sign(value) < 0 & prod_score == 0) %>% mutate(stats = "decrease.no_change")

full.validation_2.3 <- full.validation_2 %>% 
  filter(sign(value) == 0 & prod_score > 0) %>% mutate(stats = "no_change.increase")

full.validation_2.4 <- full.validation_2 %>% 
  filter(sign(value) == 0 & prod_score < 0) %>% mutate(stats = "no_change.decrease")

full.validation_2 <- rbind(full.validation_2.1,full.validation_2.2,full.validation_2.3,full.validation_2.4)

full.validation_3 <- full.validation %>% 
  filter((sign(value) > 0 & sign(prod_score) < 0) |(sign(value) < 0 & sign(prod_score) > 0) ) %>% mutate(agreement = -1)

full.validation_3.1 <- full.validation_3 %>% 
  filter(sign(value) > 0) %>% mutate(stats = "increase.decrease")

full.validation_3.2 <- full.validation_3 %>% 
  filter(sign(value) < 0) %>% mutate(stats = "decrease.increase")

full.validation_3 <- rbind(full.validation_3.1,full.validation_3.2)

new_vald <- rbind(full.validation_1,full.validation_2,full.validation_3)
new_vald <- new_vald %>% mutate(correct = (agreement >0))

# Use GG Plot to make TIMBR validation heatmap 
x.triangle.size = 0.45

# Create position vectors 
drug_hindex <-  1:length(unique(new_vald$drug_id.y))
met_hindex <- 1:length(unique(new_vald$met))

#d1 <- data_frame(x_met = c(-1, 1, 1), y_met = -c( 1, 1,-1),x_tim = c(-1,-1, 1), y_tim = -c(-1, 1,-1))
position_xvector <- data.frame(p1 = -1, p2 = 1, p3 = 1, p4 = -1, p5 = -1, p6 = 1)  
position_yvector <- data.frame(p1 = -1, p2 = -1, p3 = 1, p4 = 1, p5 = -1, p6 = 1)

drug_pos <- data.frame(matrix(0,ncol = 6,nrow = 8))
for (i in 1:length(drug_hindex)) {
  if (i == 1) {
    drug_pos[i,] <-  position_xvector
  } else {
    drug_pos[i,] <-  drug_hindex[i-1]  + (2*(i-1)) + position_xvector
  }
}

colnames(drug_pos) <- c("xmet.1","xmet.2","xmet.3","xtim.1","xtim.2","xtim.3")
met_pos <- data.frame(matrix(0,ncol = 6,nrow = 39)) # Change 16 for unique metabolites 
for (i in 1:length(met_hindex)) {
  if (i == 1) {
    met_pos[i,] <-  position_yvector
  } else {
    met_pos[i,] <-  met_hindex[i-1]  + (2*(i-1)) + position_yvector
  }
}
colnames(met_pos) <- c("ymet.1","ymet.2","ymet.3","ytim.1","ytim.2","ytim.3")
rownames(met_pos) <- unique(sort(new_vald$met, decreasing = TRUE))

full_combo.1 <- cbind(met_pos,drug_pos[1,])
full_combo.2 <- cbind(met_pos,drug_pos[2,])
full_combo.3 <- cbind(met_pos,drug_pos[3,])
full_combo.4 <- cbind(met_pos,drug_pos[4,])
full_combo.5 <- cbind(met_pos,drug_pos[5,])
full_combo.6 <- cbind(met_pos,drug_pos[6,])
full_combo.7 <- cbind(met_pos,drug_pos[7,])
full_combo.8 <- cbind(met_pos,drug_pos[8,])


full_combo <- rbind(full_combo.1,full_combo.2,full_combo.3,full_combo.4,full_combo.5,full_combo.6,full_combo.7,full_combo.8)
condition_vector <- c("acetaminophen_t1","acetaminophen_t2","carbontetrachloride_t1","carbontetrachloride_t2",
                      "TCDD_t1","TCDD_t2","trichloroethylene_t1","trichloroethylene_t2") 
newdrugs <- data.frame(drug_id.x = rep(condition_vector,each = 39)) # change to number of unique metabolites
full_combo <- cbind(full_combo,newdrugs)
newmets <- unique(sort(new_vald$met, decreasing = TRUE))
met_vectors <- data.frame(met = rep(newmets,8))
full_combo <- cbind(full_combo,met_vectors) %>% dplyr::rename(drug_id.y = "drug_id.x")

full_combination <- inner_join(full_combo,new_vald, by = "drug_id.y") %>% filter(met.x == met.y) %>% 
  dplyr::rename(exper_val = "value") 

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

full_combination <- full_combination %>% dplyr::rename(drug_id = "drug_id.y", met = "met.x")

values <- data.frame(full_combination,id = ids) %>% dplyr::select(drug_id,met,prod_score,id) %>% dplyr::rename(filler = "prod_score")
values2 <- data.frame(full_combination,id = ids2) %>% dplyr::select(drug_id,met,exper_val,id) %>% dplyr::rename(filler = "exper_val")

value_join <- rbind(values,values2)

positions <- data.frame(id = rep(ids,each = 3),x = c(t(x1)), y = c(t(y1)))
positions2 <- data.frame(id = rep(ids2,each = 3), x = c(t(x2)), y = c(t(y2)))

mega_values <- rbind(values,values2)
mega_position <- rbind(positions,positions2)

metaframe <- merge(value_join,mega_position, by = c("id"))


metaframe$filler[metaframe$filler > 1] <- 1
metaframe$filler[metaframe$filler < -1] <- -1
metaframe$filler <- as.factor(metaframe$filler)
# Heatmap figure

f5.1 <- ggplot(metaframe,aes(x = x, y = y)) +
  geom_point(alpha = 0) + 
  scale_x_continuous(breaks = seq(0,21,3), labels = c(acetaminophen_t1 = expression(paste("APAP 6hr")), acetaminophen_t2 = expression(paste("APAP 24hr")),
                                                      carbontetrachloride_t1 = expression(paste("CCl"[4]*" 6hr")), carbontetrachloride_t2 = expression(paste("CCl"[4]*" 24hr")), 
                                                      TCDD_t1 = expression(paste("TCDD 6hr")), TCDD_t2 = expression(paste("TCDD 24hr")),
                                                      trichloroethylene_t1 = expression(paste("TCE 6hr")),trichloroethylene_t2 = expression(paste("TCE 24hr"))), expand = c(0,0)) + 
  scale_y_continuous(limits = c(-1,119),breaks = seq(0,116,3), labels = sort(unique(new_vald$met), decreasing = TRUE),expand = c(0,0)) + 
  geom_polygon(color = "gray1",aes(fill = metaframe$filler, group = metaframe$id), alpha = 1) +
  theme_bw() +
  theme(panel.border = element_blank()) + 
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14))+ theme(axis.text = element_text(face = "bold", size = 10, family = "Arial")) +
  theme(legend.text = element_text(size = 12, family = "Arial")) + labs(x = "", y = "") +
  scale_fill_manual(labels = c("-1" = "Reduced Production","0" = "No change","1" = "Elevated Production"), values = c("-1" = "blue","0" = "white","1" = "red"))



# Use ggplot to create bar charts on statistics 
f5.2 <- new_vald %>% ggplot(aes(x = stats)) +
  geom_bar(width = 0.5, aes(fill = correct)) + 
  geom_text(stat='count', aes(label = ..count..), hjust = -0.5) + ylim(0,205) +
  scale_x_discrete(limits = c("decrease.increase","increase.decrease","decrease.no_change", 
                              "increase.no_change","no_change.increase","no_change.decrease",
                              "increase.increase", "decrease.decrease","no_change.no_change"), 
                   labels = c(no_change.no_change = expression(paste("E: NC\nP: NC")),
                              decrease.decrease = expression(paste("E:  ↓  \nP: ↓ ")),
                              increase.increase = expression(paste("E:  ↑  \nP: ↑ ")),
                              no_change.decrease = expression(paste("E: NC\nP: ↓ ")),
                              no_change.increase = expression(paste("E: NC\nP: ↑ ")),
                              increase.no_change = expression(paste("E: ↑  \nP: NC")),
                              decrease.no_change = expression(paste("E: ↓  \nP: NC")),
                              increase.decrease = expression(paste("E:  ↑  \nP: ↓ ")),
                              decrease.increase = expression(paste("E:  ↓  \nP: ↑ ")))) + 
  coord_flip() + theme_bw() + 
  scale_fill_manual(limits=c("TRUE","FALSE"),labels = c("TRUE" = "Agreement","FALSE" = "Disagreement"), 
                    values = c("TRUE" = "darkgreen","FALSE" = "darkorchid")) + 
  labs(x = "", y = "") + 
  ggtitle("Summary of TIMBR Predictions") + 
  theme(legend.title = element_blank()) + 
  theme(legend.text = element_text(face = "bold",family = "Arial", size = 12)) + 
  theme(legend.position = c(0.8,0.2)) + 
  theme(axis.text.y = element_text(family = "Arial",face = "bold",size = 10))

ggarrange(f5.1,f5.2, labels = c("A","B"), ncol = 1, nrow = 2, heights = c(6.5,3.5))


