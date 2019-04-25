# Kristopher Rawls
# 04/25/19
# This code was adapted from blais_timbr_weight.R available at github.com/csbl/ratcon1 published at
# http://www.nature.com/articles/ncomms14250
# This code takes in gene expression data and list of Reactions from iMAT model and creates TIMBR Reaction Weights
# Uses Supplementary Data 1 and 6 

#options(java.parameters = "-Xmx4g")
options(stringsAsFactors = FALSE)

# Source needed files from Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("S4Vectors")
biocLite("limma")

# Load necessary libraries
library(Biobase)
library(reshape2)
library(readr)
library(dplyr)
library(ggplot2)
library(openxlsx)

# Source helper file from Blais et al located at github.com/csbl/ratcon1
source("ncomm_helper.R")

# Load Matlab Directory for output
path.matlab.directory = "C:/Users/kr2up/Documents/MATLAB/"

# Read in Reaction information from annotation table in Supplementary Data 3 of the manuscript
rxn.info.load = readWorkbook("C:/Users/kr2up/Documents/ncomm_blais_submission_113016/ncomm_blais_supplementary_data3.xlsx",
                             startRow = 2) %>% as.tbl


# Convert number of human and rat genes to numeric values instead of whatever it was previously
rxn.info = rxn.info.load %>% 
  mutate(n_hsa = as.numeric(as.character(n_gene_hsa)),
         n_rno = as.numeric(as.character(n_gene_rno)))

# Pick a pruining method if desired
subset_rxns <- readWorkbook("C:/Users/kr2up/Documents/MATLAB/Rawls_Supplementary_data6.xlsx",
       startRow = 1, colNames = FALSE) %>% as.tbl() %>% dplyr::rename("rxn_id" = X1)

rxn.info <- inner_join(subset_rxns,rxn.info, by = "rxn_id")

# Transform GPR rules into non Date/Time format
rxn.gene =  rxn.info %>% filter(enabled) %>% 
  select(rxn_id, hsa = gpr_hsa, rno = gpr_rno) %>%
  reshape2::melt(c("rxn_id")) %>% ef_df %>% 
  # Excel transforms GPR rules into date/time values unfortunately
  mutate(value = gsub("193.5625","4645:30",value,fixed = T)) %>%
  mutate(value = gsub("[\\(\\)\\;\\:]+",";",value)) %>%
  ef_split_value(";") %>% ef_df %>% 
  rename(organism_id = variable, gene_id = value) %>% 
  filter(nchar(gene_id) > 0) %>% distinct 

# all gene_id values should be integers
rxn.gene %>% count(grepl("^[0-9]+$",gene_id))

# Load Expression Data
apap <- read.xlsx("C:/Users/kr2up/Documents/MATLAB/Rawls_Supplementary_data1.xlsx",sheet = 2) %>% 
  as.tbl() %>% dplyr::select(etz_gene,log2FoldChange,pvalue,padj,fdr,baseMean) %>% 
  mutate(organism_id = "rno",dose_id = "d1", drug_id = "apap",time_id = "6hr") %>% 
  dplyr::rename("logfc" = log2FoldChange) %>% dplyr::rename("pval" = pvalue, "ave" = baseMean) %>%
  mutate(gene_id = as.character(etz_gene)) %>%
  mutate(efit_root = paste(drug_id,time_id, sep = "_"))

ccl4 <- read.xlsx("C:/Users/kr2up/Documents/MATLAB/Rawls_Supplementary_data1.xlsx",sheet = 3) %>% 
  as.tbl() %>% dplyr::select(etz_gene,log2FoldChange,pvalue,padj,fdr,baseMean) %>% 
  mutate(organism_id = "rno",dose_id = "d1", drug_id = "ccl4",time_id = "6hr") %>% 
  dplyr::rename("logfc" = log2FoldChange) %>% dplyr::rename("pval" = pvalue, "ave" = baseMean) %>%
  mutate(gene_id = as.character(etz_gene)) %>%
  mutate(efit_root = paste(drug_id,time_id, sep = "_"))

tcdd <- read.xlsx("C:/Users/kr2up/Documents/MATLAB/Rawls_Supplementary_data1.xlsx",sheet = 4) %>% 
  as.tbl() %>% dplyr::select(etz_gene,log2FoldChange,pvalue,padj,fdr,baseMean) %>% 
  mutate(organism_id = "rno",dose_id = "d1", drug_id = "tcdd",time_id = "6hr") %>% 
  dplyr::rename("logfc" = log2FoldChange) %>% dplyr::rename("pval" = pvalue, "ave" = baseMean) %>%
  mutate(gene_id = as.character(etz_gene)) %>%
  mutate(efit_root = paste(drug_id,time_id, sep = "_"))


tce <- read.xlsx("C:/Users/kr2up/Documents/MATLAB/Rawls_Supplementary_data1.xlsx",sheet = 5) %>% 
  as.tbl() %>% dplyr::select(etz_gene,log2FoldChange,pvalue,padj,fdr,baseMean) %>% 
  mutate(organism_id = "rno",dose_id = "d1", drug_id = "tce",time_id = "6hr") %>% 
  dplyr::rename("logfc" = log2FoldChange) %>% dplyr::rename("pval" = pvalue, "ave" = baseMean) %>%
  mutate(gene_id = as.character(etz_gene)) %>%
  mutate(efit_root = paste(drug_id,time_id, sep = "_"))


differential.expression.load <- rbind(apap,ccl4,tcdd,tce)

# map to these two expression datasets
# 2121 rat genes and 0 human genes
differential.expression.load %>% 
  select(organism_id, gene_id) %>% distinct %>%
  count(organism_id, gene_id %in% c(rxn.gene[["gene_id"]]))

# Establish FDR Cut off
fdr.cutoff = 0.1

# Removed efit_id from original code because gene expression data isn't needed
metabolic.differential.expression = differential.expression.load %>%
  semi_join(rxn.gene %>% select(gene_id, organism_id) %>% distinct) %>% 
  mutate(significant = fdr < fdr.cutoff,
         direction = ifelse(significant, sign(logfc),0)) %>%
  group_by(efit_root, organism_id, drug_id, time_id, dose_id) %>%
  mutate(metabolic_gene_count = n(),
         n_up = sum(direction > 0), 
         n_dn = sum(direction < 0), 
         n_significant = sum(direction != 0)) %>% ungroup %>%
  mutate(pct_significant = 100 * n_significant / metabolic_gene_count) %>%
  mutate(efit_ok = pct_significant > 1)

# Filter to look at summary of up and down genes
metabolic.differential.expression %>% 
  select(efit_root, n_up, n_dn, n_significant, pct_significant, efit_ok) %>% 
  distinct %>% data.frame

rxn.pubmed = rxn.info %>% filter(enabled) %>% select(rxn_id, pubmed_id) %>% 
  reshape2::melt(c("rxn_id")) %>% ef_df %>% 
  mutate(value = gsub("\\-","",value)) %>% 
  mutate(value = gsub("PMID[\\:\\-]*",";PMID:",value)) %>% 
  ef_split_value(";") %>% ef_df %>% 
  filter(grepl("PMID|DOI|UNIPROT", value)) %>%
  filter(grepl("[0-9]+",value)) %>% distinct


# Editing TIMBR Code to remove hsa from all examples since we only have rat expression data
timbr.weights.default = rxn.info %>% filter(enabled) %>% 
  select(rxn_id,rxn_class) %>% distinct %>% 
  left_join(rxn.pubmed %>% count(rxn_id) %>% ungroup %>% rename(pubmed_count = n)) %>% 
  mutate(pubmed_count = ifelse(!is.na(pubmed_count), pubmed_count, 0)) %>% 
  left_join(bind_rows(list(
    rxn.gene %>%
      mutate(variable = paste0(organism_id, "_count_all")) %>% 
      group_by(rxn_id, variable) %>% 
      summarize(value = length(unique(setdiff(gene_id,"0")))) %>% ungroup,
    rxn.gene %>%
      semi_join(metabolic.differential.expression %>% select(organism_id, gene_id) %>% distinct) %>% 
      mutate(variable = paste0(organism_id, "_count")) %>% 
      group_by(rxn_id, variable) %>% 
      summarize(value = length(unique(setdiff(gene_id,"0")))) %>% ungroup)) %>%
      reshape2::dcast(rxn_id ~ variable, value.var = "value", fill = 0) %>% ef_df) %>% 
  mutate(rno_count = ifelse(!is.na(rno_count), rno_count, 0),
         rno_count_all = ifelse(!is.na(rno_count_all), rno_count_all, 0)) %>% 
  mutate(rxn_enzymatic =  rno_count > 0,
         rxn_referenced = pubmed_count > 0 ) %>%
  mutate(timbr_weight_enzymatic = ifelse(rxn_enzymatic,1,2),
         timbr_weight_referenced = ifelse(rxn_referenced,1,2),
         timbr_weight_class = ifelse(rxn_class == "boundary",2, ifelse(rxn_class == "transport",2,1)),
         timbr_weight_default = timbr_weight_enzymatic * timbr_weight_referenced * timbr_weight_class)



timbr.rxn.setup = timbr.weights.default %>% select(rxn_id, timbr_weight_default) %>%
  left_join(rxn.info %>% select(rxn_id, rno = gpr_rno, hsa = gpr_hsa)) %>%
  reshape2::melt(c("rxn_id","timbr_weight_default")) %>% ef_df %>%
  mutate(value = gsub("193.5625","4645:30",value, fixed = T))

timbr.expression.setup = metabolic.differential.expression %>% 
  mutate(limma_id = paste0(organism_id, "_", drug_id, "_", time_id, "_", dose_id),
         limma_ok = efit_ok) %>% filter(limma_ok) %>% 
  ef_df_slice("limma_id")


timbr.weights.list = timbr.expression.setup %>% lapply(ef_timbr_weights,timbr.rxn.setup,0,0)
timbr.weights = timbr.weights.list %>% bind_rows

timbr.weights %>% with(qplot(timbr_weight_ctl, timbr_weight_trt))

rxn.irreversible = bind_rows(list(
  timbr.weights %>% select(rxn_id) %>% distinct %>% mutate(irxn_id = paste0(rxn_id,"_f")),
  timbr.weights %>% select(rxn_id) %>% distinct %>% mutate(irxn_id = paste0(rxn_id,"_r"))))

timbr.weights.irreversible = bind_rows(list(
  timbr.weights %>% mutate(rxn_irreversible = paste0(rxn_id,"_f")),
  timbr.weights %>% mutate(rxn_irreversible = paste0(rxn_id,"_r"))))

rno.timbr.weights = bind_rows(list(
  timbr.weights.irreversible  %>%
    filter(organism_id == "rno") %>% 
    mutate(timbr_id = paste0(limma_id, "_ctl")) %>%
    select(timbr_id, organism_id, rxn_id, rxn_irreversible, rxn_weight = timbr_weight_ctl),
  timbr.weights.irreversible  %>%
    filter(organism_id == "rno") %>% 
    mutate(timbr_id = paste0(limma_id, "_trt")) %>%
    select(timbr_id, organism_id, rxn_id, rxn_irreversible, rxn_weight = timbr_weight_trt))) %>%
  reshape2::dcast(organism_id + rxn_id + rxn_irreversible ~ timbr_id, value.var = "rxn_weight") %>% ef_df

# Change file name to write to
write.table(rno.timbr.weights, file = paste0(path.matlab.directory,"kdr_hep_expression_042519_imat_timbr_weights_rno.txt"), 
            sep = "\t", quote = F, row.names = F)


