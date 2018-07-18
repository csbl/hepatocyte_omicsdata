# Rawls TIMBR Weights  
# This script is involved in generating biomarker predictions 
# based on gene expression changes using the algorithm,
# TIMBR (Transcriptionally-Inferred Biomarker Response).
#
# This script anaylzes the gene expression data generated from
# the rat hepatocytes treated with various 
# pharmaceutical compounds and environmental toxicants, and 
# creates weights for each reaction in the model
# based on the expression data.
#
# Code Adapted from blais_timbr_weight.R which can be found at 
# https://github.com/csbl/ratcon1


# Source helper files for analysis 
source("https://bioconductor.org/biocLite.R")
source("ncomm_helper.R")

# Load libraries
biocLite("S4Vectors")
library(limma)
library(Biobase)
library(reshape2)
library(readr)
library(ggplot2)
library(xlsx)
library(openxlsx)
library(dplyr)

# Expression Directory load 
path.expression.directory = "HepatocyteExpression/"

# Read in Reaction information from annotation table in Supplementary Data 3 of the manuscript
rxn.info.load = readWorkbook("ncomm_blais_supplementary_data3.xlsx",startRow = 2) %>% as.tbl

# Convert number of human and rat genes to numeric values instead of whatever it was previously
rxn.info = rxn.info.load %>% 
  mutate(n_hsa = as.numeric(as.character(n_gene_hsa)),
         n_rno = as.numeric(as.character(n_gene_rno)))

# Transform GPR rules into non Date/Time format
rxn.gene =  rxn.info %>% filter(enabled) %>% 
  select(rxn_id, hsa = gpr_hsa, rno = gpr_rno) %>%
  melt(c("rxn_id")) %>% ef_df %>% 
  # Excel annoyingly transforms GPR rules into date/time values
  mutate(value = gsub("193.5625","4645:30",value,fixed = T)) %>%
  mutate(value = gsub("[\\(\\)\\;\\:]+",";",value)) %>%
  ef_split_value(";") %>% ef_df %>% 
  dplyr::rename(organism_id = variable, gene_id = value) %>% 
  filter(nchar(gene_id) > 0) %>% distinct 

# all gene_id values should be integers
rxn.gene %>% dplyr::count(grepl("^[0-9]+$",gene_id))

# Create directory for expression data 
differential.expression.info = data_frame(
  efit_root = c("apap_6hr_TIMBR","apap_24hr_TIMBR",
                "ccl4_6hr_TIMBR","ccl4_24hr_TIMBR",
                "tcdd_6hr_TIMBR","tcdd_24hr_TIMBR",
                "tce_6hr_TIMBR","tce_24hr_TIMBR")) %>% 
  mutate(efit_file = paste0(path.expression.directory,efit_root,".csv")) %>%
  mutate(organism_id = "rno",
         dose_id = "d1",
         drug_id = gsub("_.*","",efit_root),
         time_id = gsub(".*_","",gsub("_TIMBR","",efit_root))) # Edits the previous and ending portions of the code to grab the time index

# Loads expression data from the directory created above - must use csv files not xlsx files 
differential.expression.load = differential.expression.info %>% 
  with(setNames(efit_file, efit_root)) %>% 
  lapply(read_csv,col_names=T) %>% 
  bind_rows(.id = "efit_root") %>% as.tbl %>% 
  mutate(gene_id = as.character(etz_gene)) %>%
  left_join(differential.expression.info)

# map to the expression dataset
differential.expression.load %>% 
  select(organism_id, gene_id) %>% distinct %>%
  dplyr::count(organism_id, gene_id %in% c(rxn.gene[["gene_id"]]))

# Establish FDR Cut off
fdr.cutoff = 0.1

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
  melt(c("rxn_id")) %>% ef_df %>% 
  mutate(value = gsub("\\-","",value)) %>% 
  mutate(value = gsub("PMID[\\:\\-]*",";PMID:",value)) %>% 
  ef_split_value(";") %>% ef_df %>% 
  filter(grepl("PMID|DOI|UNIPROT", value)) %>%
  filter(grepl("[0-9]+",value)) %>% distinct

# Set up TIMBR Weight calculation
timbr.weights.default = rxn.info %>% filter(enabled) %>% 
  select(rxn_id,rxn_class) %>% distinct %>% 
  left_join(rxn.pubmed %>% dplyr::count(rxn_id) %>% ungroup %>% dplyr::rename(pubmed_count = n)) %>% 
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
      dcast(rxn_id ~ variable, value.var = "value", fill = 0) %>% ef_df) %>% 
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
  melt(c("rxn_id","timbr_weight_default")) %>% ef_df %>%
  mutate(value = gsub("193.5625","4645:30",value, fixed = T))

timbr.expression.setup = metabolic.differential.expression %>% 
  mutate(limma_id = paste0(organism_id, "_", drug_id, "_", time_id, "_", dose_id),
         limma_ok = efit_ok) %>% filter(limma_ok) %>% 
  ef_df_slice("limma_id")

# Create TIMBR Weights
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
  dcast(organism_id + rxn_id + rxn_irreversible ~ timbr_id, value.var = "rxn_weight") %>% ef_df

# Write out timbr weights for analysis in MATLAB
write.table(rno.timbr.weights, file = paste0("Rawls_supplement_timbrweights.txt"), sep = "\t", quote = F, row.names = F)
