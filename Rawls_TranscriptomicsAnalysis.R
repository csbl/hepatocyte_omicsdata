# Rawls Transcritpomics Analysis
# Imports Results from Kallisto alignment
# Output DEGs
# Uses Tximport + DESeq2 to analyze results

# Load Necessary Packages 
library("readr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("ggdendro")
library("RColorBrewer")
library("gplots")
library("SummarizedExperiment")
library("genefilter")
library("openxlsx")
library("xlsx")
library("tibble")
library("grid")
library("ggpubr")
library("gridExtra")
library("cowplot")

#Source files from bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("tximport") 
biocLite("biomaRt")
biocLite("DESeq2")
biocLite('vsn')

# Read Results from Kallisto 
# Create file names to get abundance files from Kallisto 
dir <- "/Users/kdrawls/Documents/CSBL/Experimental_Results/RPTEC_transcriptomics"

# List files in each directory
list.files(dir)

# Read all samples in the file folder
samples <- read.table(file.path("/Users/kdrawls/Documents/CSBL/Experimental_Results","papin_RPTEC_genewiz_rnaseq_info.txt"), header = TRUE)

# Specify the individual files in a directory
files <- file.path(dir,samples$run_accession,"abundance.tsv")
files2 <- file.path(dir2,samples2$run_accession,"abundance.tsv")
files3 <- file.path(dir3,samples3$run_accession,"abundance.tsv")
files6 <- file.path(dir6,samples6$run_accession,"abundance.tsv")

# Label the files as sample numbers
names(files) <- paste0("sample",1:dim(samples)[1])
names(files2) <- paste0("sample",1:dim(samples2)[1])
names(files3) <- paste0("sample",1:dim(samples3)[1])
names(files6) <- paste0("sample",1:dim(samples6)[1])

# Gather Transcript and Gene IDs from biomaRt
mart  <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "rnorvegicus_gene_ensembl")
t2g6 <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id","external_gene_name", "entrezgene", "transcript_version","description","phenotype_description"), mart = mart)
t2g6 <- within(t2g6,target_id <- paste(ensembl_transcript_id, transcript_version,sep = "."))
t2g6 <- dplyr::rename(t2g6, transcript_id = ensembl_transcript_id,ens_gene = ensembl_gene_id, ext_gene = external_gene_name, etz_gene = entrezgene, descript = description, phenotype = phenotype_description)
tx2gene <- t2g6[,c(8,2)]
dds2model <- t2g6[,c(2,3,4)]
dds2model <- dds2model[!duplicated(dds2model),]

# Manually filter out duplicate or misrepresented genes for metabolic model purposes
dds2model <- dds2model %>% 
  filter(ens_gene != "ENSRNOG00000006589" & etz_gene != 81683) %>% 
  filter(ens_gene != "ENSRNOG00000007758" & etz_gene != 286964) %>% 
  filter(ens_gene != "ENSRNOG00000015267" & etz_gene != 108348083) %>%
  filter(ens_gene != "ENSRNOG00000018736" & etz_gene != 100911615) %>%
  filter(ens_gene != "ENSRNOG00000020025" & etz_gene != 108348052) %>%
  filter(ens_gene != "ENSRNOG00000046053" & etz_gene != 680248) %>% 
  filter(ens_gene != "ENSRNOG00000048573" & etz_gene != 367747) %>% 
  filter(ens_gene != "ENSRNOG00000048668" & etz_gene != 100911305) %>%
  filter(ens_gene != "ENSRNOG00000048949" & etz_gene != 171402) %>%
  filter(ens_gene != "ENSRNOG00000049895" & etz_gene != 305626) %>%
  filter(ens_gene != "ENSRNOG00000049944" & etz_gene != 309111) %>% 
  filter(ens_gene != "ENSRNOG00000049964" & etz_gene != 103690044) %>%
  filter(ens_gene != "ENSRNOG00000050201" & etz_gene != 65194) %>%
  filter(ens_gene != "ENSRNOG00000050885" & etz_gene != 171072) %>%
  filter(ens_gene != "ENSRNOG00000055277" & etz_gene != 298098) %>%
  filter(ens_gene != "ENSRNOG00000056076" & etz_gene != 103694877) %>% 
  filter(ens_gene != "ENSRNOG00000043249" & etz_gene != 290646) %>% 
  filter(ens_gene != "ENSRNOG00000029861" & etz_gene != 494499 ) %>% 
  filter(ens_gene != "ENSRNOG00000061714" & etz_gene != 54302)

# Use TxImport to count transcripts 
txi <- tximport::tximport(files,type = "kallisto", txOut = FALSE, tx2gene = tx2gene, importer = read_tsv)
txi2 <- tximport::tximport(files2,type = "kallisto", txOut = FALSE, tx2gene = tx2gene, importer = read_tsv)
txi3 <- tximport::tximport(files3,type = "kallisto", txOut = FALSE, tx2gene = tx2gene, importer = read_tsv)
txi6 <- tximport::tximport(files6,type = "kallisto", txOut = FALSE, tx2gene = tx2gene, importer = read_tsv)

# Label all conditions in the matrix in order of appearance in the each of the files variables
condition = c("APAP","APAP","APAP","APAP","CCL4","CCL4","CCL4","Ctrl","Ctrl","Ctrl",
              "TCDD","TCDD","TCDD","TCDD","TCE","TCE","TCE","TCE")
condition2 = c("Ctrl","Ctrl","Ctrl","TCDD","TCDD","TCDD","TCE","TCE","TCE","TCE")
condition3 = c("APAP","APAP","APAP","APAP","Ctrl","Ctrl","Ctrl")
condition6 = c("Ctrl","Ctrl","Ctrl","Ctrl","CCL4","CCL4","CCL4","CCL4")

# Change condition vector from character vector to factor vector
condition = factor(condition)
condition2 = factor(condition2)
condition3 = factor(condition3)
condition6 = factor(condition6)

# Convert conditions to data frame for maniupaltion 
sampleTable <- data.frame(condition)
sampleTable2 <- data.frame(condition2)
sampleTable3 <- data.frame(condition3)
sampleTable6 <- data.frame(condition6)

# Create row names based on the genes in the same order
rownames(sampleTable) <- colnames(txi$counts)
rownames(sampleTable2) <- colnames(txi2$counts)
rownames(sampleTable3) <- colnames(txi3$counts)
rownames(sampleTable6) <- colnames(txi6$counts)

# Import TxImport Results into DESeq2 structure
dds <- DESeq2::DESeqDataSetFromTximport(txi,sampleTable,~condition)
dds2 <- DESeq2::DESeqDataSetFromTximport(txi2,sampleTable2,~condition2)
dds3 <- DESeq2::DESeqDataSetFromTximport(txi3,sampleTable3,~condition3)
dds6 <- DESeq2::DESeqDataSetFromTximport(txi6,sampleTable6,~condition6)

# Pre-filtering of genes that have at least 1 count in a condition
dds <- dds[ rowSums(DESeq2::counts(dds)) > 1, ]
dds2 <- dds2[ rowSums(DESeq2::counts(dds2)) > 1, ]
dds3 <- dds3[ rowSums(DESeq2::counts(dds3)) > 1, ]
dds6 <- dds6[ rowSums(DESeq2::counts(dds6)) > 1, ]

# Create Gene level analysis in DESeq2 format
dds <- DESeq2::DESeq(dds, minReplicatesForReplace = Inf)
dds2 <- DESeq2::DESeq(dds2, minReplicatesForReplace = Inf)
dds3 <- DESeq2::DESeq(dds3, minReplicatesForReplace = Inf)
dds6 <- DESeq2::DESeq(dds6, minReplicatesForReplace = Inf)

# Acetaminophen 6hr DEGs (we want treated vs. untreated below)
res05_apap <- DESeq2::results(dds,contrast = c("condition", "APAP", "Ctrl") ,alpha = 0.05, pAdjustMethod = "BH", cooksCutoff = FALSE,independentFiltering = FALSE)
res05_apap <- res05_apap[order(res05_apap$padj),]

# CCl4 6hr DEGs
res05_ccl4 <- DESeq2::results(dds, contrast = c("condition", "CCL4","Ctrl"), alpha = 0.05, pAdjustMethod = "BH", cooksCutoff = FALSE,independentFiltering = FALSE)
res05_ccl4 <- res05_ccl4[order(res05_ccl4$padj),]

# TCDD 6hr DEGs
res05_tcdd <- DESeq2::results(dds, contrast = c("condition","TCDD", "Ctrl"), alpha = 0.05,pAdjustMethod = "BH", cooksCutoff = FALSE,independentFiltering = FALSE)
res05_tcdd <- res05_tcdd[order(res05_tcdd$padj),]

# TCE 6hr DEGs
res05_tce <- DESeq2::results(dds, contrast = c("condition", "TCE","Ctrl"), alpha = 0.05,pAdjustMethod = "BH", cooksCutoff = FALSE,independentFiltering = FALSE)
res05_tce <- res05_tce[order(res05_tce$padj),]

# Acetaminophen 24hr DEGs
res05_apap24 <- DESeq2::results(dds3, contrast = c("condition3", "APAP","Ctrl"), alpha = 0.05, pAdjustMethod = "BH", cooksCutoff = FALSE,independentFiltering = FALSE)
res05_apap24 <- res05_apap24[order(res05_apap24$padj),]

# CCL4 24hr DEGs
res05_ccl4.24 <- DESeq2::results(dds6, contrast = c("condition6", "CCL4","Ctrl"), alpha = 0.05, pAdjustMethod = "BH", cooksCutoff = FALSE,independentFiltering = FALSE)
res05_ccl4.24 <- res05_ccl4.24[order(res05_ccl4.24$padj),]

# TCDD 24hr DEGs
res05_tcdd24 <- DESeq2::results(dds2, contrast = c("condition2", "TCDD","Ctrl"), alpha = 0.05,pAdjustMethod = "BH", cooksCutoff = FALSE,independentFiltering = FALSE)
res05_tcdd24 <- res05_tcdd24[order(res05_tcdd24$padj),]

# TCE 24hr DEGs
res05_tce24 <- DESeq2::results(dds2, contrast = c("condition2", "TCE","Ctrl"), alpha = 0.05,pAdjustMethod = "BH", cooksCutoff = FALSE,independentFiltering = FALSE)
res05_tce24 <- res05_tce24[order(res05_tce24$padj),]

# Create results data frame
apap6hr_DEGs <- data.frame(DESeq2::results(dds, contrast = c("condition", "APAP", "Ctrl"), alpha = 0.05, pAdjustMethod = "BH", cooksCutoff = FALSE,independentFiltering = FALSE)) 
apap6hr_DEGs$Names <- rownames(apap6hr_DEGs)
apap6hr_expression <- inner_join(apap6hr_DEGs,dds2model,by=c("Names"="ens_gene"))

ccl46hr_DEGs <- data.frame(DESeq2::results(dds, contrast = c("condition" , "CCL4", "Ctrl"), alpha = 0.05, pAdjustMethod = "BH", cooksCutoff = FALSE,independentFiltering = FALSE)) 
ccl46hr_DEGs$Names <- rownames(ccl46hr_DEGs)
ccl46hr_expression <- inner_join(ccl46hr_DEGs,dds2model,by=c("Names"="ens_gene"))
  
tcdd6hr_DEGs <- data.frame(DESeq2::results(dds, contrast = c("condition" , "TCDD", "Ctrl"), alpha = 0.05, pAdjustMethod = "BH", cooksCutoff = FALSE,independentFiltering = FALSE)) 
tcdd6hr_DEGs$Names <- rownames(tcdd6hr_DEGs)
tcdd6hr_expression <- inner_join(tcdd6hr_DEGs,dds2model,by=c("Names"="ens_gene"))

tce6hr_DEGs <- data.frame(DESeq2::results(dds, contrast = c("condition" , "TCE", "Ctrl"), alpha = 0.05, pAdjustMethod = "BH", cooksCutoff = FALSE,independentFiltering = FALSE)) 
tce6hr_DEGs$Names <- rownames(tce6hr_DEGs)
tce6hr_expression <- inner_join(tce6hr_DEGs,dds2model,by=c("Names"="ens_gene"))

apap24hr_DEGs <- data.frame(DESeq2::results(dds3, contrast = c("condition3" , "APAP", "Ctrl"), alpha = 0.05, pAdjustMethod = "BH", cooksCutoff = FALSE,independentFiltering = FALSE)) 
apap24hr_DEGs$Names <- rownames(apap24hr_DEGs)
apap24hr_expression <- inner_join(apap24hr_DEGs,dds2model,by=c("Names"="ens_gene"))

ccl424hr_DEGs <- data.frame(DESeq2::results(dds6, contrast = c("condition6" , "CCL4", "Ctrl"), alpha = 0.05, pAdjustMethod = "BH", cooksCutoff = FALSE,independentFiltering = FALSE)) 
ccl424hr_DEGs$Names <- rownames(ccl424hr_DEGs)
ccl424hr_expression <- inner_join(ccl424hr_DEGs,dds2model,by=c("Names"="ens_gene"))

tcdd24hr_DEGs <- data.frame(DESeq2::results(dds2, contrast = c("condition2" , "TCDD", "Ctrl"), alpha = 0.05, pAdjustMethod = "BH", cooksCutoff = FALSE,independentFiltering = FALSE)) 
tcdd24hr_DEGs$Names <- rownames(tcdd24hr_DEGs)
tcdd24hr_expression <- inner_join(tcdd24hr_DEGs,dds2model,by=c("Names"="ens_gene"))

tce24hr_DEGs <- data.frame(DESeq2::results(dds2, contrast = c("condition2" , "TCE", "Ctrl"), alpha = 0.05, pAdjustMethod = "BH", cooksCutoff = FALSE,independentFiltering = FALSE)) 
tce24hr_DEGs$Names <- rownames(tce24hr_DEGs)
tce24hr_expression <- inner_join(tce24hr_DEGs,dds2model,by=c("Names"="ens_gene"))

# Write csv files for all differentially expressed genes 
apap6hr_expression %>%  filter(padj < 0.1) %>% 
  select(etz_gene,log2FoldChange,padj) %>% unique() %>% 
  write_csv("~/Desktop/050118_DavidAnalysis/apap6hr_allDEGs.csv")

ccl46hr_expression %>%  filter(padj < 0.1) %>% 
  select(etz_gene,log2FoldChange,padj) %>% unique() %>% 
  write_csv("~/Desktop/050118_DavidAnalysis/ccl46hr_allDEGs.csv")

tcdd6hr_expression %>%  filter(padj < 0.1) %>% 
  select(etz_gene,log2FoldChange,padj) %>% unique() %>% 
  write_csv("~/Desktop/050118_DavidAnalysis/tcdd6hr_allDEGs.csv")

tce6hr_expression %>%  filter(padj < 0.1) %>% 
  select(etz_gene,log2FoldChange,padj) %>% unique() %>% 
  write_csv("~/Desktop/050118_DavidAnalysis/tce6hr_allDEGs.csv")

apap24hr_expression %>%  filter(padj < 0.1) %>% 
  select(etz_gene,log2FoldChange,padj) %>% unique() %>% 
  write_csv("~/Desktop/050118_DavidAnalysis/apap24hr_allDEGs.csv")

ccl424hr_expression %>%  filter(padj < 0.1) %>% 
  select(etz_gene,log2FoldChange,padj) %>% unique() %>% 
  write_csv("~/Desktop/050118_DavidAnalysis/ccl424hr_allDEGs.csv")

tcdd24hr_expression %>%  filter(padj < 0.1) %>% 
  select(etz_gene,log2FoldChange,padj) %>% unique() %>% 
  write_csv("~/Desktop/050118_DavidAnalysis/tcdd24hr_allDEGs.csv")

tce24hr_expression %>%  filter(padj < 0.1) %>% 
  select(etz_gene,log2FoldChange,padj) %>% unique() %>% 
  write_csv("~/Desktop/050118_DavidAnalysis/tce24hr_allDEGs.csv")


#Pull out only metabolic genes 
# Find Genes that Map to the Metabolic Model (Locate all metabolic genes)
rno.rxn.gene = readWorkbook("ncomm_blais_data_rno_raven.xlsx",sheet = 5, startRow = 1) %>% as.tbl
rno.rxn.gene <- as.data.frame(rno.rxn.gene)
names(rno.rxn.gene)[5] <- "GENE_NAME"
rno.rxn.gene$etz_gene <- as.integer(rno.rxn.gene$GENE.NAME)

# Right join to see shared genes
apap6hr_genes <- semi_join(apap6hr_expression,rno.rxn.gene, by="etz_gene")
apap6hr_genes <- apap6hr_genes %>% dplyr::select(baseMean,log2FoldChange,pvalue,padj,etz_gene) %>% drop_na()

ccl46hr_genes <- semi_join(ccl46hr_expression,rno.rxn.gene, by="etz_gene")
ccl46hr_genes <- ccl46hr_genes %>% dplyr::select(baseMean,log2FoldChange,pvalue,padj,etz_gene) %>% drop_na()

tcdd6hr_genes <- semi_join(tcdd6hr_expression,rno.rxn.gene, by="etz_gene")
tcdd6hr_genes <- tcdd6hr_genes %>% dplyr::select(baseMean,log2FoldChange,pvalue,padj,etz_gene) %>% drop_na()

tce6hr_genes <- semi_join(tce6hr_expression,rno.rxn.gene, by="etz_gene")
tce6hr_genes <- tce6hr_genes %>% dplyr::select(baseMean,log2FoldChange,pvalue,padj,etz_gene) %>% drop_na()

apap24hr_genes <- semi_join(apap24hr_expression,rno.rxn.gene, by="etz_gene")
apap24hr_genes <- apap24hr_genes %>% dplyr::select(baseMean,log2FoldChange,pvalue,padj,etz_gene) %>% drop_na()

ccl424hr_genes <- semi_join(ccl424hr_expression,rno.rxn.gene, by="etz_gene")
ccl424hr_genes <- ccl424hr_genes %>% dplyr::select(baseMean,log2FoldChange,pvalue,padj,etz_gene) %>% drop_na()

tcdd24hr_genes <- semi_join(tcdd24hr_expression,rno.rxn.gene, by="etz_gene")
tcdd24hr_genes <- tcdd24hr_genes %>% dplyr::select(baseMean,log2FoldChange,pvalue,padj,etz_gene) %>% drop_na()

tce24hr_genes <- semi_join(tce24hr_expression,rno.rxn.gene, by="etz_gene")
tce24hr_genes <- tce24hr_genes %>% dplyr::select(baseMean,log2FoldChange,pvalue,padj,etz_gene) %>% drop_na()

# Add Column for significance 
apap6hr_genes$sig <- as.numeric(apap6hr_genes$padj < 0.1)
apap24hr_genes$sig <- as.numeric(apap24hr_genes$padj < 0.1)
ccl46hr_genes$sig <- as.numeric(ccl46hr_genes$padj < 0.1)
ccl424hr_genes$sig <- as.numeric(ccl424hr_genes$padj < 0.1)
tcdd6hr_genes$sig <- as.numeric(tcdd6hr_genes$padj < 0.1)
tcdd24hr_genes$sig <- as.numeric(tcdd24hr_genes$padj < 0.1)
tce6hr_genes$sig <- as.numeric(tce6hr_genes$padj < 0.1)
tce24hr_genes$sig <- as.numeric(tce24hr_genes$padj < 0.1)

merg_vec1 <- inner_join(apap24hr_genes,apap6hr_genes,by="etz_gene") %>% select(etz_gene,sig.x,sig.y) %>% rename(sig.x = "apap24hr_sig", sig.y = "apap6hr_sig") %>% unique()
merg_vec2 <- inner_join(ccl424hr_genes,ccl46hr_genes,by="etz_gene") %>% select(etz_gene,sig.x,sig.y) %>% rename(sig.x = "ccl424hr_sig", sig.y = "ccl46hr_sig") %>% unique()
merg_vec3 <- inner_join(tcdd24hr_genes,tcdd6hr_genes,by="etz_gene") %>% select(etz_gene,sig.x,sig.y) %>% rename(sig.x = "tcdd24hr_sig", sig.y = "tcdd6hr_sig") %>% unique()
merg_vec4 <- inner_join(tce24hr_genes,tce6hr_genes,by="etz_gene") %>% select(etz_gene,sig.x,sig.y) %>% rename(sig.x = "tce24hr_sig", sig.y = "tce6hr_sig") %>% unique()

merg_vec5 <- inner_join(merg_vec1,merg_vec2, by="etz_gene") %>% unique()
merg_vec6 <- inner_join(merg_vec3,merg_vec4, by="etz_gene") %>% unique()
merg_vec7 <- inner_join(merg_vec6,merg_vec5, by="etz_gene") %>% unique()

gene_significance <- merg_vec7
gene_significance$summary <- apply(gene_significance[,c(2:9)],1,sum) 
gene_sigs <- gene_significance %>% filter(summary > 0) %>% select(etz_gene)

# Heat map of Expression Changes 
tce_metgenes <- inner_join(tce24hr_genes,tce6hr_genes, by="etz_gene")
tce_metgenes <- tce_metgenes[!duplicated(tce_metgenes),] %>% 
  select(etz_gene,log2FoldChange.x,log2FoldChange.y) %>% 
  dplyr::rename("tce24hr" = log2FoldChange.x,"tce6hr" = log2FoldChange.y)

tcdd_metgenes <- inner_join(tcdd24hr_genes, tcdd6hr_genes, by="etz_gene")
tcdd_metgenes <- tcdd_metgenes[!duplicated(tcdd_metgenes),]%>% 
  select(etz_gene,log2FoldChange.x,log2FoldChange.y) %>% 
  dplyr::rename("tcdd24hr" = log2FoldChange.x,"tcdd6hr" = log2FoldChange.y)

ccl4_metgenes <- inner_join(ccl424hr_genes, ccl46hr_genes, by="etz_gene")
ccl4_metgenes <- ccl4_metgenes[!duplicated(ccl4_metgenes),] %>% 
  select(etz_gene,log2FoldChange.x,log2FoldChange.y) %>% 
  dplyr::rename("ccl424hr" = log2FoldChange.x,"ccl46hr" = log2FoldChange.y)

apap_metgenes <- inner_join(apap24hr_genes, apap6hr_genes, by="etz_gene")
apap_metgenes <- apap_metgenes[!duplicated(apap_metgenes),]%>% 
  select(etz_gene,log2FoldChange.x,log2FoldChange.y) %>% 
  dplyr::rename("apap24hr" = log2FoldChange.x,"apap6hr" = log2FoldChange.y)

joint_tbl1 <- inner_join(apap_metgenes,ccl4_metgenes, by="etz_gene")
joint_tbl1 <- joint_tbl1[!duplicated(joint_tbl1),]
joint_tbl2 <- inner_join(tcdd_metgenes,tce_metgenes, by="etz_gene")
joint_tbl2 <- joint_tbl2[!duplicated(joint_tbl2),]

all_metgenes <- inner_join(joint_tbl1,joint_tbl2, by="etz_gene")
all_metgenes <- all_metgenes[!duplicated(all_metgenes),]

all_sig_metgenes <- left_join(gene_sigs,all_metgenes,by="etz_gene")

metgenes_hmap <- all_sig_metgenes[,c(2:9)]
trans_metgenes.hmap <- t(metgenes_hmap)
conds <- c("apap6hr","apap24hr","ccl46hr","ccl424hr","tcdd6hr","tcdd24hr",
           "tce6hr","tce24hr")
trans_metgenes.hmap <- trans_metgenes.hmap[conds,]

# Color Palatte
my_palette <- colorRampPalette(c("blue","white","red"))
x_break <- seq(-0.75,0.75,0.01)

metgenes_hmap <- all_sig_metgenes[,c(2:9)]
rownames(metgenes_hmap) <- all_sig_metgenes$etz_gene

# Vertical HEaetmap 
hep_rnaseq <- heatmap.2(as.matrix(metgenes_hmap), trace = "none", density.info = "none",dendrogram = "both", 
  breaks = x_break,col=my_palette,keysize = 1, margins = c(8,8),ColSideColors = c(rep("green",2),rep("yellow",2),
  rep("blue",2),rep("red",2)),sepcolor = "black", cexCol = 2,labCol = colnames(metgenes_hmap), labRow = all_sig_metgenes$etz_gene )


# Prepare to write out a table of differentially expressed genes
output_matrix <- metgenes_hmap[rev(hep_rnaseq$rowInd), hep_rnaseq$colInd]
t2g6$etz_gene <- as.integer(t2g6$etz_gene)
output_matrix <- rownames_to_column(output_matrix, var = "etz_gene") 
output_matrix$etz_gene <- as.integer(output_matrix$etz_gene)
output_heatmap <- output_matrix %>% left_join(t2g6,by="etz_gene") %>% 
  select(apap6hr,tcdd6hr,apap24hr,ccl46hr,tce6hr,tcdd24hr,tce24hr,ccl424hr,etz_gene,descript,phenotype)

# Write out differentially expressed genes 
write.xlsx(output_heatmap,file = "Rawls_supplement_transcriptomics.xlsx",row.names = TRUE)

# Recreate Figure 2 
transcriptomics <- read.xlsx("Rawls_supplement_transcriptomics.xlsx",sheet= 1) %>% 
  select(etz_gene,apap6hr,tcdd6hr,apap24hr,ccl46hr,tce6hr,tcdd24hr,tce24hr,ccl424hr) %>% unique()

entrez_genes <- transcriptomics$etz_gene
transcriptomics <- transcriptomics[,c(2:9)]
transcriptomics.bin <- transcriptomics
transcriptomics.bin[transcriptomics.bin > 0.75] <- 0.75
transcriptomics.bin[transcriptomics.bin < -0.75] <- -0.75

transcriptomics.mat <- transcriptomics
rownames(transcriptomics.mat) <- entrez_genes

# Condition and gene clustering 
transcriptomics.clust <- hclust(dist(transcriptomics.mat, method = "euclidean"), method = "complete")
condition.clust <- hclust(dist(t(transcriptomics.mat), method = "euclidean"), method = "complete")
transcriptomics.bin$etz_gene <- entrez_genes

# Reshape data into plottable format for ggplot2
transcriptomics.melt <- melt(transcriptomics.bin, id.vars = "etz_gene",variable.name = "Condition", value.name = "FoldChange")
transcriptomics.melt$Condition <- factor(transcriptomics.melt$Condition, 
                                         levels = c("ccl424hr","apap6hr","tcdd6hr","apap24hr","tcdd24hr","tce24hr","ccl46hr","tce6hr"))
transcriptomics.melt$etz_gene <- factor(transcriptomics.melt$etz_gene,
                                        levels = entrez_genes[transcriptomics.clust$order],ordered = TRUE)

f2.1 <- ggplot(transcriptomics.melt, aes(x = Condition,y = etz_gene)) + geom_tile(aes(fill = as.numeric(FoldChange)),size = 0.01) +
  scale_fill_gradient2(low = "blue",mid = "white",high = "red",limits = c(-0.75,.75), labels = c("Downregulated","Upregulated"), breaks = c(-0.75,0.75)) + 
  scale_x_discrete(labels = c(ccl424hr = expression(paste("CCl"[4]*" 24hr")), 
                              apap6hr = expression(paste("APAP 6hr")), tcdd6hr = expression(paste("TCDD 6hr")), apap24hr = expression(paste("APAP 24hr")),
                              tcdd24hr = expression(paste("TCDD 24hr")), tce24hr = expression(paste("TCE 24hr")),
                              ccl46hr = expression(paste("CCl"[4]*" 6hr")), tce6hr = expression(paste("TCE 6hr")))) + 
  theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank()) +
  theme(axis.text.y = element_blank(),axis.text.x = element_text(family = "Times",face = "bold",size = 10)) +
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.margin = unit(c(0.5,-0.5,0.5,0.5),"cm")) +
  theme(legend.position = "left") + 
  guides(fill=guide_legend(title = expression(paste("Metabolic log"[2]*"Fold Changes")))) + 
  labs(x = "", y = "")

# Create GG Plot for dendrogram of genes
f2.2 <- ggdendrogram(transcriptomics.clust, rotate = TRUE, theme_dendro = TRUE) + 
  scale_x_reverse() + 
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) + 
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank())

f2.3 <- ggdendrogram(condition.clust, rotate = FALSE, theme_dendro = TRUE) +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggdraw() + 
  draw_plot(f2.3, x = 0.34, y = 0.765, width = 0.53, height = 0.2) + 
  draw_plot(f2.2, x = 0.855, y = 0.053, width = 0.15, height = 0.795) + 
  draw_plot(f2.1, x = 0, y = 0, width =  0.85, height = 0.83)

# Output Genes for TIMBR 
apap6hr_genes$fdr <- p.adjust(apap6hr_genes$pvalue, method = "fdr")
apap24hr_genes$fdr <- p.adjust(apap24hr_genes$pvalue, method = "fdr")
ccl46hr_genes$fdr <- p.adjust(ccl46hr_genes$pvalue, method = "fdr")
ccl424hr_genes$fdr <- p.adjust(ccl424hr_genes$pvalue, method = "fdr")
tcdd6hr_genes$fdr <- p.adjust(tcdd6hr_genes$pvalue, method = "fdr")
tcdd24hr_genes$fdr <- p.adjust(tcdd24hr_genes$pvalue, method = "fdr")
tce6hr_genes$fdr <- p.adjust(tce6hr_genes$pvalue, method = "fdr")
tce24hr_genes$fdr <- p.adjust(tce24hr_genes$pvalue, method = "fdr")

# Rename columns and write to csv files
apap6hr_TIMBR <- apap6hr_genes %>% 
  mutate(ave = baseMean, logfc = log2FoldChange, pval = pvalue) %>% 
  select(etz_gene,pval,padj,fdr,ave,logfc) %>% write_csv("apap_6hr_TIMBR.csv")

apap24hr_TIMBR <- apap24hr_genes %>% 
  mutate(ave = baseMean, logfc = log2FoldChange, pval = pvalue) %>% 
  select(etz_gene,pval,padj,fdr,ave,logfc) %>% write_csv("apap_24hr_TIMBR.csv")

ccl46hr_TIMBR <- ccl46hr_genes %>% 
  mutate(ave = baseMean, logfc = log2FoldChange, pval = pvalue) %>% 
  select(etz_gene,pval,padj,fdr,ave,logfc) %>% write_csv("ccl4_6hr_TIMBR.csv")

ccl424hr_TIMBR <- ccl424hr_genes %>% 
  mutate(ave = baseMean, logfc = log2FoldChange, pval = pvalue) %>% 
  select(etz_gene,pval,padj,fdr,ave,logfc) %>% write_csv("ccl4_24hr_TIMBR.csv")

tcdd6hr_TIMBR <- tcdd6hr_genes %>% 
  mutate(ave = baseMean, logfc = log2FoldChange, pval = pvalue) %>% 
  select(etz_gene,pval,padj,fdr,ave,logfc) %>% write_csv("tcdd_6hr_TIMBR.csv")

tcdd24hr_TIMBR <- tcdd24hr_genes %>% 
  mutate(ave = baseMean, logfc = log2FoldChange, pval = pvalue) %>% 
  select(etz_gene,pval,padj,fdr,ave,logfc) %>% write_csv("tcdd_24hr_TIMBR.csv")

tce6hr_TIMBR <- tce6hr_genes %>% 
  mutate(ave = baseMean, logfc = log2FoldChange, pval = pvalue) %>% 
  select(etz_gene,pval,padj,fdr,ave,logfc) %>% write_csv("tce_6hr_TIMBR.csv")

tce24hr_TIMBR <- tce24hr_genes %>% 
  mutate(ave = baseMean, logfc = log2FoldChange, pval = pvalue) %>% 
  select(etz_gene,pval,padj,fdr,ave,logfc) %>% write_csv("tce_24hr_TIMBR.csv")
