# Rawls Transcritpomics Data Analysis for Liver Data
# Uses Tximport + DESeq2 to analyze results
# Analyze Primary Hepatocytes Transcritopmic Data
# Outputs DEGs spreadsheets 
# Imports Gene Enrichment Data 
# Outputs Figure 2. 
# 04/17/19

# Load Pacakges 
library("tidyverse")
library("openxlsx")
library("matrixStats")
library("gplots")
library("reshape2")
library("ggdendro")
library("grid")
library("gridExtra")
library("ggpubr")

#Source files from bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("tximport") 
biocLite("biomaRt")
biocLite("DESeq2")
biocLite('vsn')


# Create file names to get abundance files from Kallisto 
dir <- "/Users/kdrawls/Documents/CSBL/Experimental_Results/011717kallisto/6hr_data"

# Check for all files in the directory to ensure they are present 
list.files(dir)

# Read in meta data table
samples <- read.table(file.path("/Users/kdrawls/Documents/CSBL/Experimental_Results/011717kallisto/","papin_genewiz6hr_rnaseq_info.txt"), header = TRUE)

# Read in each file based off the format, and label sample numbers
files <- file.path(dir,samples$run_accession,"abundance.tsv")
names(files) <- paste0("sample",1:dim(samples)[1])

# Use biomArt to gather transcript and gene IDs
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
  filter(ens_gene != "ENSRNOG00000061714" & etz_gene != 54302)


###### Use TxImport to summarize transcript changes to the gene level #####################################
txi <- tximport::tximport(files,type = "kallisto", txOut = FALSE, tx2gene = tx2gene, importer = read_tsv)

# Label all conditions in the matrix in order of appearance in the each of the files variables
condition = c("APAP","APAP","APAP","APAP","CCL4","CCL4","CCL4","Ctrl","Ctrl","Ctrl",
              "TCDD","TCDD","TCDD","TCDD","TCE","TCE","TCE","TCE")

# Change condition vector from character vector to factor vector
condition = factor(condition)

# Convert conditions to data frame for maniupaltion 
sampleTable <- data.frame(condition)

# Create row names based on the genes in the same order
rownames(sampleTable) <- colnames(txi$counts)

###################################### Use DESeq2 to analyze differential gene expression ########################################
# Import TxImport Results into DESeq2 structure
dds <- DESeq2::DESeqDataSetFromTximport(txi,sampleTable,~condition)

# Pre-filtering of genes that have at least 1 count in a condition
dds <- dds[ rowSums(DESeq2::counts(dds)) > 1, ]

# Create Gene level analysis in DESeq2 format
dds <- DESeq2::DESeq(dds, minReplicatesForReplace = Inf)

# Create results data frames
apap6hr_DEGs <- data.frame(DESeq2::results(dds, contrast = c("condition", "APAP", "Ctrl"), alpha = 0.05, pAdjustMethod = "BH", cooksCutoff = FALSE,independentFiltering = FALSE)) 
apap6hr_DEGs$Names <- rownames(apap6hr_DEGs)
apap6hr_DEGs$fdr <- p.adjust(apap6hr_DEGs$pvalue, method = "fdr")
apap6hr_expression <- inner_join(apap6hr_DEGs,dds2model,by=c("Names"="ens_gene"))

ccl46hr_DEGs <- data.frame(DESeq2::results(dds, contrast = c("condition" , "CCL4", "Ctrl"), alpha = 0.05, pAdjustMethod = "BH", cooksCutoff = FALSE,independentFiltering = FALSE)) 
ccl46hr_DEGs$Names <- rownames(ccl46hr_DEGs)
ccl46hr_DEGs$fdr <- p.adjust(ccl46hr_DEGs$pvalue, method = "fdr")
ccl46hr_expression <- inner_join(ccl46hr_DEGs,dds2model,by=c("Names"="ens_gene"))

tcdd6hr_DEGs <- data.frame(DESeq2::results(dds, contrast = c("condition" , "TCDD", "Ctrl"), alpha = 0.05, pAdjustMethod = "BH", cooksCutoff = FALSE,independentFiltering = FALSE)) 
tcdd6hr_DEGs$Names <- rownames(tcdd6hr_DEGs)
tcdd6hr_DEGs$fdr <- p.adjust(tcdd6hr_DEGs$pvalue, method = "fdr")
tcdd6hr_expression <- inner_join(tcdd6hr_DEGs,dds2model,by=c("Names"="ens_gene"))

tce6hr_DEGs <- data.frame(DESeq2::results(dds, contrast = c("condition" , "TCE", "Ctrl"), alpha = 0.05, pAdjustMethod = "BH", cooksCutoff = FALSE,independentFiltering = FALSE)) 
tce6hr_DEGs$Names <- rownames(tce6hr_DEGs)
tce6hr_DEGs$fdr <- p.adjust(tce6hr_DEGs$pvalue, method = "fdr")
tce6hr_expression <- inner_join(tce6hr_DEGs,dds2model,by=c("Names"="ens_gene"))

## Save workbook which is supplemental file 1 -------------

# Create workbook 
transcriptomics <- createWorkbook()

# Add blank sheets to the workbook 
addWorksheet(transcriptomics, "Overview") # Add a sheet to describe the supplemental file 
addWorksheet(transcriptomics, "apap6hr_allgenes") # Add sheet for APAP
addWorksheet(transcriptomics, "ccl46hr_allgenes") # Add sheet for CCl4
addWorksheet(transcriptomics, "tcdd6hr_allgenes") # Add sheet for TCDD
addWorksheet(transcriptomics,  "tce6hr_allgenes") # Add sheet for TCE

# Write data for each sheet 
writeData(transcriptomics, sheet = "apap6hr_allgenes", x = apap6hr_expression)
writeData(transcriptomics, sheet = "ccl46hr_allgenes", x = ccl46hr_expression)
writeData(transcriptomics, sheet = "tcdd6hr_allgenes", x = tcdd6hr_expression)
writeData(transcriptomics, sheet =  "tce6hr_allgenes", x =  tce6hr_expression)

# Export Data File
saveWorkbook(transcriptomics, "Rawls_Supplementary_Data1.xlsx")

## Create Heatmap from differentially expressed genes -------------

# Filter for metabolic genes
# Read in workbook from "Reconciled rat and human metabolic networks for comparative toxicogenomics analyses and biomarker predictions" by Blais et al.
# Article and data files are published at https://www.nature.com/articles/ncomms14250
rno.rxn.gene = readWorkbook("ncomm_blais_data_rno_raven.xlsx",sheet = 5, startRow = 1) %>% as.tbl
rno.rxn.gene <- as.data.frame(rno.rxn.gene)
names(rno.rxn.gene)[5] <- "GENE_NAME"
rno.rxn.gene$etz_gene <- as.integer(rno.rxn.gene$GENE.NAME)

apap6hr_genes <- inner_join(apap6hr_expression,rno.rxn.gene, by="etz_gene")
apap6hr_genes <- apap6hr_genes[!duplicated(apap6hr_genes),]
apap6hr_genes <- apap6hr_genes %>% dplyr::select(Names,baseMean,log2FoldChange,pvalue,padj,etz_gene,fdr) %>% drop_na()

ccl46hr_genes <- inner_join(ccl46hr_expression,rno.rxn.gene, by="etz_gene")
ccl46hr_genes <- ccl46hr_genes[!duplicated(ccl46hr_genes),]
ccl46hr_genes <- ccl46hr_genes %>% dplyr::select(Names,baseMean,log2FoldChange,pvalue,padj,etz_gene,fdr) %>% drop_na()

tcdd6hr_genes <- inner_join(tcdd6hr_expression,rno.rxn.gene, by="etz_gene")
tcdd6hr_genes <- tcdd6hr_genes[!duplicated(tcdd6hr_genes),]
tcdd6hr_genes <- tcdd6hr_genes %>% dplyr::select(Names,baseMean,log2FoldChange,pvalue,padj,etz_gene,fdr) %>% drop_na()

tce6hr_genes <- inner_join(tce6hr_expression,rno.rxn.gene, by="etz_gene")
tce6hr_genes <- tce6hr_genes[!duplicated(tce6hr_genes),]
tce6hr_genes <- tce6hr_genes %>% dplyr::select(Names,baseMean,log2FoldChange,pvalue,padj,etz_gene,fdr) %>% drop_na()


clust1_metgenes <- inner_join(tce6hr_genes,tcdd6hr_genes, by="Names")
clust1_metgenes <- clust1_metgenes[!duplicated(clust1_metgenes),] %>% 
  select(Names,log2FoldChange.x,log2FoldChange.y) %>% 
  dplyr::rename("tce6hr" = log2FoldChange.x,"tcdd6hr" = log2FoldChange.y)


clust2_metgenes <- inner_join(ccl46hr_genes, apap6hr_genes, by="Names")
clust2_metgenes <- clust2_metgenes[!duplicated(clust2_metgenes),] %>% 
  select(Names,log2FoldChange.x,log2FoldChange.y) %>% 
  dplyr::rename("ccl46hr" = log2FoldChange.x,"apap6hr" = log2FoldChange.y)


all_metgenes <- inner_join(clust1_metgenes,clust2_metgenes, by="Names")
all_metgenes <- all_metgenes[!duplicated(all_metgenes),]

data_hmap <- all_metgenes[,c(2:5)]


# Normal ggplot2 Heatmap 
rownames(data_hmap) <- all_metgenes$Names
hr <- hclust(dist(data_hmap, method = "euclidean"), method = "complete")
hc <- hclust(as.dist(1-cor(data_hmap, method = "spearman")), method = "complete")

melted_metgenes <- melt(all_metgenes, id.vars = "Names", 
                        variable.name = "condition",value.name = "FoldChange") %>% 
  rename(ens_gene = "Names")

gene_order <- as.list(all_metgenes[hr$order,1]) %>% unlist()
condition_unorder <- colnames(all_metgenes[,c(2:5)])
condition_order <- condition_unorder[hc$order]
melted_metgenes$ens_gene <- factor(melted_metgenes$ens_gene,levels = gene_order, ordered = TRUE)
melted_metgenes$condition <- factor(melted_metgenes$condition,levels = condition_order, ordered = TRUE)

melted_metgenes$FoldChange[melted_metgenes$FoldChange < -1] <- -1
melted_metgenes$FoldChange[melted_metgenes$FoldChange > 1] <- 1

f1.1 <- ggplot(melted_metgenes, aes(x = condition,y = ens_gene)) + geom_tile(aes(fill = as.numeric(FoldChange)),size = 0.01) +
  scale_fill_gradient2(low = "blue",mid = "white",high = "red",limits = c(-1,1), labels = c("Downregulated","Upregulated"), breaks = c(-1,1)) + 
  scale_x_discrete(position = "bottom", labels = c(apap6hr = expression(paste("APAP 6hr")), tcdd6hr = expression(paste("TCDD 6hr")),
      ccl46hr = expression(paste("CCl"[4]*" 6hr")), tce6hr = expression(paste("TCE 6hr")))) + 
  theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank()) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(family = "Arial",face = "bold",size = 16, 
                                   angle = 30, vjust = 0.8, hjust = 1)) +
  theme(axis.ticks.y = element_blank()) +
  theme(plot.margin = unit(c(0.5,-0.5,0.5,0.5),"cm")) +
  theme(legend.position = "left") + theme(legend.title=element_text(size=20), legend.text=element_text(size=16)) +
  guides(fill=guide_legend(title = expression(paste("Metabolic log"[2]*"Fold Changes")))) + 
  labs(x = "", y = "")


# Create GG Plot for dendrogram of genes
f1.2 <- ggdendrogram(hr, rotate = TRUE, theme_dendro = TRUE) + 
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) + 
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank())

f1.3 <- ggdendrogram(hc, rotate = FALSE, theme_dendro = TRUE) +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())


## Create Transcriptomics Figure from Paper (Figure 2) -------------

# Read in enrichment results 
apap6_enrich <- read.xlsx("Rawls_Supplementary_data2.xlsx", sheet = 2) %>% as.tbl %>% 
  select(Term,Count,Genes,Fold.Enrichment,FDR) %>% 
  mutate(Term = gsub("^.*:","",Term)) %>% 
  rename(Enrichment = "Fold.Enrichment") %>% 
  filter(FDR < 0.1)

ccl46_enrich <-  read.xlsx("Rawls_Supplementary_data2.xlsx", sheet = 3) %>% as.tbl %>% 
  select(Term,Count,Genes,Fold.Enrichment,FDR) %>% 
  mutate(Term = gsub("^.*:","",Term)) %>% 
  rename(Enrichment = "Fold.Enrichment") %>% 
  filter(FDR < 0.1)

tcdd6_enrich <- read.xlsx("Rawls_Supplementary_data2.xlsx", sheet = 4) %>% as.tbl %>% 
  select(Term,Count,Genes,Fold.Enrichment,FDR) %>% 
  mutate(Term = gsub("^.*:","",Term)) %>% 
  rename(Enrichment = "Fold.Enrichment") %>% 
  filter(FDR < 0.1)

tce6_enrich <- read.xlsx("Rawls_Supplementary_data2.xlsx", sheet = 5) %>% as.tbl %>% 
  select(Term,Count,Genes,Fold.Enrichment,FDR) %>% 
  mutate(Term = gsub("^.*:","",Term)) %>% 
  rename(Enrichment = "Fold.Enrichment") %>% 
  filter(FDR < 0.1)


## Shorten and Rename Enrichment Categories 
apap6_enrich$Term[3] <- "NAFLD"
apap6_enrich$Term[10] <- "Gly, Ser, and Thr Metabolism"
apap6_enrich$Term[12] <- "Ubiquitin med proteolysis"
apap6_enrich$Term[14] <- "Comp and coag cascades"
apap6_enrich$Term[9] <- "Protein processing in ER"
ccl46_enrich$Term[1] <- "Protein processing in ER"
tcdd6_enrich$Term[3] <- "Protein processing in ER"
tce6_enrich$Term[2]  <- "Protein processing in ER"

apap6_enrich$Drug <- rep("APAP",length(apap6_enrich$Term))
ccl46_enrich$Drug <- rep("CCl4",length(ccl46_enrich$Term))
tcdd6_enrich$Drug <- rep("TCDD",length(tcdd6_enrich$Term))
tce6_enrich$Drug <- rep("TCE",length(tce6_enrich$Term))

mega_enrich <- rbind(apap6_enrich,ccl46_enrich,tcdd6_enrich,tce6_enrich)

f2.1 <- mega_enrich %>% 
  ggplot(aes(x = reorder(Term,Count), y = Count)) + 
  geom_bar(stat = "identity",color = "black",aes(fill = Drug),position = "dodge") + 
  theme_classic() + coord_flip() +
  ggtitle("Enrichment of KEGG Pathways") + 
  ylab("Number of Genes") + xlab("KEGG Pathways") + 
  scale_fill_manual("legend",values = 
      c("APAP" = "grey 90","CCl4" = "grey 60","TCDD" = "grey 30","TCE" = "black"),
      labels = c("APAP",expression(paste("CCl"[4]*"")),"TCDD","TCE")) +
  theme(legend.text.align = 0) +
  theme(legend.position = c(0.7,0.4)) +
  theme(legend.text = element_text(family = "Arial", size = 30)) +
  theme(legend.title = element_blank()) + 
  theme(plot.title = element_text(family = "Arial",face = "bold",size = 32)) +
  theme(axis.text = element_text(family = "Arial",face = "bold",size = 30)) + 
  theme(axis.title = element_text(family = "Arial",face = "bold",size = 32)) 


figure2 <- ggarrange(f2.1,f1.1, nrow = 1, ncol = 2, labels = c("A","B"), font.label = list(size = 26))
figure2
ggsave(figure2,filename = "Figure2_Transcriptomicsfigure.png",width = 28,height = 20)

