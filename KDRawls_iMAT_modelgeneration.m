% Kristopher Rawls
% 4/25/19
% Code to create iMAT model to use for TIMBR predictions. 

% Start fresh 
clear all
clc

% Read in cobra model from .mat file. Can also read in cobra model file
% from excel
% Model is apated form the iRno model from excel worksheet available
% at https://github.com/csbl/ratcon1 and 
% http://www.nature.com/articles/ncomms14250 from the published paper
% "Reconciled rat and human metabolic networks for comparative 
% toxicogenomics analyses and biomarker predictions" by Blais et al.
% This model has genes specified for exchange reactions to include them for
% TIMBR predictions for each metabolite. 

load('rno_biomass_kdr.mat');

% Change objective function to biomass
rno_biomass_kdr = changeObjective(rno_biomass_kdr,'RCR99999');
[ex_a,ex_b,ex_c] = xlsread('Rawls_Supplementary_data6.xlsx','extra_genes_kdr');

% Read in gene expression data
apap_allgenes = xlsread('Rawls_Supplementary_data1.xlsx',2);
ccl4_allgenes = xlsread('Rawls_Supplementary_data1.xlsx',3);
tcdd_allgenes = xlsread('Rawls_Supplementary_data1.xlsx',4);
tce_allgenes =  xlsread('Rawls_Supplementary_data1.xlsx',5);

% Convert gene ID lists to cell strings instead of numbers
apap6_genes = cellstr(num2str(apap_allgenes(:,11)));
ccl46_genes = cellstr(num2str(ccl4_allgenes(:,11)));
tcdd6_genes = cellstr(num2str(tcdd_allgenes(:,11)));
tce6_genes = cellstr(num2str(tce_allgenes(:,11)));
exc_genes = ex_b;

% Pull out False Discovery Rate values
apap6_fdr = apap_allgenes(:,9);
ccl46_fdr = ccl4_allgenes(:,9);
tcdd6_fdr = tcdd_allgenes(:,9);
tce6_fdr = tce_allgenes(:,9);

% Grab indicies of genes who are differentially expressed. 
apap6_index = find(apap6_fdr <= 0.1);
ccl46_index = find(ccl46_fdr <= 0.1);
tcdd6_index = find(tcdd6_fdr <= 0.1);
tce6_index = find(tce6_fdr <= 0.1);

% Identify differentially expressed genes
apap6_degs = apap6_genes(apap6_index);
ccl46_degs = ccl46_genes(ccl46_index);
tcdd6_degs = tcdd6_genes(tcdd6_index);
tce6_degs = tce6_genes(tce6_index);

% Combine differentially expressed genes into one vector
all_degs = [apap6_degs; ccl46_degs; tcdd6_degs; tce6_degs;];

% Trim leading and trailing white space from gene ID strings
all_degs = strtrim(unique(all_degs));

% Create a list of all genes found in the experiment and trim the white
% space
all_genes = [apap6_genes; ccl46_genes; tcdd6_genes; tce6_genes];
all_genes = strtrim(unique(all_genes));

% Find all metabolic genes that are in the transcriptomics data
[model_genes, model_idx] = intersect(all_genes,rno_biomass_kdr.genes,'stable');

% Create a logic vector to indicate which genes in the model found in the
% experiment are differentially expressed or not.
imat_data = ismember(model_genes,all_degs);

% Update files to ensure exchagne reaction genes and data is present. 
exc_data = ex_a(:,3);
model_genes = [model_genes; exc_genes];
imat_data = [imat_data; exc_data];


% Create the iMAT data structure
imat_struct = [];
imat_struct.Locus = model_genes;
imat_struct.Data = double(imat_data);

% Initialize the cobratoolbox and set the solver (Uses gurobi6). 
initCobraToolbox

% Run the iMAT algorithm 
[hepatocyte_imat, imat_Rxns] = createTissueSpecificModel(rno_biomass_kdr,imat_struct,1,1,[],'iMAT',[],0);

% Save the model and list of reactions
save('hepatocyte_imat_042519.mat','hepatocyte_imat')
xlswrite('Rawls_Supplementary_data6.xlsx',hepatocyte_imat.rxns,'iMAT_Rxns');