% Kristopher Rawls
% 4/25/19
% Code to create Raw TIMBR Production Scores.
% Uses TIMBR Weights from KDRawls_TIMBR_Weights.R
% Outputs raw scores to go into Rawls_Supplementary_Data4.xlsx
% This code was adapted from blais_timbr_prediction.m 
% Script published with manuscript at https://www.nature.com/articles/ncomms14250
% Original script available at www.github.com/csbl/ratcon1


clear all
close all
clc

% TIMBR software requirements: COBRA toolbox and linear programming solver

% Specify path to COBRA model files in .xls or .sbml format
path_model_directory = 'C:/Users/kr2up/Documents/ncomm_blais_submission_113016/ncomm_blais_supplementary_data_models/';

% Specify output directory for TIMBR predictions
path_timbr_directory = ['C:/Users/kr2up/Documents/MATLAB/cobratoolbox-master/timbr_directory/iMAT_042519/'];

% initCobraToolbox % this only needs to be run once after installation
changeCobraSolver('gurobi6'); % or glpk

% Load rat-specific metabolic models
load('hepatocyte_imat_042519.mat')
rno_cobra_load = hepatocyte_imat;

% % Each model contains the same sets reactions initially.
% % Reactions conisdered "off" have lb and ub equal to 0
 off_in_rno = rno_cobra_load.rxns(...
     rno_cobra_load.lb == 0 & rno_cobra_load.ub == 0);

% % Remove reactions disabled in each model and set the objective to biomass
 rno_cobra = removeRxns(rno_cobra_load,off_in_rno);

% Convert COBRA models into irreversible format
% Specify default production requirements for all TIMBR simulations
% Originally, we required that both models produced biomass at a
% growth rate of 1 / week with an ATP maintenance requirement.
% However, we found that it was too difficult to resolve 
% how much TIMBR production scores reflected the global demand for 
% biomass production vs biomarker production.


% Set unconstrained reaction bounds to +/- 10^6 because the largest
% consumption rate is fairly close to the default 10^3 reaction bound.
default_bound = 100000; % Changed from 10^6 to 10^5 KDR
rno_irrev_base = ncomm_blais_model2irrev(rno_cobra);

default_rxns_rno = rno_irrev_base.rxns(rno_irrev_base.ub == 1000);
rno_irrev_default = changeRxnBounds(rno_irrev_base, ...
    default_rxns_rno, default_bound, 'u');

% Didn't use biomass or ATP requirement
% obj_value_required = 1 / 7 / 24;
% atp_value_required = 5000;
% Set unconstrained reaction bounds to +/- 10^6 because the largest
% consumption rate is fairly close to the default 10^3 reaction bound.
default_bound = 100000; % Changed from 10^6 to 10^5
rno_irrev_base = ncomm_blais_model2irrev(rno_cobra);

% rno_irrev_atp = changeRxnBounds(rno_irrev_default,...
%     'RCR11017_f',atp_value_required,'l');
% 
% rno_irrev_biomass = changeRxnBounds(rno_irrev_atp,...
%     'RCR99999_f',obj_value_required,'l');


% Note: as mentioned above, TIMBR predictions were not performed while 
% requiring biomass production.
rno_irrev = rno_irrev_default;


% Add metabolomics data constraints, to consume metabolites being consumed - KDR
consumed_mets = {'RCR30005_f','RCR30015_f','RCR30016_f','RCR30031_f','RCR30055_f','RCR30095_f',...
   'RCR30117_f','RCR30118_f','RCR30145_f','RCR30184_f','RCR30185_f','RCR30233_f','RCR30327_f',...
   'RCR30528_f','RCR30529_f','RCR30531_f','RCR30533_f','RCR30555_f'};
    
% List of consumed metabolites in the model: 
% Arginine, Glutamate, Citrulline, Tyrosine, Cystine, Cysteine
% Proline, Choline, Phenylalanine, Valine, Nicotinamide, Glycerol
% Isoleucine, Leucine, Methionine, Phosphate, Succinate
rno_irrev = changeRxnBounds(rno_irrev,consumed_mets,0,'b'); 

% Determine maxmimum flux through each exchange reaction for iRno
[rno_exchange, rno_demand] = findExcRxns(rno_irrev);
rno_production = rno_irrev.rxns((cellfun(@length, ...
    regexp(rno_irrev.rxns,'_r')) == 0) & rno_exchange);
rno_fva = zeros(length(rno_production),1);

for exchange_index = 1:length(rno_production)
    rno_consumption = setdiff(intersect(...
        regexprep(rno_production(exchange_index),'_f','_r'),...
        rno_irrev.rxns),...
        rno_production(exchange_index));
    if (~isempty(rno_consumption))
        rno_production_fba = optimizeCbModel(changeObjective(...
            changeRxnBounds(changeRxnBounds(...
            rno_irrev,rno_production,default_bound,'u'),...
            rno_consumption,0,'b'),...
            rno_production(exchange_index)));
        rno_fva(exchange_index,1) = rno_production_fba.f;
    else 
        rno_production_fba = optimizeCbModel(changeObjective(...
            changeRxnBounds(...
            rno_irrev,rno_production,default_bound,'u'),...
            rno_production(exchange_index)));
        rno_fva(exchange_index,1) = rno_production_fba.f;
    end
    disp(rno_production_fba)
end

% Specify minimimum required flux for each exchange reaction
fba_threshold = 1e-4;

rno_production_ok = rno_fva > fba_threshold;
rno_production_id = rno_production(rno_production_ok,1);
rno_production_requirement = min(0.90 * rno_fva(rno_production_ok,1), 100);
rno_production_count = length(rno_production_id);

% Load TIMBR reaction weights
% These files are generated by the R script:
% KDRawls_TIMBR_Weights.R

rno_timbr_weights_load = readtable([...
    'kdr_hep_expression_042519_imat_timbr_weights_rno.txt'],'Delimiter','\t');

[rno_model_has_weight, rno_model2weight_index] = ismember(...
    rno_irrev.rxns, rno_timbr_weights_load.rxn_irreversible);
[rno_weight_in_model, rno_weight2model_index] = ismember(...
    rno_timbr_weights_load.rxn_irreversible, rno_irrev.rxns);

if (all(rno_model_has_weight))
    % first 3 columns are organism_id, rxn_id, and rxn_irrerversible
    rno_timbr_weights = rno_timbr_weights_load(rno_model2weight_index,4:end);
    rno_timbr_id = rno_timbr_weights.Properties.VariableNames';
    rno_timbr_count = length(rno_timbr_id);

else
    warning('reaction weights not specified for all reactions in rno')
end


% Run TIMBR algorithm and save results as .txt files
% Estimate the global network demand of producing each exchange metabolite
% under treatment and control conditions for various compounds. 
% Saved outputs will be analyzed in the R/Bioconductor script:
% KDRawls_TIMBR_Analysis.R

for rno_timbr_index = 1:rno_timbr_count
    
    rno_timbr_network_demand = zeros(rno_timbr_count, 1);
    disp(rno_timbr_id{rno_timbr_index});
    for rno_production_index = 1:rno_production_count
        disp(rno_production_id(rno_production_index,1));
        rno_timbr_network_demand(rno_production_index,1) = ...
            timbr(rno_irrev, ...
            rno_production_id(rno_production_index,1), ...
            rno_production_requirement(rno_production_index,1), ...
            rno_timbr_weights(:,rno_timbr_index));
    end
    rno_timbr_file_name = [path_timbr_directory ...
        'ncomm_blais_timbr_' rno_timbr_id{rno_timbr_index} '.txt'];
    rno_timbr_table =  table(...
        rno_production_id,...
        rno_timbr_network_demand,...
        'VariableNames',{'timbr_id' 'timbr_value'});
    writetable(rno_timbr_table,rno_timbr_file_name,'Delimiter','\t');
end



