clear
clc

% initCobraToolbox();

modelName1 = 'TRRUST'
% modelName1 = 'Dorothea'
% modelName1 = 'Omnipath'

% nLayers = 0;
nLayers = 1;
% nLayers = 2;
% nLayers = 3;
name = 'Human-GEM-1.14.0'

signalingNetwork = readtable(['../DATA/' modelName1 '.txt']);

modelName = [modelName1 '_' name]

load(['../0_Metabolic_model_generation/simplified/' name '.mat']);

[G, G_ind, related, n_genes_KO, G_time, G_rxns, iMRmodel_info] = buildGmatrix_iMRmodel(modelName, model, signalingNetwork, nLayers);
% save(['G_matrices_non_simplified/G_' model '_Human-GEM-1.14.0_' num2str(layer) '_layers.mat'], 'G', 'G_ind', 'G_rxns', 'related', 'n_genes_KO', 'G_time');

maxlength = 5;
[G, G_ind, related, n_genes_KO] = filterGmatrix(G, G_ind, maxlength, G_time, G_rxns);
save(['G_matrices_simplified/G_' model '_Human-GEM-1.14.0_' num2str(layer) '_layers.mat'], 'G', 'G_ind', 'G_rxns', 'related', 'n_genes_KO', 'G_time');

[gmcs, gmcs_time, gmcs_onlynew] = calculateGeneMCS(modelName, model, 20000, 7, 'timelimit', 300, 'numWorkers', numWorkers, 'forceLength', 0);

