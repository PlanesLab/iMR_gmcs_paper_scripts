% function calculate_gMCSs_animalGEM_function_hpc(name)
clear; clc; close all;

addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1210\cplex\matlab\x64_win64\')
addpath(genpath('C:\Users\nbarrenac\Documents\cobratoolbox'))
addpath(genpath('C:\Users\nbarrenac\Documents\RAVEN'))
numWorkers = 0
name_array = {'Human-GEM-1.10.0', 'Human-GEM-1.11.0', 'Human-GEM-1.12.0', 'Human-GEM-1.13.0', 'Human-GEM-1.14.0'};

TT = table(name_array', ...
    zeros(length(name_array),1), zeros(length(name_array),1), ...
    zeros(length(name_array),1), zeros(length(name_array),1), ...
    'VariableNames', {'name', 'f_raw', 'f_simp', 'n_raw', 'n_simp'});

TT_1 = TT;
TT_2 = TT;

for ii = 1:length(name_array)
    ii
    %     for ii = 3:3
    name = name_array{ii};
    
    % Load HumanGem and rebuild rules (one is bad written)
    model2 = readCbModel(['.' filesep 'non_simplified' filesep name '.mat']);
    if ~isfield(model2, 'rev')
        model2.rev = (model2.lb<0)*1;
    end
    
    [model3, deletedDeadEndRxns] = simplifyModel(model2, true, true, true, true, false, false);
    size(model3.S)
    model3 = removeReactions(model2,deletedDeadEndRxns,true,false);
    size(model3.S)
    
    model = model3;

    f2 = FBA(model2, 'maximize', 0);
    f = FBA(model, 'maximize', 0);
    
    TT_1{ii,1} = {name};
    TT_1{ii,2} = f2.objval;
    TT_1{ii,3} = f.objval;
    TT_1{ii,4} = length(model2.rxns);
    TT_1{ii,5} = length(model.rxns);
    
    model = generateRules(model);
    model = buildRxnGeneMat(model);

    model.csense = repelem('E',size(model.S,1),1);
    writeCbModel(model, ['simplified' filesep name '.mat'])

end

