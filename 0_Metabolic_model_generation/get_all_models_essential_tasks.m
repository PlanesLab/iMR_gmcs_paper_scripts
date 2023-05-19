close all
clear
clc

addpath(genpath('C:\Users\nbarrenac\Documents\cobratoolbox'))
addpath(genpath('C:\Users\nbarrenac\Documents\RAVEN'))

changeCobraSolver('ibm_cplex','all')
setRavenSolver('cobra')

% checkInstallation

%% Load and prepare the reference GEM

% name = 'Human-GEM-1.12.0'
% name_array = {'Human-GEM-1.10.0', 'Human-GEM-1.11.0', 'Human-GEM-1.12.0', 'Human-GEM-1.13.0', 'Human-GEM-1.14.0'};

name_array = {'Human-GEM-1.14.0'}

for ii = 1:length(name_array)
    
    name = name_array{ii};
    
    addpath(genpath(['.\' name ]))

    load(['.\' name '\model\Human-GEM.mat']);  % loads model as a structure named "ihuman"

    model = ihuman;

    ihuman2 = addBoundaryMets(ihuman);

    %% Detect all the metabolites requiered

    essentialTasks = parseTaskList('./metabolicTasks_Different_growth_media_v3.txt')
    checkTasks(ihuman2, [], true, false, false, essentialTasks);


    %% generate_model_essential_task

    model = generateRules(model);
    model = buildRxnGeneMat(model);

    model_raw = model;

    f_array = zeros(1, length(essentialTasks));
    % TAB = cell(length(essentialTasks), 1);
    taskReport = cell(length(essentialTasks), 1);

    error_array = []

    for i = 1:length(essentialTasks)
        [{i} {essentialTasks(i).id} {essentialTasks(i).description}]
        [model, taskReport{i}] = calculate_model_for_task(model_raw, essentialTasks(i));

        save(['non_simplified' filesep name '.mat'], 'model')
%         f3 = optimizeCbModel(model);
%         f_array(i) = f3.f;
%         if f3.f < 1e-2
%             f3
%         else
%             f3.f
%         end
    end

    tabulate(f_array)


end
