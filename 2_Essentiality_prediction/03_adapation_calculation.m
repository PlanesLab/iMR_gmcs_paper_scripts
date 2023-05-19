clear
clc

% initCobraToolbox()
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1210\cplex\matlab\x64_win64')
addpath(genpath('C:\Users\nbarrenac\Documents\cobratoolbox'))
addpath(genpath('G:\My Drive\PhD - Signaling_gMCSs\Code\2022-12-07_Functions_to_GitHub'))  
% modelName_array = {'TRRUST', 'Omnipath', 'Dorothea'};

modelName_array = {'TRRUST', 'Omnipath', 'Dorothea'};
modelName_array_2 = {'Human1-T1', 'Human1-T2', 'Human1-T3',...
    'Human1-O1', 'Human1-O2', 'Human1-O3', ... 
    'Human1-D1', 'Human1-D2', 'Human1-D3'};

% modelName_array = {'O_D_union', 'T_D_union','T_O_D_union', 'T_O_union', 'O_D_intersection', 'T_D_intersection','T_O_D_intersection', 'T_O_intersection'};
% modelName_array_2 = {'Human1-D1UO1', 'Human1-D1UT1','Human1-D1UO1UT1', 'Human1-O1UT1', 'Human1-D1?O1', 'Human1-D1?T1','Human1-D1?O1?T1', 'Human1-O1?T1'};
cont = 0;
maxlength = 5;
for count_model = 1:length(modelName_array)
    modelName = modelName_array{count_model}
    signalingNetwork = readtable(['C:\Users\nbarrenac\OneDrive - Tecnun\2023-03-24_Signaling_Major_revision\Update_reg_networks\' modelName '.txt']);

    load('C:\Users\nbarrenac\OneDrive - Tecnun\2023-03-24_Signaling_Major_revision\G_matrix_comparison_Human_models\1-Generate_models\simplified_v2\Human-GEM-1.14.0.mat')
    for nLayers = 1:1
%     for nLayers = 1:3
        cont = cont+1;
        [reason, layersPerReaction, layerGenes, time, AllBooleanRules] = addSingnalingLayers(model, nLayers, signalingNetwork, modelName);
        gMCSs_data = readtable("../calculated_gMCSs_of_each_model/gMCSs_list.xlsx", "Sheet","Human1-T1");
        [nf, nc] = size(gMCSs_data);
        gmcs = {};
        for j = 1:nf
            data = table2cell(gMCSs_data(j,:));
            data = data(find(~cellfun(@isempty,data)));
            gmcs{j,1} = data';
        end
        lengths = cell2mat(cellfun(@length, gmcs, 'UniformOutput', false));
        position = find((lengths > 1 & lengths < (maxlength+1)));

        gMCSs_all = gmcs(position,:);
        AllGenes = unique(vertcat(layerGenes{:}));
        target = strseq('target', 1:length(AllBooleanRules));
        AllGenes = vertcat(AllGenes, target);

        AllBooleanRules = AllBooleanRules(~cellfun(@isempty, AllBooleanRules));


        for j = 1:length(AllBooleanRules)
            [RxnFormulas, auxKO_Nodes, inputs] = generateFormulas(AllBooleanRules{j}, j); 
            RxnFormulas_all{j} = RxnFormulas;
            auxKO_Nodes_all{j} = auxKO_Nodes;
            inputs_all{j} = inputs;

        end


        RxnFormulas_cat = vertcat(RxnFormulas_all{:});
        auxKO_Nodes_cat = unique(vertcat(auxKO_Nodes_all{:}));
        inputs_cat = unique(vertcat(inputs_all{:}));

        all_nodes = cellfun(@unique, regexp(RxnFormulas_cat, '([^\(\)\!\+\|\&\=\s]+)', 'match'), 'Uniform', 0);
        all_nodes = unique([all_nodes{~cellfun(@isempty,all_nodes)}])';
        all_nodes(ismember(all_nodes,'->')) = [];

        original_nodes = all_nodes(~contains(all_nodes, '_auxKO'));

        [~, working_pos] = ismember(strtrim(extractBefore(auxKO_Nodes_cat, '_auxKO')), original_nodes);

        all_nodes = [original_nodes; auxKO_Nodes_cat];

        [~, working_pos2] = ismember(auxKO_Nodes_cat, all_nodes);
        [~, input_pos] = ismember(strcat(inputs_cat, '_input'), all_nodes);

        a = 1;
        for pos = 1:length(AllBooleanRules)
            RxnFormulas = RxnFormulas_all{pos};
            [~, index_sort] = sort(extractAfter(RxnFormulas, '-> '));
            RxnFormulas = RxnFormulas(index_sort);
            gene = extractAfter(RxnFormulas{1}, '-> ');
            genes_equal = extractAfter(RxnFormulas, '-> ');

            gene_unique = unique(strtrim(genes_equal));
            times = 1;

            or_nodes = cell(length(RxnFormulas),2);
            for i=2:length(RxnFormulas)
                if strcmp(extractAfter(RxnFormulas{i}, '-> '), gene)
                    times = times + 1;
                    or_nodes{i-1,1} = gene;
                    if i == length(RxnFormulas)
                        or_nodes{i,2} = times;
                        or_nodes {i,1} = gene;
                    end
                else
                    or_nodes{i-1,2} = times;
                    or_nodes{i-1,1} = gene;
                    times = 1;
                    gene = extractAfter(RxnFormulas{i},'-> ');
                    if i == length(RxnFormulas)
                        or_nodes{i,2} = times;
                        or_nodes{i,1} = gene;
                    end
                end
            end

            or_nodes(:,1) = strtrim(or_nodes(:,1));

            and_nodes = cellfun(@unique, regexp(extractBefore(RxnFormulas,' -> '), '([^\(\)\+\|\&\=\s]+)', 'match'), 'Uniform', 0);

            for i=1:length(RxnFormulas)
                nodes = cell(1);
                if ~isempty(or_nodes{i,2})
                    GPR_node = extractAfter(RxnFormulas{i,1}, '-> ');
                    GPR_node = regexprep(GPR_node, ' ' ,'');
                    Aux_nodes{a,1} = find(strcmp(GPR_node, all_nodes));
                    GPR_and = extractBefore(RxnFormulas{i,1}, ' -> ');
                    nodes = unique(regexp(GPR_and, '([^\(\)\+\|\&\=\s]+)', 'match'));
                    if ~isempty(find(contains(RxnFormulas{i,1}, '+'), 1))
                        Aux_nodes{a,2} = 'AND';
                        for j=1:length(nodes)
                            if strfind(nodes{j}, '!')
                                nodes{j} = regexprep(nodes{j}, '!', '');
                                Aux_nodes{a,4}{j} = 1;
                            else
                                Aux_nodes{a,4}{j} = 0;
                            end
                            Aux_nodes{a,3}{j} = find(strcmp(nodes{j}, all_nodes));
                        end
                    else
                        Aux_nodes{a,2} = 'OR';
                        if ~isempty(or_nodes{i,2})
                            for j = 1:or_nodes{i,2}
                                if contains(and_nodes{i-j+1,1}, '!')
                                    nodes{j} = regexprep(and_nodes{i-j+1}, '!', '');
                                    Aux_nodes{a,4}{j} = 1;
                                else
                                    nodes{j} = and_nodes{i-j+1};
                                    Aux_nodes{a,4}{j} = 0;
                                end
                                Aux_nodes{a,3}{j} = find(strcmp(all_nodes, nodes{j}));
                            end
                        end
                    end
                    a = a+1;
                end
            end

        end

        new_variables = length(original_nodes);

        vars.y = 1:length(all_nodes);
        vars_n = vars.y(end);

        cons = struct();
        cons.Eq1 = 1:length([Aux_nodes{:,3}, Aux_nodes{:,2}]);
        cons_n = cons.Eq1(end);

        A = spalloc(cons_n, vars_n, 5*cons_n);
        rhs = zeros(cons_n, 1);
        lhs = zeros(cons_n, 1);

        a = 1;
        for i = 1:size(Aux_nodes,1)
            if strcmp(Aux_nodes{i,2}, 'AND')
                A(a+length(Aux_nodes{i,3}),Aux_nodes{i,1}) = -1;
                rhs(a+length(Aux_nodes{i,3}),1) = length(Aux_nodes{i,3})-1;
                lhs(a+length(Aux_nodes{i,3}),1) = -inf;

                idx = a-1 + (1:length(Aux_nodes{i,3}));
                A(idx, Aux_nodes{i,1})=1;
                idx2 = [Aux_nodes{i,4}{:}]==1;
                % j = 1
                A(a+length(Aux_nodes{i,3}), [Aux_nodes{i,3}{idx2}]) = -1;
                rhs(a+length(Aux_nodes{i,3}),1) = rhs(a+length(Aux_nodes{i,3}),1) - 1*(sum(idx2));
                lhs(a+length(Aux_nodes{i,3}),1) = -inf;
                A(idx(idx2), [Aux_nodes{i,3}{idx2}]) = eye(sum(idx2));
                rhs(idx(idx2),1) = 1;
                lhs(idx(idx2),1)= -inf;
                % j = 0
                A(a+length(Aux_nodes{i,3}), [Aux_nodes{i,3}{~idx2}]) = 1;
                A(idx(~idx2), [Aux_nodes{i,3}{~idx2}]) = -eye(sum(~idx2));
                rhs(idx(~idx2),1) = 0;
                lhs(idx(~idx2),1) = -inf;

                a = a+length(Aux_nodes{i,3})+1;
            else
                A(a+length(Aux_nodes{i,3}),Aux_nodes{i,1}) = 1;
                rhs(a+length(Aux_nodes{i,3}),1) = 0;
                lhs(a+length(Aux_nodes{i,3}),1) = -inf;

                idx = a-1 + (1:length(Aux_nodes{i,3}));
                A(idx,Aux_nodes{i,1}) = -1;
                idx2 = [Aux_nodes{i,4}{:}]==1;
                % j = 1
                A(a+length(Aux_nodes{i,3}), [Aux_nodes{i,3}{idx2}]) =  1;
                rhs(a+length(Aux_nodes{i,3}),1) = rhs(a+length(Aux_nodes{i,3}),1)+1*(sum(idx2));
                lhs(a+length(Aux_nodes{i,3}),1) = -inf;
                A(idx(idx2),[Aux_nodes{i,3}{idx2}]) = -eye(sum(idx2));
                rhs(idx(idx2),1) = -1;
                lhs(idx(idx2),1)= -inf;
                % j = 0
                A(a+length(Aux_nodes{i,3}), [Aux_nodes{i,3}{~idx2}]) = - 1;
                A(idx(~idx2), [Aux_nodes{i,3}{~idx2}]) =  + eye(sum(~idx2));
                rhs(idx(~idx2),1) = 0;
                lhs(idx(~idx2),1) = -inf;

                a = a+length(Aux_nodes{i,3})+1;
            end
        end
        [cons_n, n_vars] = size(A);

        %ub
        ub = zeros(n_vars, 1);
        ub(vars.y) = 1;

        %lb
        lb = zeros(n_vars, 1);
        lb(vars.y) = 0;
        lb(input_pos) = 0;
        lb(working_pos2) = 1;

        ctype = '';
        ctype(vars.y) = 'B';

        obj = zeros(1, n_vars);
        sense = 'minimize';
        pos = 1;
        final_array = {};

        for i = 1:length(gMCSs_all)
            gMCSs = gMCSs_all{i,:};

            for j = 1:length(gMCSs)
                rest = gMCSs;
                rest(j) = [];

                obj = zeros(1, n_vars);

                ub = zeros(n_vars, 1);
                ub(vars.y) = 1;

                %lb
                lb = zeros(n_vars, 1);
                lb(vars.y) = 0;
                lb(input_pos) = 0;
                lb(working_pos2) = 1;

                ub(find(strcmp(all_nodes, gMCSs{j}))) = 0;
                lb(find(strcmp(all_nodes, gMCSs{j}))) = 0;
                for k = 1:length(rest)
                    obj(find(strcmp(all_nodes, rest{k}))) = 1;
                end

                cplex = Cplex('Prueba');
                Model = struct();

                [Model.A, Model.rhs, Model.lhs, Model.ub, Model.lb, Model.obj, Model.ctype, Model.sense] = deal(A, rhs, lhs, ub, lb, obj, ctype, sense);
                cplex.Model = Model;
                cplex.Param.output.clonelog.Cur = 0;
                cplex.solve()

                if ismember(cplex.Solution.status,[1, 101]) 
                    if cplex.Solution.bestobjval > 0
                        final_array{pos,1} = gMCSs;
                        final_array{pos,2} = rest;
                        final_array{pos,3} = {};
                        solution = cplex.Solution.x(vars.y);
                        for k = 1:length(rest)
                            final_array{pos,3} = [final_array{pos,3}, solution(find(strcmp(all_nodes, rest{k})))];
                        end
                        pos = pos + 1;
                    end
                end


            end


        end
%         save(['Results_relations_' model '_' num2str(nLayer) '.mat'], 'final_array')
        
        if length(final_array)>0
            lengths = cell2mat(cellfun(@length, final_array(:,1), 'UniformOutput', false));
            position = find((lengths > 1 & lengths <=maxlength));
            final_array_final = final_array(position,:);

         
            for i = 1:length(final_array_final)
                KO_genes{i,1} = setdiff(final_array_final{i,1}, final_array_final{i,2});
            end

            gMCSs2file = cell(length(final_array_final),length(final_array_final{length(final_array_final)}));
            gMCSs2file(:) = {''};
            rest2file = cell(1);
            for i = 1:length(final_array_final)
                gmcs = final_array_final{i,1};
                rest = final_array_final{i,2};
                pos = find(cell2mat(final_array_final{i,3}));
                rest2file{i,:} = rest{pos};
                for j = 1:length(gmcs)
                    gMCSs2file(i,j) = gmcs(j);
                end
            end

            writetable(cell2table(gMCSs2file),['adaptation_txt/Results_gMCSs_' modelName_array_2{cont} '_' num2str(maxlength) '.txt'],'FileType','text',...
                'WriteVariableNames',false,'Delimiter',';');

            writetable(cell2table(KO_genes),['adaptation_txt/Results_KO_' modelName_array_2{cont} '_' num2str(maxlength) '.txt'],'FileType','text',...
                'WriteVariableNames',false,'Delimiter',';');

            writetable(cell2table(rest2file),['adaptation_txt/Results_rest_' modelName_array_2{cont} '_' num2str(maxlength) '.txt'],'FileType','text',...
                'WriteVariableNames',false,'Delimiter',';');
        end
        
    end
end

function [reason, layersPerReaction, layerGenes, time, AllBooleanRules] = addSingnalingLayers(ihuman, n_layers, signalingNetwork, model_name)
tic

% Unique GPR rules are calculated. The sgMCSs are computed once per each
% unique GPR rule.
[uniqueGPR, ~, posGPR] = unique(ihuman.grRules, 'stable');
metGenes = cellfun(@unique, regexp(uniqueGPR, '([^\(\)\+\!\|\&\=\s\or\(and)]+)', 'match'), 'Uniform', 0);

layersPerReaction = zeros(1, length(uniqueGPR));
reason = zeros(1, length(uniqueGPR));
nerwork_final = cell(1, length(uniqueGPR));

for y = 1:length(uniqueGPR)
    disp(y)
    if uniqueGPR(y)~=""
        network_2 = {};
        reactions = {};
        abr = {};
        names = {};
        BooleanRules = {};
        BooleanRules(1,1) = strcat('target', num2str(y), ' = ' ,{' '}, uniqueGPR(y));
        if  ~isempty(find(ismember(table2cell(signalingNetwork(:,2)), metGenes{y,1})))
            for layer = 1:n_layers
                network_pre = network_2;
                if layer == 1
                    genes_sel = metGenes{y,1};
                    new_network = table2cell(signalingNetwork(ismember(table2cell(signalingNetwork(:,2)), genes_sel),:));
                else
                    genes_sel = network_2(:,1);
                    new_network = table2cell(signalingNetwork(ismember(table2cell(signalingNetwork(:,2)), genes_sel),:));
                end
                network_2 = [network_2; new_network];
                network_2_new = cell2table(network_2);
                network_2 = table2cell(unique(network_2_new, 'rows'));
                if size(network_pre,1) == size(network_2,1)
                    layersPerReaction(y) = layer - 1;
                    if size(network_2,1) == 1
                        layerGenes{y} = unique([unique([network_2(:,1) network_2(:,2)]'); metGenes{y,1}']);
                    else
                        layerGenes{y} =  unique([unique([network_2(:,1) network_2(:,2)]); metGenes{y,1}']);
                    end
                    reason(y) = 1;
                    network_2 = network_pre;
                    break
                else
                    layersPerReaction(y) = layer;
                end
                for i = 1:size(network_2,1)
                    reactions(i,1) = strcat(network_2(i,1), ' -> ', {' '}, network_2(i,2));
                    abr{i,1} = ['R' num2str(i)];
                    names{i,1} = ['Reaction_' num2str(i)];
                end

                %Check of the presence of oscilations
                submodel = createModel(abr, names, reactions, 'printLevel', 0);
                S = submodel.S;
                [nrow, ncol] = size(S);

                A = sparse(nrow+1, ncol);
                A(1:nrow,:) = S;
                A(nrow+1,:) = ones(1, ncol);

                rhs = zeros(nrow+1,1);
                lhs = zeros(nrow+1,1);
                rhs(1:nrow,1) = 0;
                lhs(1:nrow,1)= 0;
                rhs(nrow+1,1) = inf;
                lhs(nrow+1,1)= 1;

                ub = zeros(ncol, 1);
                lb = zeros(ncol, 1);

                ub(1:ncol,1) = inf;
                lb(1:ncol,1) = 0;
                ctype = '';
                ctype(1:ncol) = 'C';
                obj = ones(1, ncol);
                sense = 'minimize';
                cplex = Cplex('Prueba');
                Model = struct();
                [Model.A, Model.rhs, Model.lhs, Model.ub, Model.lb, Model.obj, Model.ctype, Model.sense] = deal(A, rhs, lhs, ub, lb, obj, ctype, sense);
                cplex.Model = Model;
                cplex.Param.threads.Cur = 12;
                cplex.DisplayFunc = [];

                cplex.solve()
                if cplex.Solution.status~=103
                    layersPerReaction(y)=layer-1;
                    if size(network_pre,1)==1
                        layerGenes{y} = unique([unique([network_pre(:,1) network_pre(:,2)]'); metGenes{y,1}']);
                    elseif size(network_pre,1)==0
                        layerGenes{y} = {};
                    else
                        layerGenes{y} =  unique([unique([network_pre(:,1) network_pre(:,2)]); metGenes{y,1}']);
                    end
                    reason(y) = 2;
                    network_2 = network_pre;
                    break
                else
                    layersPerReaction(y) = layer;
                end

            end
        else
            reason(y) = 3;
        end
        if size(network_2,1)==1
            layerGenes{y} = unique([unique([network_2(:,1) network_2(:,2)]'); metGenes{y,1}']);
        elseif size(network_2,1)==0 && ~isempty(metGenes{y,1})
            layerGenes{y} = metGenes{y,1}';
        else
            layerGenes{y} =  unique([unique([network_2(:,1) network_2(:,2)]); metGenes{y,1}']);
        end
        nerwork_final{y} = network_2;
        for i = 1:size(network_2,1)
            if cell2mat(network_2(i,3)) == 1
                BooleanRules(i+1,1) = strcat(network_2(i,2), ' = ', {' '}, network_2(i,1));
            elseif cell2mat(network_2(i,3)) == -1
                BooleanRules(i+1,1) = strcat(network_2(i,2), ' = !', network_2(i,1));
            end
        end
        AllBooleanRules{y} = BooleanRules;
    end 
        
end
time = toc;

% final_filename = [pwd filesep 'GPRs_' model_name '_' num2str(n_layers) '_layer_unique.mat'];
% save(final_filename, 'reason', 'layersPerReaction', 'layerGenes', 'time', 'AllBooleanRules', 'nerwork_final');

end

function reactions = modelParser_aux(reactions, head, nodeParent, parent)
% Name nodes according to position and layer inside the Boolean rule

node = cell(size(head.children));
if length(head.children) == 1 %If it has an only child, it is named as its parent.
    node{1} = parent;
else %If the node has more than one child, they will be named as : nodeparent_1, nodeparent_2...
    node = strcat(nodeParent,'_', strsplit(num2str(1:length(head.children))));
    if isa(head,'OrNode') %if it is an OrNode:
        for child = 1:length(head.children)
            reactions.rxns{end+1} = node{child};
            reactions.formulas{end+1} = [node{child},' -> ',parent];
        end
    else %if it is an AndNode
        reactions.rxns{end+1} = strjoin(node,' & ');
        reactions.formulas{end+1} = [strjoin(node,' + '),' -> ',parent];
    end
end

%  Expand model with every node, using OR/AND rules
for i = 1:length(head.children)
    child = head.children(i);
    % check if the child is final layer
    if isa(child,'LiteralNode')
        reactions.rxns{end+1} = [child.id, '_', node{i}];
        reactions.formulas{end+1} = [child.id, ' -> ', node{i}];
    else %If the child is not the last layer.
        reactions = modelParser_aux(reactions, child,node{i}, node{i});
    end
end

end

function [RxnFormulas, auxKO_Nodes, inputs] = generateFormulas(BooleanRules, pos)
% The Boolean rules of the network are parsed to RxnFormula format.
% Auxiliary nodes and aux nodes (KO nodes) are added to the rules.

BooleanRulesParsed = firstPreParse(BooleanRules); %Boolean rules are preparsed.

% Input nodes are found.
leftNodes = cellfun(@unique, regexp(extractAfter(BooleanRulesParsed,' = '), '([^\(\)\+\!\|\&\=\s]+)', 'match'), 'Uniform', 0);
leftNodes = unique([leftNodes{~cellfun(@isempty,leftNodes)}])';

rigthNodes = cellfun(@unique, regexp(extractBefore(BooleanRulesParsed,' = '), '([^\(\)\+\!\|\&\=\s]+)', 'match'), 'Uniform', 0);
rigthNodes = unique([rigthNodes{~cellfun(@isempty,rigthNodes)}])';

inputs = setdiff(leftNodes, rigthNodes);
outputs = setdiff(rigthNodes, leftNodes);

% Aux KO nodes, which are the ones for doing the knockout of the node,
% are added to each Boolean rule.
[BooleanRulesParsed, auxKO_Nodes] = AddingauxKO_Nodes(BooleanRulesParsed, inputs, outputs);

leftNodes = regexprep(BooleanRulesParsed, ' =.*$', ''); %genes before =

BoolRulesRight = extractAfter(BooleanRulesParsed,'= ');

fp = FormulaParser();
model = cell(length(BoolRulesRight),1);

zz = num2str(length(BoolRulesRight));
for i = 1:length(BoolRulesRight)
%     if printLevel == 1
        disp(['Parsing Boolean ON rules ' num2str(i) ' / ' zz])
%     end
    head = fp.parseFormula(BoolRulesRight{i});
    model{i} = modelParser(head,['GPR_', num2str(pos), '_', leftNodes{i,1}, '_N_1'], [leftNodes{i,1}]);
    model{i} = ReduceModel(model{i});
end

% The RxnFormulas of the model are saved and reduced in order to eliminate
% the redundant nodes.
RxnFormulas = cell(size(model));
for i = 1:length(model)
    RxnFormulas{i} = reshape(model{i},1, []);
end
RxnFormulas = reshape([RxnFormulas{:}],[],1);
RxnFormulas = unique(RxnFormulas);

end

function [preParsedBoolRules] = firstPreParse(BoolRules)
% This function preparses Boolean rules. Then, it moves all the inhibitory
% relationships into the rigth part of the rule and joins all the OR
% rules for the same node in the same line by | operator. 

auxPreParsedBoolRules = regexprep(BoolRules, '[\]\}]',')'); % replace other brackets by parenthesis
auxPreParsedBoolRules = regexprep(auxPreParsedBoolRules, '[\[\{]','('); % replace other brackets by parenthesis
auxPreParsedBoolRules = regexprep(auxPreParsedBoolRules, '*=','='); % replace *= by =
auxPreParsedBoolRules = regexprep(auxPreParsedBoolRules, '"',''); % eliminate "
auxPreParsedBoolRules = regexprep(auxPreParsedBoolRules, '/','_'); % replace / by _
auxPreParsedBoolRules = regexprep(auxPreParsedBoolRules, '-','_'); % replace - by _
auxPreParsedBoolRules = regexprep(auxPreParsedBoolRules, '([\)]\s?|\s)\s*(?i)(and)\s*?(\s?[\(]|\s)\s*', '$1&$3'); % replace all ands
auxPreParsedBoolRules = regexprep(auxPreParsedBoolRules, '([\)]\s?|\s)\s*(?i)(+)\s*?(\s?[\(]|\s)\s*', '$1&$3'); % replace all ands
auxPreParsedBoolRules = regexprep(auxPreParsedBoolRules, '([\)]\s?|\s)\s*(?i)(or)\s*?(\s?[\(]|\s)\s*', '$1|$3'); % replace all ors
auxPreParsedBoolRules = regexprep(auxPreParsedBoolRules, '([\)]\s?|\s)\s*(?i)(not)\s*?(\s?[\(]|\s)\s*', '$1!$3'); % replace all nots
auxPreParsedBoolRules = regexprep(auxPreParsedBoolRules, '[\s]?&[\s]?', ' & '); % introduce spaces around ands
auxPreParsedBoolRules = regexprep(auxPreParsedBoolRules, '[\s]?\|[\s]?', ' | '); % introduce spaces around ors
auxPreParsedBoolRules = regexprep(auxPreParsedBoolRules, '! ', '!'); % eliminate spaces after !

notGenes = find(contains(extractBefore(auxPreParsedBoolRules,' ='),'!')); % nodes in the left side that have !.
% For those nodes, the ! is moved to the other side and parenthesis are
% added to the expression, so as ! affects to all the rule and not only to a
% unique node.
% Example:              !A = B | C ==> A = !(B | C)

if ~isempty(notGenes)
    leftGenes = extractBefore(auxPreParsedBoolRules,' =');
    leftGenes = regexprep(leftGenes, '!', '');
    BoolRuleRight = extractAfter(auxPreParsedBoolRules,'= ');
    for i = 1:length(notGenes)
        BoolRuleRight{notGenes(i),1} = insertBefore(BoolRuleRight{notGenes(i),1}, 1, '!(');
        BoolRuleRight{notGenes(i),1} = insertAfter(BoolRuleRight{notGenes(i),1}, BoolRuleRight{notGenes(i),1}(end), ')');
    end
    auxPreParsedBoolRules = strcat(leftGenes, ' = ', BoolRuleRight);
    auxPreParsedBoolRules = insertAfter(auxPreParsedBoolRules, '=', ' ');
end

% The rules are sorted to see if different OR rules are written in different
% lines.
auxPreParsedBoolRules = sort(auxPreParsedBoolRules);
leftGenes = strtrim(extractBefore(auxPreParsedBoolRules, '='));
leftGenesUnique = unique(leftGenes);

preParsedBoolRules = cell(size(leftGenesUnique));
if length(leftGenesUnique)<length(leftGenes)
    ind = 1;
    i = 1;
    while i <= size(auxPreParsedBoolRules, 1) % times a gene is in the 
        % left side of the rule. 
        % if it is more than one time, all its rules are joined in one 
        % line by | operators.
        Index = find(strcmp(strtrim(extractBefore(auxPreParsedBoolRules, '=')),strtrim(extractBefore(auxPreParsedBoolRules{i,1}, '='))));
        preParsedBoolRules(ind,1) = auxPreParsedBoolRules(i,1);
        if length(Index) > 1
            for j = 2:length(Index)
                preParsedBoolRules(ind,1) = strcat(preParsedBoolRules{ind,1}, {' '}, '|', extractAfter(auxPreParsedBoolRules{Index(j),1}, '='));
                i = i + 1;
            end
        end
        i = i + 1;
        ind = ind + 1;
    end
else
    preParsedBoolRules = auxPreParsedBoolRules;
end
end

function [BooleanRules, auxKO_Node] = AddingauxKO_Nodes(BooleanRules, inputs, outputs)
% This function adds the corresponding aux node to each node. It makes
% possible the knockout of the nodes.
% In the case of the input nodes of the model, in order to facilitate the
% indexing of the aux nodes, it is added a new rule for each input.
%

if ~isempty(inputs)
    if iscell(inputs{1})
        inputs = reshape([inputs{:}],[],1);
    else
        inputs = reshape(inputs,[],1);
    end
    BooleanRulesInput = strcat(inputs, ' = ', {' '} , inputs, '_input');
%     auxKO_NodeInput = strcat(inputs, '_input');
    BooleanRules = vertcat(BooleanRules, BooleanRulesInput);
end

BooleanRules2 = BooleanRules(~strcmp(regexprep(BooleanRules, ' =.*$', ''), outputs));
BooleanRules_out = BooleanRules(strcmp(regexprep(BooleanRules, ' =.*$', ''), outputs));
auxKO_Node = strcat(regexprep(BooleanRules2, ' =.*$', ''), '_auxKO');
BooleanRules_new = strcat(regexprep(BooleanRules2 ,'= ', '= ('), ') & ', {' '}, auxKO_Node);

BooleanRules = [BooleanRules_out; BooleanRules_new];
auxKO_Node = unique(sort([auxKO_Node]));

end

function model = modelParser(head, nodeParent, parent)
% This function builds auxiliary nodes for each layer of the Boolean rules.
% For example:
%     g1 & g2 & g3 is one layer in the model,   RXN = g1 + g2 + g3
%     (g1 | g2) & g3 is a two layer model:      RXN = (g1 | g2) + g3
%                                               (g1 | g2) = g1
%                                               (g1 | g2) = g2

reactions = struct();
reactions.rxns = {'aux'; 'aux'};
reactions.formulas = {'aux'; 'aux'};
reactions = modelParser_aux(reactions, head, nodeParent, parent);

reactions.rxns = reactions.rxns(3:end);
reactions.formulas = reactions.formulas(3:end);

model = reactions.formulas;

end

function [rxnFormulasReduced] = ReduceModel(rxnFormulas)
% As the algorithm generates redundant nodes, this function is used to
% remove these nodes.

LeftHand = regexprep(rxnFormulas, ' ->.*$', '');
RightHand = regexprep(rxnFormulas, '^.* -> ', '');

idx = contains(LeftHand, '+');
LeftHandCell(~idx) = cellfun(@(x)({{x}}), LeftHand(~idx));
if sum(idx)>0
    LeftHandCell(idx) = cellfun(@strsplit, LeftHand(idx), repmat({' + '}, sum(idx),1), 'UniformOutput', 0);
end

idx = contains(RightHand, '+');
RightHandCell(~idx) = cellfun(@(x)({{x}}), RightHand(~idx));
if sum(idx)>0
    RightHandCell(idx) = cellfun(@strsplit, RightHand(idx), repmat({' + '}, sum(idx),1), 'UniformOutput', 0);
end

nodes = [reshape(LeftHandCell,1,[]), reshape(RightHandCell,1,[])];
nodes = reshape(nodes,1,[]);
nodes = reshape(unique([nodes{:}]),[],1);

intermediateNodesRaw = nodes(startsWith(nodes, 'GPR_'));
intermediateNodes = false(size(intermediateNodesRaw));

for i = 1:length(intermediateNodes)
    if sum(cellfun(@length, strfind(RightHand, intermediateNodesRaw{i}))~=0)==1
        if sum(cellfun(@length, strfind(LeftHand, intermediateNodesRaw{i}))~=0)==1
            intermediateNodes(i) = 1;
        end
    end
end
intermediateNodes = intermediateNodesRaw(intermediateNodes);

% generate synonims
intermediateNodesTranslation = cell(numel(intermediateNodes),2);
intermediateNodesTranslation(:,1) = intermediateNodes;
for i = 1:length(intermediateNodes)
    if strcmp(LeftHand,intermediateNodes{i})
        intermediateNodesTranslation{i,2} = RightHand{strcmp(LeftHand, intermediateNodes{i})};
    else
        intermediateNodesTranslation{i,2} = LeftHand{strcmp(RightHand, intermediateNodes{i})};
    end
end

for i = 1:length(intermediateNodes)
    idx = find(cellfun(@length,strfind(LeftHand, intermediateNodesTranslation{i,1}))>0);
    for j = 1:1:length(idx)
        LeftHandCell{idx(j)} = regexprep(LeftHandCell{idx(j)}, intermediateNodesTranslation{i,1}, intermediateNodesTranslation{i,2});
    end
    
    idx = find(cellfun(@length,strfind(RightHand, intermediateNodesTranslation{i,1}))>0);
    for j=1:1:length(idx)
        RightHandCell{idx(j)} = regexprep(RightHandCell{idx(j)}, intermediateNodesTranslation{i,1}, intermediateNodesTranslation{i,2});
    end
end

LeftHandCell = cellfun(@sort, LeftHandCell, 'UniformOutput', 0);
RightHandCell = cellfun(@sort, RightHandCell, 'UniformOutput', 0);

idx = cellfun(@length, LeftHandCell)>1;
LeftHand(~idx) = [LeftHandCell{~idx}];

if sum(idx)>0
    pos = find(idx>0);
    for i = 1:sum(idx)
        LeftHand(pos(i)) = cellfun(@strjoin, LeftHandCell(pos(i)), repmat({' + '}, sum(idx(pos(i))),1), 'UniformOutput', 0);
    end
end

idx = cellfun(@length, RightHandCell)>1;
RightHand(~idx) = [RightHandCell{~idx}];
if sum(idx)>0
    RightHand(idx) = cellfun(@strjoin, RightHandCell(idx), repmat({' + '}, sum(idx),1), 'UniformOutput', 0);
end

idx = ~strcmp(LeftHand, RightHand);
LeftHand = LeftHand(idx);
RightHand = RightHand(idx);

rxnFormulasReduced = strcat(LeftHand, {'  -> '}, RightHand, {' '});
end
