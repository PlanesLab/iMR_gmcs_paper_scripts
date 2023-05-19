function [tModel, taskReport] = calculate_model_for_task(model, taskStructure)


% get index of exchange reactions
unbal_rxn_inds = find(sum(model.S ~= 0) == 1)';
% obtain the names of the mets participating in these rxns
[unbal_met_inds,~] = find(model.S(:,unbal_rxn_inds) ~= 0);




taskReport_in = table(taskStructure.inputs, 'VariableNames', {'Name'});
taskReport_out = table(taskStructure.outputs, 'VariableNames', {'Name'});
taskReport_rxn = table(cellfun(@(x)([x(1:15) '...']), taskStructure.equations, 'UniformOutput', 0),...
    'VariableNames', {'Name'});

taskReport_in = taskReport_in(~strcmp(taskReport_in.Name,'ALLMETSIN[e]'),:);
taskReport_out = taskReport_out(~strcmp(taskReport_out.Name,'ALLMETSIN[e]'),:);

% correct little mistake for NEFA BLOOD POOL
taskStructure.inputs(strcmp(taskStructure.inputs, 'NEFA blood pool in[x]')) = {'NEFA blood pool in[s]'};
taskStructure.outputs(strcmp(taskStructure.outputs, 'NEFA blood pool in[x]')) = {'NEFA blood pool in[s]'};

IN_OUT_mets = upper([taskStructure.inputs; taskStructure.outputs]);
IN_OUT_mets = setdiff(IN_OUT_mets , {'ALLMETSIN[E]'});

modelMets = upper(strcat(model.metNames,'[',model.comps(model.metComps),']'));

tModel = model;

[I, J]=ismember(IN_OUT_mets,modelMets(unbal_met_inds));
if ~all(I)
    [~, J] = ismember(IN_OUT_mets(J==0),modelMets);
    tModel = addExchangeRxn(tModel, tModel.mets(J));
    
    % recalculate exchanges
    unbal_rxn_inds = find(sum(tModel.S ~= 0) == 1)';
    % obtain the names of the mets participating in these rxns
    [unbal_met_inds,~] = find(tModel.S(:,unbal_rxn_inds) ~= 0);
    unbal_coefs = diag(tModel.S(unbal_met_inds,unbal_rxn_inds));
    assert(all(unbal_coefs == -1))
end

model2 = tModel;
model2.mets = modelMets;

%Set the inputs
tModel.ub(unbal_rxn_inds) = + 0;
tModel.lb(unbal_rxn_inds) = + 0;

if ~isempty(taskStructure.inputs)
    [I, J]=ismember(upper(taskStructure.inputs),modelMets);
    J=J(I); %Only keep the ones with matches
    K=ismember(upper(taskStructure.inputs),'ALLMETS');
    L=~cellfun('isempty',strfind(upper(taskStructure.inputs),'ALLMETSIN'));
    %Check that all metabolites are either real metabolites or
    %ALLMETS/ALLMETSIN
    if ~all(I|K|L)
        error(['ERROR: Could not find all inputs in "[' taskStructure.id '] ' taskStructure.description '"\n']);
    end
    if numel(J)~=numel(unique(J))
        EM=['The constraints on some input(s) in "[' taskStructure.id '] ' taskStructure.description '" are defined more than one time'];
        error(EM);
    end
    %If all metabolites should be added
    if any(K)
        %Check if ALLMETS is the first metabolite. Otherwise print a
        %warning since it will write over any other constraints that
        %are set
        warning('K non zero')
        if K(1)==0
            EM=['ALLMETS is used as an input in "[' taskStructure.id '] ' taskStructure.description '" but it it not the first metabolite in the list. Constraints defined for the metabolites before it will be over-written'];
            error(EM,false);
        end
        %Use the first match of ALLMETS. There should only be one, but
        %still..
        tModel.b(:,1)=taskStructure.UBin(find(K,1))*-1;
    end
    %If metabolites in a specific compartment should be used
    if any(L)
        L=find(L);
        for j=1:numel(L)
            %The compartment defined
            compartment=upper(taskStructure.inputs{L(j)}(11:end-1));
            %Check if it exists in the model
            C=find(ismember(upper(model.comps),compartment));
            if any(C)
                %Match to metabolites
                tModel.b(model.metComps==C,1)=taskStructure.UBin(L(j))*-1;
            else
                EM=['The compartment defined for ALLMETSIN in "[' taskStructure.id '] ' taskStructure.description '" does not exist'];
                error(EM);
            end
        end
    end
    %Then add the normal constraints
    if any(J)
        [~,JJ] = ismember(J, unbal_met_inds); % change to index in reactions
        tModel.lb(unbal_rxn_inds(JJ)) = taskStructure.UBin(I)*-1;
        tModel.ub(unbal_rxn_inds(JJ)) = taskStructure.LBin(I)*-1;
    end
end

taskReport_in.class(:) = {'in'};
taskReport_in.rxn = tModel.rxns(unbal_rxn_inds(JJ));
taskReport_in.rxnNum = unbal_rxn_inds(JJ);
taskReport_in.rxnFormula = printRxnFormula( model2, 'rxnAbbrList', tModel.rxns(unbal_rxn_inds(JJ)), 'printFlag', 0);
taskReport_in.lb = tModel.lb(unbal_rxn_inds(JJ));
taskReport_in.ub = tModel.ub(unbal_rxn_inds(JJ));
taskReport_in.obj = zeros(size(JJ));


%Set the outputs
if ~isempty(taskStructure.outputs)
    [I, J]=ismember(upper(taskStructure.outputs),modelMets);
    J=J(I); %Only keep the ones with matches
    K=ismember(upper(taskStructure.outputs),'ALLMETS');
    L=~cellfun('isempty',strfind(upper(taskStructure.outputs),'ALLMETSIN'));
    %Check that all metabolites are either real metabolites or
    %ALLMETS/ALLMETSIN
    if ~all(I|K|L)
        error(['ERROR: Could not find all outputs in "[' taskStructure.id '] ' taskStructure.description '"\n']);
    end
    if numel(J)~=numel(unique(J))
        EM=['The constraints on some output(s) in "[' taskStructure.id '] ' taskStructure.description '" are defined more than one time'];
        error(EM);
    end
    %If all metabolites should be added
    if any(K)
        %Check if ALLMETS is the first metabolite. Otherwise print a
        %warning since it will write over any other constraints that
        %are set
        warning('K non zero')
        if K(1)==0
            EM=['ALLMETS is used as an output in "[' taskStructure.id '] ' taskStructure.description '" but it it not the first metabolite in the list. Constraints defined for the metabolites before it will be over-written'];
            error(EM,false);
        end
        %Use the first match of ALLMETS. There should only be one, but
        %still..
        tModel.b(:,2)=taskStructure.UBout(find(K,1));
    end
    %If metabolites in a specific compartment should be used
    if any(L)
        L=find(L);
        for j=1:numel(L)
            %The compartment defined
            compartment=upper(taskStructure.outputs{L(j)}(11:end-1));
            %Check if it exists in the model
            C=find(ismember(upper(model.comps),compartment));
            if any(C)
                %Match to metabolites
                tModel.ub(unbal_rxn_inds(model.metComps(unbal_met_inds)==C))=taskStructure.UBout(L(j));
            else
                EM=['The compartment defined for ALLMETSIN in "[' taskStructure.id '] ' taskStructure.description '" does not exist'];
                error(EM);
            end
        end
    end
    %Then add the normal constraints
    if any(J)
        %Verify that IN and OUT bounds are consistent. Cannot require
        %that a metabolite is simultaneously input AND output at some
        %nonzero flux.
        I = find(I);  % otherwise indexing becomes confusing
        [~,JJ] = ismember(J, unbal_met_inds); % change to index in reactions
        nonzero_LBin = tModel.lb(unbal_rxn_inds(JJ)) < 0;
        nonzero_LBout = taskStructure.LBout(I) > 0;
        % printRxnFormula(tModel,'rxnAbbrList',tModel.rxns(unbal_rxn_inds(JJ(nonzero_LBin))),'printBounds', 1);
        if any(nonzero_LBin & nonzero_LBout)
            EM=['The IN LB and OUT LB in "[' taskStructure.id '] ' taskStructure.description '" cannot be nonzero for the same metabolite'];
            error(EM);
        end
        tModel.lb(unbal_rxn_inds(JJ(nonzero_LBout))) = taskStructure.LBout(I(nonzero_LBout));
        tModel.ub(unbal_rxn_inds(JJ)) = taskStructure.UBout(I);
    end
end

taskReport_out.class(:) = {'out'};
taskReport_out.rxn = tModel.rxns(unbal_rxn_inds(JJ));
taskReport_out.rxnNum = unbal_rxn_inds(JJ);
taskReport_out.rxnFormula = printRxnFormula( model2, 'rxnAbbrList', tModel.rxns(unbal_rxn_inds(JJ)), 'printFlag', 0);
taskReport_out.lb = tModel.lb(unbal_rxn_inds(JJ));
taskReport_out.ub = tModel.ub(unbal_rxn_inds(JJ));
taskReport_out.obj = zeros(size(JJ));


%Add new rxns
if ~isempty(taskStructure.equations)
    rxn.equations=taskStructure.equations;
    rxn.lb=taskStructure.LBequ;
    rxn.ub=taskStructure.UBequ;
    rxn.rxns=strcat({'TEMPORARY_'},num2str((1:numel(taskStructure.equations))'));
    %Allow for new metabolites to be added. This is because it should
    %be possible to add, say, a whole new pathway
    tModel=addRxns(tModel,rxn,3,[],true);
end


if ~isempty(taskStructure.equations)
    model2 = tModel;
    model2.mets = modelMets;

    taskReport_rxn.class(:) = {'equation'};
    taskReport_rxn.rxn = rxn.rxns;
    [~, taskReport_rxn.rxnNum] = ismember(rxn.rxns, tModel.rxns);
    taskReport_rxn.rxnFormula = printRxnFormula( model2, 'rxnAbbrList', rxn.rxns, 'printFlag', 0);
    taskReport_rxn.lb = rxn.lb;
    taskReport_rxn.ub = rxn.ub;
    taskReport_rxn.obj = zeros(size(rxn.ub));
end

%Add changed bounds
if ~isempty(taskStructure.changed)
    tModel=setParam(tModel,'lb',taskStructure.changed,taskStructure.LBrxn);
    tModel=setParam(tModel,'ub',taskStructure.changed,taskStructure.UBrxn);
end

if ~isempty(taskStructure.equations)
    taskReport = [taskReport_in; taskReport_out; taskReport_rxn];
else
    taskReport = [taskReport_in; taskReport_out];
end



% set objective function for this task
set_lb_more_0 = taskReport.lb > 0;
set_ub_less_0 = taskReport.ub < 0;
taskReport.obj(set_lb_more_0 | set_ub_less_0) = 1;
tModel = changeObjective(tModel, taskReport.rxn(set_lb_more_0 | set_ub_less_0));

tModel = setParam(tModel,'lb',taskReport.rxn(set_lb_more_0), 0);
tModel = setParam(tModel,'ub',taskReport.rxn(set_lb_more_0), 1000);

tModel = setParam(tModel,'ub',taskReport.rxn(set_ub_less_0), 0);
tModel = setParam(tModel,'lb',taskReport.rxn(set_ub_less_0), -1000);



% add auxiliar metabolite if more than one objective function
if sum(set_ub_less_0 |set_lb_more_0)>1
    tModel = addMetabolite(tModel, {'temp_XYZ'});
    tModel = addExchangeRxn(tModel, {'temp_XYZ'}, 0, max(tModel.ub));

    tModel.S(findMetIDs(tModel, {'temp_XYZ'}), findRxnIDs(tModel, taskReport.rxn(taskReport.obj>0))) = 1;
    tModel.S(findMetIDs(tModel, {'temp_XYZ'}), findRxnIDs(tModel, taskReport.rxn(taskReport.obj>0 & strcmp(taskReport.class, 'in')))) = -1;
     
    tModel = changeObjective(tModel, {'EX_temp_XYZ'});
    
    model2 = tModel;
    model2.mets(1:(end-1)) = modelMets;
    
    taskReport.Name(end+1) = {'EX_temp_XYZ'};
    taskReport.class(end) = {'aux'};
    taskReport.rxn(end) = tModel.rxns(end);
    taskReport.rxnNum(end)= length(tModel.rxns);
    taskReport.rxnFormula = printRxnFormula( model2, 'rxnAbbrList', taskReport.rxn, 'printFlag', 0);
    taskReport.lb(end) = tModel.lb(end);
    taskReport.ub(end) = tModel.ub(end);
    taskReport.obj(end) = 1;
    
elseif sum(set_ub_less_0 |set_lb_more_0)<1
    error('no objective function found')
end

 
% save previous boundaries
taskReport.lb_prev = taskReport.lb;
taskReport.ub_prev = taskReport.ub;
taskReport.obj_prev = taskReport.obj;
taskReport.lb = tModel.lb(findRxnIDs(tModel, taskReport.rxn));
taskReport.ub = tModel.ub(findRxnIDs(tModel, taskReport.rxn));
taskReport.obj = tModel.c(findRxnIDs(tModel, taskReport.rxn));

if strcmp(taskReport.class(taskReport.obj>0), 'in')
    
    rxnID = findRxnIDs(tModel, taskReport.rxn(taskReport.obj>0));
    
    tModel.S(:, rxnID) = - tModel.S(:, rxnID);
    
    tModel.lb(rxnID) = 0;
    tModel.ub(rxnID) = max(tModel.ub);
        
    model2 = tModel;
    model2.mets(1:length(modelMets)) = modelMets;
        
    taskReport.rxnFormula = printRxnFormula( model2, 'rxnAbbrList', taskReport.rxn, 'printFlag', 0);
    taskReport.lb = tModel.lb(findRxnIDs(tModel, taskReport.rxn));
    taskReport.ub = tModel.ub(findRxnIDs(tModel, taskReport.rxn));
    taskReport.obj = tModel.c(findRxnIDs(tModel, taskReport.rxn));
end  

assert(sum(taskReport.obj)==1)

taskReport



%Solve and print
% tic
% tModel2 = simplifyModel(tModel,true,false,true,true,true);
% toc

end


