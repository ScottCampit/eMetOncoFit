%% PRIME algorithm
  % PRIME consists of the following functions:
    % Find the maximal bounds
    % Find the minimum and maximum values for the normalization
    % Running PRIME

%% First, we need to find all reversible reactions and split them
rev_rxns = find(model.lb<0);
[num_mets, num_rxns] = size(model.S)

model.rev_map(:,1) = [num_rxns + 1:num_rxns + length(rev_rxns)];
model.rev_map(:,2) = rev_rxns;

rev_mao = model.rev_map
model.S = [model.S - model.S(:, rev_rxns)];

orig_lb = model.lb(rev_rxns);
model.lb(rev_rxns) = 0;
model.rev(rev_rxns) = 0;
model.rev = [model.rev; zeros(length(rev_rxns), 1)];
model.lb = [model.lb; zeros(length(rev_rxns), 1)];

model.ub = [model.ub; -orig_lb];
bkwd_rxns = model.rxns(rev_rxns);
bkwd_rxnNames = model.rxns(rev_rxns);

for i=1:length(rev_rxns)
        name = model.rxns{rev_rxns(i)};
        name_long = model.rxnNames{rev_rxns(i)};
        model.rxns{rev_rxns(i)} = sprintf ('%s_fwd',name);
        bkwd_rxns{i} = sprintf ('%s_bkwd',name);
        model.rxnNames{rev_rxns(i)} = sprintf ('%s_fwd',name_long);
        bkwd_rxnNames{i} = sprintf ('%s_bkwd',name_long);
end

model.rxns = [model.rxns; bkwd_rxns];
model.rxnNames = [model.rxnNames; bkwd_rxnNames];
model.rxnGeneMat = [model.rxnGeneMat; model.rxnGeneMat(rev_rxns,:)];
model.rules = [model.rules; model.rules(rev_rxns)];
model.grRules = [model.grRules; model.grRules(rev_rxns)];
model.subSystems = [model.subSystems; model.subSystems(rev_rxns)];
model.c = [model.c; model.c(rev_rxns)];

%% Second, we need to find a new bound that does not affect cell growth by reducing the bounds on all reactions at once

%% Third, we need to identify the minimum and maximum values for the normalization range.

% Split the reaction to its forward and backwards direction

% Set the new maximum bounds

% Compute the max biomass and set the biomass reaction to 10% of its max

% Calculate the minimum values necessary to produce 10% maximum biomass

% Reset the biomass bound

% Compute the effect of bound change on growth rate

% Compute the max range by searching for the highest point that mostly affects growth


%% PRIME
% Identify the set of reactions associated with growth

% Compute the correlation between reaction expression and growth rate

% Correct for multiple hypotheses using FDR and significance level
p = sort(p(:));
V = length(p);
I = (1:V)';

cVID = 1;
cVN = sum(1./(1:V));

pID = p(max(find(p<=I/V*q/cVID)));
pN = p(max(find(p<=I/V*q/cVN)));

% Identify the set of significantly correlated reactions

% Update all vars

% Transform expression values to bounds according to normalization and direction of correlation

% Make context specific model
