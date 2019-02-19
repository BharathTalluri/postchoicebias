function pdb_Numerical_Stratified_ChoiceSelectiveModel(behdata)

% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code calculates
% fixed effects model based weights for the cognitive task after balancing
% the trial mean between consistent and inconsistent distributions, this
% avoids the need to add an additional parameter.

global psycho_fits;
% initialise some variables
Finalparams = NaN (500, 6);
FinalNlogL = NaN(500,1);
rng shuffle;
global psycho_noise psycho_bias
for iter = 1:500
    subj_dat             =  behdata;
    % use only choice trials in this paper
    choicetrials = find(subj_dat.condition == 1 & abs(subj_dat.binchoice) == 1);
    subj_dat = subj_dat(choicetrials,:);
    % since consistency cannot be defined for trials where x2 = 0, remove those
    % trials
    trls2use    = find(subj_dat.x2_relative ~= 0);
    subj_dat             = subj_dat(trls2use, :);
    % get stratified distributions- there is imbalance in trial mean
    % distributions between consistent and inconsistent trials
    cons_trls = find(sign(subj_dat.x2_relative) == subj_dat.binchoice);
    incons_trls = find(sign(subj_dat.x2_relative) ~= subj_dat.binchoice);
    cons_mean = subj_dat.xavg(cons_trls);
    incons_mean = subj_dat.xavg(incons_trls);
    [cons_indices, incons_indices] = stratified_distribution(cons_mean, incons_mean);
    stratified_cons_trials = cons_trls(cons_indices);
    stratified_incons_trials = incons_trls(incons_indices);
    stratified_indices = union(stratified_cons_trials, stratified_incons_trials);
    subj_dat = subj_dat(stratified_indices,:);
    % plugin pre-calculated psychometric values for the fixed effects
    % subject
     psycho_noise = nanmean(psycho_fits.logisticFit(:, 2));
    psycho_bias = -nanmean(psycho_fits.logisticFit(:, 1));
    starting_pt = [datasample(1:5:25, 1) datasample(1:5:25, 1) datasample(0.05:0.05:1, 1) datasample(0.05:0.05:1, 1) datasample(0.05:0.05:1, 1) datasample(0.05:0.05:1, 1)];
    [Finalparams(iter,:), FinalNlogL(iter)] = fit_model(subj_dat, starting_pt);
end
% calculate the confidence intervals for the confirmation bias and print
CI_params = prctile(Finalparams(:,5)-Finalparams(:,6), [2.5, 50, 97.5]);
fprintf('Median = %.06f; Confidence Intervals = %.06f to %.06f', CI_params(2), CI_params(1), CI_params(3));
end

function [subj_params, subj_NlogL] = fit_model(subj_dat, starting_pt)
% function evaluation params
global dat
subj_params = NaN(1,6);
options = optimset('Display', 'notify') ;
options.MaxFunEvals = 1e10; % limit the nr of func evals
options.MaxIter = 500000;
options.TolX = 0.00001; % dont make this too small, will take forever to converge
options.TolFun = 0.00001;
options.Robust = 'on';

% first fit the consistent trials
% define a random starting point for the fitting algorithm
startingpoint = starting_pt([1,3,5]);
TrlsConsistent = find(sign(subj_dat.binchoice) == sign(subj_dat.x2_relative));
dat = subj_dat(TrlsConsistent,:);
[individualparams_consistent, ~]=subplex('model_Numerical_Stratified_ChoiceSelective', startingpoint);
% optimise again, just to make sure we are at the minimum
[Finalparams_consistent,FinalNlogL_consistent] = fminsearchbnd(@(individualparams) model_Numerical_Stratified_ChoiceSelective(individualparams),individualparams_consistent,[0,0,0],[80,1000,1000],options);

% now fit the inconsistent trials
% define a random starting point for the fitting algorithm
startingpoint = starting_pt([2,4,6]);
TrlsInconsistent = find(sign(subj_dat.binchoice) ~= sign(subj_dat.x2_relative));
dat = subj_dat(TrlsInconsistent,:);
[individualparams_inconsistent, ~] = subplex('model_Numerical_Stratified_ChoiceSelective', startingpoint);
% optimise again, just to make sure we are at the minimum
[Finalparams_inconsistent,FinalNlogL_inconsistent] = fminsearchbnd(@(individualparams) model_Numerical_Stratified_ChoiceSelective(individualparams),individualparams_inconsistent,[0,0,0],[80,1000,1000],options);
subj_params([1,3,5]) = Finalparams_consistent;
subj_params([2,4,6]) = Finalparams_inconsistent;
subj_NlogL = FinalNlogL_consistent + FinalNlogL_inconsistent;
end