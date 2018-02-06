function pdb_Stratified_ModelBased_Cognitive()

% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code calculates
% fixed effects model based weights for the cognitive task after balancing
% the trial mean between consistent and inconsistent distributions, this
% avoids the need to add an additional parameter.
close all;
clc; dbstop if error;
% specify the path to the data
behdatapath = '../Data';
behdata = readtable(sprintf('%s/Task_Cognitive.csv', behdatapath));
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
    psycho_noise = 4.7643;
    psycho_bias = 0.6337;
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
[individualparams_consistent, ~]=subplex('pdb_model', startingpoint);
% optimise again, just to make sure we are at the minimum
[Finalparams_consistent,FinalNlogL_consistent] = fminsearchbnd(@(individualparams) pdb_model(individualparams),individualparams_consistent,[0,-1000,-1000],[80,1000,1000],options);

% now fit the inconsistent trials
% define a random starting point for the fitting algorithm
startingpoint = starting_pt([2,4,6]);
TrlsInconsistent = find(sign(subj_dat.binchoice) ~= sign(subj_dat.x2_relative));
dat = subj_dat(TrlsInconsistent,:);
[individualparams_inconsistent, ~] = subplex('pdb_model', startingpoint);
% optimise again, just to make sure we are at the minimum
[Finalparams_inconsistent,FinalNlogL_inconsistent] = fminsearchbnd(@(individualparams) pdb_model(individualparams),individualparams_inconsistent,[0,-1000,-1000],[80,1000,1000],options);
subj_params([1,3,5]) = Finalparams_consistent;
subj_params([2,4,6]) = Finalparams_inconsistent;
subj_NlogL = FinalNlogL_consistent + FinalNlogL_inconsistent;
end

function optimal_funcval = pdb_model(individualparams)
global dat;global psycho_noise; global psycho_bias;
% get the relevant values first
X1              =  dat.x1_relative;
X2              = dat.x2_relative;
RealDecision    = dat.binchoice;
RealEvaluation  = dat.estim;
LEvaluation = NaN(length(X1),1);
step=0.2;
X = -500:step:500; % the range of values over which we calculate the pdf
for i = 1:num_trls
    % get the pdf for the first interval
    Y1cc = normpdf(X, (X1(i) + psycho_bias)*individualparams(2), psycho_noise*abs(individualparams(2)));
    % set the pdf of the unchosen side to zero, assuming subjects base
    % their estimations only on the chosen side
    Y1cc(sign(X) ~= sign(RealDecision(i))) = 0;
    % we need the pdf, so normalise the resulting distribution such that
    % the area under the curve is 1
    Y1cc = Y1cc./trapz(X,Y1cc);
    % get the pdf of the second interval
    Y2cc = normpdf(X, (X2(i) + psycho_bias)*individualparams(3), psycho_noise*abs(individualparams(3)));
    % convolve the two pdfs corresponding to the two inervals respectively-
    % this is similar to adding two random variables drawn from the
    % distributions
    N1 = conv(Y1cc,Y2cc,'same');
    % we need the pdf, so normalise the resulting distribution such that
    % the area under the curve is 1
    N1=N1./trapz(X,N1);
    % we will also add an additive term to account for the attraction-to-mean effect in the numbers data to a zero mean Normal distribution that represents the estimation
    % noise
    NoiseDist = normpdf(X,0,individualparams(1));
    % add this noise to the estimations from the two intervals
    NN1 = conv(N1,NoiseDist,'same');
    % we need the pdf, so normalise the resulting distribution such that
    % the area under the curve is 1; this will be the final pdf
    % corresponding to the estimations of this trial
    NN1 = NN1./trapz(X,NN1);
    % get the likelihood of the subject's estimation from the pdf obtained
    % above
    LEvaluation(i) = interp1(X+50,NN1,RealEvaluation(i)); % we add 50 because 50 is the reference in this data, unlike 0 in the perceptual data
    % to make sure the pdfs do not go over the range of values at which
    % subjects can respond (0-100)
    if nanmean(NN1) > 50 || nanmean(NN1) < -50
        LEvaluation(i) = NaN;
    end
end
% total log likelihood, one value for each trial
PlogObservation = sum(log(LEvaluation));
% allow for nan and inf in the loglikelihood function, but set their probability to be super small
if isinf(abs(PlogObservation))
    PlogObservation = -10e100;
end
if isnan(PlogObservation)
    PlogObservation = -10e100;
end
% take the negative ll for fminsearch!
optimal_funcval = -PlogObservation;
end