function pdb_ModelComparison_Perceptual()

% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code outputs a BIC
% table comparing various models described in the paper for the perceptual
% data
close all;
clc; dbstop if error;
% specify the path to the data
behdatapath = '../Data';
behdata = readtable(sprintf('%s/Task_Perceptual.csv', behdatapath));
% initialise some variables
subjects = unique(behdata.subj)';
unbiased.Finalparams = NaN (length(subjects), 3);
choicebased.Finalparams = NaN (length(subjects), 6);
evidencebased.Finalparams = NaN (length(subjects), 6);
unbiased.FinalNlogL= NaN(length(subjects)+1,1);
choicebased.FinalNlogL= NaN(length(subjects)+1,1);
evidencebased.FinalNlogL= NaN(length(subjects)+1,1);
unbiased.BIC= NaN(length(subjects)+1,1);
choicebased.BIC= NaN(length(subjects)+1,1);
evidencebased.BIC= NaN(length(subjects)+1,1);
num_trls = NaN(length(subjects),1);
% get the noise and bias parameters from psychometric fits
psycho_fits = pdb_Behaviour('Perceptual', 0);
rng shuffle;
global psycho_noise psycho_bias
for sj = subjects
    subj_dat             =  behdata(find(behdata.subj == sj),:);
    % use only choice trials in this paper
    choicetrials = find(subj_dat.condition == 1 & abs(subj_dat.binchoice) == 1);
    subj_dat = subj_dat(choicetrials,:);
    % use the same subset of trials in all models to compare- if
    % consistency has to be defined on the basis of external evidence, then
    % we need to remove trials where x1 or x2 is zero
    trls2use    = find(subj_dat.x2 ~= 0 & subj_dat.x1~= 0);
    subj_dat             = subj_dat(trls2use, :);
    psycho_noise = psycho_fits.logisticFit(sj, 2);
    psycho_bias = -psycho_fits.logisticFit(sj, 1);
    starting_pt = [datasample(1:5:25, 1) datasample(1:5:25, 1) datasample(0.05:0.05:1, 1) datasample(0.05:0.05:1, 1) datasample(0.05:0.05:1, 1) datasample(0.05:0.05:1, 1)];
    [choicebased.Finalparams(sj,:), choicebased.FinalNlogL(sj), choicebased.BIC(sj)] = fit_choicebasedmodel(subj_dat, starting_pt);
    [evidencebased.Finalparams(sj,:), evidencebased.FinalNlogL(sj), evidencebased.BIC(sj)] = fit_evidencebasedmodel(subj_dat, starting_pt);
    [unbiased.Finalparams(sj,:), unbiased.FinalNlogL(sj), unbiased.BIC(sj)] = fit_unbiasedmodel(subj_dat, starting_pt);
    num_trls(sj) = size(subj_dat,1);
end
% get the grand BIC and Negll values across the group
choicebased.FinalNlogL(length(subjects) + 1) = sum(choicebased.FinalNlogL(1:length(subjects)));
[~, choicebased.BIC(length(subjects) + 1)] = aicbic(-choicebased.FinalNlogL(length(subjects)+1), 6*length(subjects), sum(num_trls));
evidencebased.FinalNlogL(length(subjects) + 1) = sum(evidencebased.FinalNlogL(1:length(subjects)));
[~, evidencebased.BIC(length(subjects) + 1)] = aicbic(-evidencebased.FinalNlogL(length(subjects)+1), 6*length(subjects), sum(num_trls));
unbiased.FinalNlogL(length(subjects) + 1) = sum(unbiased.FinalNlogL(1:length(subjects)));
[~, unbiased.BIC(length(subjects) + 1)] = aicbic(-unbiased.FinalNlogL(length(subjects)+1), 3*length(subjects), sum(num_trls));

fprintf('Delta BIC between choice-based and unbiased model = %d \n', choicebased.BIC(length(subjects)+1)-unbiased.BIC(length(subjects)+1));
fprintf('Delta BIC between choice-based and evidence-based model = %d \n', choicebased.BIC(length(subjects)+1)-evidencebased.BIC(length(subjects)+1));
end

function [subj_params, subj_NlogL, subj_BIC] = fit_choicebasedmodel(subj_dat, starting_pt)
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
TrlsConsistent = find(sign(subj_dat.binchoice) == sign(subj_dat.x2));
dat = subj_dat(TrlsConsistent,:);
[individualparams_consistent, ~]=subplex('pdb_model', startingpoint);
% optimise again, just to make sure we are at the minimum
[Finalparams_consistent,FinalNlogL_consistent] = fminsearchbnd(@(individualparams) pdb_model(individualparams),individualparams_consistent,[0,-1000,-1000],[80,1000,1000],options);

% now fit the inconsistent trials
% define a random starting point for the fitting algorithm
startingpoint = starting_pt([2,4,6]);
TrlsInconsistent = find(sign(subj_dat.binchoice) ~= sign(subj_dat.x2));
dat = subj_dat(TrlsInconsistent,:);
[individualparams_inconsistent, ~] = subplex('pdb_model', startingpoint);
% optimise again, just to make sure we are at the minimum
[Finalparams_inconsistent,FinalNlogL_inconsistent] = fminsearchbnd(@(individualparams) pdb_model(individualparams),individualparams_inconsistent,[0,-1000,-1000],[80,1000,1000],options);
subj_params([1,3,5]) = Finalparams_consistent;
subj_params([2,4,6]) = Finalparams_inconsistent;
subj_NlogL = FinalNlogL_consistent + FinalNlogL_inconsistent;
[~,subj_BIC] = aicbic(-subj_NlogL, 6 , size(subj_dat, 1));
end

function [subj_params, subj_NlogL, subj_BIC] = fit_evidencebasedmodel(subj_dat, starting_pt)
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
TrlsConsistent = find(sign(subj_dat.x1) == sign(subj_dat.x2));
dat = subj_dat(TrlsConsistent,:);
[individualparams_consistent, ~]=subplex('pdb_model', startingpoint);
% optimise again, just to make sure we are at the minimum
[Finalparams_consistent,FinalNlogL_consistent] = fminsearchbnd(@(individualparams) pdb_model(individualparams),individualparams_consistent,[0,-1000,-1000],[80,1000,1000],options);

% now fit the inconsistent trials
% define a random starting point for the fitting algorithm
startingpoint = starting_pt([2,4,6]);
TrlsInconsistent = find(sign(subj_dat.x1) ~= sign(subj_dat.x2));
dat = subj_dat(TrlsInconsistent,:);
[individualparams_inconsistent, ~] = subplex('pdb_model', startingpoint);
% optimise again, just to make sure we are at the minimum
[Finalparams_inconsistent,FinalNlogL_inconsistent] = fminsearchbnd(@(individualparams) pdb_model(individualparams),individualparams_inconsistent,[0,-1000,-1000],[80,1000,1000],options);
subj_params([1,3,5]) = Finalparams_consistent;
subj_params([2,4,6]) = Finalparams_inconsistent;
subj_NlogL = FinalNlogL_consistent + FinalNlogL_inconsistent;
[~,subj_BIC] = aicbic(-subj_NlogL, 6 , size(subj_dat, 1));
end

function [subj_params, subj_NlogL, subj_BIC] = fit_unbiasedmodel(subj_dat, starting_pt)
% function evaluation params
global dat
options = optimset('Display', 'notify') ;
options.MaxFunEvals = 1e10; % limit the nr of func evals
options.MaxIter = 500000;
options.TolX = 0.00001; % dont make this too small, will take forever to converge
options.TolFun = 0.00001;
options.Robust = 'on';

% fit all the trials in the unbiased model
% define a random starting point for the fitting algorithm
startingpoint = starting_pt([1,3,5]);
dat = subj_dat;
[individualparams, ~]=subplex('pdb_model', startingpoint);
% optimise again, just to make sure we are at the minimum
[subj_params,subj_NlogL] = fminsearchbnd(@(individualparams) pdb_model(individualparams),individualparams,[0,-1000,-1000],[80,1000,1000],options);

[~,subj_BIC] = aicbic(-subj_NlogL, 3 , size(subj_dat, 1));
end

function optimal_funcval = pdb_model(individualparams)
global dat;global psycho_noise; global psycho_bias;
% get the relevant values first
X1 = dat.x1;
X2 = dat.x2;
RealDecision = dat.binchoice;
RealEvaluation = dat.estim;
num_trls = length(X1);
LEvaluation = NaN(num_trls,1);
step = 0.05;
X = -180:step:180; % the range of values over which we calculate the pdf
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
    % get a zero mean Normal distribution that represents the estimation
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
    LEvaluation(i) = interp1(X,NN1,RealEvaluation(i));
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