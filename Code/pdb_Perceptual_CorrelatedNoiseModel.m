function [params, BIC] = pdb_Perceptual_CorrelatedNoiseModel(behdata)

% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code calculates
% model based weights for the baseline model.
global subjects;global psycho_fits;
% initialise some variables
Finalparams = NaN (length(subjects), 4);
FinalNlogL = NaN(length(subjects),1);
BIC = NaN(length(subjects),1);
rng shuffle;
global psycho_noise psycho_bias
for sj = subjects
    subj_dat             =  behdata(find(behdata.subj == sj),:);
    % use only choice trials in this paper
    choicetrials = find(subj_dat.condition == 1 & abs(subj_dat.binchoice) == 1);
    subj_dat = subj_dat(choicetrials,:);
    % have the same number of trials across all models for comparing BIC
    % values
    trls2use    = find(subj_dat.x2 ~= 0 & subj_dat.x1 ~= 0);
    subj_dat             = subj_dat(trls2use, :);
    psycho_noise = psycho_fits.logisticFit(find(sj==subjects), 2);
    psycho_bias = -psycho_fits.logisticFit(find(sj==subjects), 1);
    starting_pt = [datasample(1:5:25, 1) datasample(0.05:0.05:1, 1) datasample(0.05:0.05:1, 1) datasample(0:0.1:1, 1)];
    [Finalparams(find(sj==subjects),:), FinalNlogL(find(sj==subjects))] = fit_model(subj_dat, starting_pt);
    [~, BIC(find(sj==subjects))] = aicbic(-FinalNlogL(find(sj==subjects)), size(Finalparams(find(sj==subjects),:),2) , size(subj_dat, 1));
end
params = Finalparams;
end

function [subj_params, subj_NlogL] = fit_model(subj_dat, starting_pt)
% function evaluation params
global dat
options = optimset('Display', 'notify') ;
options.MaxFunEvals = 1e10; % limit the nr of func evals
options.MaxIter = 500000;
options.TolX = 0.00001; % dont make this too small, will take forever to converge
options.TolFun = 0.00001;
options.Robust = 'on';

% define a random starting point for the fitting algorithm
dat = subj_dat;
[individualparams, ~]=subplex('pdb_model', starting_pt);
% optimise again, just to make sure we are at the minimum
[subj_params,subj_NlogL] = fminsearchbnd(@(individualparams) pdb_model(individualparams),individualparams,[0,-1000,-1000, 0],[80,1000,1000, 1],options);
end

function optimal_funcval = pdb_model(individualparams)
global dat;global psycho_noise; global psycho_bias;
niter = 1e4;
% get the relevant values first
X1 = dat.x1;
X2 = dat.x2;
RealDecision = dat.binchoice;
RealEvaluation = dat.estim;
num_trls = length(X1);
LEvaluation = NaN(num_trls,1);
X = -180:180; % the range of values over which we calculate the pdf
for i = 1:num_trls
    % generate a distribution of the correlated noise
X1_noise = psycho_noise.*randn(niter,1);
NoisyX1 = X1(i) + psycho_bias + X1_noise;
% get all the trials where the binary choice matches the simulated choice
trls2use = find(sign(NoisyX1) == RealDecision(i));
% use only these trials in the future
iter_valid = length(trls2use);
NoisyX1 = individualparams(2)*NoisyX1(trls2use);
% now generate the noise distribution for second interval
X2_noise = psycho_noise.*randn(iter_valid,1);
% combine this noise with correlated noise and evidence presented in
% interval 2
NoisyX2 = individualparams(3)*(X2(i) + psycho_bias + individualparams(4)*X1_noise(trls2use) + (1-individualparams(4))*X2_noise);
% get the final evaluation corrupted with estimation noise
Evaluation = NoisyX1 + NoisyX2 + randn(iter_valid,1).*individualparams(1);
Evaluation = [Evaluation; X'];
pd = fitdist(Evaluation, 'Kernel', 'Kernel','epanechnikov');
LEvaluation(i) = pdf(pd,RealEvaluation(i));
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