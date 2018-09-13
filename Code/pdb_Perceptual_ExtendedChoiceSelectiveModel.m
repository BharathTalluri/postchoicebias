function pdb_Perceptual_ExtendedChoiceSelectiveModel(behdata, selective_params)

% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code calculates
% residual shift parameter for the extended choice-based selective gain model.
global subjects;global psycho_fits;global fixed_params;
% initialise some variables
Finalparams = NaN (length(subjects), 1);
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
    starting_pt = datasample(1:5:25, 1);
    fixed_params = selective_params(find(sj==subjects),:);
    Finalparams(find(sj==subjects),1) = fit_model(subj_dat, starting_pt);
end
plot_residual_shift(Finalparams);
end

function subj_params = fit_model(subj_dat, starting_pt)
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
subj_params = fminsearchbnd(@(individualparams) pdb_model(individualparams),individualparams,-1000,1000,options);
end

function optimal_funcval = pdb_model(individualparams)
global dat;global psycho_noise; global psycho_bias;global fixed_params;
% get the relevant values first
X1 = dat.x1;
X2 = dat.x2;
RealDecision = dat.binchoice;
RealEvaluation = dat.estim;
num_trls = length(X1);
LEvaluation = NaN(num_trls,1);
step = 0.05;
X = -180:step:180; % the range of values over which we calculate the pdf
% identify the indices of consistent and inconsistent trials
consistent_trls = find(RealDecision == sign(X2));
inconsistent_trls = find(RealDecision ~= sign(X2));

% do the fitting separately for consistent and inconsistent trials
% consistent trls
for i = consistent_trls'
    % get the pdf for the first interval
    Y1cc = normpdf(X, (X1(i) + psycho_bias)*fixed_params(3), psycho_noise*abs(fixed_params(3)));
    % set the pdf of the unchosen side to zero, assuming subjects base
    % their estimations only on the chosen side
    Y1cc(sign(X) ~= sign(RealDecision(i))) = 0;
    % we need the pdf, so normalise the resulting distribution such that
    % the area under the curve is 1
    Y1cc = Y1cc./trapz(X,Y1cc);
    % get the pdf of the second interval
    Y2cc = normpdf(X, (X2(i) + psycho_bias)*fixed_params(5) + RealDecision(i).*individualparams, psycho_noise*abs(fixed_params(5)));
    % convolve the two pdfs corresponding to the two inervals respectively-
    % this is similar to adding two random variables drawn from the
    % distributions
    N1 = conv(Y1cc,Y2cc,'same');
    % we need the pdf, so normalise the resulting distribution such that
    % the area under the curve is 1
    N1=N1./trapz(X,N1);
    % get a zero mean Normal distribution that represents the estimation
    % noise
    NoiseDist = normpdf(X,0,fixed_params(1));
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
% inconsistent trls
for i = inconsistent_trls'
    % get the pdf for the first interval
    Y1cc = normpdf(X, (X1(i) + psycho_bias)*fixed_params(4), psycho_noise*abs(fixed_params(4)));
    % set the pdf of the unchosen side to zero, assuming subjects base
    % their estimations only on the chosen side
    Y1cc(sign(X) ~= sign(RealDecision(i))) = 0;
    % we need the pdf, so normalise the resulting distribution such that
    % the area under the curve is 1
    Y1cc = Y1cc./trapz(X,Y1cc);
    % get the pdf of the second interval
    Y2cc = normpdf(X, (X2(i) + psycho_bias)*fixed_params(6) + RealDecision(i).*individualparams, psycho_noise*abs(fixed_params(6)));
    % convolve the two pdfs corresponding to the two inervals respectively-
    % this is similar to adding two random variables drawn from the
    % distributions
    N1 = conv(Y1cc,Y2cc,'same');
    % we need the pdf, so normalise the resulting distribution such that
    % the area under the curve is 1
    N1=N1./trapz(X,N1);
    % get a zero mean Normal distribution that represents the estimation
    % noise
    NoiseDist = normpdf(X,0,fixed_params(2));
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

function plot_residual_shift(params)
% specify the color map and figure properties
cols = linspecer(10, 'qualitative');
colormap(linspecer);
figure;
subplot(4,4,1); hold on;
% polish the figure
set(gca, 'XLim', [0. 0.5], 'XTick', [], 'ylim',[-1.5 4.5], 'ytick', -1.5:1.5:4.5);
plot([0. 0.5], [0 0], 'k', 'LineWidth', 0.25);
myscatter(0.25*ones(1,length(params)), params, [], 75, cols, cols, cols);
[pval] = permtest(zeros(1,length(params)), params, 0, 100000); % permutation test
ylabel('Residual Shift parameter');
offsetAxes;
title(sprintf('p = %.4f', pval));
end