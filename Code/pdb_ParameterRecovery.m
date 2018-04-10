function pdb_ParameterRecovery()

% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code reproduces
% figure S2 showing that the model succesfully recovers parameters from
% simulated data
close all;
clc; dbstop if error;
% specify the path to the data
behdatapath = '../Data';
behdata = readtable(sprintf('%s/Task_Perceptual.csv', behdatapath));
% initialise some variables
subjects = unique(behdata.subj)';
recovery.actualparams = NaN (500, 3);
recovery.recoveredparams = NaN (500, 3);
recovery.recoveredNlogL = NaN(500,1);
% get the noise and bias parameters from psychometric fits
psycho_fits = pdb_Behaviour('Perceptual', 0);
rng shuffle;
global psycho_noise
for iter = 1:500
    % pick a random subject to simulate
    sj = datasample(subjects,1);
    subj_dat             =  behdata(find(behdata.subj == sj),:);
    % use only choice trials in this paper
    choicetrials = find(subj_dat.condition == 1 & abs(subj_dat.binchoice) == 1);
    subj_dat = subj_dat(choicetrials,:);
    psycho_noise = psycho_fits.logisticFit(sj, 2);
    % simulate binary choice and estimation behaviour using random
    % parameters
    params2simulate = [datasample(5:5:20,1), datasample(0.05:0.05:1,1) datasample(0.05:0.05:1, 1)];
    est_noise = params2simulate(1);
    w1 = params2simulate(2);
    w2 = params2simulate(3);
    % start simulating estimations and decisions
    x1 = subj_dat.x1;
    x2 = subj_dat.x2;
    % generate a distribution of noisy stimulus representations
    NoisyX1 = bsxfun(@plus, x1, psycho_noise .* randn(length(x1),1)); % add noise
    NoisyX2 = bsxfun(@plus, x2, psycho_noise .* randn(length(x2),1)); % add noise
    Noise_est = est_noise .*randn(length(x1),1);
    Binchoice = sign(NoisyX1);
    % get evaluation
    Evaluation = bsxfun(@times, NoisyX1, w1) ...
        + bsxfun(@times, NoisyX2, w2) + Noise_est;
    subj_dat.estim = Evaluation;
    subj_dat.binchoice = Binchoice;
    % end of simulated estimations
    starting_pt = [psycho_noise 0.5 0.5];
    [recovery.recoveredparams(iter,:), recovery.recoveredNlogL(iter)] = fit_model(subj_dat, starting_pt);
    recovery.actualparams(iter,:) = params2simulate;
end
plot_params(recovery);
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

% fit all the trials in the model
dat = subj_dat;
[individualparams, ~]=subplex('pdb_model', starting_pt);
% optimise again, just to make sure we are at the minimum
[subj_params,subj_NlogL] = fminsearchbnd(@(individualparams) pdb_model(individualparams),individualparams,[0,0,0],[20,1,1],options);
end

function optimal_funcval = pdb_model(individualparams)
global dat;global psycho_noise; 
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
    Y1cc = normpdf(X, X1(i)*individualparams(2), psycho_noise*abs(individualparams(2)));
    % set the pdf of the unchosen side to zero, assuming subjects base
    % their estimations only on the chosen side
    Y1cc(sign(X) ~= sign(RealDecision(i))) = 0;
    % we need the pdf, so normalise the resulting distribution such that
    % the area under the curve is 1
    Y1cc = Y1cc./trapz(X,Y1cc);
    % get the pdf of the second interval
    Y2cc = normpdf(X, X2(i)*individualparams(3), psycho_noise*abs(individualparams(3)));
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

function plot_params(params)
% specify the color map and figure properties
cols = linspecer(9, 'qualitative');
myfigureprops;
figure;
% just in case there are points where the fitting algorithm did not give a
% valid solution, i.e., infinite or nan loglikelihood, remove those points
trls = find(params.recoveredNlogL < 1e98);

% first the parameter recovery
subplot(4,4,1);hold on;
scatter(params.actualparams(trls,1),params.recoveredparams(trls,1),30,'filled','MarkerEdgeColor',[1,1,1], 'MarkerFaceColor', cols(8,:));
set(gca, 'ylim',[0 40], 'ytick', 0:20:40, 'xlim', [0 40], 'xtick', 0:20:40);
xlabel('Estimation Noise');
axis square;EquateAxis;
offsetAxes;

subplot(4,4,2);hold on;
scatter(params.actualparams(trls,2),params.recoveredparams(trls,2),30,'filled','MarkerEdgeColor',[1,1,1], 'MarkerFaceColor', cols(8,:));
set(gca, 'ylim',[0 1], 'ytick', 0:0.5:1, 'xlim', [0 1], 'xtick', 0:0.5:1);
xlabel('W1');
axis square;EquateAxis;
offsetAxes;

subplot(4,4,3);hold on;
scatter(params.actualparams(trls,3),params.recoveredparams(trls,3),30,'filled','MarkerEdgeColor',[1,1,1], 'MarkerFaceColor', cols(8,:));
set(gca, 'ylim',[0 1], 'ytick', 0:0.5:1, 'xlim', [0 1], 'xtick', 0:0.5:1);
xlabel('W2');
axis square;EquateAxis;
offsetAxes;
suplabel('Actual parameters','x');
suplabel('Recovered parameters','y');

% next we see if the model introduced any spurious correlations in the
% recovered params
subplot(4,4,5);hold on;
scatter(params.actualparams(trls,2),params.actualparams(trls,3),30,'filled','MarkerEdgeColor',[1,1,1], 'MarkerFaceColor', cols(8,:));
set(gca, 'ylim',[0 1], 'ytick', 0:0.5:1, 'xlim', [0 1], 'xtick', 0:0.5:1);
axis square;EquateAxis;
offsetAxes;
xlabel('W1');
ylabel('W2');

subplot(4,4,6);hold on;
scatter(params.actualparams(trls,1),params.actualparams(trls,2),30,'filled','MarkerEdgeColor',[1,1,1], 'MarkerFaceColor', cols(8,:));
set(gca, 'xlim',[0 40], 'xtick', 0:20:40, 'ylim', [0 1], 'ytick', 0:0.5:1);
axis square;
offsetAxes;
xlabel('Estimation Noise');
ylabel('W1');

subplot(4,4,7);hold on;
scatter(params.actualparams(trls,1),params.actualparams(trls,3),30,'filled','MarkerEdgeColor',[1,1,1], 'MarkerFaceColor', cols(8,:));
set(gca, 'xlim',[0 40], 'xtick', 0:20:40, 'ylim', [0 1], 'ytick', 0:0.5:1);
axis square;
offsetAxes;
xlabel('Estimation Noise');
ylabel('W2');

subplot(4,4,9);hold on;
scatter(params.recoveredparams(trls,2),params.recoveredparams(trls,3),30,'filled','MarkerEdgeColor',[1,1,1], 'MarkerFaceColor', cols(8,:));
set(gca, 'ylim',[0 1], 'ytick', 0:0.5:1, 'xlim', [0 1], 'xtick', 0:0.5:1);
axis square;EquateAxis;
offsetAxes;
xlabel('W1');
ylabel('W2');

subplot(4,4,10);hold on;
scatter(params.recoveredparams(trls,1),params.recoveredparams(trls,2),30,'filled','MarkerEdgeColor',[1,1,1], 'MarkerFaceColor', cols(8,:));
set(gca, 'xlim',[0 40], 'xtick', 0:20:40, 'ylim', [0 1], 'ytick', 0:0.5:1);
axis square;
offsetAxes;
xlabel('Estimation Noise');
ylabel('W1');

subplot(4,4,11);hold on;
scatter(params.recoveredparams(trls,1),params.recoveredparams(trls,3),30,'filled','MarkerEdgeColor',[1,1,1], 'MarkerFaceColor', cols(8,:));
set(gca, 'xlim',[0 40], 'xtick', 0:20:40, 'ylim', [0 1], 'ytick', 0:0.5:1);
axis square;
offsetAxes;
xlabel('Estimation Noise');
ylabel('W2');
end