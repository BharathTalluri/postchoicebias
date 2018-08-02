function [params, BIC] = pdb_Perceptual_StimulusSelectiveModel(behdata, isplot)

% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code calculates
% model based weights for the stimulus-based selective gain model.
global subjects;global psycho_fits;
% initialise some variables
Finalparams = NaN (length(subjects), 6);
FinalNlogL = NaN(length(subjects),1);
BIC = NaN(length(subjects),1);
anova.x = [];
anova.sj = [];
anova.f1 = [];
anova.f2 = [];
rng shuffle;
global psycho_noise psycho_bias
for sj = subjects
    subj_dat             =  behdata(find(behdata.subj == sj),:);
    % use only choice trials in this paper
    choicetrials = find(subj_dat.condition == 1 & abs(subj_dat.binchoice) == 1);
    subj_dat = subj_dat(choicetrials,:);
    % since stimulus consistency cannot be defined for trials where x2/x1 = 0, remove those
    % trials
    trls2use    = find(subj_dat.x2 ~= 0 & subj_dat.x1 ~= 0);
    subj_dat             = subj_dat(trls2use, :);
    psycho_noise = psycho_fits.logisticFit(find(sj==subjects), 2);
    psycho_bias = -psycho_fits.logisticFit(find(sj==subjects), 1);
    starting_pt = [datasample(1:5:25, 1) datasample(1:5:25, 1) datasample(0.05:0.05:1, 1) datasample(0.05:0.05:1, 1) datasample(0.05:0.05:1, 1) datasample(0.05:0.05:1, 1)];
    [Finalparams(find(sj==subjects),:), FinalNlogL(find(sj==subjects))] = fit_model(subj_dat, starting_pt);
    [~, BIC(find(sj==subjects))] = aicbic(-FinalNlogL(find(sj==subjects)), size(Finalparams(find(sj==subjects),:),2) , size(subj_dat, 1));
    anova.x     = [anova.x; Finalparams(find(sj == subjects),3); Finalparams(find(sj == subjects),5); Finalparams(find(sj == subjects),4); Finalparams(find(sj == subjects),6)];
    anova.sj    = [anova.sj; find(sj == subjects) * ones(4, 1)];
    anova.f1    = [anova.f1; [1 2 1 2]'];
    anova.f2    = [anova.f2; [1 1 2 2]'];
end
params = Finalparams;
% compute 2-way ANOVA measures
anov = rm_anova(anova.x, anova.sj, {anova.f1; anova.f2});
if isplot
    % to the plotting function now
    plot_params(Finalparams, anov);
end
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

function plot_params(boot_params, anov)
% specify the color map and figure properties
cols = linspecer(10, 'qualitative');
myfigureprops;
figure;
subplot(4,4,1); hold on;
% weights of second interval
dat1 = boot_params(:, 5);
dat2 = boot_params(:, 6);
myscatter(dat1, dat2, [], 75, cols, cols, cols);
% plot the group mean +/- s.e.m
plot([nanmean(dat1)-nansem(dat1) nanmean(dat1)+nansem(dat1)], [nanmean(dat2) nanmean(dat2)], 'Color', cols(1,:),'LineWidth',2);
plot([nanmean(dat1) nanmean(dat1)], [nanmean(dat2)-nansem(dat2) nanmean(dat2)+nansem(dat2)], 'Color', cols(1,:),'LineWidth',2);
% polish the figure
set(gca, 'XLim', [-0.2 1.0], 'XTick', -0.2:0.4:1.0,'ylim',[-0.2 1.0], 'ytick', -0.2:0.4:1.0);
plot([0 00], [-0.2 1.0], 'k:', 'LineWidth', 0.5);
plot([-0.2 1.0], [0 0], 'k:', 'LineWidth', 0.5);
EquateAxis;
[pval] = permtest(dat1, dat2, 0, 100000); % permutation test
axis square;
xlabel({'Weight for subsequent', 'inconsistent stimulus'});
ylabel({'Weight for subsequent', 'consistent stimulus'});
offsetAxes;
title({'Model-based results',sprintf('Consistent vs. Inconsistent: p = %.4f', pval)});

subplot(4,4,3); hold on;
plot([1 2], nanmean([boot_params(:, 3) boot_params(:, 5)], 1), 'O-', 'Color', [0 0 0], 'MarkerSize', 7.5, 'LineWidth', 1.5, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [1,1,1]);
plot([1 2], nanmean([boot_params(:, 4) boot_params(:, 6)], 1), 'O-', 'Color', [0.6 0.6 0.6], 'MarkerSize', 7.5, 'LineWidth', 1.5, 'MarkerFaceColor', [0.6 0.6 0.6], 'MarkerEdgeColor', [1,1,1]);
errbar([1 2], nanmean([boot_params(:, 3) boot_params(:, 5)], 1), nanstd([boot_params(:, 3) boot_params(:, 5)], 1) ./ sqrt(length(boot_params(:, 3))), '-','Color', [0 0 0], 'LineWidth',1);
errbar([1 2], nanmean([boot_params(:, 4) boot_params(:, 6)], 1), nanstd([boot_params(:, 4) boot_params(:, 6)], 1) ./ sqrt(length(boot_params(:, 4))), '-','Color', [0.6 0.6 0.6], 'LineWidth',1);
set(gca, 'XLim', [0.5 2.5], 'XTick', 1:2, 'XTickLabel', {'Stimulus 1', 'Stimulus 2'},  'ylim',[0 0.6], 'ytick', 0:0.3:0.6);
ylabel('Weights');
axis square;offsetAxes;
[pval] = permtest(boot_params(:, 3), boot_params(:, 5), 0, 100000); % ranksum much faster than permtest
mysigstar([1, 2], 0.4, pval, 0);
[pval] = permtest(boot_params(:, 4), boot_params(:, 6), 0, 100000); % ranksum much faster than permtest
mysigstar([1, 2], 0, pval, 0);
[pval] = permtest(boot_params(:, 3), boot_params(:, 4), 0, 100000); % ranksum much faster than permtest
mysigstar(0.8, 0.4, pval, 0);
[pval] = permtest(boot_params(:, 5), boot_params(:, 6), 0, 100000); % ranksum much faster than permtest
mysigstar(2.2, 0.4, pval, 0);
legend('Consistent trials', 'Inconsistent trials');
title(sprintf('F_{(%d,%d)} = %.2f, p = %.3f', anov.f1xf2.df, anov.f1xf2.fstats, anov.f1xf2.pvalue))


subplot(4,4,9); hold on;
% estimation noise
dat1 = boot_params(:, 1);
dat2 = boot_params(:, 2);
myscatter(dat1, dat2, [], 75, cols, cols, cols);
% plot the group mean +/- s.e.m
plot([nanmean(dat1)-nansem(dat1) nanmean(dat1)+nansem(dat1)], [nanmean(dat2) nanmean(dat2)], 'Color', cols(1,:),'LineWidth',2);
plot([nanmean(dat1) nanmean(dat1)], [nanmean(dat2)-nansem(dat2) nanmean(dat2)+nansem(dat2)], 'Color', cols(1,:),'LineWidth',2);
% polish the figure
set(gca, 'XLim', [0 20], 'XTick', 0:10:20,'ylim',[0 20], 'ytick', 0:10:20);
EquateAxis;
[pval] = permtest(dat1, dat2, 0, 100000); % permutation test
axis square;
xlabel({'Estimation noise', 'for Inconsistent trials'});
ylabel({'Estimation noise', 'for Consistent trials'});
offsetAxes;
title(sprintf('Consistent vs. Inconsistent: p = %.4f', pval));
end