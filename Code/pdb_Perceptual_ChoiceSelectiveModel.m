function [params, BIC] = pdb_Perceptual_ChoiceSelectiveModel(behdata, trials, bootstrap, compare, isplot)

% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code calculates
% model based weights and reproduces figure 2B, C,  S2A, C, D of the paper.

% specify the path to the data
global subjects;global psycho_fits;
% initialise some variables
Finalparams.actual = NaN (length(subjects), 6);
if bootstrap
    Finalparams.lowCI = NaN (length(subjects), 6);
    Finalparams.highCI = NaN (length(subjects), 6);
end
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
    if strcmp(trials, 'all')
        choicetrials = find(subj_dat.condition == 1 & abs(subj_dat.binchoice) == 1);
    elseif strcmp(trials, 'error')
        choicetrials = find(subj_dat.condition == 1 & abs(subj_dat.binchoice) == 1 & subj_dat.bincorrect == 0);
    elseif strcmp(trials, 'zero')
        choicetrials = find(subj_dat.condition == 1 & abs(subj_dat.binchoice) == 1 & subj_dat.x1 == 0);
    end
    subj_dat = subj_dat(choicetrials,:);
    % since consistency cannot be defined for trials where x2 = 0, remove those
    % trials
    if compare
        trls2use    = find(subj_dat.x2 ~= 0 & subj_dat.x1 ~= 0);
    else
        trls2use    = find(subj_dat.x2 ~= 0);
    end
    subj_dat             = subj_dat(trls2use, :);
    psycho_noise = psycho_fits.logisticFit(find(sj==subjects), 2);
    psycho_bias = -psycho_fits.logisticFit(find(sj==subjects), 1);
    starting_pt = [datasample(1:5:25, 1) datasample(1:5:25, 1) datasample(0.05:0.05:1, 1) datasample(0.05:0.05:1, 1) datasample(0.05:0.05:1, 1) datasample(0.05:0.05:1, 1)];
    [Finalparams.actual(find(sj==subjects),:), FinalNlogL(find(sj==subjects))] = fit_model(subj_dat, starting_pt);
    [~, BIC(find(sj==subjects))] = aicbic(-FinalNlogL(find(sj==subjects)), size(Finalparams.actual(find(sj==subjects),:),2) , size(subj_dat, 1));
    anova.x     = [anova.x; Finalparams.actual(find(sj == subjects),3); Finalparams.actual(find(sj == subjects),5); Finalparams.actual(find(sj == subjects),4); Finalparams.actual(find(sj == subjects),6)];
    anova.sj    = [anova.sj; find(sj == subjects) * ones(4, 1)];
    anova.f1    = [anova.f1; [1 2 1 2]'];
    anova.f2    = [anova.f2; [1 1 2 2]'];
    if bootstrap
        % now to bootstrapped parameters
        numtrls = size(subj_dat,1);
        bootstrap_subj_params = NaN(500,6);
        for k = 1:500
            % sample trials with replacement
            bootstrap_subj_dat = datasample(subj_dat,numtrls);
            % use the actual parameters as starting point for the bootstraps
            bootstrap_starting_pt = Finalparams.actual(find(sj==subjects),:);
            [bootstrap_subj_params(k, :), ~] = fit_model(bootstrap_subj_dat, bootstrap_starting_pt);
        end
        CI_params = prctile(bootstrap_subj_params, [100/6, 500/6]);
        Finalparams.lowCI(find(sj==subjects),:) = CI_params(1,:);
        Finalparams.upCI(find(sj==subjects),:) = CI_params(2,:);
    end
end
params = Finalparams.actual;
% compute 2-way ANOVA measures
anov = rm_anova(anova.x, anova.sj, {anova.f1; anova.f2});
if isplot
    % to the plotting function now
    plot_params(Finalparams, anov, bootstrap, trials);
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

function plot_params(boot_params, anov, bootstrap, trials)
% specify the color map and figure properties
cols = linspecer(10, 'qualitative');
colormap(linspecer);
figure;
subplot(4,4,1); hold on;
% weights of second interval
dat1 = boot_params.actual(:, 6);
dat2 = boot_params.actual(:, 5);
if bootstrap
    % get the 66% confidence intervals
    dat11 = [boot_params.lowCI(:,6)'; boot_params.highCI(:,6)'];
    dat22 = [boot_params.lowCI(:,5)'; boot_params.highCI(:,5)'];
else
    dat11 = []; dat22 = [];
end
% polish the figure
if strcmp(trials, 'all')
    set(gca, 'XLim', [-0.25 0.75], 'XTick', -0.25:0.25:0.75,'ylim',[-0.25 0.75], 'ytick', -0.25:0.25:0.75);
    axis square;
    plot([0 00], [-0.25 0.75], 'k', 'LineWidth', 0.25);
    plot([-0.25 0.75], [0 0], 'k', 'LineWidth', 0.25);
else
    set(gca, 'XLim', [-0.4 0.8], 'XTick', -0.4:0.4:0.8,'ylim',[-0.4 0.8], 'ytick', -0.4:0.4:0.8);
    axis square;
    plot([0 00], [-0.4 0.8], 'k', 'LineWidth', 0.25);
    plot([-0.4 0.8], [0 0], 'k', 'LineWidth', 0.25);
end
EquateAxis;
myscatter(dat1, dat2, [dat11 dat22], 75, cols, cols, cols);
% plot the group mean +/- s.e.m
plot([nanmean(dat1)-nansem(dat1) nanmean(dat1)+nansem(dat1)], [nanmean(dat2) nanmean(dat2)], 'Color', cols(1,:),'LineWidth',2);
plot([nanmean(dat1) nanmean(dat1)], [nanmean(dat2)-nansem(dat2) nanmean(dat2)+nansem(dat2)], 'Color', cols(1,:),'LineWidth',2);
[pval] = permtest(dat1, dat2, 0, 100000); % permutation test
xlabel({'Weight for subsequent', 'inconsistent stimulus'});
ylabel({'Weight for subsequent', 'consistent stimulus'});
offsetAxes;
if strcmp(trials, 'all')
title({'Figure 2B',sprintf('Consistent vs. Inconsistent: p = %.4f', pval)});
elseif strcmp(trials, 'error')
title({'Figure S2C; Error trials',sprintf('Consistent vs. Inconsistent: p = %.4f', pval)});
elseif strcmp(trials, 'zero')
title({'Figure S2D; Zero-degree trials',sprintf('Consistent vs. Inconsistent: p = %.4f', pval)});
end

if strcmp(trials, 'all')
    subplot(4,4,3); hold on;
    plot([1 2], nanmean([boot_params.actual(:, 3) boot_params.actual(:, 5)], 1), 'O-', 'Color', [0 0 0], 'MarkerSize', 7.5, 'LineWidth', 1.5, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [1,1,1]);
    plot([1 2], nanmean([boot_params.actual(:, 4) boot_params.actual(:, 6)], 1), 'O-', 'Color', [0.6 0.6 0.6], 'MarkerSize', 7.5, 'LineWidth', 1.5, 'MarkerFaceColor', [0.6 0.6 0.6], 'MarkerEdgeColor', [1,1,1]);
    errbar([1 2], nanmean([boot_params.actual(:, 3) boot_params.actual(:, 5)], 1), nanstd([boot_params.actual(:, 3) boot_params.actual(:, 5)], 1) ./ sqrt(length(boot_params.actual(:, 3))), '-','Color', [0 0 0], 'LineWidth',1);
    errbar([1 2], nanmean([boot_params.actual(:, 4) boot_params.actual(:, 6)], 1), nanstd([boot_params.actual(:, 4) boot_params.actual(:, 6)], 1) ./ sqrt(length(boot_params.actual(:, 4))), '-','Color', [0.6 0.6 0.6], 'LineWidth',1);
    set(gca, 'XLim', [0.5 2.5], 'XTick', 1:2, 'XTickLabel', {'Stimulus 1', 'Stimulus 2'},  'ylim',[0 0.6], 'ytick', 0:0.3:0.6);
    ylabel('Weights');
    axis square;offsetAxes;
    [pval] = permtest(boot_params.actual(:, 3), boot_params.actual(:, 5), 0, 100000); % ranksum much faster than permtest
    mysigstar([1, 2], 0.4, pval, 0);
    [pval] = permtest(boot_params.actual(:, 4), boot_params.actual(:, 6), 0, 100000); % ranksum much faster than permtest
    mysigstar([1, 2], 0, pval, 0);
    [pval] = permtest(boot_params.actual(:, 3), boot_params.actual(:, 4), 0, 100000); % ranksum much faster than permtest
    mysigstar(0.8, 0.4, pval, 0);
    [pval] = permtest(boot_params.actual(:, 5), boot_params.actual(:, 6), 0, 100000); % ranksum much faster than permtest
    mysigstar(2.2, 0.4, pval, 0);
    legend('Consistent trials', 'Inconsistent trials');
    title({'Figure 2C', sprintf('F_{(%d,%d)} = %.2f, p = %.3f', anov.f1xf2.df, anov.f1xf2.fstats, anov.f1xf2.pvalue)});
    
    
    subplot(4,4,9); hold on;
    % estimation noise
    dat1 = boot_params.actual(:, 1);
    dat2 = boot_params.actual(:, 2);
    if bootstrap
        % get the 66% confidence intervals
        dat11 = [boot_params.lowCI(:,1)'; boot_params.highCI(:,1)'];
        dat22 = [boot_params.lowCI(:,2)'; boot_params.highCI(:,2)'];
    else
        dat11 = []; dat22 = [];
    end
    % polish the figure
    set(gca, 'XLim', [0 20], 'XTick', 0:10:20,'ylim',[0 20], 'ytick', 0:10:20);
    axis square;
    EquateAxis;
    myscatter(dat1, dat2, [dat11 dat22], 75, cols, cols, cols);
    % plot the group mean +/- s.e.m
    plot([nanmean(dat1)-nansem(dat1) nanmean(dat1)+nansem(dat1)], [nanmean(dat2) nanmean(dat2)], 'Color', cols(1,:),'LineWidth',2);
    plot([nanmean(dat1) nanmean(dat1)], [nanmean(dat2)-nansem(dat2) nanmean(dat2)+nansem(dat2)], 'Color', cols(1,:),'LineWidth',2);
    [pval] = permtest(dat1, dat2, 0, 100000); % permutation test
    xlabel({'Estimation noise', 'for Inconsistent trials'});
    ylabel({'Estimation noise', 'for Consistent trials'});
    offsetAxes;
    title({'Figure S2A', sprintf('Consistent vs. Inconsistent: p = %.4f', pval)});
end
end