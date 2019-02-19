function pdb_DistanceMatched_ModelBased_Perceptual(behdata, isplot)

% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code calculates
% model based weights for consistent and inconsistent trials matched in
% terms of the deviation between first and second intervals and and
% reproduces figure S2E  of the paper.
% specify the path to the data

global subjects;global psycho_fits;
% initialise some variables
Finalparams = NaN (length(subjects), 6);
FinalNlogL = NaN(length(subjects),1);
rng shuffle;
global psycho_noise psycho_bias
for sj = subjects
    subj_dat             =  behdata(find(behdata.subj == sj),:);
    % use only choice trials in this paper
    choicetrials = find(subj_dat.condition == 1 & abs(subj_dat.binchoice) == 1);
    subj_dat = subj_dat(choicetrials,:);
    % since consistency cannot be defined for trials where x2 = 0, remove those
    % trials
    trls2use    = find(subj_dat.x2 ~= 0);
    subj_dat             = subj_dat(trls2use, :);
    % now get trials that are distance matched for consistent and
    % inconsistent trials by stratification of consistent and inconsistent
    % trial distributions
    abs_dist = abs(subj_dat.x2 - subj_dat.x1);
    choice_consistent = find(sign(subj_dat.binchoice) == sign(subj_dat.x2));
    choice_inconsistent = find(sign(subj_dat.binchoice) ~= sign(subj_dat.x2));
    matched_trls = [];
    for i = unique(abs_dist)'
        this_dist = find(abs_dist == i);
        this_dist_cons = intersect(this_dist, choice_consistent);
        this_dist_incons = intersect(this_dist, choice_inconsistent);
        [this_dist_cons_indices, this_dist_incons_indices] = stratified_distribution(this_dist_cons, this_dist_incons);
        matched_trls = [matched_trls;this_dist_cons(this_dist_cons_indices);this_dist_incons(this_dist_incons_indices)];
    end
    subj_dat = subj_dat(matched_trls,:);
    % now to fitting, as usual
    psycho_noise = psycho_fits.logisticFit(find(sj==subjects), 2);
    psycho_bias = -psycho_fits.logisticFit(find(sj==subjects), 1);
    starting_pt = [datasample(1:5:25, 1) datasample(1:5:25, 1) datasample(0.05:0.05:1, 1) datasample(0.05:0.05:1, 1) datasample(0.05:0.05:1, 1) datasample(0.05:0.05:1, 1)];
    [Finalparams(find(sj==subjects),:), FinalNlogL(find(sj==subjects))] = fit_model(subj_dat, starting_pt);
end
if isplot
    % to the plotting function now
    plot_params(Finalparams);
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
[individualparams_consistent, ~]=subplex('model_Perceptual_ChoiceSelective', startingpoint);
% optimise again, just to make sure we are at the minimum
[Finalparams_consistent,FinalNlogL_consistent] = fminsearchbnd(@(individualparams) model_Perceptual_ChoiceSelective(individualparams),individualparams_consistent,[0,-1000,-1000],[80,1000,1000],options);

% now fit the inconsistent trials
% define a random starting point for the fitting algorithm
startingpoint = starting_pt([2,4,6]);
TrlsInconsistent = find(sign(subj_dat.binchoice) ~= sign(subj_dat.x2));
dat = subj_dat(TrlsInconsistent,:);
[individualparams_inconsistent, ~] = subplex('model_Perceptual_ChoiceSelective', startingpoint);
% optimise again, just to make sure we are at the minimum
[Finalparams_inconsistent,FinalNlogL_inconsistent] = fminsearchbnd(@(individualparams) model_Perceptual_ChoiceSelective(individualparams),individualparams_inconsistent,[0,-1000,-1000],[80,1000,1000],options);
subj_params([1,3,5]) = Finalparams_consistent;
subj_params([2,4,6]) = Finalparams_inconsistent;
subj_NlogL = FinalNlogL_consistent + FinalNlogL_inconsistent;
end

function plot_params(params)
% specify the color map and figure properties
cols = linspecer(10, 'qualitative');
colormap(linspecer);
figure;
subplot(4,4,1); hold on;
% weights of second interval
dat1 = params(:, 5);
dat2 = params(:, 6);
dat11 = []; dat22 = [];
% polish the figure
set(gca, 'XLim', [-0.2 0.8], 'XTick', -0.2:0.2:0.8,'ylim',[-0.2 0.8], 'ytick', -0.2:0.2:0.8);
axis square;
plot([0 00], [-0.2 0.8], 'k', 'LineWidth', 0.25);
plot([-0.2 0.8], [0 0], 'k', 'LineWidth', 0.25);
EquateAxis;
myscatter(dat1, dat2, [dat11 dat22], 75, cols, cols, cols);
% plot the group mean +/- s.e.m
plot([nanmean(dat1)-nansem(dat1) nanmean(dat1)+nansem(dat1)], [nanmean(dat2) nanmean(dat2)], 'Color', cols(1,:),'LineWidth',2);
plot([nanmean(dat1) nanmean(dat1)], [nanmean(dat2)-nansem(dat2) nanmean(dat2)+nansem(dat2)], 'Color', cols(1,:),'LineWidth',2);
[pval] = permtest(dat1, dat2, 0, 100000); % permutation test
xlabel({'Weight for subsequent', 'inconsistent stimulus'});
ylabel({'Weight for subsequent', 'consistent stimulus'});
offsetAxes;
title({'Figure S2E',sprintf('Consistent vs. Inconsistent: p = %.4f', pval)});
end