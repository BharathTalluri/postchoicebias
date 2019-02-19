function pdb_Numerical_ChoiceSelectiveModel(behdata, isplot)

% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code calculates
% model based weights for the cognitive task and reproduces figure 3B.
global subjects;global psycho_fits;
% initialise some variables
Finalparams.actual = NaN (length(subjects), 8);
FinalNlogL = NaN(length(subjects),1);
rng shuffle;
anova.x = [];
anova.sj = [];
anova.f1 = [];
anova.f2 = [];
global psycho_noise psycho_bias
for sj = subjects
    subj_dat             =  behdata(find(behdata.subj == sj),:);
    % use only choice trials in this paper
    choicetrials = find(subj_dat.condition == 1 & abs(subj_dat.binchoice) == 1);
    subj_dat = subj_dat(choicetrials,:);
    % since consistency cannot be defined for trials where x2 = 0, remove those
    % trials
    trls2use    = find(subj_dat.x2_relative ~= 0);
    subj_dat             = subj_dat(trls2use, :);
    psycho_noise = psycho_fits.logisticFit(sj, 2);
    psycho_bias = -psycho_fits.logisticFit(sj, 1);
    starting_pt = [datasample(1:5:25, 1) datasample(1:5:25, 1) datasample(0.05:0.05:1, 1) datasample(0.05:0.05:1, 1) datasample(0.05:0.05:1, 1) datasample(0.05:0.05:1, 1) datasample(0:0.25:5,1) datasample(0:0.25:5,1)];
    [Finalparams.actual(sj,:), FinalNlogL(sj)] = fit_model(subj_dat, starting_pt);
    anova.x     = [anova.x; Finalparams.actual(find(sj == subjects),3); Finalparams.actual(find(sj == subjects),5); Finalparams.actual(find(sj == subjects),4); Finalparams.actual(find(sj == subjects),6)];
    anova.sj    = [anova.sj; find(sj == subjects) * ones(4, 1)];
    anova.f1    = [anova.f1; [1 2 1 2]'];
    anova.f2    = [anova.f2; [1 1 2 2]'];
end
% compute 2-way ANOVA measures
anov = rm_anova(anova.x, anova.sj, {anova.f1; anova.f2});
% to the plotting function now
if isplot
plot_params(Finalparams, anov);
end
end

function [subj_params, subj_NlogL] = fit_model(subj_dat, starting_pt)
% function evaluation params
global dat
subj_params = NaN(1,8);
options = optimset('Display', 'notify') ;
options.MaxFunEvals = 1e10; % limit the nr of func evals
options.MaxIter = 500000;
options.TolX = 0.00001; % dont make this too small, will take forever to converge
options.TolFun = 0.00001;
options.Robust = 'on';

% first fit the consistent trials
% define a random starting point for the fitting algorithm
startingpoint = starting_pt([1,3,5,7]);
TrlsConsistent = find(sign(subj_dat.binchoice) == sign(subj_dat.x2_relative));
dat = subj_dat(TrlsConsistent,:);
[individualparams_consistent, ~]=subplex('model_Numerical_ChoiceSelective', startingpoint);
% optimise again, just to make sure we are at the minimum
[Finalparams_consistent,FinalNlogL_consistent] = fminsearchbnd(@(individualparams) model_Numerical_ChoiceSelective(individualparams),individualparams_consistent,[0,0,0,-10],[80,1000,1000,10],options);

% now fit the inconsistent trials
% define a random starting point for the fitting algorithm
startingpoint = starting_pt([2,4,6,8]);
TrlsInconsistent = find(sign(subj_dat.binchoice) ~= sign(subj_dat.x2_relative));
dat = subj_dat(TrlsInconsistent,:);
[individualparams_inconsistent, ~] = subplex('model_Numerical_ChoiceSelective', startingpoint);
% optimise again, just to make sure we are at the minimum
[Finalparams_inconsistent,FinalNlogL_inconsistent] = fminsearchbnd(@(individualparams) model_Numerical_ChoiceSelective(individualparams),individualparams_inconsistent,[0,0,0,-10],[80,1000,1000,10],options);
subj_params([1,3,5,7]) = Finalparams_consistent;
subj_params([2,4,6,8]) = Finalparams_inconsistent;
subj_NlogL = FinalNlogL_consistent + FinalNlogL_inconsistent;
end

function plot_params(final_params, anov)
% specify the color map and figure properties
myfigureprops;
figure;
colormap('gray');
cols = gray(25);
cols = cols(4:24,:);
subplot(4,4,1); hold on;
% weights of second interval
dat1 = final_params.actual(:, 5);
dat2 = final_params.actual(:, 6);
myscatter(dat1, dat2,[], mrksize, cols,cols,cols);
lim = [0 1.2];
limtick = 0:0.4:1.2;
% plot the group mean +/- s.e.m
plot([nanmean(dat1)-nansem(dat1) nanmean(dat1)+nansem(dat1)], [nanmean(dat2) nanmean(dat2)], 'Color', [0 0 0],'LineWidth',2);
plot([nanmean(dat1) nanmean(dat1)], [nanmean(dat2)-nansem(dat2) nanmean(dat2)+nansem(dat2)], 'Color', [0 0 0],'LineWidth',2);
% polish the figure
set(gca, 'XLim', lim, 'XTick', limtick,'ylim',lim, 'ytick', limtick);
EquateAxis;
[pval] = permtest(dat1, dat2, 0, 100000); % permutation test
axis square;
xlabel({'Weight for subsequent', 'Inconsistent information'});
ylabel({'Weight for subsequent', 'Consistent information'});
title({'Figure 3B',sprintf('Consistent vs. Inconsistent: p = %.4f', pval)});

% supplementary figure
figure;
colormap('gray');
cols = gray(25);
cols = cols(4:24,:);
subplot(5,5,1);hold on;
dat1 = final_params.actual(:, 3);
dat2 = final_params.actual(:, 4);
xlab = {'Weight for first interval', 'in inconsistent trials'};
ylab = {'Weight for first interval', 'in consistent trials'};
lim = [0 1.5];
limtick = 0:0.5:1.5;
set(gca, 'XLim', lim, 'XTick', limtick,'ylim',lim, 'ytick', limtick );
axis square;
EquateAxis;
myscatter(dat1, dat2,[], mrksize, cols,cols,cols);
% plot the group mean +/- s.e.m
plot([nanmean(dat1)-nansem(dat1) nanmean(dat1)+nansem(dat1)], [nanmean(dat2) nanmean(dat2)], 'Color', [0 0 0],'LineWidth',2);
plot([nanmean(dat1) nanmean(dat1)], [nanmean(dat2)-nansem(dat2) nanmean(dat2)+nansem(dat2)], 'Color', [0 0 0],'LineWidth',2);
[pval] = permtest(dat1, dat2, 0, 100000); % permutation test
xlabel(xlab);
ylabel(ylab);
offsetAxes;
title({'Figure S3A' ,sprintf('Consistent vs. Inconsistent: p = %.4f',pval)});

subplot(5,5,2);hold on;
dat1 = final_params.actual(:, 1);
dat2 = final_params.actual(:, 2);
xlab = {'Estimation noise for', 'inconsistent trials'};
ylab = {'Estimation noise for', 'consistent trials'};
lim = [0 8];
limtick = 0:4:8;
set(gca, 'XLim', lim, 'XTick', limtick,'ylim',lim, 'ytick', limtick );
axis square;
EquateAxis;
myscatter(dat1, dat2,[], mrksize, cols,cols,cols);
% scatter(dat1, dat2, 50, 'filled', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', cols(4,:));
% plot the group mean +/- s.e.m
plot([nanmean(dat1)-nansem(dat1) nanmean(dat1)+nansem(dat1)], [nanmean(dat2) nanmean(dat2)], 'Color', [0 0 0],'LineWidth',2);
plot([nanmean(dat1) nanmean(dat1)], [nanmean(dat2)-nansem(dat2) nanmean(dat2)+nansem(dat2)], 'Color', [0 0 0],'LineWidth',2);
[pval] = permtest(dat1, dat2, 0, 100000); % permutation test
% [pval] = signrank(dat1, dat2, 'method','exact');
xlabel(xlab);
ylabel(ylab);
offsetAxes;
title({'Figure S3B' ,sprintf('Consistent vs. Inconsistent: p = %.4f',pval)});

subplot(5,5,3);hold on;
dat1 = final_params.actual(:, 7);
dat2 = final_params.actual(:, 8);
xlab = {'Theta parameter', 'in inconsistent trials'};
ylab = {'Theta parameter', 'in consistent trials'};
lim = [-10 10];
limtick = -10:5:10;
set(gca, 'XLim', lim, 'XTick', limtick,'ylim',lim, 'ytick', limtick );
axis square;
EquateAxis;
plot([0,0], lim, 'k', 'LineWidth', 0.25);
plot(lim, [0,0], 'k', 'LineWidth', 0.25);
myscatter(dat1, dat2,[], mrksize, cols,cols,cols);
% plot the group mean +/- s.e.m
plot([nanmean(dat1)-nansem(dat1) nanmean(dat1)+nansem(dat1)], [nanmean(dat2) nanmean(dat2)], 'Color', [0 0 0],'LineWidth',2);
plot([nanmean(dat1) nanmean(dat1)], [nanmean(dat2)-nansem(dat2) nanmean(dat2)+nansem(dat2)], 'Color', [0 0 0],'LineWidth',2);
[pval] = permtest(dat1, dat2, 0, 100000); % permutation test
xlabel(xlab);
ylabel(ylab);
offsetAxes;
title({'Figure S3C' ,sprintf('Consistent vs. Inconsistent: p = %.4f',pval)});

subplot(5,5,5); hold on;
plot([1 2], nanmean([final_params.actual(:, 4) final_params.actual(:, 6)], 1), 'O-', 'Color', [0 0 0], 'MarkerSize', 7.5, 'LineWidth', 1.5, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [1,1,1]);
errbar([1 2], nanmean([final_params.actual(:, 4) final_params.actual(:, 6)], 1), nanstd([final_params.actual(:, 4) final_params.actual(:, 6)], 1) ./ sqrt(length(final_params.actual(:, 1))), '-','Color', [0 0 0], 'LineWidth',1);
plot([1 2], nanmean([final_params.actual(:, 3) final_params.actual(:, 5)], 1), 'O-', 'Color', cols(14,:), 'MarkerSize', 7.5, 'LineWidth', 1.5, 'MarkerFaceColor', cols(14,:), 'MarkerEdgeColor', [1,1,1]);
errbar([1 2], nanmean([final_params.actual(:, 3) final_params.actual(:, 5)], 1), nanstd([final_params.actual(:, 3) final_params.actual(:, 5)], 1) ./ sqrt(length(final_params(:,1))), '-','Color', cols(14,:), 'LineWidth',1);
set(gca, 'XLim', [0.5 2.5], 'XTick', 1:2, 'XTickLabel', {'Stimulus 1', 'Stimulus 2'},  'ylim',[0.2 0.6], 'ytick', 0.2:0.2:0.6);
ylabel('Weights');
xlabel('Interval');
axis square;offsetAxes;
[pval] = permtest(params.w1c, params.w2c, 0, 100000); % ranksum much faster than permtest
mysigstar([1, 2], 0.6, pval, 0);
[pval] = permtest(params.w1i, params.w2i, 0, 100000); % ranksum much faster than permtest
mysigstar([1, 2], 0.1, pval, 0);
[pval] = permtest(params.w1c, params.w1i, 0, 100000); % ranksum much faster than permtest
mysigstar(0.8, 0.6, pval, 0);
[pval] = permtest(params.w2c, params.w2i, 0, 100000); % ranksum much faster than permtest
mysigstar(2.2, 0.6, pval, 0);
title({'Figure S3D', sprintf('F_{(%d,%d)} = %.2f, p = %.3f', anov.f1xf2.df, anov.f1xf2.fstats, anov.f1xf2.pvalue)});
end