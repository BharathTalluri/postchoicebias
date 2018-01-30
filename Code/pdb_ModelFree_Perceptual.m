function roc_index = pdb_ModelFree_Perceptual()
% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code calculates model free ROC indices
% and reproduces figure 2B of the paper.
close all;
clc; dbstop if error;
% specify the path to the data
behdatapath = '../Data';
behdata = readtable(sprintf('%s/Task_Perceptual.csv', behdatapath));
% initialise some variables
subjects = unique(behdata.subj)';
for sj = subjects
    dat             =  behdata(find(behdata.subj == sj),:);
    % use only choice trials in this paper
    trls2use = find(dat.condition == 1 & abs(dat.binchoice) == 1);
    dat = dat(trls2use,:);
    % get the roc indeces for consistent and inconsistent trials
    roc_index.consistent.actual(sj) = runchoiceselectivegain(dat, 1, 0);
    roc_index.inconsistent.actual(sj) = runchoiceselectivegain(dat, 0, 0);
    % bootstrap roc indices to obtain confidence intervals
    for k = 1:500
        roc_index.consistent.bootstrap(k,sj) = runchoiceselectivegain(dat, 1, 1);
        roc_index.inconsistent.bootstrap(k,sj) = runchoiceselectivegain(dat, 0, 1);
    end
end
plot_roc(roc_index);
end

function roc = runchoiceselectivegain(dat, consistent, bootstrap)
% get the roc indices for choice based selective gain mechanism

dat.idx = transpose(1:height(dat));
% since consistency cannot be defined for trials where x2 = 0, remove those
% trials
choicetrls = find(dat.x2 ~= 0);
dat = dat(choicetrls,:);
if bootstrap
    % resample trials with replacement
    numtrls = size(dat,1);
    dat = datasample(dat,numtrls);
end
% project out the effect of X1 on the estimation
dat.estim = dat.estim-0.5*dat.x1;
switch consistent % split by consistency effect
    case 1
        trls2use    = find(sign(dat.binchoice) == sign(dat.x2));
    case 0
        trls2use    = find(sign(dat.binchoice) ~= sign(dat.x2));
end
dat = dat(trls2use,:);
estimation = dat.estim;
% loop over different x2-pairs
for x2s = 1:2
    % collect the estimation distributions for different values of x2
    switch x2s
        case 1
            trlsP = estimation(find(dat.x2 == 20));
            trlsM = estimation(find(dat.x2 == 10));
        case 2
            trlsP = estimation(find(dat.x2 == -10));
            trlsM = estimation(find(dat.x2 == -20));
    end
    % calculate the roc index for the two distributions
    tmpRoc(x2s)     = rocAnalysis(trlsM, trlsP, 0, 1);
end
% average roc
roc = nanmean([tmpRoc(:).i]);
end

function plot_roc(roc_index)
% specify the color map and figure properties
cols = linspecer(9, 'qualitative');
myfigureprops;
figure;
subplot(4,4,1); hold on;
dat1 = roc_index.inconsistent.actual;
dat2 = roc_index.consistent.actual;
% get the 95% confidence intervals
dat11 = prctile(roc_index.inconsistent.bootstrap, [2.5, 97.5]);
dat22 = prctile(roc_index.consistent.bootstrap, [2.5, 97.5]);
% plot the confidence intervals for each subject
for i = 1:length(dat1)
    plot([dat11(1,i) dat11(2,i)], [dat2(i) dat2(i)], 'Color', cols(2,:),'LineWidth',0.25);
    plot([dat1(i) dat1(i)], [dat22(1,i)  dat22(2,i)], 'Color', cols(2,:),'LineWidth',0.25);
end
scatter(dat1, dat2, 50, 'filled', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', cols(2,:));
% plot the group mean +/- s.e.m
plot([nanmean(dat1)-nansem(dat1) nanmean(dat1)+nansem(dat1)], [nanmean(dat2) nanmean(dat2)], 'Color', cols(1,:),'LineWidth',2);
plot([nanmean(dat1) nanmean(dat1)], [nanmean(dat2)-nansem(dat2) nanmean(dat2)+nansem(dat2)], 'Color', cols(1,:),'LineWidth',2);
% polish the figure
set(gca, 'XLim', [0.3 0.7], 'XTick', 0.3:0.2:0.7,'ylim',[0.3 0.7], 'ytick', 0.3:0.2:0.7);
plot([0.5 0.5], [0.3 0.7], 'k:', 'LineWidth', 0.5);
plot([0.3 0.7], [0.5 0.5], 'k:', 'LineWidth', 0.5);
EquateAxis;
[pval] = permtest(dat1, dat2, 0, 100000); % permutation test
axis square;
xlabel({'Sensitivity for', 'Inconsistent second evidence'});
ylabel({'Sensitivity for', 'Consistent second evidence'});
offsetAxes;
title({'Model-free results',sprintf('Consistent vs. Inconsistent: p = %.4f', pval)});
end