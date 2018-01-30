function roc_index = pdb_ModelFree_Cognitive()
% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code calculates model free ROC indices
% and reproduces figure 3C of the paper.
close all;
clc; dbstop if error;
% specify the path to the data
behdatapath = '../Data';
behdata = readtable(sprintf('%s/Task_Cognitive.csv', behdatapath));
% initialise some variables
subjects = unique(behdata.subj)';
for sj = subjects
    dat             =  behdata(find(behdata.subj == sj),:);
    % some preprocessing of the data is required for the analysis here
    dat.x1_relative = dat.x1 - 50;
    dat.x2_relative = dat.x2 - 50;
    dat.binchoice(dat.binchoice == 1) = -1;
    dat.binchoice(dat.binchoice == 2) = 1;
    estimtrls = find(~isnan(dat.estim));
    dat = dat(estimtrls,:);
    % remove outliers
    estim = dat.estim;
    temp1 = iqr(estim); % get the inter quartile range
    temp2 = quantile(estim,3);%get quartiles
    usetrls = find((estim < temp2(3) + 1.5*temp1) & (estim > temp2(1) - 1.5*temp1)); % remove outliers
    dat = dat(usetrls,:);
    % use only choice trials in this paper
    trls2use = find(dat.condition == 1 & abs(dat.binchoice) == 1);
    dat = dat(trls2use,:);
    % get the roc indeces for consistent and inconsistent trials
    roc_index.consistent.actual(sj) = runchoiceselectivegain(dat, 1);
    roc_index.inconsistent.actual(sj) = runchoiceselectivegain(dat, 0);
end
plot_roc(roc_index);
end

function roc = runchoiceselectivegain(dat, consistent)
% get the roc indices for choice based selective gain mechanism
dat.idx = transpose(1:height(dat));
% since consistency cannot be defined for trials where x2 = 0, remove those
% trials
choicetrls = find(dat.x2_relative ~= 0);
dat = dat(choicetrls,:);
% project out the effect of X1 on the estimation
dat.estim = dat.estim-0.5*dat.x1;
switch consistent % split by consistency effect
    case 1
        trls2use    = find(sign(dat.binchoice) == sign(dat.x2_relative));
    case 0
        trls2use    = find(sign(dat.binchoice) ~= sign(dat.x2_relative));
end
dat = dat(trls2use,:);
estimation = dat.estim;
% loop over different x2-pairs
for x2s = 1:2
    % collect the estimation distributions for different values of x2
    switch x2s
        case 1
            trlsM = estimation(find(dat.x2_relative < -7 & dat.x2_relative >= -13));
            trlsP = estimation(find(dat.x2_relative < -1 & dat.x2_relative >= -6));
        case 2
            trlsM = estimation(find(dat.x2_relative < 6 & dat.x2_relative >= 0));
            trlsP = estimation(find(dat.x2_relative <= 13 & dat.x2_relative >= 7));
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
scatter(dat1, dat2, 50, 'filled', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', cols(2,:));
% plot the group mean +/- s.e.m
plot([nanmean(dat1)-nansem(dat1) nanmean(dat1)+nansem(dat1)], [nanmean(dat2) nanmean(dat2)], 'Color', cols(1,:),'LineWidth',2);
plot([nanmean(dat1) nanmean(dat1)], [nanmean(dat2)-nansem(dat2) nanmean(dat2)+nansem(dat2)], 'Color', cols(1,:),'LineWidth',2);
% polish the figure
set(gca, 'XLim', [0.5 0.8], 'XTick', 0.5:0.1:0.8,'ylim',[0.5 0.8], 'ytick', 0.5:0.1:0.8);
EquateAxis;
[pval] = permtest(dat1, dat2, 0, 100000); % permutation test
axis square;
xlabel({'Sensitivity for', 'Inconsistent second evidence'});
ylabel({'Sensitivity for', 'Consistent second evidence'});
offsetAxes;
title({'Model-free results',sprintf('Consistent vs. Inconsistent: p = %.4f', pval)});
end