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
rejectsj = [2,6,12,13];
subjects(rejectsj) = [];
for sj = subjects
    dat             =  behdata(find(behdata.subj == sj),:);
    % use only choice trials in this paper
    trls2use = find(dat.condition == 1 & abs(dat.binchoice) == 1);
    dat = dat(trls2use,:);
    % get the roc indeces for consistent and inconsistent trials
    roc_index.consistent.actual(find(sj==subjects)) = runchoiceselectivegain(dat, 1, 0);
    roc_index.inconsistent.actual(find(sj==subjects)) = runchoiceselectivegain(dat, 0, 0);
    % bootstrap roc indices to obtain confidence intervals
    for k = 1:500
        roc_index.consistent.bootstrap(k,find(sj==subjects)) = runchoiceselectivegain(dat, 1, 1);
        roc_index.inconsistent.bootstrap(k,find(sj==subjects)) = runchoiceselectivegain(dat, 0, 1);
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
switch consistent % split by consistency effect
    case 1
        trls2use    = find(sign(dat.binchoice) == sign(dat.x2));
    case 0
        trls2use    = find(sign(dat.binchoice) ~= sign(dat.x2));
end
dat = dat(trls2use,:);
estimation = dat.estim;
% loop over different x1/x2 pairs
unique_x1 = unique(dat.x1)';
tmpRoc_all = nan(1,length(unique_x1));
num_trls_all = nan(1,length(unique_x1));
for x1 = unique_x1
    num_trls_x1 = NaN(1,2);
    tmpRoc_x1 = NaN(1,2);
    for x2s = 1:2
        % collect the estimation distributions for different values of x2
        switch x2s
            case 1
                trlsP = estimation(find(dat.x2 == 20 & dat.x1 == x1));
                trlsM = estimation(find(dat.x2 == 10 & dat.x1 == x1));
            case 2
                trlsP = estimation(find(dat.x2 == -10 & dat.x1 == x1));
                trlsM = estimation(find(dat.x2 == -20 & dat.x1 == x1));
        end
        % exclude distributions with no or few trials- we cannot get robust roc
        % estimates. This was done after manually inspecting the trial
        % distributions
        if isempty(trlsP) || isempty(trlsM)
            continue
        end
%         if x1 == -10 && consistent == 1 && x2s == 1
%             continue
%         end
%         if x1 == -20 && consistent == 0 && x2s == 2
%             continue
%         end
%          if x1 == 20 && consistent == 0 && x2s == 1
%             continue
%         end
%         if x1 == -10 && consistent == 0 && x2s == 2
%             continue
%         end
%         if x1 == 10 && consistent == 1 && x2s == 2
%             continue
%         end
%         if x1 == 10 && consistent == 0 && x2s == 1
%             continue
%         end
        % calculate the roc index for the two distributions
        tmpRoc     = rocAnalysis(trlsM, trlsP, 0, 1);
        num_trls_x1(x2s) = length(trlsP) + length(trlsM);
        tmpRoc_x1(x2s) = num_trls_x1(x2s)*(tmpRoc.i');
    end
    tmpRoc_all(find(x1 == unique_x1)) = nansum(tmpRoc_x1)/nansum(num_trls_x1);
    num_trls_all(find(x1 == unique_x1)) = nansum(num_trls_x1);
end
% weighted average roc
tmpRoc_all(isnan(tmpRoc_all)) = 0;
roc = (tmpRoc_all*num_trls_all')/sum(num_trls_all);
end

function plot_roc(roc_index)
% specify the color map and figure properties

myfigureprops;
figure;
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2)/4 pos(3:4)*2])
subplot(5,5,1); hold on;
cols = linspecer(10, 'qualitative');
colormap(linspecer);
dat1 = roc_index.inconsistent.actual;
dat2 = roc_index.consistent.actual;
% get the 95% confidence intervals
dat11 = prctile(roc_index.inconsistent.bootstrap, [100/6, 500/6]);
dat22 = prctile(roc_index.consistent.bootstrap, [100/6, 500/6]);
% polish the figure
set(gca, 'XLim', [0.4 0.7], 'XTick', 0.4:0.1:0.7,'ylim',[0.4 0.7], 'ytick', 0.4:0.1:0.7);
axis square;
EquateAxis;
plot([0.5 0.5], [0.45 0.75], 'k', 'LineWidth', 0.25);
plot([0.45 0.75], [0.5 0.5], 'k', 'LineWidth', 0.25);
% plot the confidence intervals for each subject
for i = 1:length(dat1)
    plot([dat11(1,i) dat11(2,i)], [dat2(i) dat2(i)], 'Color', cols(i,:),'LineWidth',0.25);
    plot([dat1(i) dat1(i)], [dat22(1,i)  dat22(2,i)], 'Color', cols(i,:),'LineWidth',0.25);
end
scatter(dat1, dat2, 50, cols, 'filled', 'MarkerEdgeColor', [1 1 1]);
% plot the group mean +/- s.e.m
plot([nanmean(dat1)-nansem(dat1) nanmean(dat1)+nansem(dat1)], [nanmean(dat2) nanmean(dat2)], 'Color', [0 0 0],'LineWidth',2);
plot([nanmean(dat1) nanmean(dat1)], [nanmean(dat2)-nansem(dat2) nanmean(dat2)+nansem(dat2)], 'Color', [0 0 0],'LineWidth',2);
[pval] = permtest(dat1, dat2, 0, 100000); % permutation test
axis square;
xlabel({'Sensitivity for subsequent', 'inconsistent stimulus'});
ylabel({'Sensitivity for subsequent', 'consistent stimulus'});
offsetAxes;
title({'Model-free results',sprintf('Consistent vs. Inconsistent: p = %.4f', pval)});
end