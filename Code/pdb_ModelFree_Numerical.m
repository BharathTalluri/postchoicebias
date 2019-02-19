function roc_index = pdb_ModelFree_Numerical(behdata, isplot)
% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code calculates model free ROC indices
% and reproduces figure 3C of the paper.
% close all;
global subjects;
for sj = subjects
    dat             =  behdata(find(behdata.subj == sj),:);
    % use only choice trials in this paper
    trls2use = find(dat.condition == 1 & abs(dat.binchoice) == 1);
    dat = dat(trls2use,:);
    % get the roc indeces for consistent and inconsistent trials
    roc_index.consistent.actual(sj) = roc_choiceselectivegain(dat, 1);
    roc_index.inconsistent.actual(sj) = roc_choiceselectivegain(dat, 0);
end
if isplot
    plot_roc(roc_index);
end
end

function roc = roc_choiceselectivegain(dat, consistent)
% get the roc indices for choice based selective gain mechanism
dat.idx = transpose(1:height(dat));
% since consistency cannot be defined for trials where x2 = 0, remove those
% trials
choicetrls = find(dat.x2_relative ~= 0);
dat = dat(choicetrls,:);
switch consistent % split by consistency effect
    case 1
        trls2use    = find(sign(dat.binchoice) == sign(dat.x2_relative));
    case 0
        trls2use    = find(sign(dat.binchoice) ~= sign(dat.x2_relative));
end
%binarise x1
dat.x1_bin = NaN(size(dat.x1_relative));
dat.x1_bin(find(dat.x1 < 50)) = -1;
dat.x1_bin(find(dat.x1 > 50)) = 1;
dat = dat(trls2use,:);
estimation = dat.estim;
% loop over different x1/x2 pairs
unique_x1 = unique(dat.x1_bin)';
tmpRoc_all = nan(1,length(unique_x1));
num_trls_all = nan(1,length(unique_x1));
for x1 = unique_x1
    num_trls_x1 = NaN(1,2);
    tmpRoc_x1 = NaN(1,2);
    for x2s = 1:2
        % collect the estimation distributions for different values of x2
        switch x2s
            case 1
                trlsM = estimation(find(dat.x2_relative < -7 & dat.x2_relative >= -13 & dat.x1_bin == x1));
                trlsP = estimation(find(dat.x2_relative < -1 & dat.x2_relative >= -6 & dat.x1_bin == x1));
            case 2
                trlsM = estimation(find(dat.x2_relative < 6 & dat.x2_relative >= 0 & dat.x1_bin == x1));
                trlsP = estimation(find(dat.x2_relative <= 13 & dat.x2_relative >= 7 & dat.x1_bin == x1));
        end
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
colormap('gray');
cols = gray(25);
cols = cols(4:24,:);
myfigureprops;
figure;
subplot(4,4,1); hold on;
dat1 = roc_index.inconsistent.actual;
dat2 = roc_index.consistent.actual;
% polish the figure
set(gca, 'XLim', [0.4 0.8], 'XTick', 0.4:0.1:0.8,'ylim',[0.4 0.8], 'ytick', 0.4:0.1:0.8);
axis square;
EquateAxis;
plot([0.5 0.5], [0.4 0.8], 'k', 'LineWidth', 0.25);
plot([0.4 0.8], [0.5 0.5], 'k', 'LineWidth', 0.25);
scatter(dat1, dat2, 50, cols, 'filled', 'MarkerEdgeColor', [1 1 1]);
% plot the group mean +/- s.e.m
plot([nanmean(dat1)-nansem(dat1) nanmean(dat1)+nansem(dat1)], [nanmean(dat2) nanmean(dat2)], 'Color', [0 0 0],'LineWidth',2);
plot([nanmean(dat1) nanmean(dat1)], [nanmean(dat2)-nansem(dat2) nanmean(dat2)+nansem(dat2)], 'Color', [0 0 0],'LineWidth',2);
[pval] = permtest(dat1, dat2, 0, 100000); % permutation test
xlabel({'ROC-index for subsequent', 'inconsistent information'});
ylabel({'Roc-index for subsequent', 'consistent information'});
title({'Figure 3C',sprintf('Consistent vs. Inconsistent: p = %.4f', pval)});
end