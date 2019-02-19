function roc_index = pdb_ModelFree_Perceptual(behdata, bootstrap, isplot)
% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code calculates model free ROC indices
% and reproduces figure 2B of the paper.
global subjects;
for sj = subjects
    dat             =  behdata(find(behdata.subj == sj),:);
    % use only choice trials in this paper
    trls2use = find(dat.condition == 1 & abs(dat.binchoice) == 1);
    dat = dat(trls2use,:);
    % get the roc indeces for consistent and inconsistent trials
    roc_index.consistent.actual(find(sj==subjects)) = get_roc(dat, 1);
    roc_index.inconsistent.actual(find(sj==subjects)) = get_roc(dat, 0);
    % bootstrap roc indices to obtain confidence intervals
    if bootstrap
        for k = 1:500
            % resample trials with replacement
            numtrls = size(dat,1);
            dat = datasample(dat,numtrls);
            roc_index.consistent.bootstrap(k,find(sj==subjects)) = get_roc(dat, 1);
            roc_index.inconsistent.bootstrap(k,find(sj==subjects)) = get_roc(dat, 0);
        end
    end
end
if isplot
    plot_roc(roc_index, bootstrap);
end
end

function plot_roc(roc_index, bootstrap)
% specify the color map and figure properties
cols = linspecer(10, 'qualitative');
colormap(linspecer);
figure;
subplot(4,4,1); hold on;
dat1 = roc_index.inconsistent.actual;
dat2 = roc_index.consistent.actual;
if bootstrap
    % get the 66% confidence intervals
    dat11 = prctile(roc_index.inconsistent.bootstrap, [100/6, 500/6]);
    dat22 = prctile(roc_index.consistent.bootstrap, [100/6, 500/6]);
else
    dat11 = []; dat22 = [];
end
% polish the figure
set(gca, 'XLim', [0.4 0.7], 'XTick', 0.4:0.1:0.7,'ylim',[0.4 0.7], 'ytick', 0.4:0.1:0.7);
axis square;
plot([0.5 0.5], [0.4 0.7], 'k', 'LineWidth', 0.25);
plot([0.4 0.7], [0.5 0.5], 'k', 'LineWidth', 0.25);
EquateAxis;
myscatter(dat1, dat2, [dat11 dat22], 75, cols, cols, cols);
% plot the group mean +/- s.e.m
plot([nanmean(dat1)-nansem(dat1) nanmean(dat1)+nansem(dat1)], [nanmean(dat2) nanmean(dat2)], 'Color', cols(1,:),'LineWidth',2);
plot([nanmean(dat1) nanmean(dat1)], [nanmean(dat2)-nansem(dat2) nanmean(dat2)+nansem(dat2)], 'Color', cols(1,:),'LineWidth',2);
[pval] = permtest(dat1, dat2, 0, 100000); % permutation test
xlabel({'Sensitivity for subsequent', 'inconsistent stimulus'});
ylabel({'Sensitivity for subsequent', 'consistent stimulus'});
title({'Figure 2F',sprintf('Consistent vs. Inconsistent: p = %.4f', pval)});
end