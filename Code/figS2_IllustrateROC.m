function figS2_IllustrateROC()
% Bharath Talluri
% code accompanying the post-decision bias paper. This reproduces
% figure S2 of the paper
% illustrate ROC index calculation using data from subject 09
close all;
clc; dbstop if error;
% specify the color map and figure properties
cols = linspecer(9, 'qualitative');
myfigureprops;
figure;
% specify the path to the data
behdatapath = '../Data';
behdata = readtable(sprintf('%s/Task_Perceptual.csv', behdatapath));
dat = behdata(find(behdata.subj == 9),:);
% get all estimation trials
trls2use = find(dat.condition == 1 & abs(dat.binchoice) == 1);
dat = dat(trls2use,:);
dat.estim = dat.estim - 0.5*dat.x1;
estimation = dat.estim;
% calculate the ROC value for a subset of consistent trials
trlsM = estimation(find(dat.binchoice == 1 & dat.x2 == 10));
trlsP = estimation(find(dat.binchoice == 1 & dat.x2 == 20));
tmproc = rocAnalysis(trlsM, trlsP, 0, 1); % ROC between two distributions
% plot
subplot(4,4,1);hold on;
histogram(trlsM,20);
histogram(trlsP,20);
set(gca, 'XLim', [0 20], 'XTick', 0:10:20, 'ylim',[0 20], 'ytick', 0.0:10:20,'box','off');
xlabel('Estimations (degrees)');
ylabel('Trial Count')
legend('X2 = 10','X2 = 20');
axis square;
subplot(4,4,2);hold on;
plot(tmproc.roc(:,1),tmproc.roc(:,2),'-', 'color', cols(2,:),'LineWidth',2);
for i = 2:length(tmproc.roc(:,1))
    area([tmproc.roc(i-1,1) tmproc.roc(i,1)],[tmproc.roc(i-1,2) tmproc.roc(i,2)],'EdgeColor',[1,1,1],'FaceColor',[0.5 0.5 0.5],'EdgeAlpha',0,'FaceAlpha',0.2);
end
plot([0 1],[0 1],'--', 'color', [0,0,0],'LineWidth',2);
set(gca, 'XLim', [0 1], 'XTick', 0:0.5:1, 'ylim',[0 1], 'ytick', 0.0:0.5:1,'box','off');
xlabel('X2 = 10');
ylabel('X2 = 20');
title(sprintf('ROC index = %.4f', tmproc.i));
axis square;