function pdb_individual_Behaviour(behdata, isplot)
% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code reproduces
% figure S1A, B of the paper.

% specify the color map and figure properties
cols = linspecer(9, 'qualitative');
figure;
% initialise some variables
subjects = unique(behdata.subj)';
x = -20:0.1:20;
xavg = [-20 -15 -10 -5 0 5 10 15 20];
subplotnum_all = [1 2 3 4 5 11 12 13 14 15 21 22 23 24 25];
for sj = subjects
    subplotnum = subplotnum_all(sj);
    dat             =  behdata(find(behdata.subj == sj),:);
    % we get all behavioural measures for the perceptual task
    % use only choice trials in this paper
    trls2use = find(abs(dat.condition) == 1 & abs(dat.binchoice) == 1);
    dat = dat(trls2use,:);
    % get the indices for Choice trials
    estimtrls = find((dat.condition == 1 & abs(dat.binchoice) == 1));
    % Psychometric function, we define the psychometric function below
    [bias, slope] = fitcumnormal((dat.x1), dat.binchoice>0);
    logisticFit   = [bias, slope];
    Noise(sj)   = slope;
    dat2 = dat.binchoice;
    dat2(dat2<0) = 0;
    logisticPoints = splitapply(@nanmean, dat2, findgroups(dat.x1));
    % Estimations
    dat = dat(estimtrls,:);
    estimations = splitapply(@nanmedian, dat.estim, findgroups(dat.xavg));
    mdl = LinearModel.fit(dat.xavg, dat.estim, 'linear', 'RobustOpts', 'on');
    Beta(sj) = mdl.Coefficients.Estimate(2:end);
    if isplot
        % psychometric function
        subplot(6,5,subplotnum); hold on;
        plot([0 0],[0 1],'-','color',cols(8,:),'LineWidth',0.5);
        plot([-20 20],[0.5 0.5],'--','color',cols(8,:),'LineWidth',0.25);
        plot([-20 -10 0 10 20], logisticPoints,'o','MarkerFaceColor', cols(8,:), 'MarkerEdgeColor', cols(8,:));
        plot(x, cumnormal(logisticFit, x), '-', 'color', [0 0 0],'MarkerSize',3);
        set(gca, 'XLim', [-20 20], 'XTick', -20:10:20, 'ylim',[0 1], 'ytick', 0.0:0.5:1.0);
        xlabel('Direction in Interval 1'); ylabel('Proportion CW choices');
        title(sprintf('Noise = %1.2f', Noise(sj)));
        axis square;
        
        % Estimations as a function of mean evidence
        subplot(6,5,subplotnum+5); hold on;
        plot([0 0],[-20 20],'-','color',cols(8,:),'LineWidth',0.25);
        plot([-20 20],[0 0],'-','color',cols(8,:),'LineWidth',0.25);
        plot(xavg,estimations,'o-','MarkerFaceColor', cols(8,:), 'LineWidth',1, 'MarkerEdgeColor', cols(8,:));
        set(gca, 'XLim', [-20 20], 'XTick', -20:20:20, 'ylim',[-20 20], 'ytick', -20:20:20);
        ylabel('Estimation  (degrees)');
        xlabel('Mean Direction across interval 1&2  (degrees)');
        title(sprintf('Beta = %1.2f', Beta(sj)));
        axis square;EquateAxis;
    end
end
%% SHOW THE DATA
subplot(6,5,30); hold on;
cols = linspecer(9, 'qualitative');
colormap(linspecer);
y = ones(1,length(Beta)) + 0.25*rand(1, length(Beta));
x =  Beta;
scatter(x,y,60,cols(8,:), 'filled', 'MarkerEdgeColor', [1 1 1]);
plot([0.3 0.3], [0 1.5], 'k', 'LineWidth', 0.25);
set(gca, 'YLim', [0.5 1.5], 'XLim', [0 1.2], 'YTick', 1.0, 'YTickLabel', {'Beta'}, 'XTick', 0:0.3:1.2, 'XTickLabel', 0:0.3:1.2);
xlabel('Beta');
title('Slope of best fitting line');
end

%% Functions used above
function [bias, slope] = fitcumnormal(x,y)
% find the best fitting parameters for the psychometric function by
% maximising the loglikelihood across 10 random starting points for bias
% and slope parameters
niter = 10;
bias_start = datasample(-20:2:20,niter);
slope_start = datasample(0:5:50,niter);
iter_params = NaN(niter,2);
NegLLiter = NaN(niter,1);
for i = 1:niter
    iter_params(i,:) = fminsearchbnd(@(p) cumnormal_LL(p, ...
        x, y), [bias_start(i) slope_start(i)], [-100 0], [100 100],optimset('MaxFunEvals',1000000,'MaxIter',100000));
    NegLLiter(i) = cumnormal_LL(iter_params(i,:),x,y);
end
[~,idx] = min(NegLLiter);
bias        = iter_params(idx,1);
slope       = iter_params(idx,2);
end

function NegLL = cumnormal_LL(p, inputs, responses)
% see http://courses.washington.edu/matlab1/Lesson_5.html#1
% compute the vector of responses for each level of intensity
w   = cumnormal(p, inputs);
% negative loglikelihood, to be minimised
NegLL = -sum(responses.*log(w) + (1-responses).*log(1-w));
end


function y = cumnormal(p, x)
% Parameters: p(1) bias
%             p(2) slope
%             x   inputs.
y = normcdf(x,p(1),p(2));
end