function pdb_Perceptual_ExtendedChoiceSelectiveModel(behdata, selective_params)

% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code calculates
% residual shift parameter for the extended choice-based selective gain model.
global subjects;global psycho_fits;global fixed_params;
% initialise some variables
Finalparams = NaN (length(subjects), 1);
rng shuffle;
global psycho_noise psycho_bias
for sj = subjects
    subj_dat             =  behdata(find(behdata.subj == sj),:);
    % use only choice trials in this paper
    choicetrials = find(subj_dat.condition == 1 & abs(subj_dat.binchoice) == 1);
    subj_dat = subj_dat(choicetrials,:);
    % have the same number of trials across all models for comparing BIC
    % values
    trls2use    = find(subj_dat.x2 ~= 0 & subj_dat.x1 ~= 0);
    subj_dat             = subj_dat(trls2use, :);
    psycho_noise = psycho_fits.logisticFit(find(sj==subjects), 2);
    psycho_bias = -psycho_fits.logisticFit(find(sj==subjects), 1);
    starting_pt = datasample(1:5:25, 1);
    fixed_params = selective_params(find(sj==subjects),:);
    Finalparams(find(sj==subjects),1) = fit_model(subj_dat, starting_pt);
end
plot_residual_shift(Finalparams);
end

function subj_params = fit_model(subj_dat, starting_pt)
% function evaluation params
global dat
options = optimset('Display', 'notify') ;
options.MaxFunEvals = 1e10; % limit the nr of func evals
options.MaxIter = 500000;
options.TolX = 0.00001; % dont make this too small, will take forever to converge
options.TolFun = 0.00001;
options.Robust = 'on';

% define a random starting point for the fitting algorithm
dat = subj_dat;
[individualparams, ~]=subplex('model_Perceptual_ExtendedChoiceSelective', starting_pt);
% optimise again, just to make sure we are at the minimum
subj_params = fminsearchbnd(@(individualparams) model_Perceptual_ExtendedChoiceSelective(individualparams),individualparams,-1000,1000,options);
end

function plot_residual_shift(params)
% specify the color map and figure properties
cols = linspecer(10, 'qualitative');
colormap(linspecer);
figure;
subplot(4,4,1); hold on;
% polish the figure
set(gca, 'XLim', [0. 0.5], 'XTick', [], 'ylim',[-1.5 4.5], 'ytick', -1.5:1.5:4.5);
plot([0. 0.5], [0 0], 'k', 'LineWidth', 0.25);
myscatter(0.25*ones(1,length(params)), params, [], 75, cols, cols, cols);
[pval] = permtest(zeros(1,length(params)), params, 0, 100000); % permutation test
ylabel('Residual Shift parameter');
title(sprintf('p = %.4f', pval));
end