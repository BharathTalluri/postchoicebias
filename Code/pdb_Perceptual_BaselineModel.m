function [params, BIC] = pdb_Perceptual_BaselineModel(behdata)

% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code calculates
% model based weights for the baseline model.
global subjects;global psycho_fits;
% initialise some variables
Finalparams = NaN (length(subjects), 3);
FinalNlogL = NaN(length(subjects),1);
BIC = NaN(length(subjects),1);
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
    starting_pt = [datasample(1:5:25, 1) datasample(0.05:0.05:1, 1) datasample(0.05:0.05:1, 1)];
    [Finalparams(find(sj==subjects),:), FinalNlogL(find(sj==subjects))] = fit_model(subj_dat, starting_pt);
    [~, BIC(find(sj==subjects))] = aicbic(-FinalNlogL(find(sj==subjects)), size(Finalparams(find(sj==subjects),:),2) , size(subj_dat, 1));
end
params = Finalparams;
end

function [subj_params, subj_NlogL] = fit_model(subj_dat, starting_pt)
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
[individualparams, ~]=subplex('model_Perceptual_Baseline', starting_pt);
% optimise again, just to make sure we are at the minimum
[subj_params,subj_NlogL] = fminsearchbnd(@(individualparams) model_Perceptual_Baseline(individualparams),individualparams,[0,-1000,-1000],[80,1000,1000],options);
end