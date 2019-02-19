function HighLevel()
% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code reproduces all
% the results and figures in the paper

% do some polishing and initialise variables
close all;clear;
clc; dbstop if error;
myfigureprops;
% specify the path to the data
global behdatapath; global subjects;global psycho_fits;
behdatapath = '../Data';
toolpath = 'Tools/';
addpath(genpath(sprintf('%s', toolpath)));

%% PERCEPTUAL TASK %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
behdata = readtable(sprintf('%s/Task_Perceptual.csv', behdatapath));
% plot individual behavioural data for all subjects, and identify the
% outliers. Reproduce Figure S1A, B.
pdb_individual_Behaviour(behdata, 1);
subjects = unique(behdata.subj)';
rejectsj = [2,6,12,13];
subjects(rejectsj) = [];
isplot = 1; % plot the results
% obtain the behavioural measures and generate Figure 1B, C.
psycho_fits = pdb_Behaviour('Perceptual', isplot);
% obtain the trial distribution and generate Figure S1C, D
pdb_permutationHistogram(behdata);
% Choice-based selective gain model
bootstrap = 0; % this can be set to 1 if 66% confidence intervals for each subjects' parameter estimates are to be plotted.
compare = 0; % we will use a subset of trials fitted here for model comparison. Set this to 0 now.
% get model-based measures and generate Figures 2B, C & S2A.
[params.choice_selective, ~] = pdb_Perceptual_ChoiceSelectiveModel(behdata, 'all', bootstrap, compare, isplot);
% get model-free measures and generate Figure 2F
roc_index.choice_selective = pdb_ModelFree_Perceptual(behdata, bootstrap, isplot); 
% get the measures of consistency to correlate, for Figure 3D
consistency.modelbased.perceptual = params.choice_selective(:,5) - params.choice_selective(:,6);
consistency.modelfree.perceptual = roc_index.choice_selective.consistent.actual - roc_index.choice_selective.inconsistent.actual;
% Model comparisons
isplot = 0; % set this to 0 for model comparisons, and to 1 to inspect parameters of different models
compare = 1;
bootstrap = 0;
[params.stimulus_selective, BIC.stimulus_selective] = pdb_Perceptual_StimulusSelectiveModel(behdata, isplot);
[~, BIC.choice_selective] = pdb_Perceptual_ChoiceSelectiveModel(behdata, 'all', bootstrap, compare, isplot);
[params.shift, BIC.shift] = pdb_Perceptual_ShiftModel(behdata);
[params.baseline, BIC.baseline] = pdb_Perceptual_BaselineModel(behdata);
[params.correlated_noise, BIC.correlated_noise] = pdb_Perceptual_CorrelatedNoiseModel(behdata);
% plot the delta BIC values and generate Figure 2A.
pdb_model_comparison(BIC); 
% compare consistency effect from the two selective models, and generate Figure 2D
pdb_compare_selective_consistency(params.choice_selective(:,5) - params.choice_selective(:,6), params.stimulus_selective(:,5) - params.stimulus_selective(:,6));
% simulate data with best fitting parameters of each model, and investigate
% ROC effect in the simulated data, and generate Figure 2E
pdb_simulate_models_bestparams(behdata, params);
% get model-based measures for error trials, and zero-degree trials, and generate Figures S2C, D.
compare = 0;
isplot = 1;
pdb_Perceptual_ChoiceSelectiveModel(behdata, 'error', bootstrap, compare, isplot);
pdb_Perceptual_ChoiceSelectiveModel(behdata, 'zero', bootstrap, compare, isplot);
% generate consistency heatmaps for simulated data with a range of
% parameters in various models, and generate Figure S2B.
pdb_simulate_models_allparams(behdata, params);
% fit the choice-based selective gain model on data with distance-matched consistent and inconsistent trial distributions, and generate Figures S2E.
pdb_DistanceMatched_ModelBased_Perceptual(behdata, isplot)
% fit theresidual shift model on data, and generate Figures S2F.
pdb_Perceptual_ExtendedChoiceSelectiveModel(behdata, params.choice_selective);

%% NUMERICAL TASK %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We do similar analysis on the data from numerical integration task from
% Bronfman et al. 2015.
behdata = readtable(sprintf('%s/Task_Numerical.csv', behdatapath));
subjects = unique(behdata.subj)';
isplot = 0;
% obtain the behavioural measures
psycho_fits = pdb_Behaviour('Numerical', isplot);
isplot = 1;
% Choice-based selective gain model
% get model-based measures and generate Figures 3B, S3A-D.
[params.choice_selective, ~] = pdb_Numerical_ChoiceSelectiveModel(behdata, 'all', isplot);
% get model-free measures and generate Figure 3C
roc_index.choice_selective = pdb_ModelFree_Numerical(behdata, isplot); 
% get the measures of consistency to correlate, for Figure 3D
consistency.modelbased.numerical = params.choice_selective(:,5) - params.choice_selective(:,6);
consistency.modelfree.numerical = roc_index.choice_selective.consistent.actual - roc_index.choice_selective.inconsistent.actual;
% correlation between the roc indices and weights; reproduce figure 3D
pdb_ROC_vs_Weights(consistency)
% finally, use a stratified model by matching the means of consistent and
% inconsistent distributions, to show that the results hold in a fixed-effects subjects even without
% the theta parameter
pdb_Numerical_Stratified_ChoiceSelectiveModel(behdata)