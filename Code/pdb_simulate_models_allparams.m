function pdb_simulate_models_allparams(behdata,params)

% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code estimates ROC
% indices for simulated data, and generates figure S2B.
global psycho_fits;
rng shuffle;
global psycho_noise psycho_bias
numiter = 100;
psycho_noise = nanmean(psycho_fits.logisticFit(:, 2));
psycho_bias = -nanmean(psycho_fits.logisticFit(:, 1));
dat             =  behdata;
% use those trials which were used to fit the model
trls2use = find(dat.condition == 1 & abs(dat.binchoice) == 1 & dat.x2 ~= 0);
dat = dat(trls2use,:);
% simulate data from choice-based selective gain model
simulated_consistency.selectivegain = simulate_allParams_choice_selective(dat, params.choice_selective, numiter);
% simulate data from shift model
simulated_consistency.shift = simulate_allParams_shift(dat, params.shift, numiter);
% simulate data from baseline model
simulated_consistency.baseline = simulate_allParams_baseline(dat, params.baseline, numiter);
% simulate data from correlated noise model
simulated_consistency.corrNoise = simulate_allParams_corrNoise(dat, params.correlated_noise, numiter);
% simulate data from extended conditioned perception model (Luu & Stocker, 2018, eLife)
simulated_consistency.conditioned_Perception = simulate_allParams_conditioned_perception(dat, numiter);
% plot heatmaps
plot_consistency(simulated_consistency);
end

function consistency = simulate_allParams_choice_selective(dat, params, numiter)
cons_W1 = 0:0.2:1;
cons_W2 = 0:0.2:1;
incons_W1 = 0:0.2:1;
incons_W2 = 0:0.2:1;
k = 1;
for w1_cons = cons_W1
    for w1_incons = incons_W1
        for w2_cons = cons_W2
            for w2_incons = incons_W2
                % fix the noise parameters at the mean across the group
                cons_noise = nanmean(params(:,1));
                incons_noise = nanmean(params(:,2));
                parameters = [cons_noise incons_noise w1_cons w1_incons w2_cons w2_incons];
                for iter = 1:numiter
                    data_simulated = simulate_choice_selective(dat, parameters);
                    % get the roc indices for consistent and inconsistent trials
                    roc_consistent(iter) = get_roc(data_simulated, 1);
                    roc_inconsistent(iter) = get_roc(data_simulated, 0);
                end
                consistency(k) = nanmean(roc_consistent-roc_inconsistent);
                k = k + 1;
            end
        end
    end
end
clear roc_consistent roc_inconsistent
consistency = reshape(consistency, length(incons_W2), length(cons_W2), length(incons_W1), length(cons_W1));
% marginalise across irrelevant parameters; we are interested in how the
% consistency landscape changes with weights in the interval 2 only. We can
% also marginalise across other parameters by changing the order of permute
consistency = nanmean(permute(consistency, [2,1,4,3]),4);
consistency = nanmean(consistency,3); % consistency is now a 2D matrix, with cons_W2 along the rows and incons_W2 along columns.
end

function data_simulated = simulate_choice_selective(data_actual, parameters)
global psycho_noise psycho_bias
% simulate the selective model
rng shuffle
% collect all the variables to work with
data_simulated = data_actual;
% simulate data
X1              = data_actual.x1;
RealDecision    = data_actual.binchoice;
All_NoisyX1 = NaN(length(X1),1);
% generate a distribution of noisy stimulus representations, and select
% those noisy representations whose sign matches the actual binary choice
% made by the subject
for i = 1:length(X1)
    while 1
        All_NoisyX1(i) = bsxfun(@plus, (X1(i) + psycho_bias), psycho_noise* randn);
        if sign(All_NoisyX1(i)) == RealDecision(i)
            break
        end
    end
end
simulated_decision = sign(All_NoisyX1);
data_simulated.binchoice = simulated_decision;

% consistent trials
cons_est_noise = parameters(1);
cons_w1 = parameters(3);
cons_w2 = parameters(5);
TrlsConsistent = find(sign(data_simulated.binchoice) == sign(data_simulated.x2));
X2              = data_simulated.x2(TrlsConsistent);
RealEvaluation  = data_actual.estim(TrlsConsistent);
% generate weighted noisy stimulus representations for both intervals
NoisyX1 = cons_w1*All_NoisyX1(TrlsConsistent);
NoisyX2   = bsxfun(@plus, (X2 + psycho_bias)*cons_w2, abs(cons_w2)*psycho_noise.* randn(length(X2), 1));
% generate the simulated estimations
simulated_estimation = NoisyX1 + NoisyX2 + cons_est_noise.*randn(length(RealEvaluation),1);
data_simulated.estim(TrlsConsistent) = simulated_estimation;

% inconsistent trials
incons_est_noise = parameters(2);
incons_w1 = parameters(4);
incons_w2 = parameters(6);
TrlsInconsistent = find(sign(data_simulated.binchoice) ~= sign(data_simulated.x2));
X2              = data_simulated.x2(TrlsInconsistent);
RealEvaluation  = data_actual.estim(TrlsInconsistent);
% generate weighted noisy stimulus representations for both intervals
NoisyX1 = incons_w1*All_NoisyX1(TrlsInconsistent);
NoisyX2   = bsxfun(@plus, (X2 + psycho_bias)*incons_w2, abs(incons_w2)*psycho_noise.* randn(length(X2), 1));
% generate the simulated estimations
simulated_estimation = NoisyX1 + NoisyX2 + incons_est_noise.*randn(length(RealEvaluation),1);
data_simulated.estim(TrlsInconsistent) = simulated_estimation;
end

function consistency = simulate_allParams_shift(dat, params, numiter)
W1 = 0:0.2:1;
W2 = 0:0.2:1;
Shift = -5:2:5;
k = 1;
for w1 = W1
    for w2 = W2
        for s = Shift
            % fix the noise parameters at the mean across the group
            est_noise = nanmean(params(:,1));
            parameters = [est_noise w1 w2 s];
            for iter = 1:numiter
                data_simulated = simulate_shift(dat, parameters);
                % get the roc indices for consistent and inconsistent trials
                roc_consistent(iter) = get_roc(data_simulated, 1);
                roc_inconsistent(iter) = get_roc(data_simulated, 0);
            end
            consistency(k) = nanmean(roc_consistent-roc_inconsistent);
            k = k + 1;
        end
    end
end
clear roc_consistent roc_inconsistent
consistency = reshape(consistency, length(Shift), length(W2), length(W1));
% marginalise across irrelevant parameters; we are interested in how the
% consistency landscape changes with weights in the interval 2 and shift only. We can
% also marginalise across other parameters by changing the order of permute
consistency = nanmean(permute(consistency, [1,2,3]),3); % consistency is now a 2D matrix, with Shift along the rows and W2 along columns.
end

function data_simulated = simulate_shift(data_actual, parameters)
global psycho_noise psycho_bias
% simulate the shift model
rng shuffle
% collect all the variables to work with
data_simulated = data_actual;
est_noise = parameters(1);
w1 = parameters(2);
w2 = parameters(3);
shift = parameters(4);
% simulate data
X1              = data_actual.x1;
X2              = data_actual.x2;
RealDecision    = data_actual.binchoice;
RealEvaluation  = data_actual.estim;
NoisyX1 = NaN(length(X1),1);
% generate a distribution of noisy stimulus representations, and select
% those noisy representations whose sign matches the actual binary choice
% made by the subject
for i = 1:length(X1)
    while 1
        NoisyX1(i) = bsxfun(@plus, (X1(i) + psycho_bias), psycho_noise* randn);
        if sign(NoisyX1(i)) == RealDecision(i)
            break
        end
    end
end
simulated_decision = sign(NoisyX1);
data_simulated.binchoice = simulated_decision;
% generate weighted noisy stimulus representations for both intervals
NoisyX1 = w1*NoisyX1;
NoisyX2   = bsxfun(@plus, (X2 + psycho_bias)*w2 , abs(w2)*psycho_noise* randn(length(X2), 1));
% generate the simulated estimations
simulated_estimation = NoisyX1 + NoisyX2 + shift.*simulated_decision + est_noise*randn(length(RealEvaluation),1);
data_simulated.estim = simulated_estimation;
end

function consistency = simulate_allParams_baseline(dat, params, numiter)
W1 = 0:0.2:1;
W2 = 0:0.2:1;
k = 1;
for w1 = W1
    for w2 = W2
        % fix the noise parameters at the mean across the group
        est_noise = nanmean(params(:,1));
        parameters = [est_noise w1 w2];
        for iter = 1:numiter
            data_simulated = simulate_baseline(dat, parameters);
            % get the roc indices for consistent and inconsistent trials
            roc_consistent(iter) = get_roc(data_simulated, 1);
            roc_inconsistent(iter) = get_roc(data_simulated, 0);
        end
        consistency(k) = nanmean(roc_consistent-roc_inconsistent);
        k = k + 1;
    end
end
clear roc_consistent roc_inconsistent
consistency = reshape(consistency, length(W2), length(W1));
% no need to marginalise because we have only 2 free parameters
% consistency is now a 2D matrix, with W2 along the rows and W1 along columns.
end

function data_simulated = simulate_baseline(data_actual, parameters)
global psycho_noise psycho_bias
% simulate the baseline model
rng shuffle
% collect all the variables to work with
data_simulated = data_actual;
est_noise = parameters(1);
w1 = parameters(2);
w2 = parameters(3);
% simulate data
X1              = data_actual.x1;
X2              = data_actual.x2;
RealDecision    = data_actual.binchoice;
RealEvaluation  = data_actual.estim;
NoisyX1 = NaN(length(X1),1);
% generate a distribution of noisy stimulus representations, and select
% those noisy representations whose sign matches the actual binary choice
% made by the subject
for i = 1:length(X1)
    while 1
        NoisyX1(i) = bsxfun(@plus, (X1(i) + psycho_bias), psycho_noise* randn);
        if sign(NoisyX1(i)) == RealDecision(i)
            break
        end
    end
end
simulated_decision = sign(NoisyX1);
data_simulated.binchoice = simulated_decision;
% generate weighted noisy stimulus representations for both intervals
NoisyX1 = w1*NoisyX1;
NoisyX2   = bsxfun(@plus, (X2 + psycho_bias)*w2 , abs(w2)*psycho_noise* randn(length(X2), 1));
% generate the simulated estimations
simulated_estimation = NoisyX1 + NoisyX2 + est_noise*randn(length(RealEvaluation),1);
data_simulated.estim = simulated_estimation;
end

function consistency = simulate_allParams_corrNoise(dat, params, numiter)
W1 = 0:0.2:1;
W2 = 0:0.2:1;
Rho = 0:0.2:1;
k = 1;
for w1 = W1
    for w2 = W2
        for rho = Rho
            % fix the noise parameters at the mean across the group
            est_noise = nanmean(params(:,1));
            parameters = [est_noise w1 w2 rho];
            for iter = 1:numiter
                data_simulated = simulate_corrNoise(dat, parameters);
                % get the roc indices for consistent and inconsistent trials
                roc_consistent(iter) = get_roc(data_simulated, 1);
                roc_inconsistent(iter) = get_roc(data_simulated, 0);
            end
            consistency(k) = nanmean(roc_consistent-roc_inconsistent);
            k = k + 1;
        end
    end
end
clear roc_consistent roc_inconsistent
consistency = reshape(consistency, length(Rho), length(W2), length(W1));
% marginalise across irrelevant parameters; we are interested in how the
% consistency landscape changes with weights in the interval 2 and rho only. We can
% also marginalise across other parameters by changing the order of permute
consistency = nanmean(permute(consistency, [1,2,3]),3); % consistency is now a 2D matrix, with Rho along the rows and W2 along columns.
end

function data_simulated = simulate_corrNoise(data_actual, parameters)
global psycho_noise psycho_bias
% simulate the correlated noise model
rng shuffle
% collect all the variables to work with
data_simulated = data_actual;
est_noise = parameters(1);
w1 = parameters(2);
w2 = parameters(3);
corr_param = parameters(4);
% simulate data
X1              = data_actual.x1;
X2              = data_actual.x2;
RealDecision    = data_actual.binchoice;
RealEvaluation  = data_actual.estim;
NoisyX1 = NaN(length(X1),1);
x1_noise = NaN(length(X1),1);
% generate a distribution of noisy stimulus representations, and select
% those noisy representations whose sign matches the actual binary choice
% made by the subject
for i = 1:length(X1)
    while 1
        x1_noise(i) = psycho_noise*randn;
        NoisyX1(i) = X1(i) + psycho_bias + x1_noise(i);
        if sign(NoisyX1(i)) == RealDecision(i)
            break
        end
    end
end
simulated_decision = sign(NoisyX1);
data_simulated.binchoice = simulated_decision;
x2_noise = psycho_noise*randn(length(X2),1);
% generate weighted noisy stimulus representations for both intervals
NoisyX1 = w1*NoisyX1;
NoisyX2 = w2*(X2 + psycho_bias + corr_param*x1_noise + (1-corr_param)*x2_noise);
% generate the simulated estimations
simulated_estimation = NoisyX1 + NoisyX2 + est_noise*randn(length(RealEvaluation),1);
data_simulated.estim = simulated_estimation;
end

function consistency = simulate_allParams_conditioned_perception(dat, numiter)
Input_Noise = [1 3 5:5:20];
Output_Noise = [1 3 5:5:20];
k = 1;
for noise_in = Input_Noise
    for noise_out = Output_Noise
        % fix the noise parameters at the mean across the group
        parameters = [noise_in noise_out rho];
        for iter = 1:numiter
            data_simulated = simulate_conditioned_perception(dat, parameters);
            % get the roc indices for consistent and inconsistent trials
            roc_consistent(iter) = get_roc(data_simulated, 1);
            roc_inconsistent(iter) = get_roc(data_simulated, 0);
        end
        consistency(k) = nanmean(roc_consistent-roc_inconsistent);
        k = k + 1;
    end
end
clear roc_consistent roc_inconsistent
consistency = reshape(consistency, length(Output_Noise), length(Input_Noise));
% no marginalisaion required since there are only 2 parameters
% consistency is now a 2D matrix, with Output Noise along the rows and Input Noise along columns.
end

function data_simulated = simulate_conditioned_perception(data_actual, parameters)
global psycho_bias
% simulate the extended conditioned perception model
% the model uses weights = 0.5 for both intervals, and conditions the pdf
% of both intervals on the choice.
rng shuffle
% collect all the variables to work with
data_simulated = data_actual;
in_noise = parameters(1);
out_noise = parameters(1);
% simulate data
X1              = data_actual.x1;
X2              = data_actual.x2;
RealDecision    = data_actual.binchoice;
% generate a distribution of noisy stimulus representations
step=0.05;
X = -180:step:180;
for i = 1:length(X1)
    % get the pdf for noisy representation of interval 1
    NoisyX1=normpdf(X,(X1(i)+psycho_bias)*0.5,in_noise*0.5);
    % condition NoisyX1 on the binary choice
    NoisyX1(sign(X) ~= sign(RealDecision(i))) = 0;
    NoisyX1 = NoisyX1./trapz(X,NoisyX1); % normalise to get the pdf
    % calculate the mean of the posterior distribution as the estimate: we can
    % also replace this with maximum aposteriori value.
    estim_X1 = sum(NoisyX1.* X) * step;
    % get the pdf for noisy representation of interval 2
    NoisyX2=normpdf(X,(X2(i)+psycho_bias)*0.5,in_noise*0.5);
    % condition NoisyX2 on the binary choice
    NoisyX2(sign(X) ~= sign(RealDecision(i))) = 0;
    NoisyX2 = NoisyX2./trapz(X,NoisyX2); % normalise to get the pdf
    % calculate the mean of the posterior distribution as the estimate: we can
    % also replace this with maximum aposteriori value.
    estim_X2 = sum(NoisyX2.* X) * step;
    % combine estimations from both intervals
    pre_estim = estim_X1 + estim_X2;
    estim = pre_estim + out_noise*randn;
    data_simulated.estim(i) = estim;
end
end

function plot_consistency(simulated_consistency)

figure;
% choice-based selective gain model
subplot(5,5,1); hold on;
imagesc(simulated_consistency.selectivegain, [-0.1, 0.1]);
% polish the figure
lims = size(simulated_consistency.selectivegain,1);
set(gca, 'XLim', [0 lims+0.5], 'XTick', 1:lims, 'XTickLabel', 0:0.2:1, 'YLim', [0 lims+0.5], 'YTick', 1:lims, 'YTickLabel', 0:0.2:1);
axis square;
xlabel('W2 Inconsistent');
ylabel('W2 Consistent');
title({'Choice-based', 'Selective Gain'});

% baseline model
subplot(5,5,2); hold on;
imagesc(simulated_consistency.baseline, [-0.1, 0.1]);
% polish the figure
lims = size(simulated_consistency.baseline,1);
set(gca, 'XLim', [0 lims+0.5], 'XTick', 1:lims, 'XTickLabel', 0:0.2:1, 'YLim', [0 lims+0.5], 'YTick', 1:lims, 'YTickLabel', 0:0.2:1);
axis square;
xlabel('W2');
ylabel('W1');
title('Baseline');

% shift model
subplot(5,5,3); hold on;
imagesc(simulated_consistency.shift, [-0.1, 0.1]);
% polish the figure
lims = size(simulated_consistency.shift,1);
set(gca, 'XLim', [0 lims+0.5], 'XTick', 1:lims, 'XTickLabel', 0:0.2:1, 'YLim', [0 lims+0.5], 'YTick', 1:lims, 'YTickLabel', -5:2:5);
axis square;
xlabel('W2');
ylabel('Shift');
title('Shift');
% draw the colorbar
pos=get(gca,'pos');
colorbar('location','southoutside','position',[pos(1) pos(2)-0.05 pos(3) 0.01]);

% correlated noise model
subplot(5,5,4); hold on;
imagesc(simulated_consistency.corrNoise, [-0.1, 0.1]);
% polish the figure
lims = size(simulated_consistency.corrNoise,1);
set(gca, 'XLim', [0 lims+0.5], 'XTick', 1:lims, 'XTickLabel', 0:0.2:1, 'YLim', [0 lims+0.5], 'YTick', 1:lims, 'YTickLabel', 0:.2:1);
axis square;
xlabel('W2');
ylabel('Rho');
title('Correlated Noise');

% extended conditioned perception model
subplot(5,5,5); hold on;
imagesc(simulated_consistency.conditioned_Perception, [-0.1, 0.1]);
% polish the figure
lims = size(simulated_consistency.conditioned_Perception,1);
set(gca, 'XLim', [0 lims+0.5], 'XTick', 1:lims, 'XTickLabel', [1 3 5:5:20], 'YLim', [0 lims+0.5], 'YTick', 1:lims, 'YTickLabel', [1 3 5:5:20]);
axis square;
xlabel('Input Noise');
ylabel('Output Noise');
title({'Conditioned', 'Perception'});

suplabel({'Figure S2B', 'Simulated Estimations with a range of parameters'}, 't');
end