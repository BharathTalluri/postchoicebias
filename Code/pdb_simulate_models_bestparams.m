function pdb_simulate_models_bestparams(behdata, params)

% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code estimates ROC
% indices for simulated data, and generates figure 2E.
global subjects;global psycho_fits;
rng shuffle;
global psycho_noise psycho_bias
numiter = 500;
for sj = subjects
    psycho_noise = psycho_fits.logisticFit(find(sj==subjects), 2);
    psycho_bias = -psycho_fits.logisticFit(find(sj==subjects), 1);
    for k = 1:numiter
        
        % simulate data from choice-based selective gain model parameters
        dat             =  behdata(find(behdata.subj == sj),:);
        % use those trials which were used to fit the model
        trls2use = find(dat.condition == 1 & abs(dat.binchoice) == 1 & dat.x2 ~= 0);
        dat = dat(trls2use,:);
        parameters = params.choice_selective;
        dat = simulate_choice_selective(dat, parameters);
        % get the roc indices for consistent and inconsistent trials
        roc_simulated.selective.consistent(k, find(sj==subjects)) = get_roc(dat, 1);
        roc_simulated.selective.inconsistent(k, find(sj==subjects)) = get_roc(dat, 0);
        clear dat parameters
        
        
        % simulate data from shift model parameters
        dat             =  behdata(find(behdata.subj == sj),:);
        % use those trials which were used to fit the model
        trls2use = find(dat.condition == 1 & abs(dat.binchoice) == 1 & dat.x2 ~= 0 & dat.x1 ~= 0);
        dat = dat(trls2use,:);
        parameters = params.shift;
        dat = simulate_shift(dat, parameters);
        % get the roc indices for consistent and inconsistent trials
        roc_simulated.shift.consistent(k, find(sj==subjects)) = get_roc(dat, 1);
        roc_simulated.shift.inconsistent(k, find(sj==subjects)) = get_roc(dat, 0);
        clear dat parameters
        
        % simulate data from correlated noise model parameters
        dat             =  behdata(find(behdata.subj == sj),:);
        % use those trials which were used to fit the model
        trls2use = find(dat.condition == 1 & abs(dat.binchoice) == 1 & dat.x2 ~= 0 & dat.x1 ~= 0);
        dat = dat(trls2use,:);
        parameters = params.correlated_noise;
        dat = simulate_corrNoise(dat, parameters);
        % get the roc indeces for consistent and inconsistent trials
        roc_simulated.corrNoise.consistent(k, find(sj==subjects)) = get_roc(dat, 1);
        roc_simulated.corrNoise.inconsistent(k, find(sj==subjects)) = get_roc(dat, 0);
        clear dat parameters
    end
end
plot_roc(roc_simulated);
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

function plot_roc(roc_simulated)
% specify the color map and figure properties
global subjects
cols = linspecer(10, 'qualitative');
colormap(linspecer);
subjidx = 1:length(subjects);
figure;
% choice-based selective gain model
subplot(4,4,1); hold on;
dat1 = prctile(roc_simulated.selective.inconsistent(:, subjidx), 50);
dat2 = prctile(roc_simulated.selective.consistent(:, subjidx), 50);
% get the 95% confidence intervals
dat11 = prctile(roc_simulated.selective.inconsistent(:, subjidx), [100/6, 500/6]);
dat22 = prctile(roc_simulated.selective.consistent(:, subjidx), [100/6, 500/6]);
% polish the figure
set(gca, 'XLim', [0.45 0.75], 'XTick', 0.45:0.15:0.75,'ylim',[0.45 0.75], 'ytick', 0.45:0.15:0.75);
axis square;
plot([0.5 0.5], [0.45 0.75], 'k', 'LineWidth', 0.25);
plot([0.45 0.75], [0.5 0.5], 'k', 'LineWidth', 0.25);
EquateAxis;
myscatter(dat1, dat2, [dat11 dat22], 75, cols, cols, cols);
% plot the group mean +/- s.e.m
plot([nanmean(dat1)-nansem(dat1) nanmean(dat1)+nansem(dat1)], [nanmean(dat2) nanmean(dat2)], 'Color', [0, 0, 0],'LineWidth',2);
plot([nanmean(dat1) nanmean(dat1)], [nanmean(dat2)-nansem(dat2) nanmean(dat2)+nansem(dat2)], 'Color', [0, 0, 0],'LineWidth',2);
% polish the figure
[pval_simul] = permtest(dat1, dat2, 0, 100000); % permutation test
title({'Choice-based Selective Gain model', sprintf('Consistent vs. Inconsistent: p = %.3f', pval_simul)});

% correlated noise model
subplot(4,4,2); hold on;
dat1 = prctile(roc_simulated.corrNoise.inconsistent(:, subjidx), 50);
dat2 = prctile(roc_simulated.corrNoise.consistent(:, subjidx), 50);
% get the 95% confidence intervals
dat11 = prctile(roc_simulated.corrNoise.inconsistent(:, subjidx), [100/6, 500/6]);
dat22 = prctile(roc_simulated.corrNoise.consistent(:, subjidx), [100/6, 500/6]);
% polish the figure
set(gca, 'XLim', [0.45 0.75], 'XTick', 0.45:0.15:0.75,'ylim',[0.45 0.75], 'ytick', 0.45:0.15:0.75);
axis square;
plot([0.5 0.5], [0.45 0.75], 'k', 'LineWidth', 0.25);
plot([0.45 0.75], [0.5 0.5], 'k', 'LineWidth', 0.25);
EquateAxis;
myscatter(dat1, dat2, [dat11 dat22], 75, cols, cols, cols);
% plot the group mean +/- s.e.m
plot([nanmean(dat1)-nansem(dat1) nanmean(dat1)+nansem(dat1)], [nanmean(dat2) nanmean(dat2)], 'Color', [0, 0, 0],'LineWidth',2);
plot([nanmean(dat1) nanmean(dat1)], [nanmean(dat2)-nansem(dat2) nanmean(dat2)+nansem(dat2)], 'Color', [0, 0, 0],'LineWidth',2);
% polish the figure
[pval_simul] = permtest(dat1, dat2, 0, 100000); % permutation test
title({'Correlated Noise model', sprintf('Consistent vs. Inconsistent: p = %.3f', pval_simul)});

% shift model
subplot(4,4,3); hold on;
dat1 = prctile(roc_simulated.shift.inconsistent(:, subjidx), 50);
dat2 = prctile(roc_simulated.shift.consistent(:, subjidx), 50);
% get the 66% confidence intervals
dat11 = prctile(roc_simulated.shift.inconsistent(:, subjidx), [100/6, 500/6]);
dat22 = prctile(roc_simulated.shift.consistent(:, subjidx), [100/6, 500/6]);
% polish the figure
set(gca, 'XLim', [0.45 0.75], 'XTick', 0.45:0.15:0.75,'ylim',[0.45 0.75], 'ytick', 0.45:0.15:0.75);
axis square;
plot([0.5 0.5], [0.45 0.75], 'k', 'LineWidth', 0.25);
plot([0.45 0.75], [0.5 0.5], 'k', 'LineWidth', 0.25);
EquateAxis;
myscatter(dat1, dat2, [dat11 dat22], 75, cols, cols, cols);
% plot the group mean +/- s.e.m
plot([nanmean(dat1)-nansem(dat1) nanmean(dat1)+nansem(dat1)], [nanmean(dat2) nanmean(dat2)], 'Color', [0, 0, 0],'LineWidth',2);
plot([nanmean(dat1) nanmean(dat1)], [nanmean(dat2)-nansem(dat2) nanmean(dat2)+nansem(dat2)], 'Color', [0, 0, 0],'LineWidth',2);
% polish the figure
[pval_simul] = permtest(dat1, dat2, 0, 100000); % permutation test
title({'Shift model', sprintf('Consistent vs. Inconsistent: p = %.3f', pval_simul)});

suplabel('ROC-index for subsequent inconsistent stimulus', 'x');
suplabel('ROC-index for subsequent consistent stimulus', 'y');
suplabel('Figure 2E; Simulated Estimations', 't');
end