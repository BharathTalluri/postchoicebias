function optimal_funcval = model_Perceptual_CorrelatedNoise(individualparams)
global dat;global psycho_noise; global psycho_bias;
niter = 1e4;
% get the relevant values first
X1 = dat.x1;
X2 = dat.x2;
RealDecision = dat.binchoice;
RealEvaluation = dat.estim;
num_trls = length(X1);
LEvaluation = NaN(num_trls,1);
X = -180:180; % the range of values over which we calculate the pdf
for i = 1:num_trls
    % generate a distribution of the correlated noise
X1_noise = psycho_noise.*randn(niter,1);
NoisyX1 = X1(i) + psycho_bias + X1_noise;
% get all the trials where the binary choice matches the simulated choice
trls2use = find(sign(NoisyX1) == RealDecision(i));
% use only these trials in the future
iter_valid = length(trls2use);
NoisyX1 = individualparams(2)*NoisyX1(trls2use);
% now generate the noise distribution for second interval
X2_noise = psycho_noise.*randn(iter_valid,1);
% combine this noise with correlated noise and evidence presented in
% interval 2
NoisyX2 = individualparams(3)*(X2(i) + psycho_bias + individualparams(4)*X1_noise(trls2use) + (1-individualparams(4))*X2_noise);
% get the final evaluation corrupted with estimation noise
Evaluation = NoisyX1 + NoisyX2 + randn(iter_valid,1).*individualparams(1);
Evaluation = [Evaluation; X'];
pd = fitdist(Evaluation, 'Kernel', 'Kernel','epanechnikov');
LEvaluation(i) = pdf(pd,RealEvaluation(i));
end
% total log likelihood, one value for each trial
PlogObservation = sum(log(LEvaluation));
% allow for nan and inf in the loglikelihood function, but set their probability to be super small
if isinf(abs(PlogObservation))
    PlogObservation = -10e100;
end
if isnan(PlogObservation)
    PlogObservation = -10e100;
end
% take the negative ll for fminsearch!
optimal_funcval = -PlogObservation;
end