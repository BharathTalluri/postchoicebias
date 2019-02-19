function optimal_funcval = model_Numerical_ChoiceSelective(individualparams)
global dat;global psycho_noise; global psycho_bias;
% get the relevant values first
X1              =  dat.x1_relative;
X2              = dat.x2_relative;
RealDecision    = dat.binchoice;
RealEvaluation  = -50 + dat.estim; % relative estimation values
LEvaluation = NaN(length(X1),1);
step=0.1;
X = -100:step:100; % the range of values over which we calculate the pdf
num_trls = length(X1);
for i = 1:num_trls
    % get the pdf for the first interval
    Y1cc = normpdf(X, (X1(i) + psycho_bias)*individualparams(2), psycho_noise*abs(individualparams(2)));
    % set the pdf of the unchosen side to zero, assuming subjects base
    % their estimations only on the chosen side
    Y1cc(sign(X) ~= sign(RealDecision(i))) = 0;
    % we need the pdf, so normalise the resulting distribution such that
    % the area under the curve is 1
    Y1cc = Y1cc./trapz(X,Y1cc);
    % get the pdf of the second interval
    Y2cc = normpdf(X, (X2(i) + psycho_bias)*individualparams(3), psycho_noise*abs(individualparams(3)));
    % convolve the two pdfs corresponding to the two inervals respectively-
    % this is similar to adding two random variables drawn from the
    % distributions
    N1 = conv(Y1cc,Y2cc,'same');
    % we need the pdf, so normalise the resulting distribution such that
    % the area under the curve is 1
    N1=N1./trapz(X,N1);
    % we will also add an additive term to account for the attraction-to-mean effect in the numbers data to a zero mean Normal distribution that represents the estimation
    % noise
    NoiseDist = normpdf(X,individualparams(4),individualparams(1));
    % add this noise to the estimations from the two intervals
    NN1 = conv(N1,NoiseDist,'same');
    % we need the pdf, so normalise the resulting distribution such that
    % the area under the curve is 1; this will be the final pdf
    % corresponding to the estimations of this trial
    NN1 = NN1./trapz(X,NN1);
    % get the likelihood of the subject's estimation from the pdf obtained
    % above
    LEvaluation(i) = interp1(X(find(X<=50 & X>=-50)),NN1(find(X<=50 & X>=-50)),RealEvaluation(i)); % in this task, subjects are restricted to the range 0-100 in their estimations
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