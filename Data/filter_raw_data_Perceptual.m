function new_data = filter_raw_data_Perceptual()
% Bharath Talluri & Anne Urai
% get raw data and pre-process it for all other analyses
dat = readtable('Rawdata_Perceptual.csv');
% remove trials where subjects did not adhere to instructions and where
% binary response was too fast
trls2use = find(sign(abs(dat.condition)) == sign(abs(dat.binchoice)) & dat.binrt > 0.2 & ~isnan(dat.binchoice));
dat = dat(trls2use,:);
subjects = unique(dat.subj)';
new_data = [];
for sj = subjects
    subj_dat = dat(find(dat.subj == sj),:);
    % now for each subject remove estimation outliers
    estim_trials = find(subj_dat.condition ~= -1);
    % take all choice-only trials
    choiceonlytrls =  find(subj_dat.condition == -1);
    new_data = [new_data;subj_dat(choiceonlytrls,:)];
    subj_dat = subj_dat(estim_trials,:);
    % check if there are any invalid estimations in estim trials
    if sum(isnan(subj_dat.estim))>0
        keyboard
    end
    % remove outliers
    for i = unique(subj_dat.xavg')
        gettrls = find(subj_dat.xavg == i);
        dat3 = subj_dat(gettrls,:);
        estim = dat3.estim;
        
        temp1 = iqr(estim); % get the inter quartile range
        temp2 = quantile(estim,3);%get quartiles
        usetrls = find((estim < temp2(3) + 1.5*temp1) & (estim > temp2(1) - 1.5*temp1)); % remove outliers
        new_data = [new_data;dat3(usetrls,:)];
    end
end
new_data = new_data(find(new_data.condition~= 0), :);