function new_data = filter_raw_data_Cognitive()
% Bharath Talluri & Anne Urai
% get raw data and pre-process it for all other analyses
dat = readtable('Rawdata_Cognitive.csv');
% remove trials where subjects did not adhere to instructions and where
% binary response was too fast
subjects = unique(dat.subj)';
new_data = [];
for sj = subjects
    subj_dat = dat(find(dat.subj == sj),:);
    subj_dat.x1_relative = subj_dat.x1 - 50;
    subj_dat.x2_relative = subj_dat.x2 - 50;
    % now for each subject remove estimation outliers
    estim_trials = find(subj_dat.condition ~= -1);
    % take all choice-only trials
    choiceonlytrls =  find(subj_dat.condition == -1);
    new_data = [new_data;subj_dat(choiceonlytrls,:)];
    subj_dat = subj_dat(estim_trials,:);
    % remove outliers
    estim = subj_dat.estim;
    temp1 = iqr(estim); % get the inter quartile range
    temp2 = quantile(estim,3);%get quartiles
    usetrls = find((estim < temp2(3) + 1.5*temp1) & (estim > temp2(1) - 1.5*temp1)); % remove outliers
    subj_dat = subj_dat(usetrls,:);
    new_data = [new_data;subj_dat];
end
new_data = new_data(find(new_data.condition~= 0), :);