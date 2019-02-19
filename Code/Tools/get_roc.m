function roc = get_roc(dat, consistent)
% get the roc indices for choice based selective gain mechanism
dat.idx = transpose(1:height(dat));
% since consistency cannot be defined for trials where x2 = 0, remove those
% trials
choicetrls = find(dat.x2 ~= 0);
dat = dat(choicetrls,:);
switch consistent % split by consistency effect
    case 1
        trls2use    = find(sign(dat.binchoice) == sign(dat.x2));
    case 0
        trls2use    = find(sign(dat.binchoice) ~= sign(dat.x2));
end
dat = dat(trls2use,:);
estimation = dat.estim;
% loop over different x1/x2 pairs
unique_x1 = unique(dat.x1)';
tmpRoc_all = nan(1,length(unique_x1));
num_trls_all = nan(1,length(unique_x1));
for x1 = unique_x1
    num_trls_x1 = NaN(1,2);
    tmpRoc_x1 = NaN(1,2);
    for x2s = 1:2
        % collect the estimation distributions for different values of x2
        switch x2s
            case 1
                trlsP = estimation(find(dat.x2 == 20 & dat.x1 == x1));
                trlsM = estimation(find(dat.x2 == 10 & dat.x1 == x1));
            case 2
                trlsP = estimation(find(dat.x2 == -10 & dat.x1 == x1));
                trlsM = estimation(find(dat.x2 == -20 & dat.x1 == x1));
        end
        % exclude distributions with no trials- we cannot get roc estimates.
        if isempty(trlsP) || isempty(trlsM)
            continue
        end
        % calculate the roc index for the two distributions
        tmpRoc     = rocAnalysis(trlsM, trlsP, 0, 1);
        num_trls_x1(x2s) = length(trlsP) + length(trlsM);
        tmpRoc_x1(x2s) = num_trls_x1(x2s)*(tmpRoc.i');
    end
    tmpRoc_all(find(x1 == unique_x1)) = nansum(tmpRoc_x1)/nansum(num_trls_x1);
    num_trls_all(find(x1 == unique_x1)) = nansum(num_trls_x1);
end
% weighted average roc
tmpRoc_all(isnan(tmpRoc_all)) = 0;
roc = (tmpRoc_all*num_trls_all')/sum(num_trls_all);
end