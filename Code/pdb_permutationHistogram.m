function pdb_permutationHistogram(behdata)
% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code reproduces
% figure S1 C, D of the paper.
global subjects;
figure;
ga_permutation_histogram = nan(numel(unique(subjects)), 5, 5);
for types = 1:3
    for sj = subjects
        dat = behdata(find(behdata.subj == sj),:);
        % define the type of trials
        switch types
            case 1
                trls2use = find(dat.condition == 1 & abs(dat.binchoice) == 1);
            case 2
                trls2use = find(dat.condition == 1 & abs(dat.binchoice) == 1 & sign(dat.binchoice) == sign(dat.x2) & dat.x2 ~= 0);
            case 3
                trls2use = find(dat.condition == 1 & abs(dat.binchoice) == 1 & sign(dat.binchoice) ~= sign(dat.x2) & dat.x2 ~= 0);
        end
        dat = dat(trls2use,:);
        permutationHistogram = nan(5,5);
        dirs = -20:10:20;
        for d = 1:length(dirs)
            x1trls = find(dat.x1 == dirs(d));
            for d2 = 1:length(dirs)
                permutationHistogram(d, d2) = length(find(dat.x2(x1trls) == dirs(d2)));
            end
        end
        ga_permutation_histogram(find(sj == subjects), :, :) = permutationHistogram;
    end
    
    %% PLOT THE HISTOGRAM SHOWING THE AVERAGE TRIAL COUNTS
    subplot(4,4,types);
    
    avg_permutation_histogram = squeeze(round(nanmean(ga_permutation_histogram)));
    colormap viridis;
    switch types
        case 1
            imagesc(dirs, dirs, avg_permutation_histogram, [38 44]);
        case 2
            imagesc(dirs, dirs, avg_permutation_histogram, [5 40]);
        case 3
            imagesc(dirs, dirs, avg_permutation_histogram, [5 40]);
    end
    
    axis square;
    set(gca, 'ydir', 'normal');
    switch types
        case 1
            title('Choice trials');
        case 2
            title('Consistent trials');
        case 3
            title('Inconsistent trials');
    end
    set(gca, 'xtick', -20:10:20, 'ytick', -20:10:20);
    hold on;
    
    plotGrid = -25:10:25;
    % show grid lines
    for k = plotGrid
        % horizontal lines
        x = [plotGrid(1) plotGrid(end)];
        y = [k k];
        plot(x,y,'Color','w','LineStyle','-');
        % vertical lines
        x = [k k];
        y = [plotGrid(1) plotGrid(end)];
        plot(x,y,'Color','w','LineStyle','-');
    end
    % coordinates
    [xlbl, ylbl] = meshgrid(dirs, dirs);
    lbl = strtrim(cellstr(num2str((avg_permutation_histogram(:)')')));
    text(xlbl(:), ylbl(:), lbl(:),'color', 'w',...
        'HorizontalAlignment','center','VerticalAlignment','middle', 'fontsize', 7);
    
    % labels
    xlabel(sprintf('Interval 2 direction (%c)', char(176)));
    ylabel(sprintf('Interval 1 direction (%c)', char(176)));
end