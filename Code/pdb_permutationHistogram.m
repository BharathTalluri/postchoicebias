function pdb_permutationHistogram()
% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code reproduces
% figure S1 of the paper.
close all;
clc; dbstop if error;
% specify the figure properties
myfigureprops;
figure;
% specify the path to the data
behdatapath = '../Data';
behdata = readtable(sprintf('%s/Task_Perceptual.csv', behdatapath));
subjects = unique(behdata.subj)';
for sj = subjects
    dat = behdata(find(behdata.subj == sj),:);
    % define the type of trials
    trls2use = find(dat.condition == 1 & abs(dat.binchoice) == 1);
    dat = dat(trls2use,:);
    permutationHistogram = nan(5,5);
    dirs = -20:10:20;
    for d = 1:length(dirs)
        x1trls = find(dat.x1 == dirs(d));
        for d2 = 1:length(dirs)
            permutationHistogram(d, d2) = length(find(dat.x2(x1trls) == dirs(d2)));
        end
    end
    % plot grid
    subplot(4,4,find(sj==subjects));
    colormap viridis;
    imagesc(dirs, dirs, permutationHistogram, [0 45]);
    axis square;
    set(gca, 'ydir', 'normal');
    title(sprintf('S%02d', sj));
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
    lbl = strtrim(cellstr(num2str((permutationHistogram(:)')')));
    text(xlbl(:), ylbl(:), lbl(:),'color', 'w',...
        'HorizontalAlignment','center','VerticalAlignment','middle', 'fontsize', 7);
end
% make the colorbar prettier, move to the side
c       = colorbar;
ax      = gca;
axpos   = ax.Position;
cpos    = c.Position;
cpos(3) = 0.3*cpos(3); % thinner
cpos(1) = cpos(1)*1.2; % to the rigth
cpos(4) = cpos(4) *0.9; % shorter
c.Position = cpos;
ax.Position = axpos;
% put a string on top of the colorbar
c.Label.String = 'Number of trials';
c.Label.Rotation = 270;
c.Label.HorizontalAlignment = 'center';
c.Box = 'off';
c.TickDirection = 'out';
clpos = c.Label.Position; clpos(1) = clpos(1)*4;
c.Label.Position = clpos;
% labels
suplabel('Interval 1 direction', 'x');
suplabel('Interval 2 direction', 'y');

