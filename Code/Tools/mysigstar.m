function h = mysigstar(xpos, ypos, pval, ns, color, whichWay)
% replaces sigstar, which doesnt work anymore in matlab 2014b

if ~exist('ns', 'var'); ns = '1'; end
if ~exist('color', 'var'); color = 'k'; end
if ~exist('whichWay', 'var'), whichWay = 'none'; end

if numel(ypos) > 1,
    assert(ypos(1) == ypos(2), 'line wont be straight!');
    ypos = ypos(1);
end

% draw line
hold on;
if numel(xpos) > 1,
    if ns
        % plot the horizontal line
        p = plot([xpos(1), xpos(2)], ...
            [ypos ypos], '-', 'LineWidth', 0.5, 'color', color);
    else
        if pval < 0.05
            p = plot([xpos(1), xpos(2)], ...
            [ypos ypos], '-', 'LineWidth', 0.5, 'color', color);
        end
    end
    
        % also add small ticks
        switch whichWay
            case 'none'
                fprintf('Not adding small vertical ticks \n');
            case 'down'
                plot([xpos(1) xpos(1)], [ypos ypos-0.05*range(get(gca, 'ylim'))], '-k', 'LineWidth', 0.5, 'color', color);
                plot([xpos(2) xpos(2)], [ypos ypos-0.05*range(get(gca, 'ylim'))], '-k', 'LineWidth', 0.5, 'color', color);
            case 'up'
                plot([xpos(1) xpos(1)], [ypos ypos+0.05*range(get(gca, 'ylim'))], '-k', 'LineWidth', 0.5, 'color', color);
                plot([xpos(2) xpos(2)], [ypos ypos+0.05*range(get(gca, 'ylim'))], '-k', 'LineWidth', 0.5, 'color', color);
        end
    
    % use white background
    txtBg = 'white';
else
    txtBg = 'none';
end

fz = 10;
fw = 'bold';
if pval < 1e-3
    txt = '***';
elseif pval < 1e-2
    txt = '**';
elseif pval < 0.05
    txt = '*';
else
    txt = 'n.s.';
    fz = 9;
    fw = 'normal';
end

% draw the stars in the bar

if ns
    h = text(mean(xpos), mean(ypos), txt, ...
        'horizontalalignment', 'center', 'backgroundcolor', txtBg, 'margin', 1, 'fontsize', fz, 'FontWeight', fw, 'color', color);
else
    if pval < 0.05
        h = text(mean(xpos), mean(ypos), txt, ...
            'horizontalalignment', 'center', 'backgroundcolor', txtBg, 'margin', 1, 'fontsize', fz, 'FontWeight', fw, 'color', color);
    end
end
end
