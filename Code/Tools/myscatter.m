function myscatter(x, y, errorbars, markersize, markercolor, xeColor, yeColor)

% Bharath Talluri
if numel(xeColor) == 3
    xeColor = repmat(xeColor, length(x),1);
end
if size(yeColor) == 3
    yeColor = repmat(yeColor, length(x),1);
end
if ~isempty(errorbars)
    xe1 = errorbars(:,1);xe2 = errorbars(:,2);ye1 = errorbars(:,3);ye2 = errorbars(:,4);
    for i = 1:length(x)
        plot([xe1(i) xe2(i)], [y(i) y(i)], 'Color', xeColor(i,:),'LineWidth',0.25);
        plot([x(i) x(i)], [ye1(i)  ye2(i)], 'Color', yeColor(i,:),'LineWidth',0.25);
    end
end
if numel(markercolor) == 3
    scatter(x, y, markersize, 'filled', 'MarkerEdgeColor',[1 1 1], 'MarkerFaceColor', markercolor);hold on;
else
    scatter(x, y, markersize, markercolor, 'filled', 'MarkerEdgeColor',[1 1 1]);hold on;
end

