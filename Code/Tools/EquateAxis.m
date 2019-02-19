%this script matches the values of the axis and then adds a black diagonal.

xl=get(gca,'xlim'); yl=get(gca,'ylim');
axlim=[min([xl(1) yl(1)]) max([xl(2) yl(2)])];
set(gca,'xlim',axlim, 'ylim', axlim);
hold on
plot(axlim, axlim,'k--','LineWidth',0.5, 'Color', [0.6 0.6 0.6]);
