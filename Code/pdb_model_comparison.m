function pdb_model_comparison(BIC)

% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code generates
% Figure 2A
global subjects;
cols = linspecer(10, 'qualitative');
colormap(linspecer);
figure;
subplot(4,4,[1,2,3]);hold on;
% draw patches to indicate regions with strong evidence in favor of one
% model or the other
x = [-225 -10 -10 -225];
y = [0 0 4.25 4.25];
patch(x,y, [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4], 'FaceAlpha', 0.2);
x = [10 625 625 10];
patch(x,y, [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4], 'FaceAlpha', 0.2);
plot([0 0], [0 4.25], 'k', 'LineWidth', 0.25);
y = 4*ones(1,length(subjects)) + 0.1*rand(1,length(subjects));
x = BIC.choice_selective - BIC.baseline;
scatter(x,y,75,cols, 'filled', 'MarkerEdgeColor', [1 1 1]);
y = 3*ones(1,length(subjects)) + 0.1*rand(1,length(subjects));
x = BIC.choice_selective - BIC.correlated_noise;
scatter(x,y,75,cols, 'filled', 'MarkerEdgeColor', [1 1 1]);
y = 2*ones(1,length(subjects)) + 0.1*rand(1,length(subjects));
x = BIC.choice_selective - BIC.shift;
scatter(x,y,75,cols, 'filled', 'MarkerEdgeColor', [1 1 1]);
y = 1*ones(1,length(subjects)) + 0.1*rand(1,length(subjects));
x = BIC.choice_selective - BIC.stimulus_selective;
scatter(x,y,75,cols, 'filled', 'MarkerEdgeColor', [1 1 1]);
set(gca, 'YLim', [1 4.25], 'XLim', [-150 650], 'YTick', 1:4, 'YTickLabel', {'Stimulus-based Selective Gain', 'Shift', 'Correlated Noise',  'Baseline'}, 'XTick', -225:75:600, 'XTickLabel', -225:75:600)
xlabel('\DeltaBIC relative to Choice-based Selective Gain model');
title({'Figure 2A', 'Model Comparison'});
offsetAxes;