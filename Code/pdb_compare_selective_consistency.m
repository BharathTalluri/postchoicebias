function pdb_compare_selective_consistency(choice_params, selective_params)
% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code reproduces
% figure 2D of the paper.

% specify the color map and figure properties
cols = linspecer(10, 'qualitative');
colormap(linspecer);
figure;
subplot(4,4,1); hold on;
% weights of second interval
dat2 = choice_params - selective_params;
dat1 = ones(1,length(dat2)) + 0.1*rand(1, length(dat2));
% polish the figure
set(gca, 'xlim', [0.75 1.25], 'xcolor', [1 1 1], 'ylim',[-0.2 0.6], 'ytick', -0.2:0.2:0.6);
plot([0.75 1.25], [0 0], 'k', 'LineWidth', 0.25);
EquateAxis;
myscatter(dat1, dat2, [], 75, cols, cols, cols);
[pval] = permtest(choice_params, selective_params, 0, 100000); % permutation test
ylabel('\DeltaW for subsequent stimulus');
offsetAxes;
title({'Figure 2D','Choice model vs. Stimulus model', sprintf('p = %.4f', pval)});