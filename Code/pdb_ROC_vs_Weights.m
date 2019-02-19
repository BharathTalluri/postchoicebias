function pdb_ROC_vs_Weights(consistency)

% Bharath Talluri & Anne Urai
% code accompanying the post-decision bias paper. This code calculates
% the correlation between roc indices and model weights, and reproduces figure 3D.

% specify the color map and figure properties
cols = linspecer(9, 'qualitative');
colormap(linspecer);
figure;
subplot(4,4,1); hold on;
numerical_samples = length(consistency.modelbased.numerical);
perceptual_samples = length(consistency.modelbased.perceptual);
dat1 = [consistency.modelfree.numerical; consistency.modelfree.perceptual];
dat2 = [consistency.modelbased.numerical; consistency.modelbased.perceptual];
dat11 = []; dat22 = [];
% polish the figure
set(gca, 'XLim', [-0.2 0.4], 'XTick', -0.2:0.2:0.4,'ylim',[-0.4 0.8], 'ytick', -0.4:0.4:0.8);
axis square;
EquateAxis;
myscatter(dat1(1:numerical_samples), dat2(1:numerical_samples), [dat11 dat22], 50, cols(2,:), cols(2,:), cols(2,:));
myscatter(dat1(numerical_samples1:end), dat2(numerical_samples + 1:end), [dat11 dat22], 50, cols(3,:), cols(3,:), cols(3,:));
[rho, pval] = corr(dat1, dat2, 'type', 'Spearman'); % permutation test
ylabel('Delta Weights');
xlabel('Delta ROC');
title({'Figure 3D',sprintf('Spearmans rho = %.4f, p = %.4f', rho, pval)});
end