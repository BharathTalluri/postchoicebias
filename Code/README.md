# Post-choice Bias


This folder contains the code to run all model-free and model-based analyses.

Note:
1. Our model-based analyses use likelihood optimisation to find the appropriate parameters. Even though this method generally gives a unique solution which is ideally the global maximum, we recommend using multiple starting points (20-25 would be sufficient) to make sure the optimisation algorithm finds the global maximum.

2. We recommend parallelising the fitting procedure - by fitting each subject, and each iteration separately. This reduces the computational time. We fit the models on a supercomputing cluster by submitting multiple jobs, one per subject and starting point. Each job took around 4 hrs to compute the best fitting parameters.

3. Most of the tools required are included in the tools subfolder. For optimisation, we use the Subplex algorithm, a generalization of the Nelder-Mead simplex method, which is well suited to optimize high dimensional noisy objective functions. The code for this algorithm is available online (https://link.springer.com/article/10.3758%2FBF03206554). please cite the corresponding paper if you use this.
For further optimisation, we use fminsearchbnd, a matlab bound constrained optimisation algorithm available here (https://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon).
