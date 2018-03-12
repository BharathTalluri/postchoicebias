# Post-choice Bias


This folder contains the code to run all model-free and model-based analyses. 

Note:
1. Our model-based analyses use likelihood optimisation to find the appropriate parameters. Even though this method generally gives a unique solution which is ideally the global maximum, we recommend using multiple starting points (20-25 would be sufficient) to make sure the optimisation algorithm finds the global maximum.

2. We recommend parallelising the fitting procedure - by fitting each subject separately. This reduces the computational time. We fit the models on a supercomputing cluster by submitting multiple jobs, one per subject and starting point. Each job took around 4 hrs to compute the best fitting parameters.