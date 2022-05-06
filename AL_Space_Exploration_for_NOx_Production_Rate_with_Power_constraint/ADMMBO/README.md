# ADMMBO

This package is an implementation of "ADMMBO, An ADMM Framework for Bayesian Optimization with Unknown Constraints'', to appear in the Journal of Machine Learning Research (JMLR), special issue on Bayesian Optimization.
See the following link the for original codes: https://github.com/SetarehAr/ADMMBO

This package uses the Gaussian process library gpml and Bayesian optimization library bayesopt. 

This package was modified by Ketong shao (May 2022), with the following modifications:
1. Reusing the history of data points. The original code will not record the evaluations and will recalculate the evaluations using the history inputs.
2. Adapt to the situations when objective and constraints can be co-evaluated.
