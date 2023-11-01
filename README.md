# RegressionSim
Data and scripts for analyses of OLS vs IRLS regression on real and simulated photometry data 

_IRLS_dFF.m_ computes an Iteratively Reweighted Least-Squares (IRLS) dF/F score, given an experimental signal, isosbestic signal and IRLS tuning constant. 

_dFF_simulation.m_ offers a sandbox environment to simulate fibre photometry data with a range of parameter options and test the consequences of different analysis options. This script:
  - simulates experimental and isosbestic signals according to input parameters;
  - computes multiple dF/F scores (and produces related peri-event plots) showing the effects of three key analytical decisions: low-pass filtering (low-pass vs no low-pass), regression type (OLS vs IRLS) and baseline correction (dF vs dF/F);
  - computes mean residuals and absolute residuals between normalized true and extracted signals.

Remaining scripts are functions used by the dFF_simulation script:

  - _dFFsim_plot_dFsignals.m_ computes the dF and dF/F, and plots them along with the experimental, isosbestic and fitted isosbestic signals
  - _dFFsim_plot_ERT.m_ plots the peri-event dF/F
  - _dFFsim_plot_linTrend.m_ plots the linear trend of the resulting dF/F
  - _double_exp_decay.m_ is used to apply the double exponential decay (photobleaching) component to the raw signals
