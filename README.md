# RegressionSim
Data and scripts for analyses of OLS vs IRLS regression on real and simulated photometry data 

_dFF_simulation_data_ contains the simulation data displayed in Figure 3 of our article, compressed into a zip file.

_IRLS_dFF.m_ computes an Iteratively Reweighted Least-Squares (IRLS) dF/F score, given an experimental signal, isosbestic signal and IRLS tuning constant. 

_dFF_simulation.m_ offers a sandbox environment to simulate fibre photometry data with various parameter options and test the consequences of different analysis options. This script:
  - simulates experimental and isosbestic signals according to input parameters;
  - computes multiple dF/F scores (and produces related peri-event plots) showing the effects of three key analytical decisions: low-pass filtering (low-pass vs no low-pass), regression type (OLS vs IRLS) and baseline correction (dF vs dF/F);
  - computes mean residuals and absolute residuals between normalized true and extracted signals;
  - computes tCI with consecutive significance threshold for extracted ERTs and identifies inferential errors compared to the true signal.

The remaining scripts are functions used by the dFF_simulation script:

  - _consec_idx.m_ is used to apply a consecutive threshold to the tCI calculated for each average ERT
  - _dFFsim_plot_dFsignals.m_ computes the dF and dF/F, and plots them along with the experimental, isosbestic, and fitted isosbestic signals
  - _dFFsim_plot_ERT.m_ plots the peri-event dF/F
  - _dFFsim_plot_linTrend.m_ plots the linear trend of the resulting dF/F
  - _double_exp_decay.m_ is used to apply the double exponential decay (photobleaching) component to the raw signals
  - _nullZ.m_ is used to z-score the data around zero

---------------------------
EDITS/ADDITIONS [2024-11-14] - Updated _dFF_simulation.m_ to include calculation and plotting of tCIs, and identification of inferential errors for ERTs of each simulated signal. Added function _consec_idx.m_ which is used in dFF simulation to apply consecutive significance threshold to the tCI.
