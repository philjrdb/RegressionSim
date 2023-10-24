function [dFF, ft_iso_signal] = IRLS_dFF(exp_signal, iso_signal, IRLS_constant)

 % IRLS_dFF - Compute Iteratively Reweighted Least Squares (IRLS) dFF
    %   [dFF] = IRLS_dFF(exp_signal, iso_signal, IRLS_constant) fits the
    %   isosbestic signal to the experimental signal using IRLS regression
    %   with a specified tuning constant and computes dFF using the following equation:
    %   dFF = (experimental signal - fitted isosbestic)/fitted isosbestic
    %
    %   Inputs:
    %     exp_signal - The experimental signal (dynamic + artifacts)
    %     iso_signal - The isosbestic signal (artifacts)
    %     IRLS_constant - A constant parameter for the robust regression. 
    %                     Smaller values equate to more aggressive downweighting of outliers.
    %                     Suggested default of 1.4
    %
    %   Outputs:
    %     dFF - The artifact-corrected dF/F score
    %     ft_iso_signal - The IRLS-fitted isosbestic signal
    %
    %   Example:
    %     dFF = IRLS_dFF(465_signal, 405_signal, 1.4);
    %
    %   See also: robustfit

% Input Validation - check is signals are same length
    
    if numel(exp_signal) ~= numel(iso_signal)
        error('exp_signal and iso_signal must have the same length.');
    end

% Compute IRLS-fitted isosbestic and dF/F
IRLS_coeffs = reshape(flipud(robustfit(iso_signal, exp_signal, 'bisquare', IRLS_constant, 'on')), [1, 2]);
              % reshape and flipud used to arrange IRLS coefficients in same manner as polyval function
ft_iso_signal = polyval(IRLS_coeffs,iso_signal);
dFF = (exp_signal-ft_iso_signal)./ft_iso_signal;

end