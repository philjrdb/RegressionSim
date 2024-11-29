% Fiber photometry data processing from Martianova et al. (2019)
%
% This code has been adapted for use in Keevers et al. (2024)
% If you want to use this code, please cite the original Jove paper: 
%   Martianova, E., Aronson, S., Proulx, C.D. Multi-Fiber Photometry to  
%   Record Neural Activity in Freely Moving Animal. J. Vis. Exp. 
%   (152), e60278, doi:10.3791/60278 (2019).

% Run section by section (Ctrl+Enter)

%% Load simulated data before running!

n_sims = 10;

for sim = 1:n_sims
    % Change the next two lines depending on your data frame
    martianova.(sim_name).simulation(sim).raw_reference = dFFsim.(sim_name).simulation(sim).iso_signal.'; 
    martianova.(sim_name).simulation(sim).raw_signal = dFFsim.(sim_name).simulation(sim).exp_signal.';
    
    % Plot raw data
    figure
    subplot(2,1,1)
    plot(martianova.(sim_name).simulation(sim).raw_reference,'m')
    subplot(2,1,2)
    plot(martianova.(sim_name).simulation(sim).raw_signal,'b')
    close all
    
    %% Use function get_zdFF.m to analyze data
    martianova.(sim_name).simulation(sim).zdFF = get_zdFF(martianova.(sim_name).simulation(sim).raw_reference,martianova.(sim_name).simulation(sim).raw_signal);
    
    % Plot z-score dF/F
    figure
    plot(martianova.(sim_name).simulation(sim).zdFF,'k')
    close all
    
    %% Analysis step by step
    % Smooth data
    martianova.parameters.smooth_win = 3;
    martianova.(sim_name).simulation(sim).smooth_reference = movmean(martianova.(sim_name).simulation(sim).raw_reference,martianova.parameters.smooth_win);
    martianova.(sim_name).simulation(sim).smooth_signal = movmean(martianova.(sim_name).simulation(sim).raw_signal,martianova.parameters.smooth_win);
    
    % Plot smoothed signals
    figure
    subplot(2,1,1)
    plot(martianova.(sim_name).simulation(sim).smooth_reference,'m')
    subplot(2,1,2)
    plot(martianova.(sim_name).simulation(sim).smooth_signal,'b')
    
    %% Remove slope using airPLS algorithm (airPLS.m)
    martianova.parameters.lambda = 5e9;
    martianova.parameters.order = 2;
    martianova.parameters.wep = 0.7;
    martianova.parameters.p = 0.5;
    martianova.parameters.itermax = 50;
    [martianova.(sim_name).simulation(sim).reference,martianova.(sim_name).simulation(sim).base_r]= airPLS(martianova.(sim_name).simulation(sim).smooth_reference,martianova.parameters.lambda,martianova.parameters.order,martianova.parameters.wep,martianova.parameters.p,martianova.parameters.itermax);
    [martianova.(sim_name).simulation(sim).signal,martianova.(sim_name).simulation(sim).base_s]= airPLS(martianova.(sim_name).simulation(sim).smooth_signal,martianova.parameters.lambda,martianova.parameters.order,martianova.parameters.wep,martianova.parameters.p,martianova.parameters.itermax);
    
    % Plot slopes
    figure
    subplot(2,1,1)
    plot(martianova.(sim_name).simulation(sim).smooth_reference,'m')
    hold on
    plot(martianova.(sim_name).simulation(sim).base_r,'k')
    hold off
    subplot(2,1,2)
    plot(martianova.(sim_name).simulation(sim).smooth_signal,'b')
    hold on
    plot(martianova.(sim_name).simulation(sim).base_s,'k')
    hold off
    
    %% Standardize signals
    martianova.(sim_name).simulation(sim).z_reference = (martianova.(sim_name).simulation(sim).reference - median(martianova.(sim_name).simulation(sim).reference)) / std(martianova.(sim_name).simulation(sim).reference);
    martianova.(sim_name).simulation(sim).z_signal = (martianova.(sim_name).simulation(sim).signal - median(martianova.(sim_name).simulation(sim).signal)) / std(martianova.(sim_name).simulation(sim).signal);
    
    % Plot signals
    figure
    subplot(2,1,1)
    plot(martianova.(sim_name).simulation(sim).z_reference,'m')
    subplot(2,1,2)
    plot(martianova.(sim_name).simulation(sim).z_signal,'b')
    
    %% Fit reference signal to calcium signal 
    % using non negative robust linear regression
    martianova.(sim_name).simulation(sim).fitdata = fit(martianova.(sim_name).simulation(sim).z_reference',martianova.(sim_name).simulation(sim).z_signal',fittype('poly1'),'Robust','on');
        
    % Plot fit
    figure
    hold on
    plot(martianova.(sim_name).simulation(sim).z_reference,martianova.(sim_name).simulation(sim).z_signal,'k.')
    plot(martianova.(sim_name).simulation(sim).fitdata,'b')
    hold off
    
    %% Align reference to signal
    martianova.(sim_name).simulation(sim).z_reference = martianova.(sim_name).simulation(sim).fitdata(martianova.(sim_name).simulation(sim).z_reference)';
   
    % Plot aligned signals
    figure
    plot(martianova.(sim_name).simulation(sim).reference,'m')
    hold on
    plot(martianova.(sim_name).simulation(sim).signal,'b')
    hold off
    
    %% Calculate z-score dF/F
    martianova.(sim_name).simulation(sim).zdFF = (martianova.(sim_name).simulation(sim).z_signal - martianova.(sim_name).simulation(sim).z_reference).';
    martianova.(sim_name).simulation(sim).nullZdFF = nullZ(martianova.(sim_name).simulation(sim).zdFF);
    
    % Plot z-score dF/F
    figure
    plot(martianova.(sim_name).simulation(sim).nullZdFF,'k')
    close all
end

%% Add results to main dFF_simulation
for sim = 1:n_sims
dFFsim.(sim_name).results(sim).nullZ_martianova_dFF = martianova.(sim_name).simulation(sim).nullZdFF;
end

%% Done
clearvars -except dFFsim savefolder sim_name n_sims martianova
toc;
