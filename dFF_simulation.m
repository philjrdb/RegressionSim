%% Parameters (PLEASE FILL OUT)
savefolder = 'G:\Current\Fibre Photometry\dFF_sim'; % Path to folder where the results to be saved
sim_name = 'final_test'; % Name of simulation for multiple iterations with different paramaters
                           % (change between iterations to avoid overriding previous plots)
n_sims = 1; % # of simulations

% Signal parameters
samp_rate = 10; % sampling rate (Hz)
n_secs = 1200; % session length (secs)
Lp_cutoff = 3; % low-pass filter Hz for experimental & isosbestic signals

noise_factor = [2]; % multiplier on gaussian noise
movement_attenuation = 50; % max % reduction in signal due to movement 

exp_base = [200]; % experimental signal baseline (starting/background fluorescence) [200 200]
exp_decay_base = [40]; % decay % from exp base [40 40]
exp_decay_rate1 = [.02]; % decay rate1 applied in double_exp_decay.m (rate1&2 of .02 & .002 works well) [.02 .02]
exp_decay_rate2 = [.002]; % decay rate2 (for double exponent decay) [.002 .002]

iso_base = [80]; % isosbestic signal baseline (starting/background fluorescence) [80 80]
iso_decay_base = [40];  % decay % from iso base [40 40]
iso_decay_rate1 = [.02]; % decay rate [.02 .05]
iso_decay_rate2 = [.002]; % decay rate [.002 .002]

IRLS_constant = [1.4 3 4.685]; % Tuning constants for iteratively-reweighted-least-squares regression
                               % Smaller values = more aggressive downweighting of outliers 
                               % Results saved as "IRLS1", "IRLS2"... "IRLS(n)" for n constants

% Neural component parameters
n_events = 100;
ERT = [0.1	0.6	0.95 1 0.8	0.5	0.3	0.2	0.1	0.05]; % event-related transient waveform 
ERT_times = linspace(10, 1100,n_events); % rand(n_events,1)*(pre_post_idx+n_secs-((ERT_length-pre_post_idx)/samp_rate));
                                         % or linspace(first time, last time,n_events)
    
ERT_peak_multiplier = 20; % multiplier on ERT height (not length)
ERT_coeff = ones(n_events,1)*3; % ERT scale multiplier (height & length [length never < base ERT])
                                % [] = randomised (norm positive distribution) 
                                % ones(n_events,1)*#, or linspace(start coeff,end coeff,n_events)
ERT_scale_variability = 3; % multiplier on ERT coeff if randomised ERT used
pre_post_idx = 5; % No. of data points plotted on peri-event graphs

% Plotting parameters
ERT_ylims = [];
ERT_Z_ylims = [-2 6];
lintrend_ylims = [-5 5];
lintrend_Z_ylims = [-.5 .5];

%% Housekeeping
tic;
clc;
close all;

% Save scripts/functions to Scripts folder
warning off;

%% Create folders
  warning off;
  mkdir(savefolder, sim_name);
  mkdir([savefolder '\' sim_name],'Scripts');
  mkdir([savefolder '\' sim_name],'Signals');
  mkdir([savefolder '\' sim_name],'ERT');
  mkdir([savefolder '\' sim_name],'Residual');
  diary([savefolder '\' sim_name '\Scripts\' sim_name ' simulation (' date ').txt']);
  copy_file = {'dFF_simulation.m','dFFsim_plot_dFsignals.m','dFFsim_plot_ERT.m','dFFsim_plot_linTrend.m','double_exp_decay.m'};
  for c = 1:length(copy_file)
     copyfile(copy_file{c},[savefolder '\' sim_name '\Scripts'])
  end
  warning on;
  
  fprintf([' --- ' sim_name ' ---\n']);
  fprintf(['Noise factor: ' num2str(noise_factor) '\n']);
  fprintf(['Experimental signal Pre-Post decay [rates]: ' num2str(exp_base) '-' num2str(exp_decay_base)...
     ' [' num2str(exp_decay_rate1) ',' num2str(exp_decay_rate2) ']\n']); 
  fprintf(['Isosbestic signal Pre-Post decay [rates]: ' num2str(iso_base) '-' num2str(iso_decay_base) ...
    ' [' num2str(iso_decay_rate1) ',' num2str(iso_decay_rate2) ']\n']);
  
  % Working parameters
  n_dtpts = n_secs * samp_rate;
  time = linspace(0,n_dtpts/samp_rate,n_dtpts)';
  ERT_idx = NaN(n_events,2);
  ERT_length = length(ERT);
  if isempty(ERT_coeff)
      ERT_coeff = abs(randn(n_events,1))*ERT_scale_variability;
  end
  
  %% Fill/pre-allocate data structure
  dFFsim.(sim_name).paramaters.samp_rate = samp_rate;
  dFFsim.(sim_name).paramaters.n_secs = n_secs;
  dFFsim.(sim_name).paramaters.Lp_cutoff = Lp_cutoff;
  dFFsim.(sim_name).paramaters.noise_factor = noise_factor;
  dFFsim.(sim_name).paramaters.movement_attenuation = movement_attenuation; 
  dFFsim.(sim_name).paramaters.exp_base = exp_base;
  dFFsim.(sim_name).paramaters.exp_decay_base = exp_decay_base;
  dFFsim.(sim_name).paramaters.exp_decay_rate1 = exp_decay_rate1;
  dFFsim.(sim_name).paramaters.exp_decay_rate2 = exp_decay_rate2;
  dFFsim.(sim_name).paramaters.iso_base = iso_base;
  dFFsim.(sim_name).paramaters.iso_decay_base = iso_decay_base; 
  dFFsim.(sim_name).paramaters.iso_decay_rate1 = iso_decay_rate1;
  dFFsim.(sim_name).paramaters.iso_decay_rate2 = iso_decay_rate2;
  for c = 1:length(IRLS_constant)
    dFFsim.(sim_name).paramaters.(['IRLS_constant_' num2str(c)]) = IRLS_constant(c);
  end
  dFFsim.(sim_name).paramaters.ERT_magnitude = ERT_peak_multiplier; 
  dFFsim.(sim_name).paramaters.ERT_scale_variability = ERT_scale_variability; 
  dFFsim.(sim_name).paramaters.n_events = n_events;
  dFFsim.(sim_name).paramaters.ERT = ERT;
  dFFsim.(sim_name).paramaters.ERT_length = ERT_length;
  dFFsim.(sim_name).paramaters.ERT_times = ERT_times;
  dFFsim.(sim_name).paramaters.ERT_coeff = ERT_coeff;
  dFFsim.(sim_name).paramaters.pre_post_idx = pre_post_idx;
  
  dFFsim.(sim_name).simulation = struct('ERT_idx',cell(1, n_sims), 'neural_component',cell(1, n_sims),...
    'exp_decay_component',cell(1, n_sims), 'iso_decay_component',cell(1, n_sims), ...
    'exp_noise_component',cell(1, n_sims), 'iso_noise_component',cell(1, n_sims), ...
    'exp_signal',cell(1, n_sims), 'iso_signal',cell(1, n_sims), ...
    'Lp_exp_signal',cell(1, n_sims), 'Lp_iso_signal',cell(1, n_sims));
  dFFsim.(sim_name).results = struct('OLS_ft405',cell(1, n_sims), ...
    'OLS_dF',cell(1, n_sims), 'OLS_dFF',cell(1, n_sims), ...
    'Lp_OLS_ft405',cell(1, n_sims), 'Lp_OLS_dF',cell(1, n_sims), 'Lp_OLS_dFF',cell(1, n_sims), ...
    'true_ev',cell(1, n_sims), 'true_bsl',cell(1, n_sims), ...
    'OLS_dF_ev',cell(1, n_sims), 'OLS_dF_bsl',cell(1, n_sims), ...
    'OLS_dFF_ev',cell(1, n_sims), 'OLS_dFF_bsl',cell(1, n_sims), 'Lp_OLS_dF_ev',cell(1, n_sims), ...
    'Lp_OLS_dF_bsl',cell(1, n_sims), 'Lp_OLS_dFF_ev',cell(1, n_sims), 'Lp_OLS_dFF_bsl',cell(1, n_sims));
  
  %% Run simulations
  for sim = 1:n_sims
    fprintf(['Simulation #' int2str(sim) '... ']);
    %% Generate dataset 
    % Event component
    neural_component = zeros(n_dtpts,1);
    for e = 1:length(ERT_times)
       [~,ERT_idx(e,1)] = min(abs(ERT_times(e)-time));
  
       if ERT_coeff(e) > 1
           scaled_idx = round(ERT_length*ERT_coeff(e));
           scaled_ERT = interp1(linspace(1,scaled_idx,ERT_length),...
               ERT*ERT_coeff(e)*ERT_peak_multiplier,1:scaled_idx);
           ERT_idx(e,2) = ERT_idx(e,1)+scaled_idx-1; 
           neural_component(ERT_idx(e,1):ERT_idx(e,2)) = ...
               neural_component(ERT_idx(e,1):ERT_idx(e,2)) + scaled_ERT';
       else
           scaled_ERT = ERT*ERT_coeff(e)*ERT_peak_multiplier;
           ERT_idx(e,2) = ERT_idx(e,1)+ERT_length-1;
           neural_component(ERT_idx(e,1):ERT_idx(e,2)) = ...
               neural_component(ERT_idx(e,1):ERT_idx(e,2)) + scaled_ERT';
       end
    end
  
    movement_component = 1-(lowpass(rand(n_dtpts,1),.1,samp_rate)...
        *(movement_attenuation/100));
  
    exp_decay_component = (double_exp_decay(exp_decay_rate1,exp_decay_rate2,time))...
        *(exp_decay_base/100)+(1-exp_decay_base/100); 
    iso_decay_component = (double_exp_decay(iso_decay_rate1,iso_decay_rate2,time))...
      *(iso_decay_base/100)+(1-iso_decay_base/100);
  
    exp_noise_component = randn(n_dtpts,1)*noise_factor;
    iso_noise_component = randn(n_dtpts,1)*noise_factor;
  
    exp_signal = ((neural_component + exp_base).*exp_decay_component)...
        .*movement_component + exp_noise_component;
    iso_signal = ((zeros(n_dtpts,1) + iso_base).*iso_decay_component)...
        .*movement_component + iso_noise_component;
  
    Lp_exp_signal = lowpass(exp_signal,Lp_cutoff,samp_rate,ImpulseResponse="iir",Steepness=0.95);
    Lp_exp_signal = Lp_exp_signal+(mean(exp_signal)-mean(Lp_exp_signal)); % correct for y-axis shift
    Lp_iso_signal = lowpass(iso_signal,Lp_cutoff,samp_rate,ImpulseResponse="iir",Steepness=0.95);
    Lp_iso_signal = Lp_iso_signal+(mean(iso_signal)-mean(Lp_iso_signal)); % correct for y-axis shift

    %% Plot signals
    % Experimental signals
    figure;
    subplot(5,1,1); hold on
    plot(time,neural_component);
    plot(xlim,[0 0],'k:');
    title('Experimental event component');
  
    subplot(5,1,2); hold on
    plot(time,exp_base*exp_decay_component);
    title('Baseline w/ bleaching component');
  
    subplot(5,1,3); hold on
    plot(time,movement_component);
    title('Movement component');
  
    subplot(5,1,4); hold on
    plot(time,exp_noise_component);
    plot(xlim,[0 0],'k:');
    title('Noise component');
  
    subplot(5,1,5); hold on
    plot(time,exp_signal,'g-');
    title('Experimental signal');
    
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[savefolder '\' sim_name '\Signals\Exp signals (sim' int2str(sim) ').png']);
    close all
              
    % Isosbestic signals
    figure;
    subplot(5,1,1); hold on
    plot(time,zeros(n_dtpts,1));
    plot(xlim,[0 0],'k:');
    title('Isosbestic event component');
  
    subplot(5,1,2); hold on
    plot(time,iso_base*iso_decay_component);
    title('Baseline w/ bleaching component');
  
    subplot(5,1,3); hold on
    plot(time,movement_component);
    title('Movement component');
  
    subplot(5,1,4); hold on
    plot(time,iso_noise_component);
    plot(xlim,[0 0],'k:');
    title('Noise component');
  
    subplot(5,1,5); hold on
    plot(time,iso_signal,'m-');
    title('Isosbestic signal');
    
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[savefolder '\' sim_name '\Signals\Iso signals (sim' int2str(sim) ').png']);
    close all
  
    %% Test regressions
    % Ordinary least squares (OLS)
    OLS_coeffs = polyfit(iso_signal,exp_signal,1);
    OLS_ft405 = polyval(OLS_coeffs,iso_signal);
    [~, OLS_dF, OLS_dFF] = ...
      dFFsim_plot_dFsignals(exp_signal,iso_signal,OLS_ft405,neural_component,time,'OLS');
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[savefolder '\' sim_name '\Signals\Regression - OLS (sim' int2str(sim) ').png']);
    close all
    
    % Low passed OLS
    Lp_OLS_coeffs = polyfit(Lp_iso_signal,Lp_exp_signal,1);
    Lp_OLS_ft405 = polyval(Lp_OLS_coeffs,Lp_iso_signal);
    [~, Lp_OLS_dF, Lp_OLS_dFF] = ...
      dFFsim_plot_dFsignals(Lp_exp_signal,Lp_iso_signal,Lp_OLS_ft405,neural_component,time,'Lp OLS');
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[savefolder '\' sim_name '\Signals\Regression - Lp OLS (sim' int2str(sim) ').png']);
    close all
    
    % IRLS regression [c = IRLS_constant]
    IRLS_ft405 = cell(1, length(IRLS_constant));
    IRLS_dF = cell(1, length(IRLS_constant));
    IRLS_dFF = cell(1, length(IRLS_constant));

    Lp_IRLS_ft405 = cell(1, length(IRLS_constant));
    Lp_IRLS_dF = cell(1, length(IRLS_constant));
    Lp_IRLS_dFF = cell(1, length(IRLS_constant));

    for c = 1:length(IRLS_constant)
       IRLS_coeffs = reshape(flipud(robustfit(iso_signal, exp_signal, 'bisquare', IRLS_constant(c), 'on')), [1, 2]);
       IRLS_ft405{c} = polyval(IRLS_coeffs,iso_signal);
       [~, IRLS_dF{c}, IRLS_dFF{c}] = ...
         dFFsim_plot_dFsignals(exp_signal,iso_signal,IRLS_ft405{c},neural_component,time,['IRLS (c =' num2str(IRLS_constant(c)) ')']);
       set(gcf,'Position',get(0,'Screensize'));
       saveas(figure(1),[savefolder '\' sim_name '\Signals\Regression - IRLS (c = ' num2str(IRLS_constant(c)) ') (sim' int2str(sim) ').png']);
       close all
   
       Lp_IRLS_coeffs = reshape(flipud(robustfit(Lp_iso_signal,Lp_exp_signal,'bisquare', IRLS_constant(c), 'on')), [1, 2]);
       Lp_IRLS_ft405{c} = polyval(Lp_IRLS_coeffs,Lp_iso_signal);
       [~, Lp_IRLS_dF{c}, Lp_IRLS_dFF{c}] = ...
         dFFsim_plot_dFsignals(Lp_exp_signal,Lp_iso_signal,Lp_IRLS_ft405{c},neural_component,time,['Lp IRLS (c =' num2str(IRLS_constant(c)) ')']);
       set(gcf,'Position',get(0,'Screensize'));
       saveas(figure(1),[savefolder '\' sim_name '\Signals\Regression - Lp IRLS (c = ' num2str(IRLS_constant(c)) ') (sim' int2str(sim) ').png']);
       close all
    end

    %% nullZ activity

    nullZ_neural_component = nullZ(neural_component);
    nullZ_OLS_dF = nullZ(OLS_dF);
    nullZ_OLS_dFF = nullZ(OLS_dFF);
    nullZ_Lp_OLS_dF = nullZ(Lp_OLS_dF);
    nullZ_Lp_OLS_dFF = nullZ(Lp_OLS_dFF);

    nullZ_IRLS_dF = cell(1, length(IRLS_constant));
    nullZ_IRLS_dFF = cell(1, length(IRLS_constant));
    nullZ_Lp_IRLS_dF = cell(1, length(IRLS_constant));
    nullZ_Lp_IRLS_dFF = cell(1, length(IRLS_constant));

    for c = 1:length(IRLS_constant)
      nullZ_IRLS_dF{c} = nullZ(IRLS_dF{c});
      nullZ_IRLS_dFF{c} = nullZ(IRLS_dFF{c});
      nullZ_Lp_IRLS_dF{c} = nullZ(Lp_IRLS_dF{c});
      nullZ_Lp_IRLS_dFF{c} = nullZ(Lp_IRLS_dFF{c});
    end

    %% Save simulations
    dFFsim.(sim_name).simulation(sim).ERT_idx = ERT_idx;
    dFFsim.(sim_name).simulation(sim).neural_component = neural_component;
    dFFsim.(sim_name).simulation(sim).nullZ_neural_component = nullZ_neural_component;
    dFFsim.(sim_name).simulation(sim).exp_decay_component = exp_decay_component;
    dFFsim.(sim_name).simulation(sim).iso_decay_component = iso_decay_component;
    dFFsim.(sim_name).simulation(sim).exp_noise_component = exp_noise_component;
    dFFsim.(sim_name).simulation(sim).iso_noise_component = iso_noise_component;
    dFFsim.(sim_name).simulation(sim).exp_signal = exp_signal;
    dFFsim.(sim_name).simulation(sim).iso_signal = iso_signal;
    dFFsim.(sim_name).simulation(sim).Lp_exp_signal = Lp_exp_signal;
    dFFsim.(sim_name).simulation(sim).Lp_iso_signal = Lp_iso_signal;

    %% Peri-event activity (raw)
    max_ev_window = max(ERT_idx(:,2)-ERT_idx(:,1));
    window_time = ...
      linspace(-pre_post_idx/samp_rate,(max_ev_window+pre_post_idx)/samp_rate,...
      max_ev_window+2*pre_post_idx+1);             
    c = turbo;    
    grad_col = interp1(c, linspace(1, size(c,1), n_events));
    linear_contrast = linspace(-1,1,n_events)';
    
    % True event
    [~, true_ev, true_bsl] = dFFsim_plot_ERT(neural_component,...
        ERT_idx,window_time,max_ev_window,pre_post_idx,'True',grad_col, ERT_ylims);
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[savefolder '\' sim_name '\ERT\true (sim' int2str(sim) ').png']);
    close all
    dFFsim_plot_linTrend(true_ev, true_bsl, linear_contrast, window_time, lintrend_ylims,...
        'True peri-event linear trend');
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[savefolder '\' sim_name '\ERT\true trend (sim' int2str(sim) ').png']);
    close all
  
    % OLS dF
    [~, OLS_dF_ev, OLS_dF_bsl] = dFFsim_plot_ERT(OLS_dF,...
        ERT_idx,window_time,max_ev_window,pre_post_idx,'OLS dF',grad_col, ERT_ylims);
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[savefolder '\' sim_name '\ERT\OLS dF (sim' int2str(sim) ').png']);
    close all
    dFFsim_plot_linTrend(OLS_dF_ev, OLS_dF_bsl, linear_contrast, window_time, lintrend_ylims,...
        'OLS dF peri-event linear trend');
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[savefolder '\' sim_name '\ERT\OLS dF trend (sim' int2str(sim) ').png']);
    close all
  
    % OLS dFF
    [~, OLS_dFF_ev, OLS_dFF_bsl] = dFFsim_plot_ERT(OLS_dFF,...
        ERT_idx,window_time,max_ev_window,pre_post_idx,'OLS dFF',grad_col, ERT_ylims);
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[savefolder '\' sim_name '\ERT\OLS dFF (sim' int2str(sim) ').png']);
    close all
    dFFsim_plot_linTrend(OLS_dFF_ev, OLS_dFF_bsl, linear_contrast, window_time, lintrend_ylims,...
        'OLS dFF peri-event linear trend');
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[savefolder '\' sim_name '\ERT\OLS dFF trend (sim' int2str(sim) ').png']);
    close all
  
    % Low-passed OLS dF
    [~, Lp_OLS_dF_ev, Lp_OLS_dF_bsl] = dFFsim_plot_ERT(Lp_OLS_dF,...
        ERT_idx,window_time,max_ev_window,pre_post_idx,['Lp' num2str(Lp_cutoff) ' OSL dF'],grad_col, ERT_ylims);
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[savefolder '\' sim_name '\ERT\Lp OLS dF(sim' int2str(sim) ').png']);
    close all
    dFFsim_plot_linTrend(Lp_OLS_dF_ev, Lp_OLS_dF_bsl, linear_contrast, window_time, lintrend_ylims,...
        ['Lp' num2str(Lp_cutoff) ' OLS dF peri-event linear trend']);
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[savefolder '\' sim_name '\ERT\Lp OLS dF trend(sim' int2str(sim) ').png']);
    close all
  
    % Low-passed OLS dFF
    [~, Lp_OLS_dFF_ev, Lp_OLS_dFF_bsl] = dFFsim_plot_ERT(Lp_OLS_dFF,...
        ERT_idx,window_time,max_ev_window,pre_post_idx,['Lp' num2str(Lp_cutoff) ' OLS dFF'],grad_col, ERT_ylims);
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[savefolder '\' sim_name '\ERT\Lp OLS dFF(sim' int2str(sim) ').png']);
    close all
    dFFsim_plot_linTrend(Lp_OLS_dFF_ev, Lp_OLS_dFF_bsl, linear_contrast, window_time, lintrend_ylims,...
        ['Lp' num2str(Lp_cutoff) ' OLS dFF peri-event linear trend']);
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[savefolder '\' sim_name '\ERT\Lp OLS dFF trend(sim' int2str(sim) ').png']);
    close all
  
    % IRLS dF [c = IRLS_constant]
    
    IRLS_dF_ev = cell(1, length(IRLS_constant));
    IRLS_dF_bsl = cell(1, length(IRLS_constant));

    for c = 1:length(IRLS_constant)
      [~, IRLS_dF_ev{c}, IRLS_dF_bsl{c}] = dFFsim_plot_ERT(IRLS_dF{c},...
          ERT_idx,window_time,max_ev_window,pre_post_idx,['IRLS (c = ' num2str(IRLS_constant(c)) ') dF'],grad_col, ERT_ylims);
      set(gcf,'Position',get(0,'Screensize'));
      saveas(figure(1),[savefolder '\' sim_name '\ERT\IRLS (c = ' num2str(IRLS_constant(c)) ') dF (sim' int2str(sim) ').png']);
      close all
      dFFsim_plot_linTrend(IRLS_dF_ev{c}, IRLS_dF_bsl{c}, linear_contrast, window_time, lintrend_ylims,...
        ['IRLS (c = ' num2str(IRLS_constant(c)) ') dF peri-event linear trend']);
      set(gcf,'Position',get(0,'Screensize'));
      saveas(figure(1),[savefolder '\' sim_name '\ERT\IRLS (c = ' num2str(IRLS_constant(c)) ') dF trend (sim' int2str(sim) ').png']);
      close all
    end

    % IRLS dFF [c = IRLS_constant]
    
    IRLS_dFF_ev = cell(1, length(IRLS_constant));
    IRLS_dFF_bsl = cell(1, length(IRLS_constant));

    for c = 1:length(IRLS_constant)
      [~, IRLS_dFF_ev{c}, IRLS_dFF_bsl{c}] = dFFsim_plot_ERT(IRLS_dFF{c},...
          ERT_idx,window_time,max_ev_window,pre_post_idx,['IRLS (c = ' num2str(IRLS_constant(c)) ') dFF'],grad_col, ERT_ylims);
      set(gcf,'Position',get(0,'Screensize'));
      saveas(figure(1),[savefolder '\' sim_name '\ERT\IRLS (c = ' num2str(IRLS_constant(c)) ') dFF (sim' int2str(sim) ').png']);
      close all
      dFFsim_plot_linTrend(IRLS_dFF_ev{c}, IRLS_dFF_bsl{c}, linear_contrast, window_time, lintrend_ylims,...
        ['IRLS (c = ' num2str(IRLS_constant(c)) ') dFF peri-event linear trend']);
      set(gcf,'Position',get(0,'Screensize'));
      saveas(figure(1),[savefolder '\' sim_name '\ERT\IRLS (c = ' num2str(IRLS_constant(c)) ') dFF trend (sim' int2str(sim) ').png']);
      close all
    end

    % Low-passed IRLS dF [c = IRLS_constant]
    
    Lp_IRLS_dF_ev = cell(1, length(IRLS_constant));
    Lp_IRLS_dF_bsl = cell(1, length(IRLS_constant));

    for c = 1:length(IRLS_constant)
      [~, Lp_IRLS_dF_ev{c}, Lp_IRLS_dF_bsl{c}] = dFFsim_plot_ERT(Lp_IRLS_dF{c},...
          ERT_idx,window_time,max_ev_window,pre_post_idx,['Lp IRLS (c = ' num2str(IRLS_constant(c)) ') dF'],grad_col, ERT_ylims);
      set(gcf,'Position',get(0,'Screensize'));
      saveas(figure(1),[savefolder '\' sim_name '\ERT\Lp IRLS (c = ' num2str(IRLS_constant(c)) ') dF (sim' int2str(sim) ').png']);
      close all
      dFFsim_plot_linTrend(Lp_IRLS_dF_ev{c}, Lp_IRLS_dF_bsl{c}, linear_contrast, window_time, lintrend_ylims,...
        ['Lp IRLS (c = ' num2str(IRLS_constant(c)) ') dF peri-event linear trend']);
      set(gcf,'Position',get(0,'Screensize'));
      saveas(figure(1),[savefolder '\' sim_name '\ERT\Lp IRLS (c = ' num2str(IRLS_constant(c)) ') dF trend (sim' int2str(sim) ').png']);
      close all
    end

    % Low-passed IRLS dFF [c = IRLS_constant]
    
    Lp_IRLS_dFF_ev = cell(1, length(IRLS_constant));
    Lp_IRLS_dFF_bsl = cell(1, length(IRLS_constant));

    for c = 1:length(IRLS_constant)
      [~, Lp_IRLS_dFF_ev{c}, Lp_IRLS_dFF_bsl{c}] = dFFsim_plot_ERT(Lp_IRLS_dFF{c},...
          ERT_idx,window_time,max_ev_window,pre_post_idx,['Lp IRLS (c = ' num2str(IRLS_constant(c)) ') dFF'],grad_col, ERT_ylims);
      set(gcf,'Position',get(0,'Screensize'));
      saveas(figure(1),[savefolder '\' sim_name '\ERT\Lp IRLS (c = ' num2str(IRLS_constant(c)) ') dFF (sim' int2str(sim) ').png']);
      close all
      dFFsim_plot_linTrend(Lp_IRLS_dFF_ev{c}, Lp_IRLS_dFF_bsl{c}, linear_contrast, window_time, lintrend_ylims,...
        ['Lp IRLS (c = ' num2str(IRLS_constant(c)) ') dFF peri-event linear trend']);
      set(gcf,'Position',get(0,'Screensize'));
      saveas(figure(1),[savefolder '\' sim_name '\ERT\Lp IRLS (c = ' num2str(IRLS_constant(c)) ') dFF trend (sim' int2str(sim) ').png']);
      close all
    end
  
    %% Peri-event activity (nullZ)
    % True event
    [~, nullZ_true_ev, nullZ_true_bsl] = dFFsim_plot_ERT(nullZ_neural_component,...
        ERT_idx,window_time,max_ev_window,pre_post_idx,'nullZ True',grad_col, ERT_Z_ylims);
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[savefolder '\' sim_name '\ERT\nullZ true (sim' int2str(sim) ').png']);
    close all
    dFFsim_plot_linTrend(nullZ_true_ev, nullZ_true_bsl, linear_contrast, window_time, lintrend_Z_ylims,...
        'nullZ True peri-event linear trend');
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[savefolder '\' sim_name '\ERT\nullZ true trend (sim' int2str(sim) ').png']);
    close all
  
    % OLS dF
    [~, nullZ_OLS_dF_ev, nullZ_OLS_dF_bsl] = dFFsim_plot_ERT(nullZ_OLS_dF,...
        ERT_idx,window_time,max_ev_window,pre_post_idx,'nullZ OLS dF',grad_col, ERT_Z_ylims);
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[savefolder '\' sim_name '\ERT\nullZ OLS dF (sim' int2str(sim) ').png']);
    close all
    dFFsim_plot_linTrend(nullZ_OLS_dF_ev, nullZ_OLS_dF_bsl, linear_contrast, window_time, lintrend_Z_ylims,...
        'nullZ OLS dF peri-event linear trend');
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[savefolder '\' sim_name '\ERT\nullZ OLS dF trend (sim' int2str(sim) ').png']);
    close all
  
    % OLS dFF
    [~, nullZ_OLS_dFF_ev, nullZ_OLS_dFF_bsl] = dFFsim_plot_ERT(nullZ_OLS_dFF,...
        ERT_idx,window_time,max_ev_window,pre_post_idx,'nullZ OLS dFF',grad_col, ERT_Z_ylims);
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[savefolder '\' sim_name '\ERT\nullZ OLS dFF (sim' int2str(sim) ').png']);
    close all
    dFFsim_plot_linTrend(nullZ_OLS_dFF_ev, nullZ_OLS_dFF_bsl, linear_contrast, window_time, lintrend_Z_ylims,...
        'nullZ OLS dFF peri-event linear trend');
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[savefolder '\' sim_name '\ERT\nullZ OLS dFF trend (sim' int2str(sim) ').png']);
    close all
  
    % Low-passed OLS dF
    [~, nullZ_Lp_OLS_dF_ev, nullZ_Lp_OLS_dF_bsl] = dFFsim_plot_ERT(nullZ_Lp_OLS_dF,...
        ERT_idx,window_time,max_ev_window,pre_post_idx,['nullZ Lp' num2str(Lp_cutoff) ' OLS dF'],grad_col, ERT_Z_ylims);
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[savefolder '\' sim_name '\ERT\nullZ Lp OLS dF(sim' int2str(sim) ').png']);
    close all
    dFFsim_plot_linTrend(nullZ_Lp_OLS_dF_ev, nullZ_Lp_OLS_dF_bsl, linear_contrast, window_time, lintrend_Z_ylims,...
        ['nullZ Lp' num2str(Lp_cutoff) ' OLS dF peri-event linear trend']);
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[savefolder '\' sim_name '\ERT\nullZ Lp OLS dF trend(sim' int2str(sim) ').png']);
    close all
  
    % Low-passed OLS dFF
   
    [~, nullZ_Lp_OLS_dFF_ev, nullZ_Lp_OLS_dFF_bsl] = dFFsim_plot_ERT(nullZ_Lp_OLS_dFF,...
        ERT_idx,window_time,max_ev_window,pre_post_idx,['nullZ Lp' num2str(Lp_cutoff) ' OLS dFF'],grad_col, ERT_Z_ylims);
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[savefolder '\' sim_name '\ERT\nullZ Lp OLS dFF(sim' int2str(sim) ').png']);
    close all
    dFFsim_plot_linTrend(nullZ_Lp_OLS_dFF_ev, nullZ_Lp_OLS_dFF_bsl, linear_contrast, window_time, lintrend_Z_ylims,...
        ['nullZ Lp' num2str(Lp_cutoff) ' OLS dFF peri-event linear trend']);
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[savefolder '\' sim_name '\ERT\nullZ Lp OLS dFF trend(sim' int2str(sim) ').png']);
    close all
  
    % IRLS dF [c = IRLS_constant]

    nullZ_IRLS_dF_ev = cell(1, length(IRLS_constant));
    nullZ_IRLS_dF_bsl = cell(1, length(IRLS_constant));

    for c = 1:length(IRLS_constant)
      [~, nullZ_IRLS_dF_ev{c}, nullZ_IRLS_dF_bsl{c}] = dFFsim_plot_ERT(nullZ_IRLS_dF{c},...
          ERT_idx,window_time,max_ev_window,pre_post_idx,['nullZ IRLS (c = ' num2str(IRLS_constant(c)) ') dF'],grad_col, ERT_Z_ylims);
      set(gcf,'Position',get(0,'Screensize'));
      saveas(figure(1),[savefolder '\' sim_name '\ERT\nullZ IRLS (c = ' num2str(IRLS_constant(c)) ') dF (sim' int2str(sim) ').png']);
      close all
      dFFsim_plot_linTrend(nullZ_IRLS_dF_ev{c}, nullZ_IRLS_dF_bsl{c}, linear_contrast, window_time, lintrend_Z_ylims,...
          ['nullZ IRLS (c = ' num2str(IRLS_constant(c)) ') dF peri-event linear trend']);
      set(gcf,'Position',get(0,'Screensize'));
      saveas(figure(1),[savefolder '\' sim_name '\ERT\nullZ IRLS (c = ' num2str(IRLS_constant(c)) ') dF trend (sim' int2str(sim) ').png']);
      close all
    end

    % IRLS dFF [c = IRLS_constant]

    nullZ_IRLS_dFF_ev = cell(1, length(IRLS_constant));
    nullZ_IRLS_dFF_bsl = cell(1, length(IRLS_constant));

    for c = 1:length(IRLS_constant)
      [~, nullZ_IRLS_dFF_ev{c}, nullZ_IRLS_dFF_bsl{c}] = dFFsim_plot_ERT(nullZ_IRLS_dFF{c},...
          ERT_idx,window_time,max_ev_window,pre_post_idx,['nullZ IRLS (c = ' num2str(IRLS_constant(c)) ') dFF'],grad_col, ERT_Z_ylims);
      set(gcf,'Position',get(0,'Screensize'));
      saveas(figure(1),[savefolder '\' sim_name '\ERT\nullZ IRLS (c = ' num2str(IRLS_constant(c)) ') dFF (sim' int2str(sim) ').png']);
      close all
      dFFsim_plot_linTrend(nullZ_IRLS_dFF_ev{c}, nullZ_IRLS_dFF_bsl{c}, linear_contrast, window_time, lintrend_Z_ylims,...
          ['nullZ IRLS (c = ' num2str(IRLS_constant(c)) ') dFF peri-event linear trend']);
      set(gcf,'Position',get(0,'Screensize'));
      saveas(figure(1),[savefolder '\' sim_name '\ERT\nullZ IRLS (c = ' num2str(IRLS_constant(c)) ') dFF trend (sim' int2str(sim) ').png']);
      close all
    end

    % Low-passed IRLS dF [c = IRLS_constant]

    nullZ_Lp_IRLS_dF_ev = cell(1, length(IRLS_constant));
    nullZ_Lp_IRLS_dF_bsl = cell(1, length(IRLS_constant));

    for c = 1:length(IRLS_constant)
      [~, nullZ_Lp_IRLS_dF_ev{c}, nullZ_Lp_IRLS_dF_bsl{c}] = dFFsim_plot_ERT(nullZ_Lp_IRLS_dF{c},...
          ERT_idx,window_time,max_ev_window,pre_post_idx,['nullZ Lp IRLS (c = ' num2str(IRLS_constant(c)) ') dF'],grad_col, ERT_Z_ylims);
      set(gcf,'Position',get(0,'Screensize'));
      saveas(figure(1),[savefolder '\' sim_name '\ERT\nullZ Lp IRLS (c = ' num2str(IRLS_constant(c)) ') dF (sim' int2str(sim) ').png']);
      close all
      dFFsim_plot_linTrend(nullZ_Lp_IRLS_dF_ev{c}, nullZ_Lp_IRLS_dF_bsl{c}, linear_contrast, window_time, lintrend_Z_ylims,...
          ['nullZ Lp IRLS (c = ' num2str(IRLS_constant(c)) ') dF peri-event linear trend']);
      set(gcf,'Position',get(0,'Screensize'));
      saveas(figure(1),[savefolder '\' sim_name '\ERT\nullZ Lp IRLS (c = ' num2str(IRLS_constant(c)) ') dF trend (sim' int2str(sim) ').png']);
      close all
    end

    % Low-passed IRLS dFF [c = IRLS_constant]

    nullZ_Lp_IRLS_dFF_ev = cell(1, length(IRLS_constant));
    nullZ_Lp_IRLS_dFF_bsl = cell(1, length(IRLS_constant));

    for c = 1:length(IRLS_constant)
      [~, nullZ_Lp_IRLS_dFF_ev{c}, nullZ_Lp_IRLS_dFF_bsl{c}] = dFFsim_plot_ERT(nullZ_Lp_IRLS_dFF{c},...
          ERT_idx,window_time,max_ev_window,pre_post_idx,['nullZ Lp IRLS (c = ' num2str(IRLS_constant(c)) ') dFF'],grad_col, ERT_Z_ylims);
      set(gcf,'Position',get(0,'Screensize'));
      saveas(figure(1),[savefolder '\' sim_name '\ERT\nullZ Lp IRLS (c = ' num2str(IRLS_constant(c)) ') dFF (sim' int2str(sim) ').png']);
      close all
      dFFsim_plot_linTrend(nullZ_Lp_IRLS_dFF_ev{c}, nullZ_Lp_IRLS_dFF_bsl{c}, linear_contrast, window_time, lintrend_Z_ylims,...
          ['nullZ Lp IRLS (c = ' num2str(IRLS_constant(c)) ') dFF peri-event linear trend']);
      set(gcf,'Position',get(0,'Screensize'));
      saveas(figure(1),[savefolder '\' sim_name '\ERT\nullZ Lp IRLS (c = ' num2str(IRLS_constant(c)) ') dFF trend (sim' int2str(sim) ').png']);
      close all
    end

    %% Save results
    % Session signals
    dFFsim.(sim_name).results(sim).OLS_ft405 = OLS_ft405;
    dFFsim.(sim_name).results(sim).OLS_dF = OLS_dF;
    dFFsim.(sim_name).results(sim).OLS_dFF = OLS_dFF;
    dFFsim.(sim_name).results(sim).nullZ_OLS_dF = nullZ_OLS_dF;
    dFFsim.(sim_name).results(sim).nullZ_OLS_dFF = nullZ_OLS_dFF;
  
    dFFsim.(sim_name).results(sim).Lp_OLS_ft405 = Lp_OLS_ft405;
    dFFsim.(sim_name).results(sim).Lp_OLS_dF = Lp_OLS_dF;
    dFFsim.(sim_name).results(sim).Lp_OLS_dFF = Lp_OLS_dFF;
    dFFsim.(sim_name).results(sim).nullZ_Lp_OLS_dF = nullZ_Lp_OLS_dF;
    dFFsim.(sim_name).results(sim).nullZ_Lp_OLS_dFF = nullZ_Lp_OLS_dFF;
    
      % Results for IRLS (per constant)
      for c = 1:length(IRLS_constant)
        dFFsim.(sim_name).results(sim).(['IRLS' num2str(c) '_ft405']) = IRLS_ft405{c};
        dFFsim.(sim_name).results(sim).(['IRLS' num2str(c) '_dF']) = IRLS_dF{c};
        dFFsim.(sim_name).results(sim).(['IRLS' num2str(c) '_dFF']) = IRLS_dFF{c};
        dFFsim.(sim_name).results(sim).(['nullZ_IRLS' num2str(c) '_dF']) = nullZ_IRLS_dF{c};
        dFFsim.(sim_name).results(sim).(['nullZ_IRLS' num2str(c) '_dFF']) = nullZ_IRLS_dFF{c};
        dFFsim.(sim_name).results(sim).(['Lp_IRLS' num2str(c) '_ft405']) = Lp_IRLS_ft405{c};
        dFFsim.(sim_name).results(sim).(['Lp_IRLS' num2str(c) '_dF']) = Lp_IRLS_dF{c};
        dFFsim.(sim_name).results(sim).(['Lp_IRLS' num2str(c) '_dFF']) = Lp_IRLS_dFF{c};
        dFFsim.(sim_name).results(sim).(['nullZ_Lp_IRLS' num2str(c) '_dF']) = nullZ_Lp_IRLS_dF{c};
        dFFsim.(sim_name).results(sim).(['nullZ_Lp_IRLS' num2str(c) '_dFF']) = nullZ_Lp_IRLS_dFF{c};
      end

    % ERT
    dFFsim.(sim_name).results(sim).true_ev = true_ev;
    dFFsim.(sim_name).results(sim).true_bsl = true_bsl;
    dFFsim.(sim_name).results(sim).nullZ_true_ev = nullZ_true_ev;
    dFFsim.(sim_name).results(sim).nullZ_true_bsl = nullZ_true_bsl;
    
    dFFsim.(sim_name).results(sim).OLS_dF_ev = OLS_dF_ev; 
    dFFsim.(sim_name).results(sim).OLS_dF_bsl = OLS_dF_bsl;
    dFFsim.(sim_name).results(sim).OLS_dFF_ev = OLS_dFF_ev; 
    dFFsim.(sim_name).results(sim).OLS_dFF_bsl = OLS_dFF_bsl;
    dFFsim.(sim_name).results(sim).nullZ_OLS_dF_ev = nullZ_OLS_dF_ev; 
    dFFsim.(sim_name).results(sim).nullZ_OLS_dF_bsl = nullZ_OLS_dF_bsl;
    dFFsim.(sim_name).results(sim).nullZ_OLS_dFF_ev = nullZ_OLS_dFF_ev; 
    dFFsim.(sim_name).results(sim).nullZ_OLS_dFF_bsl = nullZ_OLS_dFF_bsl;
    
    dFFsim.(sim_name).results(sim).Lp_OLS_dF_ev = Lp_OLS_dF_ev; 
    dFFsim.(sim_name).results(sim).Lp_OLS_dF_bsl = Lp_OLS_dF_bsl;
    dFFsim.(sim_name).results(sim).Lp_OLS_dFF_ev = Lp_OLS_dFF_ev; 
    dFFsim.(sim_name).results(sim).Lp_OLS_dFF_bsl = Lp_OLS_dFF_bsl;
    dFFsim.(sim_name).results(sim).nullZ_Lp_OLS_dF_ev = nullZ_Lp_OLS_dF_ev; 
    dFFsim.(sim_name).results(sim).nullZ_Lp_OLS_dF_bsl = nullZ_Lp_OLS_dF_bsl;
    dFFsim.(sim_name).results(sim).nullZ_Lp_OLS_dFF_ev = nullZ_Lp_OLS_dFF_ev; 
    dFFsim.(sim_name).results(sim).nullZ_Lp_OLS_dFF_bsl = nullZ_Lp_OLS_dFF_bsl;

      % Results for IRLS (per constant)
      for c = 1:length(IRLS_constant)
        dFFsim.(sim_name).results(sim).(['IRLS' num2str(c) '_dF_ev']) = IRLS_dF_ev{c};
        dFFsim.(sim_name).results(sim).(['IRLS' num2str(c) '_dF_bsl']) = IRLS_dF_bsl{c};
        dFFsim.(sim_name).results(sim).(['IRLS' num2str(c) '_dFF_ev']) = IRLS_dFF_ev{c};
        dFFsim.(sim_name).results(sim).(['IRLS' num2str(c) '_dFF_bsl']) = IRLS_dFF_bsl{c};
        dFFsim.(sim_name).results(sim).(['nullZ_IRLS' num2str(c) '_dF_ev']) = nullZ_IRLS_dF_ev{c};
        dFFsim.(sim_name).results(sim).(['nullZ_IRLS' num2str(c) '_dF_bsl']) = nullZ_IRLS_dF_bsl{c};
        dFFsim.(sim_name).results(sim).(['nullZ_IRLS' num2str(c) '_dFF_ev']) = nullZ_IRLS_dFF_ev{c};
        dFFsim.(sim_name).results(sim).(['nullZ_IRLS' num2str(c) '_dFF_bsl']) = nullZ_IRLS_dFF_bsl{c};
        dFFsim.(sim_name).results(sim).(['Lp_IRLS' num2str(c) '_dF_ev']) = Lp_IRLS_dF_ev{c};
        dFFsim.(sim_name).results(sim).(['Lp_IRLS' num2str(c) '_dF_bsl']) = Lp_IRLS_dF_bsl{c};
        dFFsim.(sim_name).results(sim).(['Lp_IRLS' num2str(c) '_dFF_ev']) = Lp_IRLS_dFF_ev{c};
        dFFsim.(sim_name).results(sim).(['Lp_IRLS' num2str(c) '_dFF_bsl']) = Lp_IRLS_dFF_bsl{c};
        dFFsim.(sim_name).results(sim).(['nullZ_Lp_IRLS' num2str(c) '_dF_ev']) = nullZ_Lp_IRLS_dF_ev{c};
        dFFsim.(sim_name).results(sim).(['nullZ_Lp_IRLS' num2str(c) '_dF_bsl']) = nullZ_Lp_IRLS_dF_bsl{c};
        dFFsim.(sim_name).results(sim).(['nullZ_Lp_IRLS' num2str(c) '_dFF_ev']) = nullZ_Lp_IRLS_dFF_ev{c};
        dFFsim.(sim_name).results(sim).(['nullZ_Lp_IRLS' num2str(c) '_dFF_bsl']) = nullZ_Lp_IRLS_dFF_bsl{c};
      end
    
    fprintf('Done. ');
    diary OFF
    toc;
  end
  
%% Calculate residuals and absolute residuals from true signal

  result_names = fieldnames(dFFsim.(sim_name).results(sim));
  res_fields = {};
  
  for i = 1:numel(result_names)
      field_name = result_names{i};
      if startsWith(field_name, 'nullZ') && endsWith(field_name, 'F')
          res_fields{end + 1} = field_name;
      end
  end

  comparison = 'nullZ_neural_component';
  
  for f = 1:length(res_fields)
    for sim = 1:n_sims
      dFFsim.(sim_name).results(sim).([res_fields{f} '_res']) = ...
        dFFsim.(sim_name).results(sim).(res_fields{f})-dFFsim.(sim_name).simulation(sim).(comparison);
     end
  end 
      
  for f = 1:length(res_fields)
    for sim = 1:n_sims
      figure; hold on
      plot(dFFsim.(sim_name).results(sim).([res_fields{f} '_res']));
      ylim([-2.5 2.5]);
      plot(xlim,[0 0],'k--');
      title([res_fields{f} ' vs ' comparison ' residual'],'interpreter','none');
      saveas(figure(1),[savefolder '\' sim_name '\Residual\' res_fields{f} ' v ' comparison ' residual (sim' int2str(sim) ').png']);
      close all
    end
  end
  
  event_idx = dFFsim.(sim_name).simulation(1).(comparison) ~= 0;
  for f = 1:length(res_fields)
    for sim = 1:n_sims
      data = dFFsim.(sim_name).results(sim).([res_fields{f} '_res']);
      dFFsim.(sim_name).results(sim).([res_fields{f} '_mean_res']) = mean(data);
      dFFsim.(sim_name).results(sim).([res_fields{f} '_meanEv_res']) = mean(data(event_idx));
      dFFsim.(sim_name).results(sim).([res_fields{f} '_meanBsl_res']) = mean(data(~event_idx));
      dFFsim.(sim_name).results(sim).([res_fields{f} '_mean_abres']) = mean(abs(data));
      dFFsim.(sim_name).results(sim).([res_fields{f} '_meanEv_abres']) = mean(abs(data(event_idx)));
      dFFsim.(sim_name).results(sim).([res_fields{f} '_meanBsl_abres']) = mean(abs(data(~event_idx)));
    end
  end

%% Done
clearvars -except dFFsim savefolder sim_name
toc;
save([savefolder '\' sim_name '_data.mat']);
fprintf('Saved. ');
toc;