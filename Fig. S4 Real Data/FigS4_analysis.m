%% Load in data 

data = readmatrix("2_R1_LateAcq_signals.csv");
EventTimes = readmatrix("2_R1_LateAcq_EventTimes.csv");

time = data(:,1);
r465 = data(:,2);
r405 = data(:,3);

%% Identify and exclude invalid periods

fig_fullSess = figure; hold on
plot(time,r465,'g-');
plot(time,r405,'m-');
xlabel('Session time (secs)');
ylabel('mV');
title('Signals from whole session');
legend({'465nm signal' '405nm signal'});

valid_idx = (time>33) & (time<1243);
valid_time = time(valid_idx);
valid_465 = r465(valid_idx);
valid_405 = r405(valid_idx);

% Plot signals across session (valid only)
fig_validSess = figure; hold on
plot(valid_time,valid_465,'g-');
plot(valid_time,valid_405,'m-');
xlabel('Session time (secs)');
ylabel('mV');
title('Signals from valid portion of session');
legend({'465nm signal' '405nm signal'});

%% Fit isosbestic to experimental with OLS and IRLS

% OLS regression
coeffs = polyfit(valid_405,valid_465,1);
fitted_405 = polyval(coeffs,valid_405);

% Robust regression
rob_coeffs = robustfit(valid_405,valid_465,'bisquare',1.4,'on');
robust_ft405 = valid_405*rob_coeffs(2)+rob_coeffs(1);

fig_fitted405 = figure; hold on
plot(valid_time,valid_465,'g-');
plot(valid_time, fitted_405,'r-');
plot(valid_time, robust_ft405,'b-');
xlabel('Session time (secs)');
ylabel('mV');
title('OLS vs IRLS fitted control signals');

%% Calculate dF/Fs

% OLS-based dF/F
OLS_dFF = (valid_465-fitted_405)./fitted_405;
OLS_dFF = OLS_dFF*100;

% Robust-based dF/F
rob_dFF = (valid_465-robust_ft405)./robust_ft405;
rob_dFF = rob_dFF*100;

%% Peri-event dF/Fs

beer_times = EventTimes(:,1);
beer_times = beer_times(~isnan(beer_times)); % remove NaN values
beer_times = beer_times(beer_times>33 & beer_times<1243); % isolate valid times
n_beer = length(beer_times);

% Seconds around event
pre = 3;
post = 10;
baseline = 2; % pre-event period to zero signals to (assumes pre-event serves as a better baseline)

% Determine sampling rate, get peri-event window
fs = length(valid_time)/(valid_time(end)-valid_time(1));
pre_idx = round(fs*pre);
post_idx = round(fs*post);
bsl_idx = round(fs*baseline);
peri_event_time = linspace(-pre,post,pre_idx+post_idx+1);

% Pre-allocate space for arrays
beer_idx = NaN(n_beer,1);
pre_beer_bsl_OLS = NaN(n_beer,1);
pre_beer_bsl_rob = NaN(n_beer,1);
beer_dFF_OLS = NaN(n_beer,pre_idx+post_idx+1);
beer_dFF_rob = NaN(n_beer,pre_idx+post_idx+1);

% Iterate through beer delivery trials 
for t = 1:n_beer
    [~,beer_idx(t)] = min(abs(valid_time-beer_times(t)));
    % If valid portion of signal exists, collate dF/F
    if (beer_idx(t)-pre_idx)<1 || (beer_idx(t)+post_idx)>length(valid_time)
        warning(['Trial #' int2str(t) ' window extends beyond valid portions of session']);
    else 
        beer_dFF_OLS(t,:) = OLS_dFF(beer_idx(t)-pre_idx:beer_idx(t)+post_idx);
        beer_dFF_rob(t,:) = rob_dFF(beer_idx(t)-pre_idx:beer_idx(t)+post_idx);
        pre_beer_bsl_OLS(t) = mean(beer_dFF_OLS(t,1:bsl_idx));
        pre_beer_bsl_rob(t) = mean(beer_dFF_rob(t,1:bsl_idx));
    end
end

% Baselined dFF around beer delivery
bsld_beer_dFF_OLS = beer_dFF_OLS-pre_beer_bsl_OLS;
bsld_beer_dFF_rob = beer_dFF_rob-pre_beer_bsl_rob;

% CI and peri-event plot for OLS dFF
[sig,~,beer_CI_OLS] = ttest(beer_dFF_OLS);
sig_idx = find(sig);
thresh = fs/3;
beer_sig = NaN(1,length(sig));
beer_sig(sig_idx(consec_idx(sig_idx,thresh))) = 1;

fig_OLSbeerCI = figure; hold on
errorplot3(beer_CI_OLS(1,:), beer_CI_OLS(2,:), [-pre post], [.8 .8 .8], .25);
plot(peri_event_time,mean(beer_dFF_OLS,1,'omitnan'),'k','LineWidth',2);
yl = ylim;
plot(peri_event_time,beer_sig*yl(2),'LineWidth',5);
xlim([-pre post]);
plot(xlim,[0 0],'k--');
plot([0 0],ylim,'k:');
xlabel('Secs from onset');
ylabel('% dFF');
title('% dFF (OLS) around beer deliveries');
legend({'95% CI' 'Mean dFF' 'Sig diff vs 0'},'Location','best');

% CI and peri-event plot for IRLS dFF
[sig,~,beer_CI_rob] = ttest(beer_dFF_rob);
sig_idx = find(sig);
thresh = fs/3;
beer_sig = NaN(1,length(sig));
beer_sig(sig_idx(consec_idx(sig_idx,thresh))) = 1;

fig_robbeerCI = figure; hold on
errorplot3(beer_CI_rob(1,:), beer_CI_rob(2,:), [-pre post], [.8 .8 .8], .25);
plot(peri_event_time,mean(beer_dFF_rob,1,'omitnan'),'k','LineWidth',2);
yl = ylim;
plot(peri_event_time,beer_sig*yl(2),'LineWidth',5);
xlim([-pre post]);
plot(xlim,[0 0],'k--');
plot([0 0],ylim,'k:');
xlabel('Secs from onset');
ylabel('% dFF');
title('% dFF (robust) around beer deliveries');
legend({'95% CI' 'Mean dFF' 'Sig diff vs 0'},'Location','best');

%% Done
toc;
