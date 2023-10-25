function F1 = dFFsim_plot_linTrend(ev, bsl, contrast_coeff, window_time, ylims, title_tags)
ev_lin_trend = mean(ev.*contrast_coeff,1);
bsl_ev_lin_trend = mean((ev-bsl).*contrast_coeff,1);

F1 = figure;
subplot(1,2,1); hold on
plot(window_time,ev_lin_trend,'r-','LineWidth',2);
xlim([window_time(1) window_time(end)]);
ylim(ylims);
plot([0 0],ylims,'k:');
plot(xlim,[0 0],'k--');
title([title_tags ' (unbaselined)']);
subplot(1,2,2); hold on
plot(window_time,bsl_ev_lin_trend,'r-','LineWidth',2);
xlim([window_time(1) window_time(end)]);
ylim(ylims);
plot([0 0],ylims,'k:');
plot(xlim,[0 0],'k--');
title([title_tags ' (baselined)']);