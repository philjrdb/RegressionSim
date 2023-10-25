function [Fig, dF, dFF] = dFFsim_plot_dFsignals(exp_signal,iso_signal,fitted_iso,event_component,time,title_tag)

dF = exp_signal-fitted_iso;
dFF = (exp_signal-fitted_iso)./fitted_iso;
Fig = figure;
subplot(4,1,1); hold on
plot(time,exp_signal,'g-');
plot(time,iso_signal,'m-');
plot(time,fitted_iso,'b-');
title([title_tag ' signals']);
subplot(4,1,2); hold on
plot(time,event_component);
plot(xlim,[0 0],'k:');
title('Event component');
subplot(4,1,3); hold on
plot(time,dF);
plot(xlim,[0 0],'k:');
title([title_tag ' dF']);
subplot(4,1,4); hold on
plot(time,dFF);
plot(xlim,[0 0],'k:');
title([title_tag ' dFF']);