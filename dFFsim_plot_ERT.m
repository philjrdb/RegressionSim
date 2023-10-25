function [F1, ev, bsl] = dFFsim_plot_ERT(signal,ERT_idx,window_time,max_ev_window,pre_post_idx,title_tag,grad_col,ylims)

n_events = size(ERT_idx,1);
win_size = max_ev_window+pre_post_idx*2+1;

ev = NaN(n_events,win_size);
bsl = NaN(n_events,1);

F1 = figure; 
for e = 1:n_events
  ev(e,:) = signal(ERT_idx(e,1)-pre_post_idx:...
    ERT_idx(e,1)+max_ev_window+pre_post_idx);
  bsl(e) = mean(ev(e,1:pre_post_idx));

  subplot(1,2,1); hold on
  plot(window_time,ev(e,:),'Color',grad_col(e,:));

  subplot(1,2,2); hold on
  plot(window_time,ev(e,:)-bsl(e),'Color',grad_col(e,:));
end

subplot(1,2,1); hold on
plot(window_time,mean(ev,1),'b-','LineWidth',2);
xlim([window_time(1) window_time(end)]);
if ~isempty(ylims)
  ylim(ylims)
end
plot([0 0],ylim,'k:');
plot(xlim,[0 0],'k--');
title([title_tag ' peri-event transient (unbaselined)']);
hold on
dummy_line1 = plot([0,0], [0,0], 'Color', grad_col(1,:), 'LineWidth', 2, 'DisplayName', 'Start Color');
dummy_line2 = plot([0,0], [0,0], 'Color', grad_col(ceil(n_events/2),:), 'LineWidth', 2, 'DisplayName', 'Middle Color');
dummy_line3 = plot([0,0], [0,0], 'Color', grad_col(end,:), 'LineWidth', 2, 'DisplayName', 'End Color');
legend([dummy_line1, dummy_line2, dummy_line3]);

subplot(1,2,2); hold on
plot(window_time,mean(ev-bsl,1),'b-','LineWidth',2);
xlim([window_time(1) window_time(end)]);
if ~isempty(ylims)
  ylim(ylims)
end
plot([0 0],ylim,'k:');
plot(xlim,[0 0],'k--');
title([title_tag ' peri-event transient (baselined)']);
hold on
dummy_line1 = plot([0,0], [0,0], 'Color', grad_col(1,:), 'LineWidth', 2, 'DisplayName', 'Start Color');
dummy_line2 = plot([0,0], [0,0], 'Color', grad_col(ceil(n_events/2),:), 'LineWidth', 2, 'DisplayName', 'Middle Color');
dummy_line3 = plot([0,0], [0,0], 'Color', grad_col(end,:), 'LineWidth', 2, 'DisplayName', 'End Color');
legend([dummy_line1, dummy_line2, dummy_line3]);