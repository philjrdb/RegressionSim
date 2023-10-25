function [decay_rate] = double_exp_decay(rate1,rate2,time)
% Double exponential decay rate

decay_rate = (((1-rate1).^time)+((1-rate2).^time))/2;