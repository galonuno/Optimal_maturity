function [average_duration, average_term_to_maturity] = Duration(steady, parameters)
% Duration
% This codes compute the duration of a portfolio
parameter_unpack;
time_preallocate;

r_bar = parameters.r_bar_ss;
delta = parameters.delta;

psi1     = delta/r_bar*(1-exp(-r_bar*tau'));
psi2     = exp(-r_bar*tau');

duration_inst = 1./(psi1 + psi2) .* ((1-delta/r_bar) .*tau'.*exp(-r_bar*tau') + psi1 /r_bar);
average_duration = sum(steady.f_ss.*duration_inst*dt)/sum(steady.f_ss*dt);
average_term_to_maturity =  sum(steady.f_ss.*tau'*dt)/sum(steady.f_ss*dt);