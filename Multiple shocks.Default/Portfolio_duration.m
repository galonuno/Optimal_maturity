function [average_duration, average_term_to_maturity] = Portfolio_duration(results, parameters)
% Duration
% This codes compute the duration of a portfolio
parameter_unpack;
time_preallocate;
delta = parameters.delta;
r_n = results.r_n;
f_n = results.f_n;

duration_bond   = zeros(N_tau,N_t);
PV              = zeros(N_tau,N_t);
w_PV            = zeros(N_tau,N_t);
average_term_to_maturity = zeros(1,N_t);
average_duration = zeros(1,N_t);

for i=1:N_t
    average_term_to_maturity(i) = sum(f_n(:,i).*tau'*dt)/sum(f_n(:,i)*dt);
    if i > N_t - N_tau - 1
            average_duration(i) = average_duration(i-1); % close to the steady state
    else             
        discount  = exp(-r_n(i:i+N_tau-1).*tau);
        PV(1,i)   = discount(1);
        w_PV(1,i) = discount(1)*tau(1);
        duration_bond(1,i) = w_PV(1,i)/PV(1,i);
        for j=2:N_tau
            PV(j,i)   = delta*sum(discount(1:j-1)) + discount(j);
            w_PV(j,i) = delta*sum(discount(1:j-1).*tau(1:j-1)) + discount(j)*tau(j);
            duration_bond(j,i) = w_PV(j,i)/PV(j,i);
        end
        average_duration(i) = sum(f_n(:,i).*duration_bond(:,i)*dt)/sum(f_n(:,i)*dt);
    end    
end
