function [indices] = issuance_indexes(parameters)

%This Function provides the indi
Issued_Maturities = parameters.Issued_Maturities;

t_max     = parameters.t_max;    % Maximum number of years (time)
tau_max   = parameters.tau_max;  % Maximum number of years (maturity)
dt        = parameters.dt;
t_max     = parameters.t_max;    % Maximum number of years (time)
tau_max   = parameters.tau_max;  % Maximum number of years (maturity)
dt        = parameters.dt;       % Monthly steps
t         = 0:dt:t_max;      % calendar time
N         = length(t);       % number of time nodes
tau       = dt:dt:tau_max;   % maturities
I         = length(tau);     % number of tau nodes
S            =  parameters.S;       %number of shocks to the interest rate
T         = tau_max;

Index_1_year = round (I/T);
indices           = round (Index_1_year*Issued_Maturities);
end

