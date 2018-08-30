%% Preallocate Time Vectors
t_max     = parameters.t_max  ; % Maximum number of years (time)
tau_max   = parameters.tau_max; % Maximum number of years (maturity)
dt        = parameters.dt     ; % Monthly steps

% Preallocate Maturity/Time Vector
tau       = dt:dt:tau_max     ; % tau vector
t         = 0:dt:t_max        ; % time vector
N_t       = length(t)         ; % number of time nodes
tau       = dt:dt:tau_max     ; % maturities
N_tau     = length(tau)       ; % number of tau nodes