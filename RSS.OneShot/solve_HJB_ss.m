function [j_ss]=solve_HJB_ss(r_ss,parameters);
% Unpacking vectors
time_preallocate;

% model parameters
delta = parameters.delta;

%% Pre-Allocation - Build HJB_ss Operator
u         = zeros (N_tau,1);   % vector of coupons 
u(:)      = delta;         % vector of coupons 
u(1)      = u(1) +1/dt;
A0        = -speye(N_tau);
aa        = ones(N_tau-1,1);
A1        = spdiags(aa,-1,N_tau,N_tau);
A         = 1/dt*(A0+A1);

% Solution 
j_ss       = (r_ss*speye(N_tau) - A)\u;
