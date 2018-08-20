function [j_rss] = solve_HJB_rss(r_rss,j_0_mat,parameters)
% This version codes Galo's solution using standard finite difference methods.
% junk....
                                
% Unpacking vectors
time_preallocate;

% model parameters
delta    = parameters.delta ; % Coupon
phi      = parameters.phi   ; % intensity of the shock
prob_vec = parameters.prob_vec;

% Pre-Allocation - Build HJB_rss Operator
u         = delta + phi*(j_0_mat*prob_vec); % vector of coupons + change in valuation
u(1)      = u(1) +1/dt                  ;
A0        = -speye(N_tau)               ;
aa        = ones(N_tau-1,1)             ;
A1        = spdiags(aa,-1,N_tau,N_tau)  ;
A         = 1/dt*(A0+A1)                ;

% Risky Steady State Solution
j_rss     = ((r_rss+phi)*speye(N_tau) - A)\u;


