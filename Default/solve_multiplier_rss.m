function mu_rss=solve_multiplier_rss(parameters,rss)
% Computes the value of the Lagrange multiplier in the RSS

% Unpacking vectors
time_preallocate;

% model parameters
rho         = parameters.rho;
r_bar_ss    = parameters.r_bar_ss;
dU          = parameters.dU;
lambda_bar  = parameters.lambda_bar;

% Variables
c_rss     = rss.c_rss   ; 
iota_rss  = rss.iota_rss;

%% Pre-Allocation - Build HJB_ss Operator
u         = dU(c_rss) .*iota_rss .*(1-0.5*lambda_bar*iota_rss)  ;   % vector of coupons 
A0        = -speye(N_tau);
aa        = ones(N_tau-1,1);
A1        = spdiags(aa,-1,N_tau,N_tau);
A         = 1/dt*(A0+A1);

% Solution 
mu_rss       = ((r_bar_ss-rho)*speye(N_tau) - A')\u;

