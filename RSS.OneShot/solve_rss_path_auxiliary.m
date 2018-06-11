function [results] = solve_rss_path(rss,steady,paths,parameters)
% Solve_RSS_path
% solves the RSS for one shock iterating over {c_rss}U{c(o)-->c(ss)}
% This version codes Galo's solution using standard finite difference methods.
% Transition path, no shocks

%% Model Parameters
parameter_unpack;
time_preallocate;

% Probability of the Shock
phi     = parameters.phi;          % Poisson arrival rate

% Recover Prices 
f0      = steady.f_ss;
c_in    = zeros(1,N_t)+steady.c_ss;
c_rss_in= steady.c_ss;

%% ----------------------------
%MAIN LOOP
tic
for iter = 1:1
    
    % Construct Deterministic path
    [path_out]=compute_path(c_n_in,f0,rss,paths,parameters);
end
if iter == max_iter
    disp('Error: maximum number of iterations reached')
end

toc
results.V_rss       = Vrss;    % Value function
results.iota_rss   = iota_rss;% Policy function
results.f_rss       = fRSS;    % Density
results.c_rss      = c_rss;   % Consumption
results.c_after    = c_after; % Consumption
results.psi_rss    = psi_rss;
results.r_rss      = rho;

%% Function that solves the consumption path
function [path_out]=compute_path(c_in,f0,rss,paths,parameters);
% Unpack Model Parameters
parameter_unpack;
time_preallocate;

% Unpack Parameters
rho   = parameters.rho  ;
sigma = parameters.gamma;
delta = parameters.delta;
lambda_bar= parameters.lambda_bar;
U_p_ratio = parameters.U_p_ratio;

% Unpack Exogenous Paths
psi_n   = paths.psi_n; 
y_n     = paths.y_n  ;

% Unpack Risky Steady State
psi_rss   = rss.psi_rss; 

%% Step 1: Construct r_n
% [!] Warning, I am assuming that path reaches a proper STEADY STATE
% [!!] doing this one at a time
% [!!!] not sure what to do at non-convergent solution
c_rss       = c_in(1);
c_0         = c_in(2); 
c_in(2)     = c_in(2:N_t+1);
r_n(1:N_t-1)= rho+sigma*(c_in(2:N_t)./c_in(1:N_t-1)-1)*dt;
r_n(N_t,1)  = rho;

%% Step 2: Value at Deterministic Path
[v_n]= solve_HJB_path(r_n,rho,parameters);
v_0=v_n(1);

%% Step 1: Solving for Risky Steady State
v_0_mat=v_0*U_p_ratio(c_0,c_rss)             ;
[v_rss]=solve_HJB_rss(rho,v_0_mat,parameters);

%% Compute paths

% Optimal policies
iota_n = 1/lambda_bar*(1-v_n./psi_n); % Optimal policies

% Construct - KFE Operator
A0        = -speye(N_tau);
aa        = ones(N_tau-1,1);
A1        = spdiags(aa,-1,N_tau,N_tau);
A         = 1/dt*(A0+A1);
D         = speye(N_tau) - dt *A';

% KFE 
f_n(:,1)   = f0;
for n=2:N_t
    h      = iota_n(:,n) *dt + f_n(:,n-1);
    f_aux  = D\h;
    f_n(:,n) = f_aux;
end

% Consumption and interest rates
c_n = y_n  - f_n(1,:) + sum((psi_n-1/2*psi_n*lambda_bar*iota_n).*iota_n-delta*f_n) * dt;
    


