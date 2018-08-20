% Code to compute the Risky Steady State
%(c) Bigio, Nuno & Passadore model
% This version codes Galo's solution using standard finite difference methods.
% The Code Solves for the Risky Steady State with default
close all; % clear;
plotit=1;

%% Code parameters
% Tolerance
parameters.tol_path   = 1e-3;  % tolerance parameter optimal path
parameters.tol        = 1e-5;  % tolerance parameter
parameters.max_iter   = 5000;  % maximum number of iterations
parameters.relax      = 0.1 ;  % relaxation coefficient
parameters.relax_path = 0.01;  % relaxation coefficient optimal path

% Numerical parameters
parameters.t_max     = 60    ; % maximum number of years (time)
parameters.tau_max   = 20    ; % maximum number of years (maturity)
parameters.dt        = 1/12  ; % monthly steps

% Method to Solve the model
parameters.solve_method='iter'; % iter for iterative, solver for solver

% Number of Shocks
N_r=   1;
N_y=   1;

%% Model parameters
% Time Period annual
% Preference Parameters
parameters.gamma   = 2.00     ; % risk aversion
parameters.rho     = 0.06   ; % discount factor - Model is Quarterly

% Output parameters
parameters.y_ss         = 1.0    ; % output at steady state
parameters.y_0(1:N_y)   = 1.0    ; % initial value, output after a shock
parameters.rho_y(1:N_y) = 0.05    ; % persistence output after shock - half life 1/rho months

%Interest Rate Shock Parameters
parameters.r_bar_ss           = 0.04                ; % steady-state short rate
parameters.delta              = parameters.r_bar_ss ; % coupon = r_ss so we have bonds at par
parameters.r_bar_0(1:N_r,1)   = 0.04;%0.041               ; % initial value
parameters.rho_r_bar(1:N_r,1) = 0.05                ; % persistence (in years)

%Risky Steady State Parameters
parameters.phi  = 0.1       ; % Poisson arrival rate

% adjustment cost
parameters.lambda_bar = 60  ; %

% Default parameters - Uniform distribution
parameters.V_low = -9;
parameters.V_high = 1;

% there's only one state
parameters.prob_vec=1         ;

% Price of Risk
parameters.U         =  @(c) (c.^(1-parameters.gamma)-1)/(1-parameters.gamma);
parameters.U_p_ratio =  @(c_a,c_b) (c_a/c_b).^(-parameters.gamma);
parameters.dU        =  @(c) (c.^(-parameters.gamma));

%% Pre-Allocation
time_preallocate;

%% Paths
% Path for the Short Rate and output
exo_path_construct;

%% Save Exogenous Paths
paths.r_bar_n = r_bar_n;
paths.y_n     = y_n    ;

%% Step 1: Bond Prices in Steady State
psi_n          = solve_HJB_path(r_bar_n,r_bar_ss,parameters) ;
psi_0_mat      = psi_n(:,1)                                  ;
psi_ss         = solve_HJB_ss(r_bar_ss,parameters)           ;
rss.psi_0_mat  = psi_0_mat;
% results
results.psi_ss  = psi_ss ;

%% Step 2: Value functions in Steady State 
steady=solve_steady(psi_ss,parameters);
% Save Steady States
steady.psi_ss  = psi_ss;
paths.psi_n    = psi_n   ;

%% Step 2: Risky-Steady State
path_out        = solve_rss_path(rss,steady,paths,parameters);
rss.psi_rss     = path_out.psi_rss  ;
results.psi_rss = path_out.psi_rss;
paths.psi_rss   = path_out.psi_rss             ;
results.yield_rss = parameters.delta-log(results.psi_rss)./(tau');


%% Step 3: Pack Results
results.f_n         = path_out.f_n          ;
results.f_rss       = path_out.f_out(:,1)   ;
results.iota_n      = path_out.iota_n       ;
results.iota_rss    = path_out.iota_out(:,1);
results.v_n         = path_out.v_n          ; 
results.v_rss       = path_out.v_out(:,1)   ; 
results.r_n         = path_out.r_n          ; 
results.c_rss       = path_out.c_out(:,1)   ; 
results.c_n         = path_out.c_n          ;

%% Step 4: Plotting
plotting

