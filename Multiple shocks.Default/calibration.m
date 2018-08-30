% Code to compute the calibration of the model
%(c) Bigio, Nuno & Passadore model

%%

% function error = calibration(guess)
% call fminsearch(@(x) calibration(x), [10,0.042])
% guess_lambda_bar = guess(1);
% guess_rho        = guess(2);
guess_lambda_bar = 7.0240;    
guess_rho        = 0.0416;

%% Calibration parameters
% guess_lambda_bar   = 10;
% guess_rho          = 0.042;
Target_debt_to_GDP   = 0.46;  
Rate_depleted        = 0.6;  % Amount sold per period 
eta  = 0.015;           % Bank spread

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

%% Model parameters
% Time Period annual
% Preference Parameters
parameters.gamma   = 2.00     ; % risk aversion
parameters.rho     = guess_rho   ; % discount factor - Model is Quarterly

% Output parameters
parameters.y_ss         = 1.0    ; % output at steady state

%Interest Rate Shock Parameters
parameters.r_bar_ss           = 0.04                ; % steady-state short rate
parameters.delta              = parameters.r_bar_ss ; % coupon = r_ss so we have bonds at par

% adjustment cost
parameters.lambda_bar = guess_lambda_bar  ; %

%% Pre-Allocation
time_preallocate;

%% Step 1: Bond Prices in Steady State
psi_ss         = solve_HJB_ss(parameters.r_bar_ss,parameters)           ;

%% Step 2: Compute the Steady State solution
steady = solve_steady(psi_ss,parameters);

%% Step 3: Compute targets
Total_debt_to_GDP = sum(steady.f_ss*dt)/parameters.y_ss; % Debt to GDP
error1            = abs(Target_debt_to_GDP/Total_debt_to_GDP -1);
Total_issuances   = sum(steady.iota_ss*dt);
mu                = Rate_depleted * Total_issuances / parameters.tau_max;
lambda_bar        = eta/2 / mu / parameters.y_ss;
error2            = abs(lambda_bar/guess_lambda_bar -1);
error             = error1 + error2;
