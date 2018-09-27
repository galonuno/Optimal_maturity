% Replication code for Optima Debt Maturity Management
%(c) Bigio, Nuno & Passadore model
close all;  clear;
mkdir figures % Creates folder to save figures

%% 0 PARAMETERS
% Code parameters
% Tolerance
parameters.tol_path   = 1e-3;  % tolerance parameter optimal path
parameters.tol        = 1e-5;  % tolerance parameter
parameters.max_iter   = 5000;  % maximum number of iterations
parameters.relax      = 0.1 ;  % relaxation coefficient
parameters.relax_path = 0.01;  % relaxation coefficient optimal path

% Numerical parameters
parameters.t_max     = 100    ; % maximum number of years (time)
parameters.tau_max   = 20    ; % maximum number of years (maturity)
parameters.dt        = 1/12  ; % monthly steps

% Method to Solve the model
parameters.solve_method='iter'; % iter for iterative, solver for solver

% Model parameters
% Time Period annual
% Preference Parameters
parameters.gamma        = 2.00     ; % risk aversion
parameters.rho          = 0.0416   ; % discount factor - Model is Quarterly

% adjustment cost
parameters.lambda_bar   = 7.0240  ; %

% Output 
parameters.y_ss         = 1.0 ; % output at steady state
 
%Interest Rate 
parameters.r_bar_ss     = 0.04 ; % steady-state short rate
parameters.delta        = parameters.r_bar_ss ; % coupon = r_ss so we have bonds at par

% Number of Shocks
N_r=   1;
N_y=   1;
parameters.N_r=N_r;
parameters.N_y=N_y;

% Update to a matrix of joint probabilities - discrete matrix of marginals
% - sum of row and columns must add up to 1
parameters.prob_mat=1; % [s*]

% Output parameters
parameters.y_0(1:N_y,1)       = 1;                 % [s*] initial value, output after a shock
parameters.rho_y(1:N_y,1)     = (1-0.95)*4;        % [s*] persistence output after shock - 0.95 in quarterly data
parameters.r_bar_0(1:N_r,1)   = 0.04           ;   % [s*] initial value volatility of 0.048 quarterly
parameters.rho_r_bar(1:N_r,1) = (1-0.95)*4      ;  % [s*] persistence (in years)


%% 1 PERFECT FORESIGHT 
% parameters.mode = 'PF'; % Perfect Foresight
% 
% %Risky Steady State Parameters - No RSS
% parameters.phi  = 0.00       ; % Poisson arrival rate
% 
% % Default parameters           - No default
% parameters.defmodel ='logit'     ; % default model
% parameters.muprobit =10;
% parameters.penalty  = 1;        % Autarky penalty

% 1.1 Output shock
% parameters.y_0(1:N_y,1)   = 0.95; % [s*] initial value, output after a shock
% execute_codes % % 
% close all; 
% cd figures; 
% save PF_Y; %figure(3); saveas(gcf,'PF_Y','pdf');saveas(gcf,'PF_Y','fig');close all;  
% cd ..
% parameters.y_0(1:N_y,1)   = 1; % [s*] initial value, output after a shock

% % 1.2 Interest rate shock
% parameters.r_bar_0(1:N_r,1)   = 0.05           ; % [s*] initial value
% execute_codes;
% close all;
% cd figures; 
% save PF_R; %figure(2);   orient landscape; saveas(gcf,'DSS','pdf'); saveas(gcf,'DSS','fig'); figure(3);  saveas(gcf,'PF_IR','pdf'); saveas(gcf,'PF_IR','fig');close all;
% cd ..
% parameters.r_bar_0(1:N_r,1)   = 0.04        ; % [s*] initial value 

% 1.3 Interest rate shock IES = 0
% parameters.r_bar_0(1:N_r,1)   = 0.05           ; % [s*] initial value
% parameters.gamma        = 0.00     ; % risk aversion
% execute_codes;
% close all;
% cd figures; 
% save PF_R_IES_0; %figure(2);   orient landscape; saveas(gcf,'DSS','pdf'); saveas(gcf,'DSS','fig'); figure(3);  saveas(gcf,'PF_IR','pdf'); saveas(gcf,'PF_IR','fig');close all;
% cd ..
% parameters.r_bar_0(1:N_r,1)   = 0.04        ; % [s*] initial value 
% parameters.gamma        = 2.00     ; % risk aversion

%% 2 RSS
parameters.mode = 'RS'; % RSS

%Risky Steady State Parameters - No RSS
% parameters.phi  = 0.02          ; % Poisson arrival rate
parameters.phi  = 0.05          ; % Poisson arrival rate

% Default parameters           - No default
parameters.defmodel ='logit'     ; % default model
parameters.muprobit =10;
parameters.penalty  = 1;        % Autarky penalty

% %2.1 Output shock
% parameters.y_0(1:N_y,1)   = 0.95; % [s*] initial value, output after a shock
% execute_codes % % 
% cd figures; save RSS_Y; cd ..
% close all; 
% parameters.y_0(1:N_y,1)   = 1; % [s*] initial value, output after a shock
% 
% 
% % 2.2 Interest rate shock
% parameters.r_bar_0(1:N_r,1)   = 0.05           ; % [s*] initial value volatility of 0.048 quarterly
% execute_codes;
% cd figures; save RSS_R; cd ..
% close all; 
% parameters.r_bar_0(1:N_r,1)   = 0.04        ; % [s*] initial value 

% 2.3 Interest rate shock IES = 0
parameters.gamma        = 0.00     ; % risk aversion

parameters.r_bar_0(1:N_r,1)   = 0.05           ; % [s*] initial value volatility of 0.048 quarterly
execute_codes;
% cd figures; save RSS_R_IES_0; cd ..
close all; 
parameters.r_bar_0(1:N_r,1)   = 0.04        ; % [s*] initial value 
parameters.gamma        = 2.00     ; % risk aversion

% %% 3 Default
% parameters.mode = 'DF'; % RSS with default
% 
% %Risky Steady State Parameters - No RSS
% parameters.phi  = 0.02          ; % Poisson arrival rate
% 
% % Default parameters           - No default
% parameters.defmodel ='logit'     ; % default model
% parameters.muprobit = 100;
% parameters.penalty  = 0.015;        % Autarky penalty
% 
% % %3.1 Output shock
% parameters.y_0(1:N_y,1)   = 0.95; % [s*] initial value, output after a shock
% execute_codes % % 
% cd figures; save DEF; cd ..
% close all; 
% parameters.y_0(1:N_y,1)   = 1; % [s*] initial value, output after a shock
% 
% % %3.2 Output shock IES = 0
% parameters.gamma          = 0.00     ; % risk aversion
% parameters.y_0(1:N_y,1)   = 0.95; % [s*] initial value, output after a shock
% execute_codes % % 
% cd figures; save DEF_IES_0; cd ..
% close all; 
% parameters.y_0(1:N_y,1)   = 1; % [s*] initial value, output after a shock
% parameters.gamma         = 2.00     ; % risk aversion

%% PLOTTING
% plotting_RSS
 plotting
% plotting_DEF
