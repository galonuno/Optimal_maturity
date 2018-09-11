% Code to compute the Risky Steady State of the model with default
%(c) Bigio, Nuno & Passadore model
% This version codes Galo's solution using standard finite difference methods.
% The Code Solves for the Risky Steady State with default
% close all; % clear;
plotit=0;

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
N_r=   1;%2;
N_y=   1;%3;
parameters.N_r=N_r;
parameters.N_y=N_y;

%% Model parameters
% Time Period annual
% Preference Parameters
parameters.gamma   = 2.00     ; % risk aversion
parameters.rho     = 0.0416   ; % discount factor - Model is Quarterly

% Output parameters
parameters.y_ss         = 1.0    ; % output at steady state
parameters.y_0(1:N_y,1)   = 0.95;%[0.99; 1.0; 0.95 ]; % [s*] initial value, output after a shock
parameters.rho_y(1:N_y,1) = 0.22;%[0.25; 0.45; 0.22]; % [s*] persistence output after shock - half life 1/rho months
 
%Interest Rate Shock Parameters
parameters.r_bar_ss           = 0.04                ; % steady-state short rate
parameters.delta              = parameters.r_bar_ss ; % coupon = r_ss so we have bonds at par
parameters.r_bar_0(1:N_r,1)   = 0.08;%[0.06; 0.08]        ; % [s*] initial value [0.042 ; 0.043] 
parameters.rho_r_bar(1:N_r,1) = 0.05;%[0.25; 0.35]        ; % [s*] persistence (in years)

%Risky Steady State Parameters
parameters.phi  = 0.05         ; % Poisson arrival rate

% adjustment cost
parameters.lambda_bar = 7.0240  ; %

% Default parameters - Uniform distribution
parameters.defmodel='logit'     ; % default model
parameters.muprobit=10          ; % >> Why?
parameters.penalty  = 0.99;        % Autarky penalty

% Update to a matrix of joint probabilities - discrete matrix of marginals
% - sum of row and columns must add up to 1
parameters.prob_mat=zeros(N_r,N_y); % [s*]
parameters.prob_mat=1;%[0.1 0.25 0.25; 0.3 0.05 0.05]; % [s*] % joint probability distribution, correlation is built in
if sum(sum(parameters.prob_mat))~=1;
    error('Prob matrix won''t add to one');
end

% Price of Risk
parameters.U         =  @(c) (c.^(1-parameters.gamma)-1)/(1-parameters.gamma);
parameters.U_p_ratio =  @(c_a,c_b) (c_a/c_b).^(-parameters.gamma);
parameters.dU        =  @(c) (c.^(-parameters.gamma));

%% Pre-Allocation
time_preallocate;

%% Paths
phi_vec    = [(0.0001:0.001:parameters.phi) (parameters.phi+0.001:0.001:parameters.phi+0.01)];
f_rss_mat  = zeros(N_tau,length(phi_vec));
iota0_rss_mat  = zeros(N_tau,length(phi_vec));
iota1_rss_mat  = zeros(N_tau,length(phi_vec));
for nn=1:length(phi_vec);
    % pickup value
    parameters.phi      = phi_vec(nn);
    
    % Path for the Short Rate and output
    exo_path_construct;

    % Autarky Values
    autarkyV_construct;
    parameters.V_a=kron(ones(N_r,1),V_a)     ;

    %% Save Exogenous Paths
    paths.r_bar_n = r_bar_n;
    paths.y_n     = y_n    ;

    %% Step 1: Bond Prices in Steady State - index of r shock in columns
    psi_n=zeros(N_tau,N_t,N_r); % [s*]
    psi_0_mat=zeros(N_tau,N_r); % [s*]
    for rr=1:N_r
        psi_n(:,:,rr)          = solve_HJB_path(r_bar_n(rr,:),r_bar_ss,parameters) ; % [s*]
        psi_0_mat(:,rr)        = psi_n(:,1,rr)                                     ; % [s*]
    end
    psi_ss         = solve_HJB_ss(r_bar_ss,parameters)                             ; % [s*]
    rss.psi_0_mat  = psi_0_mat                                                     ; % [s*]
    % results
    results.psi_ss  = psi_ss                                                       ; % [s*]

    %% Step 2: Value functions in Steady State 
    steady=solve_steady(psi_ss,parameters)      ;
    % Save Steady States
    steady.psi_ss  = psi_ss                     ;
    paths.psi_n    = psi_n                      ;

    %% Step 2: Risky-Steady State
    path_out        = solve_rss_path(rss,steady,paths,parameters); % [s*<-]
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
    
    %% Step 4: Save results for Compstat
    f_rss_mat(:,nn)   = path_out.f_n(:,1);
    iota0_rss_mat(:,nn)= path_out.iota_n(:,1);
    iota1_rss_mat(:,nn)= path_out.iota_n(:,2);
end
% close all;
close all;

figure
surf([phi_vec],tau,[f_rss_mat],'FaceAlpha',0.5,'EdgeColor','None'); hold on;
surf([phi_vec],tau,[iota0_rss_mat]*0,'FaceAlpha',0.1,'EdgeColor','w'); 
title('$f_{rss}$','interpreter','latex');
xlabel('$\bar{\lambda}$','interpreter','latex');
ylabel('$\tau$','interpreter','latex'); axis tight;

figure
surf([phi_vec],tau,[iota0_rss_mat],'FaceAlpha',0.5,'EdgeColor','None'); hold on;
surf([phi_vec],tau,[iota0_rss_mat]*0,'FaceAlpha',0.1,'EdgeColor','w'); hold on;
% surf([0 phi_vec],tau,[path_out.iota_n(:,end) iota1_rss_mat],'FaceAlpha',0.5,'EdgeColor','None'); hold on;
title('$\iota_{rss}$','interpreter','latex');
xlabel('$\bar{\lambda}$','interpreter','latex');
ylabel('$\tau$','interpreter','latex'); axis tight;

figure
surf([phi_vec],tau,[iota0_rss_mat],'FaceAlpha',0.5,'EdgeColor','None'); hold on;
surf([phi_vec],tau,[iota1_rss_mat],'FaceAlpha',0.5,'EdgeColor','None'); hold on;
surf([phi_vec],tau,[iota0_rss_mat]*0,'FaceAlpha',0.1,'EdgeColor','w'); hold on;
% surf([0 phi_vec],tau,[path_out.iota_n(:,end) iota1_rss_mat],'FaceAlpha',0.5,'EdgeColor','None'); hold on;
title('$\iota(\tau)_{rss}$ vs. $\iota(\tau,1)$','interpreter','latex');
xlabel('$\bar{\lambda}$','interpreter','latex');
ylabel('$\tau$','interpreter','latex'); axis tight; 
