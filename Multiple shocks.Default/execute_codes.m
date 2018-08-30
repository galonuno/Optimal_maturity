% Code to compute the Deterministic and Risky Steady State of the model
%(c) Bigio, Nuno & Passadore model
% This version codes Galo's solution using standard finite difference methods.
% The Code Solves for the Risky Steady State with default
plotit=1;

% Price of Risk
parameters.U         =  @(c) (c.^(1-parameters.gamma)-1)/(1-parameters.gamma);
parameters.U_p_ratio =  @(c_a,c_b) (c_a/c_b).^(-parameters.gamma);
parameters.dU        =  @(c) (c.^(-parameters.gamma));

%% Pre-Allocation
time_preallocate;

%% Paths
% Path for the Short Rate and output
exo_path_construct;

% Autarky Values
autarkyV_construct;
parameters.V_a = kron(ones(N_r,1),V_a)     ;

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
results.psi_n = psi_n;
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

