%% Construct Exogenous Paths
% exo_path_construct
dt            = parameters.dt;
r_bar_ss      = parameters.r_bar_ss;
r_bar_0       = parameters.r_bar_0(1:N_r,1); % initial value
rho_r_bar     = parameters.rho_r_bar(1:N_r,1); % persistence (in months)
r_bar_n       = zeros(N_r,N_t)  ;
r_bar_n(:,1)  = r_bar_0         ;
for tt=2:N_t   
    r_bar_n(:,tt)= r_bar_ss*(rho_r_bar(:)*dt) + (1-rho_r_bar(:)*dt).*r_bar_n(:,tt-1); 
end

% Path for output
y_ss     = parameters.y_ss       ;
y_0      = parameters.y_0(1:N_y,1)     ; % initial value
rho_y     = parameters.rho_y(1:N_y,1) ; % persistence (in months)
y_n       = zeros(N_y,N_t)      ;
y_n(:,1)  = y_0                 ;
for tt=2:N_t   
    y_n(:,tt)= (y_ss*(rho_y(:)*dt)+ (1-rho_y(:)*dt).*y_n(:,tt-1));
end