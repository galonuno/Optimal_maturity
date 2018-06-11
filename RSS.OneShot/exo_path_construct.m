%% Construct Exogenous Paths
% exo_path_construct
dt       = parameters.dt;
r_bar_ss = parameters.r_bar_ss;
r_bar_0  = parameters.r_bar_0(1:N_r,1); % initial value
rho_r_bar= parameters.rho_r_bar(1:N_r,1); % persistence (in months)
r_bar_n       = zeros(N_r,N_t)  ;
r_bar_n(N_r,1)= r_bar_0         ;
for ii=2:N_t   
    % r_bar_n(:,ii)= (r_bar_ss(:)*(1-rho_r_bar(:)*dt)+ rho_r_bar(:)*dt*r_bar_n(:,ii-1)); 
    r_bar_n(:,ii)= (r_bar_ss(:)*(1-rho_r_bar(:))+ rho_r_bar(:)*r_bar_n(:,ii-1)); 
end

% Path for output
y_ss    = parameters.y_ss       ;
y_0=parameters.y_0(1:N_y,1)     ; % initial value
rho_y=parameters.rho_y(1:N_y,1) ; % persistence (in months)
y_n       = zeros(N_y,N_t)      ;
y_n(N_y,1)= y_0                 ;
for ii=2:N_t   
    y_n(:,ii)= (y_ss(:)*(1-rho_y(:)*dt)+ rho_y(:)*dt*y_n(:,ii-1)); 
    y_n(:,ii)= (y_ss(:)*(1-rho_y(:))+ rho_y(:)*y_n(:,ii-1));
end