function [rss]=solve_rss(psi_rss,v_0_mat,parameters);
% Numerical parameters
parameter_unpack;
time_preallocate;

% solve solution
rho       = parameters.rho;
lambda_bar= parameters.lambda_bar;
delta     = parameters.delta;
y_ss      = parameters.y_ss;

%% STEADY STATE
v_rss    = solve_HJB_rss(rho,v_0_mat,parameters)   ; % HJB: Solution at Steady State
iota_rss = 1/lambda_bar * (1-v_rss./psi_rss); % Optimal policies
f_rss     = flip(cumsum(flip(iota_rss)))*dt      ; % KFE
c_rss = y_ss  - f_rss(1) + sum((psi_rss.*(1-1/2*lambda_bar.*iota_rss).*iota_rss)-delta*f_rss)*dt; % Consumption
r_rss = rho;

if c_rss<=0 % The case lambda <lambda_min
    warning('Negative Consumption at RSS');
    
%     r_ss = (r_bar + rho) /2; % First guess
%     for iter = 1:max_iter
%         v_ss     = -delta/r_ss *(1-exp(-r_ss*tau'))-exp(-r_ss*tau');
%         iota_ss  = 1./lambda_bar .*  (1-v_ss./psi_ss);% Optimal policies
%         f_rss     = flip(cumsum(flip(iota_ss)))*dt;% KFE
%         c_ss = y_ss  - f_rss(1) + sum((psi_ss - 1/2 *lambda_bar .* iota_ss).*iota_ss-delta*f_rss) * dt;% Consumption
%         if  abs(c_ss) < tol_path
%             break  % Consumption should be zero
%         else
%             r_ss = r_ss + relax_path*c_ss; 
%         end
%     end
%     if iter == max_iter
%         disp('Error: maximum number of iterations reached at steady state loop')
%     end
end
rss.c_rss=c_rss;
rss.iota_rss=iota_rss;
rss.v_rss=v_rss;
rss.f_rss=f_rss;
rss.r_rss=r_rss;