function [steady]=solve_steady(psi_ss,parameters)
% Numerical parameters
parameter_unpack;
time_preallocate;

% solve solution
rho       = parameters.rho;
lambda_bar= parameters.lambda_bar;
delta     = parameters.delta;
y_ss      = parameters.y_ss ;
r_bar_ss  = parameters.r_bar_ss ;

%% STEADY STATE
% v_ss    = solve_HJB_ss(rho,parameters)   ; % HJB: Solution at Steady State
v_ss      = delta/rho*(1-exp(-rho*tau'))+exp(-rho*tau');
iota_ss  = 1/lambda_bar*(1-v_ss./psi_ss); % Optimal policies
f_ss     = flip(cumsum(flip(iota_ss)))*dt      ; % KFE
c_ss     = y_ss-f_ss(1)+sum(psi_ss.*(1-1/2*lambda_bar.*iota_ss).*iota_ss-delta*f_ss)*dt; % Consumption
r_ss     = rho;
steady.type = 'ss';
if c_ss<=0 % The case lambda <lambda_min
    r_ss = (r_bar_ss + rho) /2; % First guess
    for iter = 1:max_iter
        v_ss     = delta/r_ss*(1-exp(-r_ss*tau'))+exp(-r_ss*tau');
        iota_ss  = 1./lambda_bar.*(1-v_ss./psi_ss); % Optimal policies
        f_ss     = flip(cumsum(flip(iota_ss)))*dt;  % KFE
        c_ss = y_ss  - f_ss(1) + sum(psi_ss.*(1-1/2 *lambda_bar.*iota_ss).*iota_ss-delta*f_ss) * dt;% Consumption
        if  abs(c_ss) < tol_path
            steady.type='asymptotic';
            c_ss=0                  ;
            disp('Message: asymptotic steady state found');
            disp(['rh0=' num2str(rho) ' whereas r_ss=' num2str(r_ss)]);
            break  % Consumption should be zero
        else
            r_ss = r_ss + relax_path*c_ss; 
        end
    end
    % Break 
    if iter == max_iter
        disp('Error: maximum number of iterations reached trying to find asymptotic value');
    end
else
    disp('Message: found a proper steady state');
    disp(['C_ss=' num2str(c_ss)]);
end
steady.c_ss=c_ss        ;
steady.iota_ss=iota_ss  ;
steady.v_ss=v_ss        ;
steady.f_ss=f_ss        ;
steady.r_ss=r_ss        ;