%% Unpack Parameters
% internal routine to unpack parameters

%% Code Parameters
% code parameters 
tol        = parameters.tol;
tol_path   = parameters.tol_path;
max_iter   = parameters.max_iter;  % Maximum number of iterations
relax      = parameters.relax;     % Relaxation coefficient
relax_path = parameters.relax_path;

% vectorization parameters
t_max     = parameters.t_max;    % Maximum number of years (time)
tau_max   = parameters.tau_max;  % Maximum number of years (maturity)
dt        = parameters.dt;       % Monthly steps

% %% Model Parameters
% % preference Parameters
% rho     = parameters.rho;       % Discount factor - Model is Quarterly
% delta   = parameters.delta;     % Coupon portion
% gamma   = parameters.gamma;     % risk aversion
% 
% % frictions parameters
% lambda_bar = parameters.lambda_bar;