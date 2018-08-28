function results = risky_SS(parameters)
% This version codes Galo's solution using standard finite difference methods.
% Transition path, no shocks

%% Model parameters
% Model is annual

% Tolerance
tol      = parameters.tol;
max_iter = parameters.max_iter;  % Maximum number of iterations
relax    = parameters.relax;     % Relaxation coefficient

% Numerical parameters
t_max     = parameters.t_max;    % Maximum number of years (time)
tau_max   = parameters.tau_max;  % Maximum number of years (maturity)
dt        = parameters.dt;       % Monthly steps

% Preference Parameters
rho     = parameters.rho;       % Discount factor - Model is Quarterly
delta   = parameters.delta;     % Coupon portion
gamma   = parameters.gamma;     % risk aversion

% Frictions parameters
lambda_bar = parameters.lambda_bar;

%Yield curve 
psi     = bond_prices(parameters); %After the Shock
psi_rss = bond_prices_rss(parameters); %Before thee Shock
r_bar   = parameters.r_bar ;    % Market discount factor

% Probability of the Shock
phi     = parameters.phi;          % Poisson arrival rate

%% Variables (preallocation)
t         = 0:dt:t_max;    % calendar time
N         = length(t);     % Number of time nodes
tau       = dt:dt:tau_max; % maturities
I         = length(tau);   % Number of tau nodes
u         = zeros (I,1);
A0        = -speye(I);
aa        = ones(I-1,1);
A1        = spdiags(aa,-1,I,I);
A         = 1/dt*(A0+A1);
%%
y_bar     = parameters.y_bar ;  %Before the Shock
y_bar_0   = parameters.y_bar_0;
rho_y     = parameters.rho_y;
y_n       = zeros(1,N);
for i=2:N
    y_n(1)=y_bar_0;
    y_n(i)=y_bar*(1-rho_y)+rho_y*y_n(i-1);
end
parameters.y_n=y_n;

%% Initial guess of the distribution
tau       = parameters.dt:parameters.dt:parameters.tau_max; % maturities
I         = length(tau);   % Number of tau nodes
f0 = zeros(I,1) /dt /I;     %Unit mass uniform - initial density
% results = optimal_path(parameters,f0);
% f0      = results.f_n(:,end);
fRSS     = f0;
% c_rss      = results.c_n(1,1);
c_rss      = y_bar_0;
%% ----------------------------
%MAIN LOOP
tic
for iter = 1:max_iter
    %% Value if the shock arrives
    results          = optimal_path(parameters,fRSS);
    parameters.y_bar = y_bar;
    v_hat         = results.V_n(:,1);
    c_after      = results.c_n(1,1);
    %% HJB: Solution at Steady State, before the shock Arives
    Brss         = ((rho+phi) *speye(I) - A);
    u(:)        = -delta + phi*v_hat(:) *(c_rss/c_after)^(gamma);
    u(1)        = u(1) -1/dt;
    Vrss         = Brss\u;                   %
       
    %% Optimal policies
    iota_rss = 1/lambda_bar * (Vrss + psi_rss );
    %% KFE
      
    fRSS_new  = flipud(cumsum(flipud(iota_rss)))*results.dt; 
 
    %% Interest rate update
    error = max(abs (fRSS_new - fRSS));
    if  error < tol
        break
    else
        fRSS = relax * fRSS_new + (1-relax) * fRSS;
                  disp('MAIN ITERATION ALGORITHM')
                  disp(iter)
                  disp(error)
    end
    
       %% Consumption and interest rates
    c_rss =y_bar  - fRSS(1) + sum((psi_rss - 1/2 *lambda_bar * iota_rss).*iota_rss-delta*fRSS) * dt;
     
end
if iter == max_iter
    disp('Error: maximum number of iterations reached')
end

toc
results.V_rss       = Vrss;    % Value function
results.iota_rss   = iota_rss;% Policy function
results.f_rss       = fRSS;    % Density
results.c_rss      = c_rss;   % Consumption
results.c_after    = c_after; % Consumption
results.psi_rss    = psi_rss;
results.r_rss      = rho;

%end;
%end