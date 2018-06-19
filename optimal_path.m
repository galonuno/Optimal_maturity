function results = optimal_path(parameters,f0)
% This version codes Galo's solution using standard finite difference methods.
% Transition path, no shocks

%% Model parameters
parameter_unpack;

% psi
psi_n= parameters.psi;
psi_ss  = psi_n(:,end);

%% Variables (preallocation)
t         = 0:dt:t_max;    % calendar time
N         = length(t);     % Number of time nodes
tau       = dt:dt:tau_max; % maturities
I         = length(tau);   % Number of tau nodes
u         = zeros (I,1);
u(:)      = -delta;
u(1)      = u(1) -1/dt;
A0        = -speye(I);
aa        = ones(I-1,1);
A1        = spdiags(aa,-1,I,I);
A         = 1/dt*(A0+A1);
V_n       = zeros(I,N);   % Value function
iota_n    = zeros(I,N);   % Policy function
f_n       = zeros(I,N);   % Density
c_n       = zeros(1,N);   % Consumption
r_n_u     = zeros(1,N);   % Interest rate (update)

%Path for Income
y_bar     = parameters.y_bar ;  %Before the Shock
y_bar_0   = parameters.y_bar_0;
rho_y     = parameters.rho_y;
y_n       = zeros(1,N);
for i=2:N
y_n(1) = y_bar_0;
y_n(i) = y_bar *(1-rho_y)+rho_y*y_n(i-1);
end

%% Initial guess of interest rates
r_n      = rho * ones(1,N);


%% STEADY STATE
V    = -delta/rho *(1-exp(-rho*tau'))-exp(-rho*tau');  % HJB: Solution at Steady State
iota  = 1/lambda_bar * (V+ psi_ss );% Optimal policies
f     = flip(cumsum(flip(iota)))*dt;% KFE
c = y_bar  - f(1) + sum((psi_ss - 1/2 *lambda_bar * iota).*iota-delta*f) * dt; % Consumption
r = rho;

if c<=0 % The case lambda <lambda_min
    r = (r_bar + rho) /2; % First guess
    for iter = 1:max_iter
        V     = -delta/r *(1-exp(-r*tau'))-exp(-r*tau');
        iota  = 1./lambda_bar .* (V+ psi_ss );% Optimal policies
        f     = flip(cumsum(flip(iota)))*dt;% KFE
        c = y_bar  - f(1) + sum((psi_ss - 1/2 *lambda_bar .* iota).*iota-delta*f) * dt;% Consumption
        if  abs(c) < tol_path
            break  % Consumption should be zero
        else
            r = r + relax_path*c; 
        end
    end
    if iter == max_iter
        disp('Error: maximum number of iterations reached at steady state loop')
    end
end

%% DYNAMICS
%%----------------------------
% Initial guess of interest rates
r_n      = r * ones(1,N);
tic
for iter = 1:max_iter
    
    % HJB: Dynamic solution
    V_n(:,N)   = V; % Terminal value
    for n = N-1:-1:1
        B = (1/dt + r_n(n))*speye(I) -A;
        d = u + V_n(:,n+1) /dt;
        V_aux = B\d;
        V_n(:,n) = V_aux;
    end
    
    % Optimal policies
    iota_n = 1/lambda_bar * (V_n + psi_n );
    
    % KFE
    f_n(:,1)   = f0;
    D         = speye(I) - dt *A';
    for n=2:N
        h      = iota_n(:,n) *dt + f_n(:,n-1);
        f_aux  = D\h;
        f_n(:,n) = f_aux;
    end
    
    % Consumption and interest rates
    c_n = y_n  - f_n(1,:) + sum((psi_n - 1/2 *lambda_bar * iota_n).*iota_n-delta*f_n) * dt;
    if c>tol_path
        r_n_u(1:N-1) = rho * ones(1,N-1)  + gamma/dt * (c_n(2:N) - c_n(1:N-1)) ./ c_n(1:N-1);
        r_n_u(N)     = rho;
    else
        r_n_u(1:N-2)   = rho * ones(1,N-2)  + gamma/dt * (c_n(2:N-1) - c_n(1:N-2)) ./ c_n(1:N-2);
        r_n_u(N-1)     = r;
        r_n_u(N)       = r;
    end
    
    % Interest rate update
    error = max(abs (r_n - r_n_u));
    if  error < tol_path
        break
    else
        r_n = relax_path * r_n_u + (1-relax_path) * r_n;
        if mod(iter,500) == 0
                   disp(iter)
                   disp(error)
        end
    end
end
if iter == max_iter
    disp('Error: maximum number of iterations reached')
end

toc
results.V_n       = V_n;   % Value function
results.iota_n    = iota_n;% Policy function
results.f_n       = f_n;   % Density
results.c_n       = c_n;   % Consumption
results.r_n       = r_n;   % Interest rate (update)
results.y_n       = y_n;
results.psi       = psi_n;
results.dt        = dt;
