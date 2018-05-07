% OPTIMAL DEBT DYNAMICS
%(c) Bigio, Nuno & Passadore model
clear; %close all;
clc

% This version codes Galo's solution using standard finite difference methods.
% Transition path, no shocks

%% Model parameters
% Model is annual

% Tolerance
tol      = 2e-4;
tolSS    = 2e-4;
max_iter = 100000;  % Maximum number of iterations
relax    = 0.0001; % Relaxation coefficient
relaxSS  = 0.0001; % Relaxation coefficient SS
% Numerical parameters
t_max     = 300     ; % Maximum number of years (time)
tau_max   = 20      ; % Maximum number of years (maturity)
dt        = 1/12    ; % Monthly steps

% Preference Parameters
rho     = 0.06;  % Discount factor - Model is Quarterly
delta   = 0.03;  % Coupon portion
y_bar   = 1   ;  % output
gamma   = 2   ;  % risk aversion

% Frictions parameters
lambda_bar = 1;

%Yield curve
r_bar   = 0.04 ; % Market discount factor

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
y_n       = y_bar * ones(1,N);

% Initial distribution
% f0 = ones(I,1) /dt /I;  %Unit mass uniform - initial density
  load g0
  f0 = g0; 

 % Bond prices
 psi    = delta/r_bar *(1-exp(-r_bar*tau')) + exp(-r_bar*tau') ;% Bond prices
 
%% STEADY STATE
V    = -delta/rho *(1-exp(-rho*tau'))-exp(-rho*tau');  % HJB: Solution at Steady State
iota  = 1/lambda_bar * (V+ psi );% Optimal policies
f     = flip(cumsum(flip(iota)))*dt;% KFE
c = y_bar  - f(1) + sum((psi - 1/2 *lambda_bar * iota).*iota-delta*f) * dt; % Consumption
r = rho;

if c<=0 % The case lambda <lambda_min
    r = (r_bar + rho) /2; % First guess
    for iter = 1:max_iter
        V     = -delta/r *(1-exp(-r*tau'))-exp(-r*tau');
        iota  = 1./lambda_bar .* (V+ psi );% Optimal policies
        f     = flip(cumsum(flip(iota)))*dt;% KFE
        c = y_bar  - f(1) + sum((psi - 1/2 *lambda_bar .* iota).*iota-delta*f) * dt;% Consumption
        if  abs(c) < tolSS
            break  % Consumption should be zero
        else
            r = r + relaxSS*c; 
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
    iota_n = 1/lambda_bar * (V_n + psi );
    
    % KFE
    f_n(:,1)   = f0;
    D         = speye(I) - dt *A';
    for n=2:N
        h      = iota_n(:,n) *dt + f_n(:,n-1);
        f_aux  = D\h;
        f_n(:,n) = f_aux;
    end
    
    % Consumption and interest rates
    c_n = y_n  - f_n(1,:) + sum((psi - 1/2 *lambda_bar * iota_n).*iota_n-delta*f_n) * dt;
    if c>tolSS
        r_n_u(1:N-1) = rho * ones(1,N-1)  + gamma/dt * (c_n(2:N) - c_n(1:N-1)) ./ c_n(1:N-1);
        r_n_u(N)     = rho;
    else
        r_n_u(1:N-2)   = rho * ones(1,N-2)  + gamma/dt * (c_n(2:N-1) - c_n(1:N-2)) ./ c_n(1:N-2);
        r_n_u(N-1)     = r;
        r_n_u(N)       = r;
    end
    
    % Interest rate update
    error = max(abs (r_n - r_n_u));
    if  error < tol
        break
    else
        r_n = relax * r_n_u + (1-relax) * r_n;
                  disp(iter)
                  disp(error)
    end
end
if iter == max_iter
    disp('Error: maximum number of iterations reached')
end

% %% SAVING THE TERMINAL DENSITY
  g0 = f_n(:,end);
  save g0 g0
%% Plotting
toc
%Steady state
figure(1)
subplot(221)
plot(tau,f_n(:,end),'Linewidth',2)
xlabel('Maturity, $\tau$','FontSize',14,'interpreter','latex')
title('Debt $f(\tau)$','FontSize',14,'interpreter','latex')
hold on

subplot(222)
plot(tau,iota_n(:,end),'Linewidth',2)
xlabel('Maturity, $\tau$','FontSize',14,'interpreter','latex')
title('Issuances $\iota(\tau)$','FontSize',14,'interpreter','latex')
hold on

subplot(223)
plot(tau,psi(:,end),'Linewidth',2)
xlabel('Maturity, $\tau$','FontSize',14,'interpreter','latex')
title('Bond prices $\psi(\tau)$','FontSize',14,'interpreter','latex')
hold on

subplot(224)
plot(tau,V_n(:,end),'Linewidth',2)
xlabel('Maturity, $\tau$','FontSize',14,'interpreter','latex')
title('Value function $v(\tau)$','FontSize',14,'interpreter','latex')
hold on
%%
% Transitional dynamics
figure(2)
subplot(121)
mesh(t,tau,f_n)
xlabel('Time, $t$','FontSize',14,'interpreter','latex')
ylabel('Maturity, $\tau$','FontSize',14,'interpreter','latex')
title('Debt $f(t,\tau)$','FontSize',14,'interpreter','latex')

subplot(122)
mesh(t,tau,iota_n)
xlabel('Time, $t$','FontSize',14,'interpreter','latex')
ylabel('Maturity, $\tau$','FontSize',14,'interpreter','latex')
title('Issuances $\iota(t,\tau)$','FontSize',14,'interpreter','latex')

figure(3)
subplot(221)
plot(t,r_n*100,'Linewidth',2)
xlabel('Time, $t$','FontSize',14,'interpreter','latex')
ylabel('$\%$','FontSize',14,'interpreter','latex')
title('Interest rate $r(t)$','FontSize',14,'interpreter','latex')
hold on
xlim([0 300])

subplot(222)
plot(t,c_n,'Linewidth',2)
xlabel('Time, $t$','FontSize',14,'interpreter','latex')
title('Consumption $c(t)$','FontSize',14,'interpreter','latex')
hold on
xlim([0 300])

subplot(223)
i_1 = sum(iota_n(1:1*12,:))*dt;
i_5 = sum(iota_n((1*12+1):5*12,:))*dt;
i_10 = sum(iota_n((5*12+1):10*12,:))*dt;
i_20 = sum(iota_n((10*12+1):end,:))*dt;
plot(t,100*[(i_1'-i_1(1)) (i_5'-i_5(1)) (i_10'-i_10(1)) (i_20'-i_20(1))],'Linewidth',2)
legend({'0-1 year','1-5 years','5-10 years','10-20 years'},'Interpreter','Latex','FontSize',14)
xlabel('Time, $t$','FontSize',14,'interpreter','latex')
ylabel('Deviations from steady state (in $\%$)','FontSize',14,'interpreter','latex')
title('Issuances $\iota(t)$','FontSize',14,'interpreter','latex')
hold on
xlim([0 300])

subplot(224)
f_1 = sum(f_n(1:1*12,:))*dt;
f_5 = sum(f_n((1*12+1):5*12,:))*dt;
f_10 = sum(f_n((5*12+1):10*12,:))*dt;
f_20 = sum(f_n((10*12+1):end,:))*dt;
plot(t,100*[(f_1'-f_1(1)) (f_5'-f_5(1)) (f_10'-f_10(1)) (f_20'-f_20(1))],'Linewidth',2)
legend({'0-1 year','1-5 years','5-10 years','10-20 years'},'Interpreter','Latex','FontSize',14)
xlabel('Time, $t$','FontSize',14,'interpreter','latex')
ylabel('Deviations from steady state (in $\%$)','FontSize',14,'interpreter','latex')
title('Debt $f(t)$','FontSize',14,'interpreter','latex')
hold on
xlim([0 300])

