function [psi_n,r_bar_n]= bond_prices(parameters)
% This version codes Galo's solution using standard finite difference methods.
% to compute bond prices for a No arbitrage General Yield Curve 
% The PDE for the bond prices r(t)psi(tau,t)=delta+dpdi/dt-dpsi/dtau
%                                 psi(0,  t)= 1
%                                 psi(tau,inf)=flat yield curve prices
% Output:  psi(tau,inf)   is a I times N Matrix with the arbitrage free
%                         bond prices

%% Model parameters
% Numerical parameters
t_max     = parameters.t_max;    % Maximum number of years (time)
tau_max   = parameters.tau_max;  % Maximum number of years (maturity)
dt        = parameters.dt;       % Monthly steps

% Preference Parameters
delta       = parameters.delta      ; % Coupon portion
r_bar       = parameters.r_bar      ; % short rate
r_bar_0     = parameters.r_bar_0    ; % initial value
rho_r_bar   = parameters.rho_r_bar  ; % persistence

%% Variables (preallocation)
t         = 0:dt:t_max;    % calendar time
N         = length(t);     % number of time nodes
tau       = dt:dt:tau_max; % maturities
I         = length(tau);   % number of tau nodes
u         = zeros (I,1);   % vector of coupons 
u(:)      = delta;         % vector of coupons 

u(1)      = u(1) +1/dt;
A0        = -speye(I);
aa        = ones(I-1,1);
A1        = spdiags(aa,-1,I,I);
A         = 1/dt*(A0+A1);
psi_n     = zeros(I,N);   % Bond Prices
r_bar_n   = zeros(N,1);

%%Path for the Short Rate
for i=2:N
    r_bar_n(1)= r_bar_0;   
    r_bar_n(i)= r_bar*(1-rho_r_bar)+ rho_r_bar*r_bar_n(i-1); 
end
      
%% HJB: Solution at Steady State
Bss         = (r_bar*speye(I) - A);
Psi_ss      = Bss\u;                   %
Psi_check   = delta/r_bar *(1-exp(-r_bar*tau))'+exp(-r_bar*tau)';

%  Check with the analytical solution: good up to 5e-5
%  error=sum((Psi_ss -Psi_check).^2);

%% HJB: Dynamic solution
psi_n(:,N)   = Psi_ss; % Terminal value
for n = N-1:-1:1
    B = (1/dt + r_bar_n(n))*speye(I) -A;
    d = u + psi_n(:,n+1) /dt;
    Psi_aux = B\d;
    psi_n(:,n) = Psi_aux;
end
    
% [psi_n,r_bar_n];   % Value function

end

