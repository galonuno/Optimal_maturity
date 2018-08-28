function [psi_n,r_bar_n]= bond_prices(parameters)
% This version codes Galo's solution using standard finite difference methods.
% to compute bond prices for a no arbitrage General Yield Curve 
% The PDE for the bond prices r(t)psi(tau,t)=delta+dpdi/dt-dpsi/dtau
%                                 psi(0,  t)= 1
%                                 psi(tau,inf)=flat yield curve prices
% Output:  psi(tau,inf)   is a I times N Matrix with the arbitrage free
%                         bond prices

%% Model parameters
% Numerical parameters
time_preallocate;

% Preference Parameters
delta       = parameters.delta      ; % Coupon portion

%% Preallocation - PDE terms
u         = zeros (N_tau,1);   % vector of coupons 
u(:)      = delta;         % vector of coupons 
u(1)      = u(1) +1/dt;
A0        = -speye(N_tau);
aa        = ones(N_tau-1,1);
A1        = spdiags(aa,-1,N_tau,N_tau);
A         = 1/dt*(A0+A1);
psi_n     = zeros(N_tau,N);   % Bond Prices
r_bar_n   = zeros(N,1);
      
%% HJB: Solution at Steady State
Bss         = (r_bar*speye(N_tau) - A);
Psi_ss      = Bss\u;                   %
% Psi_check   = delta/r_bar *(1-exp(-r_bar*tau))'+exp(-r_bar*tau)';

%  Check with the analytical solution: good up to 5e-5
%  error=sum((Psi_ss -Psi_check).^2);

%% HJB: Dynamic solution
psi_n(:,N)   = Psi_ss; % Terminal value
for n = N-1:-1:1
    B = (1/dt + r_bar_n(n))*speye(N_tau) -A;
    d = u + psi_n(:,n+1) /dt;
    Psi_aux = B\d;
    psi_n(:,n) = Psi_aux;
end
    
% [psi_n,r_bar_n];   % Value function

end

