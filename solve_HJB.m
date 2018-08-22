function [j_t]= HJB(r_vec,r_ss,parameters)
% This version codes the solution to the standard HJB equation in the paper
% to compute bond prices for a no arbitrage General Yield Curve 
% The PDE for the bond prices r(t)j_t(tau,t)=delta+dj/dt-dj/dtau
%                                 j(0,  t)= 1
%                                 j(tau,inf)=flat yield curve prices
% Output:  j(tau,inf)   is a I times N Matrix 

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
Bss       = (r_ss*speye(N_tau) - A);
j_ss      = Bss\u                  ;

%% Check Solution
checkit=1;
if checkit==1
    j_check   = delta/r_bar *(1-exp(-r_bar*tau))'+exp(-r_bar*tau)';
    %  Check with the analytical solution: good up to 5e-5
    res=sum((j_ss -j_check).^2);
end

%% HJB: Dynamic solution
j_n(:,N)   = j_ss; % Terminal value
for n = N-1:-1:1
    B = (1/dt + r_vec(n))*speye(N_tau) -A;
    d = u + j_n(:,n+1) /dt;
    j_n(:,n) = B\d;
end
    
% [psi_n,r_bar_n];   % Value function

end

