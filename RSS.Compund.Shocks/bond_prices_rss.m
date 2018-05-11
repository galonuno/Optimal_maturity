function results = bond_prices_rss(parameters,precompute)

% To compute bond prices for a No arbitrage General Yield Curve with a
% risky shock.The PDE for the bond prices after the shock is given by 
%                                 r(t)psi(tau,t)=delta+dpdi/dt-dpsi/dtau
%                                 psi(0,  t)= 1
%                                 psi(tau,inf)=flat yield curve prices
% The PDE for the bond prices before the shock is 
%                                 r_rss*psi(tau)=delta+pi*E_r'[psi(tau,0,r')-psi(tau)]
%                                    -dpsi/dtau
%                                 
% Output:  psi(tau,inf)   is a I times N Matrix with the arbitrage free
%                         bond prices


% %Step 1: Compute the Bond prices after the shock
%          OR get the  Bond prices after the shock.  

 psi = precompute.psi;

% WE OMIT THIS STEP BECAUSE NOW WE ENTER AS PARAMETERS THE PRICES.  
% psi=bond_prices(parameters);

%Step 2: Compute the Risky Steady State

% Preference Parameters
delta   = parameters.delta;  %Coupon
phi     = parameters.phi;    %Probability of the Shock
r_bar   = parameters.r_bar;  %Steady State Short rate before the shock

% Shock Parameters

S                    = parameters.S; %number of shocks to the interest rate
p_r_bar_0            = (1/S)*ones(S,1);
 

% Numerical parameters
t_max     = parameters.t_max;    % Maximum number of years (time)
tau_max   = parameters.tau_max;  % Maximum number of years (maturity)
dt        = parameters.dt;       % Monthly steps
t         = 0:dt:t_max;    % calendar time
N         = length(t);     % number of time nodes
tau       = dt:dt:tau_max; % maturities
I         = length(tau);   % number of tau nodes
u         = zeros (I,1);   % vector of coupons 



u(:)      = delta +phi*reshape (psi(:,1,:), [I S])*p_r_bar_0 ;  % vector of coupons 

u(1)      = u(1) +1/dt;
A0        = -speye(I);
aa        = ones(I-1,1);
A1        = spdiags(aa,-1,I,I);
A         = 1/dt*(A0+A1);

Bss       = ((r_bar+phi)*speye(I) - A);
psi_rss   = Bss\u;  

% Reported results
results = psi_rss;

%end