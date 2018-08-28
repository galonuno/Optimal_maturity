function y_curve= yield_curve(parameters)

t_max     = parameters.t_max;    % Maximum number of years (time)
tau_max   = parameters.tau_max;  % Maximum number of years (maturity)
dt        = parameters.dt;       % Monthly steps

% Preference Parameters
delta       = parameters.delta      ; % Coupon portion
r_bar       = parameters.r_bar      ; % short rate
r_bar_0     = parameters.r_bar_0    ; % initial value
rho_r_bar   = parameters.rho_r_bar  ; % persistence
psi         = parameters.psi;       % Input prices

%% Variables (preallocation)
t         = 0:dt:t_max;    % calendar time
N         = length(t);     % number of time nodes
tau       = dt:dt:tau_max; % maturities
I         = length(tau);
P_0       = zeros(I,N);
P         = psi;
y_curve   = zeros(I,N);
A         = tril(ones(I,I));
Id        = eye(I);
A1        = (delta*dt*A+Id);

for i=1:N;
    P_0(:,i)=A1\P(:,i);
end

for i=1:N;
   y_curve(:,i)=-log(P_0(:,i))./tau';
end









        
        
        
        


