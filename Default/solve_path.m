function [path]=solve_path(rr,yy,f_0,v_n,paths,parameters)
% Unpacking Code Parameters
parameter_unpack;
time_preallocate;

% Unpack Parameters
delta = parameters.delta;
lambda_bar= parameters.lambda_bar;

% Unpack Exogenous Paths
psi_n   = paths.psi_n(:,:,rr); 
y_n     = paths.y_n(yy,:);

% Optimal policies
iota_n = 1/lambda_bar*(1-v_n./psi_n); % Optimal policies

% Construct - KFE Operator
A0        = -speye(N_tau);
aa        = ones(N_tau-1,1);
A1        = spdiags(aa,-1,N_tau,N_tau);
A         = 1/dt*(A0+A1);
D         = speye(N_tau) - dt *A';

% KFE 
f_n(:,1)   = f_0;
for n=2:N_t
    h      = iota_n(:,n) *dt + f_n(:,n-1);
    f_aux  = D\h;
    f_n(:,n) = f_aux;
end

% Consumption and interest rates
c_n = y_n  - f_n(1,:) + sum(psi_n.*(1-1/2*lambda_bar.*iota_n).*iota_n-delta*f_n)*dt;
path.c_n=c_n      ;
path.f_n=f_n      ;
path.iota_n=iota_n;
