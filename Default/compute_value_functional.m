%% Function that computes the value functional
function V=compute_value_functional(c_in,parameters)

% Unpack Model Parameters
parameter_unpack;
time_preallocate;

U = parameters.U;
rho   = parameters.rho  ;
phi      = parameters.phi   ; % intensity of the shock

% Variables 
c_n    = c_in(:,2:end);
[N_c,~]=size(c_n);

% Compute the value after the shock
V= sum(U(c_n).*(ones(N_c,1)*exp(-rho*t)),2)*dt; %[s*:]