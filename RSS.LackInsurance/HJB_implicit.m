%% HJB_implicit
function [V]=HJB_implicit(U,s,muF,muB,sigma,rho,Delta,x,Ib,If)
% ------------------------------------------
% USES:
% ------------------------------------------
%
% Solving: r*V(s) - mu(s)*V'(s) - sigma(s)^2*V''(s)/2 = U(s)
%
% rho is the discount rate 
% s = [s(1), s(2) ... s(N)] is the state space (in increasing order)
% mu is a drift vector of length N, with MU(1) >= 0 and MU(N) <= 0
% sigma is a volatility vector of lenght N, with S(1) = S(N) = 0
% U is a payoff flow of length N
% In implicit method we are solving for:
% --------------------------------------------------------------------------
% r*V = G + MV + SV,
% Where MV is drift of value and SV is volatility
% This expression can be written as:
% A*V=G where
%
% A = I*rho + (-1 translation of diagonal) negative difference + (+1
% translation of diagonal) positive difference
% --------------------------------------------------------------------------
% Note, 
% 1. we are assuming no volatility at the right or left ----> code is
%    not general yet
% 2. It is compatible with 1/Delta = 0 in Moll's algorithm. 
% --------------------------------------------------------------------------
% Vectorization
N = length(s)         ; % Length
ds = s(2:N) - s(1:N-1); % Change in S.  

% Construction of Matrix
% Volatility Component 
S = zeros(N,1)                                       ; 
S(2:N-1) = sigma(2:N-1).^2./(ds(1:N-2) + ds(2:N-1))  ; % Note we are using 2 first differences first difference here

%% Saki's Method
% % Plus diagonal
% MU = zeros(N,1)              ; 
% MU(1:N-1) = max(muF(1:N-1),0) ;
MU  = If.*muF;
SU = zeros(N,1)              ;
SU(1:N-1) = S(1:N-1)         ;
DU = -(MU+SU)./[ds; 1]       ; % Negative from right of HJB 

% Minus Diagonal Terms
% MD = zeros(N,1)          ;
% MD(2:N) = -min(muB(2:N),0);
MD = -Ib.*muB;
SD = zeros(N,1)          ;
SD(2:N) = S(2:N)         ;
DD = -(MD+SD)./[1; ds]   ; % Negative from right of HJB

% On Diagonal Terms
D0        = rho*ones(N,1) - DU - DD;%+ (MU + MD + 2*S)/ds(1);


% Making it Compatible with Sparse Matrix Notation
DU=[0; DU(1:end-1)];
DD=[DD(2:end); 0];

A= spdiags([DD D0 DU],-1:1,N,N);
% test=spdiags(A);

V = (A+speye(N)/Delta)\(U + x/Delta);
end