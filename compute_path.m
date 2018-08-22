%% Function that solves the consumption path
function [path_out]=compute_path(c_in,rss,steady,paths,parameters)

% Unpack Model Parameters
parameter_unpack ;
time_preallocate ;
type=steady.type ;

% Unpack Parameters
rho   = parameters.rho          ;
sigma = parameters.gamma        ;
phi   = parameters.phi          ; % [s*: didn't see intensity anywhere]
U_p_ratio = parameters.U_p_ratio;
N_r     = parameters.N_r        ;
N_y     = parameters.N_y        ;
prob_mat= parameters.prob_mat   ;

% Unpack Risky Steady State
psi_rss   = rss.psi_rss; 

%% Step 1: Construct r_n [s*]
% [!] Warning, I am assuming that path reaches a proper STEADY STATE
% [!!] doing this one at a time
% [!!!] not sure what to do at non-convergent solution
c_rss         = c_in(:,1);
c_0           = c_in(:,2); 
c_in          = c_in(:,2:N_t+1);
r_n(:,1:N_t-1)= rho+sigma*(c_in(:,2:N_t)./c_in(:,1:N_t-1)-1)*dt;
if strcmp(type,'ss')
    r_n(1:N_r*N_y,N_t)    = rho;
elseif strcmp(type,'asymptotic')
    r_n(1:N_r*N_y,N_t)    = steady.r_ss;    
end

%% Step 2: Value at Deterministic Path [s*]
ii=0;
v_0_mat=zeros(N_tau,1);
v_n    =zeros(N_tau,N_t,N_r*N_y)*NaN;
for rr=1:N_r
    for yy=1:N_y
        ii=1+ii;
        v_n(:,:,ii) = solve_HJB_path(r_n(ii,:),rho,parameters);
        v_0 = v_n(:,1,ii)                                     ;

        %% Step 2.b: Solving for Values at RSS -> not sure where phi is...
        v_0_mat=v_0_mat+phi*prob_mat(rr,yy)*v_0*U_p_ratio(c_0(ii),c_rss(ii))             ; % [s*:<- important change here]
        % adapt with rss.a term
    end
end

%% Step 3: Solving for Risky Steady State
rss_out=solve_rss(psi_rss,v_0_mat,parameters);

%% Step 4: Solving for Path
f_0=rss_out.f_rss;
c_out=zeros(N_r*N_y,1+N_t);
f_out=zeros(N_tau,1+N_t,N_r*N_y);
iota_out=zeros(N_tau,1+N_t,N_r*N_y);

% new iteration
ii=0;
for rr=1:N_r
    for yy=1:N_y
        ii=1+ii;
        path_out=solve_path(rr,yy,f_0,v_n(:,:,ii),paths,parameters);
            
        %% Step 5: Packaging Output
        c_out(ii,:)=[rss_out.c_rss path_out.c_n];
        f_out(:,:,ii)=[rss_out.f_rss path_out.f_n];
        iota_out(:,:,ii)=[rss_out.iota_rss path_out.iota_n];

    end
end

% Packaging Output
ii=0;
for rr=1:N_r
    for yy=1:N_y
        ii=1+ii;
        path_out.v_out(:,:,ii) = [rss_out.v_rss v_n(:,:,ii)]    ;
    end
end
path_out.v_n     = v_n                  ;
path_out.c_out   = c_out                ;
path_out.f_out   = f_out                ;
path_out.iota_out= iota_out             ;
path_out.r_n     = r_n                  ;
path_out.rss     = rss_out              ;
