function [path_out] = solve_rss_path(rss,steady,paths,parameters)
% Solve_RSS_path
% solves the RSS for one shock iterating over {c_rss}U{c(o)-->c(ss)}
% This version codes Galo's solution using standard finite difference methods.
% Transition path, no shocks

%% Model Parameters
parameter_unpack;
time_preallocate;

% update solver method
solve_method=parameters.solve_method;

% Update Outer Parameters
gamma   = parameters.gamma   ; 
phi     = parameters.phi     ; 
rho     = parameters.rho     ;    
prob_mat= parameters.prob_mat;
N_r     = parameters.N_r     ;
N_y     = parameters.N_y     ;

% value of autarky
V_a=parameters.V_a;

%% Update Guess
% update guess for all possible paths
if strcmp(steady.type,'ss')
    c_in                  = zeros(N_r*N_y,N_t)+steady.c_ss  ;
    c_rss_in(1:N_r*N_y,1)= steady.c_ss*1.2                 ;
else
    strcmp(steady.type,'asymptotic')
    r_ss                = steady.r_ss;
    c_rss_in(1:N_r*N_y,1) = parameters.y_ss;
    c_in     =  0.5*c_rss_in.*exp(-(rho-r_ss)/gamma.*t);
end
%% ----------------------------
%MAIN LOOP
tic
figure
% solve_method='solver';
% Initial Guess!
c_in=[c_rss_in c_in]; 
Phi        = ones(N_r,N_y); % changed to cap Phi
rss.a      = 0; 

if strcmp(solve_method,'iter')
    for iter = 1:max_iter
        %% Bond prices in rss
        r_bar_ss           = parameters.r_bar_ss; 
        psi_0_mat          = rss.psi_0_mat     ;
        
        % Stack Solutions [s*]
        Epsi_0_mat=psi_0_mat*sum(Phi.*prob_mat,2);

        % Update RSS solution for pricce with default
        rss.psi_rss        = solve_HJB_rss(r_bar_ss,Epsi_0_mat,parameters); % PHI multiplies the price after default
        
        % Construct Deterministic path
        path_out  = compute_path(c_in,rss,steady,paths,parameters); % [s*:<- important changes inside, related to phi]
        c_out     = path_out.c_out; % [s*: first entry should be the same!, worried that first entry is not crss]
        
        c_in = relax_path * c_out + (1-relax_path) * c_in;
        V    = compute_value_functional(c_in,parameters) ;       
        
        residual=max(max(abs(c_in./c_out-1)));
        if iter>100
            if residual<tol_path
%                 iter=max_iter+1;
                disp('Converged in less than max iterations');
                break
            end
        end
        
        % Report
        if mod(iter,20) == 0
            if 1==0
                plot(c_out'); hold on; drawnow
             end
            disp(iter)
            disp(residual)
         end
        
         % Plot Iteration
         if rand>0.9
            figure(10)
            plot([c_out(:,1)*ones(1,10) c_out]'); hold on; drawnow
         end
        
        % Update default value
        % here goes the line recomputing PHI, PHI_prime and a
        rss_iter= path_out.rss             ; % [s*:<- DOUBLE CHECK THIS IS STEADY STATE, should be saved correctly]
        Phi     = phi_cdf(parameters,V,V_a); % [s*:<- need to fix here...need values of autarky [done]]
        Phi     = reshape(Phi,[N_y N_r])'  ;
        Phi_prime  = phi_pdf(parameters,V,V_a); % [s*:<- need to fix here...discuss with Galo [done]]
        mu_rss     = solve_multiplier_rss(parameters,rss_iter); % [s*:<- need to fix here...discuss with Galo [done]]
        EPhi=0;  EPhi_prime=0; EPhiMu=0; ii=0;
        for rr=1:N_r
            for yy=1:N_y
                ii=1+ii;
                EPhi =  EPhi   +prob_mat(rr,yy)*Phi(rr,yy);
                % EPhi_prime = EPhi_prime+prob_mat(rr,yy)*Phi_prime(ii);
                EPhiMu= EPhiMu +prob_mat(rr,yy)*Phi_prime(ii)*sum(mu_rss.*rss.psi_0_mat(:,rr)*dt);
            end
        end
        rss.a      = EPhi + EPhiMu; 
    end
    if iter == max_iter
        error('Error: maximum number of iterations reached')
    end
    
 

else
    strcmp(solve_method,'solver');  % THIS METHOD IS NOT IMPLEMENTED FOR DEFAULT
    fc2solve= @(x) x-c_path(x,rss,steady,paths,parameters);
    options=optimset('Algorithm','levenberg-marquardt');
    c_in   = fsolve(@(x) fc2solve(x),c_in,options);
    
end
path_out= compute_path(c_in,rss,steady,paths,parameters);
toc
path_out.psi_rss   = rss.psi_rss;
path_out.V         = V;
path_out.mu_rss    = mu_rss;

    


