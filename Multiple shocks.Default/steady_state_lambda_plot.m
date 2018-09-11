% STEADY STATE ANALYSIS
%(c) Bigio, Nuno & Passadore model
clear; close all;

% This code solves Galo's approach to the steady state for any value of
% lambda

%% Model parameters
% Model is annual

% Tolerance
tol      = 1e-5;
max_iter = 5000;  % Maximum number of iterations
relax    = 0.0001; % Relaxation coefficient
% Numerical parameters
tau_max   = 20      ; % Maximum number of years (maturity)
dt        = 1/12    ; % Monthly steps

% Preference Parameters
rho     = 0.06;  % Discount factor - Model is Quarterly
delta   = 0.03;  % Coupon portion
y_bar   = 1   ;  % output
gamma   = 2   ;  % risk aversion

% Frictions parameters
lambda_bar_vec = 0.01:0.01:2;
N_lambda         = length(lambda_bar_vec) ;

%Yield curve
r_bar   = 0.04 ; % Market discount factor


%% Variables (preallocation)
tau       = dt:dt:tau_max; % maturities
I         = length(tau);   % Number of tau nodes
u         = zeros (I,1);

u(:)      = -delta;
u(1)      = u(1) -1/dt;
A0        = -speye(I);
aa        = ones(I-1,1);
A1        = spdiags(aa,-1,I,I);
A         = 1/dt*(A0+A1);
psi      = delta/r_bar*(1-exp(-r_bar*tau')) + exp(-r_bar*tau') ;% Bond prices
 
%%
%%----------------------------
%MAIN LOOP
tic

for i=1:N_lambda 

lambda_bar=lambda_bar_vec(i);    
% The case lambda >lambda_min (Guess)
V     = -delta/rho *(1-exp(-rho*tau'))-exp(-rho*tau');
% Optimal policies
iota  = 1/lambda_bar * (V+ psi );
% KFE
f     = flip(cumsum(flip(iota)))*dt;
% Consumption
c = y_bar  - f(1) + sum((psi - 1/2 *lambda_bar * iota).*iota-delta*f) * dt;

if c<=0 % The case lambda <lambda_min
    r = (r_bar + rho) /2; % First guess
    for iter = 1:max_iter
        V     = -delta/r *(1-exp(-r*tau'))-exp(-r*tau');
        % Optimal policies
        iota  = 1/lambda_bar * (V+ psi );
        % KFE
        f     = flip(cumsum(flip(iota)))*dt;
        % Consumption
        c = y_bar  - f(1) + sum((psi - 1/2 *lambda_bar * iota).*iota-delta*f) * dt;
        if  abs(c) < tol
            break  % Consumption should be zero
        else
            r = r + relax*c;
            
        end
    end
    if iter == max_iter
        disp('Error: maximum number of iterations reached')
    end 
    if i>1
        if c_ss(i-1)>0
           lambda_o=lambda_bar; 
           I_lambda_o=i;
        end
    end
end


c_ss(i)=c;
r_ss(i)=r;
f_ss(i,:)=f;
iota_ss(i,:)=iota;
end
toc

% %% Plotting
% %Steady state
% figure(1)
% subplot(121)
% plot(lambda_bar_vec,c_ss,'Linewidth',2)
% xlabel('$\lambda$', 'FontSize',14,'interpreter','latex')
% title('$c_{ss}$','FontSize',14,'interpreter','latex')
% hold on
% 
% subplot(122)
% plot(lambda_bar_vec,r_ss,'Linewidth',2)
% xlabel('$\lambda$', 'FontSize',14,'interpreter','latex')
% title('$r_{ss}$','FontSize',14,'interpreter','latex')
% 
% figure(2)
% subplot(221)
% plot(lambda_bar_vec,c_ss,'Linewidth',2)
% xlabel('$\lambda$', 'FontSize',14,'interpreter','latex')
% title('$c_{ss}$','FontSize',14,'interpreter','latex');
% hold on
% 
% subplot(222)
% plot(lambda_bar_vec,r_ss,'Linewidth',2)
% xlabel('$\lambda$', 'FontSize',14,'interpreter','latex')
% title('$r_{ss}$','FontSize',14,'interpreter','latex');
% 
% subplot(223)
% surf(tau,lambda_bar_vec,f_ss,'edgecolor','none','FaceAlpha',0.5);
% xlabel('$\tau$', 'FontSize',14,'interpreter','latex');
% ylabel('$\lambda$', 'FontSize',14,'interpreter','latex');
% title('$f_{ss}$','FontSize',14,'interpreter','latex');
% 
% subplot(224)
% surf(tau,lambda_bar_vec,iota_ss,'edgecolor','none','FaceAlpha',0.5)
% xlabel('$\tau$', 'FontSize',14,'interpreter','latex');
% ylabel('$\lambda$', 'FontSize',14,'interpreter','latex');
% title('$\iota_{ss}$','FontSize',14,'interpreter','latex');
% 
% % orient landscape
% % saveas(gcf,'F_steadystate.lambdacomp','pdf')
% 
% figure(3)
% subplot(221)
% plot(lambda_bar_vec,c_ss,'Linewidth',2)
% xlabel('$\lambda$', 'FontSize',14,'interpreter','latex')
% title('$c_{ss}$','FontSize',14,'interpreter','latex');
% hold on
% 
% subplot(222)
% plot(lambda_bar_vec,r_ss,'Linewidth',2)
% xlabel('$\lambda$', 'FontSize',14,'interpreter','latex')
% title('$r_{ss}$','FontSize',14,'interpreter','latex');
% 
% subplot(223)
% contourf(tau,lambda_bar_vec,f_ss);
% xlabel('$\tau$', 'FontSize',14,'interpreter','latex');
% ylabel('$\lambda$', 'FontSize',14,'interpreter','latex');
% title('$f_{ss}$','FontSize',14,'interpreter','latex');
% 
% subplot(224)
% contour(tau,lambda_bar_vec,iota_ss)
% xlabel('$\tau$', 'FontSize',14,'interpreter','latex');
% ylabel('$\lambda$', 'FontSize',14,'interpreter','latex');
% title('$\iota_{ss}$','FontSize',14,'interpreter','latex');

figure(1)
subplot(221)
plot(lambda_bar_vec,c_ss,'Linewidth',2,'Color','k'); grid on;
[y1s]=ylim;
line([lambda_o lambda_o],[ y1s(1) y1s(2)],'Color','k','LineWidth',2,'LineStyle','--');
xlabel('Liquidity, $\bar{\lambda}$', 'FontSize',14,'interpreter','latex')
title('$c_{\infty}$','FontSize',14,'interpreter','latex');
text(lambda_o,y1s(1)-0.1*(y1s(2)-y1s(1)),'$\lambda_{o}$','interpreter','latex');

subplot(222)
plot(lambda_bar_vec,r_ss,'Linewidth',2,'Color','k')
[y1s]=ylim;
line([lambda_o lambda_o],[ y1s(1) y1s(2)],'Color','k','LineWidth',2,'LineStyle','--');
xlabel('Liquidity, $\bar{\lambda}$', 'FontSize',14,'interpreter','latex')
title('$r_{\infty}$','FontSize',14,'interpreter','latex');
text(lambda_o,y1s(1)-0.1*(y1s(2)-y1s(1)),'$\lambda_{o}$','interpreter','latex'); grid on;

subplot(223)
for ii=1:10:N_lambda 
    plot(tau,iota_ss(ii,:),'LineWidth',0.1,'Color',[ii/N_lambda ii/N_lambda ii/N_lambda]); hold on;
end
plot(tau,iota_ss(I_lambda_o,:),'LineWidth',2,'Color','r'); hold on;

    xlabel('Maturity, $\tau$', 'FontSize',14,'interpreter','latex');
    %ylabel('$\lambda$', 'FontSize',14,'interpreter','latex');
    title('$\iota_{\infty}$','FontSize',14,'interpreter','latex');  grid on;

    
subplot(224)
for ii=1:10:N_lambda 
    plot(tau,f_ss(ii,:),'LineWidth',0.1,'Color',[ii/N_lambda ii/N_lambda ii/N_lambda]); hold on;
end
plot(tau,f_ss(I_lambda_o,:),'LineWidth',2,'Color','r'); hold on;
xlabel('Maturity, $\tau$', 'FontSize',14,'interpreter','latex');
%ylabel('$\lambda$', 'FontSize',14,'interpreter','latex');
title('$f_{\infty}$','FontSize',14,'interpreter','latex');  grid on;

orient landscape;
saveas(gcf,'F_PFlambdacomp','pdf');