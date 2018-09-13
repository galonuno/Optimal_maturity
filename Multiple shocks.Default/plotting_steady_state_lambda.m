% Change in parameters

parameters.relax_path = 0.0001;  % relaxation coefficient optimal path

% Price of Risk
parameters.U         =  @(c) (c.^(1-parameters.gamma)-1)/(1-parameters.gamma);
parameters.U_p_ratio =  @(c_a,c_b) (c_a/c_b).^(-parameters.gamma);
parameters.dU        =  @(c) (c.^(-parameters.gamma));

%% Pre-Allocation
time_preallocate;
psi_ss = parameters.delta/parameters.r_bar*(1-exp(-parameters.r_bar*tau')) + exp(-parameters.r_bar*tau') ;% Bond prices

%% Main Loop
lambda_bar_vec = 0.01:0.005:1.5;
N_lambda         = length(lambda_bar_vec) ;

for i=1:N_lambda
    parameters.lambda_bar   = lambda_bar_vec(i) ; %
    steady=solve_steady(psi_ss,parameters)      ;
    
    c_ss(i)=steady.c_ss;
    r_ss(i)=steady.r_ss;
    f_ss(i,:)=steady.f_ss;
    iota_ss(i,:)=steady.iota_ss;
    if i>1
        if c_ss(i-1)==0 & c_ss(i)>0
           lambda_o=parameters.lambda_bar; 
           I_lambda_o=i;
        end
    end
end

%% Plotting
figure(1)
subplot(221)
plot(lambda_bar_vec,c_ss,'Linewidth',2,'Color','k'); grid on;
[y1s]=ylim;
line([lambda_o lambda_o],[ y1s(1) y1s(2)],'Color','k','LineWidth',2,'LineStyle','--');
xlabel('Liquidity, $\bar{\lambda}$', 'FontSize',14,'interpreter','latex')
title('Consumption, $c_{\infty}$','FontSize',14,'interpreter','latex');
text(lambda_o,y1s(1)-0.1*(y1s(2)-y1s(1)),'$\lambda_{o}$','interpreter','latex','FontSize',14);

subplot(222)
plot(lambda_bar_vec,r_ss,'Linewidth',2,'Color','k')
[y1s]=ylim;
line([lambda_o lambda_o],[ y1s(1) y1s(2)],'Color','k','LineWidth',2,'LineStyle','--');
xlabel('Liquidity, $\bar{\lambda}$', 'FontSize',14,'interpreter','latex')
title('Discount, $r_{\infty}$','FontSize',14,'interpreter','latex');
text(lambda_o,y1s(1)-0.1*(y1s(2)-y1s(1)),'$\lambda_{o}$','interpreter','latex','FontSize',14); grid on;
ylim([parameters.r_bar parameters.rho+0.0001])

subplot(223)
for ii=1:10:N_lambda 
    plot(tau,iota_ss(ii,:),'LineWidth',0.1,'Color',[ii/N_lambda ii/N_lambda ii/N_lambda]); hold on;
end
plot(tau,iota_ss(I_lambda_o,:),'LineWidth',2,'Color','r'); hold on;

    xlabel('Maturity, $\tau$', 'FontSize',14,'interpreter','latex');
    %ylabel('$\lambda$', 'FontSize',14,'interpreter','latex');
    title('Issuances, $\iota_{\infty}$','FontSize',14,'interpreter','latex');  grid on;

    
subplot(224)
for ii=1:10:N_lambda 
    plot(tau,f_ss(ii,:),'LineWidth',0.1,'Color',[ii/N_lambda ii/N_lambda ii/N_lambda]); hold on;
end
plot(tau,f_ss(I_lambda_o,:),'LineWidth',2,'Color','r'); hold on;
xlabel('Maturity, $\tau$', 'FontSize',14,'interpreter','latex');
%ylabel('$\lambda$', 'FontSize',14,'interpreter','latex');
title('Distribution, $f_{\infty}$','FontSize',14,'interpreter','latex');  grid on;

cd figures
orient landscape;
saveas(gcf,'steady_state_lambda','pdf');
saveas(gcf,'steady_state_lambda','fig');
cd ..