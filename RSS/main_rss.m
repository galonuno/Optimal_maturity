% Code to compute the Risky Steady State
%(c) Bigio, Nuno & Passadore model
% This version codes Galo's solution using standard finite difference methods.
% The Code Solves for the Risky Steady State
clear all
close all

%% Model parameters
% Model is annual
% Tolerance
parameters.tol_path = 1e-5;  % tolerance parameter optimal path
parameters.tol      = 1e-5;  % tolerance parameter
parameters.max_iter = 12000;  % maximum number of iterations
parameters.relax      = .2;%0.0001 ;  % relaxation coefficient
parameters.relax_path = .01;%0.0001 ;  % relaxation coefficient optimal path

% Numerical parameters
parameters.t_max     = 300   ; % maximum number of years (time)
parameters.tau_max   = 20    ; % maximum number of years (maturity)
parameters.dt        = 1/12  ; % monthly steps

% Preference Parameters
parameters.gamma   = 4.00 ;    % risk aversion
parameters.rho     = 0.06 ;    % discount factor - Model is Quarterly

% Output parameters
parameters.y_bar     = 1.0    ; % output at steady state
parameters.y_bar_0   = 1.0    ; % initial value, output after a shock
parameters.rho_y     = 1.0    ; % persistence output after shock - half life 1/rho months

%Interest Rate Shock Parameters
parameters.r_bar     = 0.04  ; % steady-state short rate
parameters.delta     = 0.04  ; % coupon portion
    parameters.r_bar_0   = 0.06;%0.08  ; % initial value
    parameters.rho_r_bar = 0;%0.995 ; % persistence (in months)

%Risky Steady State Parameters
parameters.phi  = 0.125;       % Poisson arrival rate

% adjustment cost
parameters.lambda_bar = 50;%50;

%% Computations
%Step 1: Bond Prices in Risky Steady State

psi    = bond_prices(parameters);
psi_rss= bond_prices_rss(parameters);

parameters.psi=psi;
parameters.psi_rss=psi_rss;



%% Step 2: Risky Steady State

results = risky_SS(parameters);

t_max     = parameters.t_max;    % Maximum number of years (time)
tau_max   = parameters.tau_max;  % Maximum number of years (maturity)
dt        = parameters.dt;       % Monthly steps
tau       = dt:dt:tau_max;

% Useful numbers
N_tau         = length(tau);                                    % number of maturity nodes

% Calculations 
parameters.yield_rss = parameters.delta-log(parameters.psi_rss)./(tau');

%% Main Plots
% Main Plot Params
pre_t=-10;
font_size=20;
line_width=4;
LineStyle_2=':';
base_color=[0.2 0.2 0.6];

% Main Plots
pre_t=-10;

%% Plots for the Transitional Dynamics
color_base=[0.1 0.1 0.6]; 
color1   =  [0.4 0.2 0.8];
color2   =  [0.2 0.6 0.8];
color3   =  [0.2 0.2 0.8];
color4   =  [0.8 0.2 0.8];
color_mat=  [color1; color2; color3; color4];
color_base2=[0.6 0.1 0.1];
color_rss=  [0.6 0.1 0.1];
color_rss_val=  [0.6 0.6 0.1];

% Ranges for debt levels
year_1=5;
year_2=10;
year_3=15;
year_4=20;
range_1  = (1:year_1*12);
range_5  = (year_1*12+1:year_2*12);
range_10 = (year_2*12+1:year_3*12);
range_20 = (year_3*12+1:N_tau);

% Time T variables
t         = 0:parameters.dt:parameters.t_max;    % calendar time

%% Plot the yield curve - RSS
figure(1)
subplot(2,2,1)
plot(tau,parameters.psi_rss,'Linewidth',line_width)
xlabel('Maturity $\tau$','FontSize',font_size,'interpreter','latex')
title('Price of Debt $\psi(\tau,rss)$','FontSize',font_size,'interpreter','latex')
hold on; axis tight;

subplot(2,2,2)
plot(tau,parameters.yield_rss ,'Linewidth',line_width)
xlabel('Maturity $\tau$','FontSize',font_size,'interpreter','latex')
title('Yield Curve $\Psi(\tau,rss)$','FontSize',font_size,'interpreter','latex')
hold on; axis tight;

% subplot(2,2,3:4)
% plot(t,parameters.r_bar_n*100,'Linewidth',3,'Color','b'); hold on;
% plot([pre_t 0],[parameters.r_bar parameters.r_bar]*100,'Linewidth',3,'Color','b'); axis tight;
% scatter(0,parameters.r_bar*100,80,'b','o');
% scatter(0,parameters.r_bar_n(1)*100,80,'b','o','filled');
% xlabel('Calendar Time $t$','FontSize',font_size,'interpreter','latex')
% title('Short Rate $\bar{r}_t$ (bps)','FontSize',font_size,'interpreter','latex')
% hold on; grid on; % ylim([0 max(parameters.r_bar_n)]);
% axis tight;

%Steady state
figure(2)
subplot(221)
plot(tau,results.f_n(:,end),'Linewidth',line_width,'Color',base_color); hold on;
plot(tau,results.f_rss(:),'Linewidth',line_width,'Color',color_rss)
xlabel('Maturity, $\tau$','FontSize',font_size,'interpreter','latex')
title('Outstanding Debt $f(\tau)$','FontSize',font_size,'interpreter','latex')
hold on; grid on;
subplot(222)
plot(tau,results.iota_n(:,end),'Linewidth',line_width,'Color',base_color); hold on;
plot(tau,results.iota_rss(:),'Linewidth',line_width,'Color',color_rss);
xlabel('Maturity, $\tau$','FontSize',font_size,'interpreter','latex')
title('Issuances $\iota(\tau)$','FontSize',font_size,'interpreter','latex')
hold on; grid on;
subplot(223)
plot(tau,results.f_n(:,end)./abs(sum(results.f_n(:,end))),'Linewidth',line_width,'Color',base_color); hold on;
plot(tau,results.f_rss(:)./abs(sum(results.f_rss(:))),'Linewidth',line_width,'Color',color_rss)      ;
xlabel('Maturity, $\tau$','FontSize',font_size,'interpreter','latex')
title('Maturity Distribution','FontSize',font_size,'interpreter','latex')
hold on; grid on;
subplot(224)
plot(tau,results.psi(:,end),'Linewidth',line_width,'Color',base_color); hold on;
plot(tau,results.psi_rss(:),'Linewidth',line_width,'Color',color_rss); hold on;
plot(tau,-results.V_n(:,end),'Linewidth',line_width,'Color',base_color,'LineStyle',LineStyle_2); 
plot(tau,-results.V_rss(:),'Linewidth',line_width,'LineStyle',LineStyle_2,'Color',color_rss); 
plot(tau,-results.V_n(:,1),'Linewidth',line_width,'Color',base_color/2,'LineStyle','-.');
xlabel('Maturity, $\tau$','FontSize',font_size,'interpreter','latex')
title('$\psi(\tau)$ vs. $-v(\tau)$','FontSize',font_size,'interpreter','latex');
hold on; grid on;
legend({'$\psi(\tau,ss)$','$\psi(\tau,rss)$','$-v(\tau,ss)$','$-v(\tau,rss)$','$-v(\tau,0)$'},...
    'Box','off','FontSize',14,'interpreter','latex');

%     'Orientation','Horizontal',...
    

%% TRANSITIONAL DYNAMICS AFTER The SHOCK

figure(3)
subplot(221)
plot(t,results.r_n*100,'Linewidth',line_width,'Color',color_base); hold on;
plot([pre_t 0],[parameters.rho parameters.rho]*100,'Linewidth',3,'Color',color_base);
scatter(0,parameters.rho*100,80,color_base,'o');
plot(0,parameters.rho*100,'*','Linewidth',line_width,'Color',color_rss_val,'MarkerSize',8);
scatter(0,results.r_n(1)*100,80,color_base,'o','filled'); 
scatter(t(end),results.r_n(end)*100,80,color_base,'>','filled'); 
xlabel('Time, $t$','FontSize',font_size,'interpreter','latex')
ylabel('$\%$','FontSize',font_size,'interpreter','latex')
title('Time Discount $r(t)$ (bps)','FontSize',font_size,'interpreter','latex')
hold on; xlim([pre_t t(end)]); grid on;

subplot(222)
plot(t,results.c_n,'Linewidth',line_width,'Color',color_base); hold on;
plot(t,results.y_n,'Linewidth',line_width,'Color',color_base2);
% scatter(0,results.c_ss,80,color_base,'o');
scatter(0,results.c_n(1),80,color_base,'o','filled'); 
scatter(t(end),results.c_n(end),80,color_base,'>','filled');
plot([pre_t 0],[results.c_rss results.c_rss],'Linewidth',3,'Color',color_base);
% plot([pre_t 0],[steady.c_ss steady.c_ss],'Linewidth',3,'Color',color_base);
% scatter(0,steady.y_ss,80,color_base2,'o');
scatter(0,results.y_n(1),80,color_base2,'o','filled'); 
scatter(t(end),results.y_n(end),80,color_base2,'>','filled');
plot([pre_t 0],[results.y_n(end) results.y_n(end)],'Linewidth',3,'Color',color_base2);
scatter(0,results.c_rss,80,color_base,'o');
plot(0,results.c_rss,'*','Linewidth',4,'Color',color_rss_val,'MarkerSize',8);
grid on; xlabel('Time, $t$','FontSize',font_size,'interpreter','latex');
title('$c(t)$ vs. $y(t)$','FontSize',font_size,'interpreter','latex');
xlim([pre_t t(end)]); legend({'$c(t)$','y(t)'},'interpreter','latex');
 
subplot(223);
i_1 =  sum(results.iota_n(range_1,:))*results.dt;
i_5 = sum(results.iota_n(range_5,:))*results.dt;
i_10 = sum(results.iota_n(range_10,:))*results.dt;
i_20 = sum(results.iota_n(range_20,:))*results.dt;
i_1_rss = sum(results.iota_rss(range_1,end))*results.dt;
i_5_rss = sum(results.iota_rss(range_5,end))*results.dt;
i_10_rss = sum(results.iota_rss(range_10,end))*results.dt;
i_20_rss = sum(results.iota_rss(range_20,end))*results.dt;
plot(t,(i_1'),'Linewidth',2,'Color',color_mat(1,:)); hold on;
plot(t,(i_5'),'Linewidth',2,'Color',color_mat(2,:));
plot(t,(i_10'),'Linewidth',2,'Color',color_mat(3,:));
plot(t,(i_20'),'Linewidth',2,'Color',color_mat(4,:));
plot([pre_t 0],[i_1_rss i_1_rss],'Linewidth',2,'Color',color_mat(1,:));
plot([pre_t 0],[i_5_rss i_5_rss],'Linewidth',2,'Color',color_mat(2,:));
plot([pre_t 0],[i_10_rss i_10_rss],'Linewidth',2,'Color',color_mat(3,:));
plot([pre_t 0],[i_20_rss i_20_rss],'Linewidth',2,'Color',color_mat(4,:));
scatter(0,i_1(1),80,color_mat(1,:),'o','filled');
scatter(0,i_5(1),80,color_mat(2,:),'o','filled'); 
scatter(0,i_10(1),80,color_mat(3,:),'o','filled'); 
scatter(0,i_20(1),80,color_mat(4,:),'o','filled');
scatter(t(end),i_1(end),80,color_mat(1,:),'>','filled');
scatter(t(end),i_5(end),80,color_mat(2,:),'>','filled'); 
scatter(t(end),i_10(end),80,color_mat(3,:),'>','filled'); 
scatter(t(end),i_20(end),80,color_mat(4,:),'>','filled');
scatter(0,i_1_rss,80,color_mat(1,:),'o');  plot(0,i_1_rss,'*','Linewidth',4,'Color',color_rss_val,'MarkerSize',8);
scatter(0,i_5_rss,80,color_mat(2,:),'o');  plot(0,i_5_rss,'*','Linewidth',4,'Color',color_rss_val,'MarkerSize',8);
scatter(0,i_10_rss,80,color_mat(3,:),'o'); plot(0,i_10_rss,'*','Linewidth',4,'Color',color_rss_val,'MarkerSize',8);
scatter(0,i_20_rss,80,color_mat(4,:),'o'); plot(0,i_20_rss,'*','Linewidth',4,'Color',color_rss_val,'MarkerSize',8);
xlim([pre_t t(end)]);
% legend({'0-5 year','5-10 years','10-15 years','15 $\geq$ years'},'Interpreter','Latex','FontSize',font_size)
xlabel('Time, $t$','FontSize',font_size,'interpreter','latex')
ylabel('$\%$ $y_{ss}$','FontSize',font_size,'interpreter','latex')
title('Issuances Rate, $\iota(t)$','FontSize',font_size,'interpreter','latex')
grid on;

subplot(224)
f_1 =  sum(results.f_n(range_1,:))*results.dt;
f_5 = sum(results.f_n(range_5,:))*results.dt;
f_10 = sum(results.f_n(range_10,:))*results.dt;
f_20 = sum(results.f_n(range_20,:))*results.dt;
f_1_rss = sum(results.f_rss(range_1,end))*results.dt;
f_5_rss = sum(results.f_rss(range_5,end))*results.dt;
f_10_rss = sum(results.f_rss(range_10,end))*results.dt;
f_20_rss = sum(results.f_rss(range_20,end))*results.dt;
plot(t,(f_1'),'Linewidth',2,'Color',color_mat(1,:)); hold on;
plot(t,(f_5'),'Linewidth',2,'Color',color_mat(2,:));
plot(t,(f_10'),'Linewidth',2,'Color',color_mat(3,:));
plot(t,(f_20'),'Linewidth',2,'Color',color_mat(4,:));
plot([pre_t 0],[f_1_rss f_1_rss],'Linewidth',2,'Color',color_mat(1,:));
plot([pre_t 0],[f_5_rss f_5_rss],'Linewidth',2,'Color',color_mat(2,:));
plot([pre_t 0],[f_10_rss f_10_rss],'Linewidth',2,'Color',color_mat(3,:));
plot([pre_t 0],[f_20_rss f_20_rss],'Linewidth',2,'Color',color_mat(4,:));
scatter(0,f_1(1),80,color_mat(1,:),'o','filled');
scatter(0,f_5(1),80,color_mat(2,:),'o','filled'); 
scatter(0,f_10(1),80,color_mat(3,:),'o','filled'); 
scatter(0,f_20(1),80,color_mat(4,:),'o','filled');
scatter(t(end),f_1(end),80,color_mat(1,:),'>','filled');
scatter(t(end),f_5(end),80,color_mat(2,:),'>','filled'); 
scatter(t(end),f_10(end),80,color_mat(3,:),'>','filled'); 
scatter(t(end),f_20(end),80,color_mat(4,:),'>','filled');
scatter(0,f_1_rss,80,color_mat(1,:),'o');  plot(0,f_1_rss,'*','Linewidth',4,'Color',color_rss_val,'MarkerSize',8);
scatter(0,f_5_rss,80,color_mat(2,:),'o');  plot(0,f_5_rss,'*','Linewidth',4,'Color',color_rss_val,'MarkerSize',8);
scatter(0,f_10_rss,80,color_mat(3,:),'o'); plot(0,f_10_rss,'*','Linewidth',4,'Color',color_rss_val,'MarkerSize',8);
scatter(0,f_20_rss,80,color_mat(4,:),'o'); plot(0,f_20_rss,'*','Linewidth',4,'Color',color_rss_val,'MarkerSize',8);
xlim([pre_t t(end)]);
  legend({'0-5 (years)','5-10 ','10-15 ','T-5 to T '},...
        'Position',[0.332239224848658 -0.00773908315260428 0.373750358484675 0.0524046428167998],'Interpreter','Latex','FontSize',16,'Location','SouthEast','Box','off','Orientation','horizontal');

xlabel('Time, $t$','FontSize',font_size,'interpreter','latex')
ylabel('$\%$ $y_{ss}$','FontSize',font_size,'interpreter','latex')
title('Debt Outstanding, $f(t)$','FontSize',font_size,'interpreter','latex')
grid on;
