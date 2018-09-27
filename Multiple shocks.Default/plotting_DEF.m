%% Main Plots YY
% Loading
file_RSS     = 'RSS_Y';
file_DEF     = 'DEF';

% Parameters
pre_t       = -10;
font_size   = 14;
line_width  =2;
base_color  = [0.2 0.2 0.6];

pre_t       = -10;
pos_t       = 60;

color_base =[0.1 0.1 0.6]; 
color1   =  [0.4 0.2 0.8];
color2   =  [0 0 0];
color3   =  [0.2 0.7 0.3];
color4   =  [0 0.4 0.1];
color_mat=  [color1; color2; color3; color4];
color_base2=[0.6 0.1 0.1];
color_rss=  [0.6 0.1 0.1];
color_rss_val=  [0.6 0.6 0.1];
markersize   = 40;
LineStyle_2 = '-.';
LineStyle_3 = '-.';
LineStyle_4 = '-.';


cd figures
load(file_RSS)
cd ..

% Variables
dt = parameters.dt;
t         = 0:parameters.dt:parameters.t_max;    % calendar time
Debt_rss = sum(results.f_rss(:,end))*dt;
Debt     = sum(results.f_n(:,:))*dt;
[average_duration, average_term_to_maturity] = Portfolio_duration(results, parameters, paths);


% %% Steady State - Comparisons
figure
subplot(321)
plot(tau,results.f_rss,'Linewidth',line_width,'Color',color2); hold on;
xlabel('Maturity, $\tau$','FontSize',font_size,'interpreter','latex')
title('(a) Debt distribution, $f(\tau)$','FontSize',font_size,'interpreter','latex')
hold on; grid on;

subplot(3,2,2)
plot(t,results.c_n,'Linewidth',line_width,'Color',color2); hold on;
scatter(0,results.c_n(1),markersize,color2,'o','filled'); 
scatter(t(end),results.c_rss,markersize,color2,'>','filled');
plot([pre_t 0],[results.c_rss results.c_rss],'Linewidth',line_width,'Color',color2);
scatter(0,results.c_rss,markersize,color2,'o');
plot(0,results.c_rss,'*','Linewidth',4,'Color',color_rss_val,'MarkerSize',2);
grid on; xlabel('Time, $t$ (years)','FontSize',font_size,'interpreter','latex');
title('(b) Consumption, $c(t)$','FontSize',font_size,'interpreter','latex');
xlim([pre_t pos_t]); 

subplot(3,2,3)
plot(t,Debt*100,'Linewidth',line_width,'Color',color2); hold on;
plot([pre_t 0],[Debt_rss Debt_rss]*100,'Linewidth',line_width,'Color',color2);
scatter(0,Debt_rss*100,markersize,color2,'o');
plot(0,Debt_rss*100,'*','Linewidth',line_width,'Color',color_rss_val,'MarkerSize',2);
scatter(0,Debt(1)*100,markersize,color2,'o','filled'); 
scatter(t(end),Debt(end)*100,markersize,color2,'>','filled'); 
xlabel('Time, $t$ (years)','FontSize',font_size,'interpreter','latex')
ylabel('$\%$ $y_{ss}$','FontSize',font_size,'interpreter','latex')
title('(c) Total debt, $b(t)$','FontSize',font_size,'interpreter','latex')
hold on; xlim([pre_t pos_t]); grid on;

subplot(3,2,4)
plot(t,average_duration,'Linewidth',line_width,'Color',color2); hold on;
plot([pre_t 0],[average_duration(1) average_duration(1)],'Linewidth',line_width,'Color',color2);
scatter(0,average_duration(1),markersize,color2,'o');
plot(0,average_duration(1),'*','Linewidth',line_width,'Color',color_rss_val,'MarkerSize',2);
scatter(0,average_duration(1),markersize,color2,'o','filled'); 
scatter(t(end),average_duration(end),markersize,color2,'>','filled'); 
xlabel('Time, $t$ (years)','FontSize',font_size,'interpreter','latex')
ylabel('Years','FontSize',font_size,'interpreter','latex')
title('(d) Average duration','FontSize',font_size,'interpreter','latex')
hold on; xlim([pre_t pos_t]); grid on;

cd figures
load(file_DEF)
cd ..
Debt_rss = sum(results.f_rss(:,end))*dt;
Debt     = sum(results.f_n(:,:))*dt;
[average_duration, average_term_to_maturity] = Portfolio_duration(results, parameters, paths);

subplot(321)
plot(tau,results.f_rss(:,end),'Linewidth',line_width,'Color',color4,'LineStyle',LineStyle_4); hold on;
legend({'RSS','Default' },...
    'Box','off','FontSize',12,'interpreter','latex');

subplot(3,2,2)
plot(t,results.c_n,'Linewidth',line_width,'Color',color4,'LineStyle',LineStyle_4); hold on;
scatter(0,results.c_n(1),markersize,color4,'o','filled'); 
scatter(t(end),results.c_n(end),markersize,color4,'>','filled');
plot([pre_t 0],[results.c_rss results.c_rss],'Linewidth',line_width,'Color',color4,'LineStyle',LineStyle_4);
scatter(0,results.c_rss,markersize,color4,'o');
plot(0,results.c_rss,'*','Linewidth',4,'Color',color_rss_val,'MarkerSize',2);

subplot(3,2,3)
plot(t,Debt*100,'Linewidth',line_width,'Color',color4,'LineStyle',LineStyle_4); hold on;
plot([pre_t 0],[Debt_rss Debt_rss]*100,'Linewidth',line_width,'Color',color4,'LineStyle',LineStyle_4);
scatter(0,Debt_rss*100,markersize,color4,'o');
plot(0,Debt_rss*100,'*','Linewidth',line_width,'Color',color_rss_val,'MarkerSize',2);
scatter(0,Debt(1)*100,markersize,color4,'o','filled'); 
scatter(t(end),Debt(end)*100,markersize,color4,'>','filled'); 

subplot(3,2,4)
plot(t,average_duration,'Linewidth',line_width,'Color',color4,'LineStyle',LineStyle_4); hold on;
plot([pre_t 0],[average_duration(1) average_duration(1)],'Linewidth',line_width,'Color',color4,'LineStyle',LineStyle_4);
scatter(0,average_duration(1),markersize,color4,'o');
plot(0,average_duration(1),'*','Linewidth',line_width,'Color',color_rss_val,'MarkerSize',2);
scatter(0,average_duration(1),markersize,color4,'o','filled'); 
scatter(t(end),average_duration(end),markersize,color4,'>','filled'); 

subplot(3,2,5)
plot(tau,steady.psi_ss(:),'Linewidth',line_width,'Color',color2,'LineStyle',LineStyle_2); hold on;
plot(tau,results.psi_rss(:),'Linewidth',line_width,'LineStyle',LineStyle_4,'Color',color4);
plot(tau,results.psi_n(:,1),'Linewidth',line_width,'Color',color3,'LineStyle',':');
xlabel('Maturity, $\tau$','FontSize',font_size,'interpreter','latex')
title('(e) Bond prices, $\psi(\tau,t)$','FontSize',font_size,'interpreter','latex');
hold on; grid on;
legend({'DSS','RSS','$\psi$($\tau$,0)'},...
    'Box','off','FontSize',12,'interpreter','latex');

subplot(3,2,6)
plot(tau,steady.v_ss(:),'Linewidth',line_width,'Color',color2,'LineStyle',LineStyle_2); hold on;
plot(tau,results.v_rss(:),'Linewidth',line_width,'LineStyle',LineStyle_4,'Color',color4);
plot(tau,results.v_n(:,1),'Linewidth',line_width,'Color',color3,'LineStyle',':');
xlabel('Maturity, $\tau$','FontSize',font_size,'interpreter','latex')
title('(f) Domestic valuations, $v(\tau,t)$','FontSize',font_size,'interpreter','latex');
hold on; grid on;
legend({'RSS','Default','$v$($\tau$,0)'},...
    'Box','off','FontSize',12,'interpreter','latex');

%% Main Plots Y IES 0
% Loading
file_RSS     = 'DEF';
file_DEF     = 'DEF_IES_0';

% Parameters
pre_t       = -10;
font_size   = 14;
line_width  =2;
base_color  = [0.2 0.2 0.6];

pre_t       = -10;
pos_t       = 60;

color_base =[0.1 0.1 0.6]; 
color1   =  [0.4 0.2 0.8];
color2   =  [0 0.4 0.1];
color3   =  [1 0 0];
color4   =  [1 0 0];
color_mat=  [color1; color2; color3; color4];
color_base2=[0.6 0.1 0.1];
color_rss=  [0.6 0.1 0.1];
color_rss_val=  [0.6 0.6 0.1];
markersize   = 40;
LineStyle_2 = '-.';
LineStyle_3 = '-.';
LineStyle_4 = '-.';


cd figures
load(file_RSS)
cd ..

% Variables
dt = parameters.dt;
t         = 0:parameters.dt:parameters.t_max;    % calendar time
Debt_rss = sum(results.f_rss(:,end))*dt;
Debt     = sum(results.f_n(:,:))*dt;
[average_duration, average_term_to_maturity] = Portfolio_duration(results, parameters, paths);


% %% Steady State - Comparisons
figure
subplot(321)
plot(tau,results.f_rss,'Linewidth',line_width,'Color',color2); hold on;
xlabel('Maturity, $\tau$','FontSize',font_size,'interpreter','latex')
title('(a) Debt distribution, $f(\tau)$','FontSize',font_size,'interpreter','latex')
hold on; grid on;

subplot(3,2,2)
plot(t,results.c_n,'Linewidth',line_width,'Color',color2); hold on;
scatter(0,results.c_n(1),markersize,color2,'o','filled'); 
scatter(t(end),results.c_rss,markersize,color2,'>','filled');
plot([pre_t 0],[results.c_rss results.c_rss],'Linewidth',line_width,'Color',color2);
scatter(0,results.c_rss,markersize,color2,'o');
plot(0,results.c_rss,'*','Linewidth',4,'Color',color_rss_val,'MarkerSize',2);
grid on; xlabel('Time, $t$ (years)','FontSize',font_size,'interpreter','latex');
title('(b) Consumption, $c(t)$','FontSize',font_size,'interpreter','latex');
xlim([pre_t pos_t]); 

subplot(3,2,3)
plot(t,Debt*100,'Linewidth',line_width,'Color',color2); hold on;
plot([pre_t 0],[Debt_rss Debt_rss]*100,'Linewidth',line_width,'Color',color2);
scatter(0,Debt_rss*100,markersize,color2,'o');
plot(0,Debt_rss*100,'*','Linewidth',line_width,'Color',color_rss_val,'MarkerSize',2);
scatter(0,Debt(1)*100,markersize,color2,'o','filled'); 
scatter(t(end),Debt(end)*100,markersize,color2,'>','filled'); 
xlabel('Time, $t$ (years)','FontSize',font_size,'interpreter','latex')
ylabel('$\%$ $y_{ss}$','FontSize',font_size,'interpreter','latex')
title('(c) Total debt, $b(t)$','FontSize',font_size,'interpreter','latex')
hold on; xlim([pre_t pos_t]); grid on;

subplot(3,2,4)
plot(t,average_duration,'Linewidth',line_width,'Color',color2); hold on;
plot([pre_t 0],[average_duration(1) average_duration(1)],'Linewidth',line_width,'Color',color2);
scatter(0,average_duration(1),markersize,color2,'o');
plot(0,average_duration(1),'*','Linewidth',line_width,'Color',color_rss_val,'MarkerSize',2);
scatter(0,average_duration(1),markersize,color2,'o','filled'); 
scatter(t(end),average_duration(end),markersize,color2,'>','filled'); 
xlabel('Time, $t$ (years)','FontSize',font_size,'interpreter','latex')
ylabel('Years','FontSize',font_size,'interpreter','latex')
title('(d) Average duration','FontSize',font_size,'interpreter','latex')
hold on; xlim([pre_t pos_t]); grid on;

cd figures
load(file_DEF)
cd ..
Debt_rss = sum(results.f_rss(:,end))*dt;
Debt     = sum(results.f_n(:,:))*dt;
[average_duration, average_term_to_maturity] = Portfolio_duration(results, parameters, paths);

subplot(321)
plot(tau,results.f_rss(:,end),'Linewidth',line_width,'Color',color4,'LineStyle',LineStyle_4); hold on;
legend({'Default','$\sigma = 0$' },...
    'Box','off','FontSize',12,'interpreter','latex');

subplot(3,2,2)
plot(t,results.c_n,'Linewidth',line_width,'Color',color4,'LineStyle',LineStyle_4); hold on;
scatter(0,results.c_n(1),markersize,color4,'o','filled'); 
scatter(t(end),results.c_n(end),markersize,color4,'>','filled');
plot([pre_t 0],[results.c_rss results.c_rss],'Linewidth',line_width,'Color',color4,'LineStyle',LineStyle_4);
scatter(0,results.c_rss,markersize,color4,'o');
plot(0,results.c_rss,'*','Linewidth',4,'Color',color_rss_val,'MarkerSize',2);

subplot(3,2,3)
plot(t,Debt*100,'Linewidth',line_width,'Color',color4,'LineStyle',LineStyle_4); hold on;
plot([pre_t 0],[Debt_rss Debt_rss]*100,'Linewidth',line_width,'Color',color4,'LineStyle',LineStyle_4);
scatter(0,Debt_rss*100,markersize,color4,'o');
plot(0,Debt_rss*100,'*','Linewidth',line_width,'Color',color_rss_val,'MarkerSize',2);
scatter(0,Debt(1)*100,markersize,color4,'o','filled'); 
scatter(t(end),Debt(end)*100,markersize,color4,'>','filled'); 

subplot(3,2,4)
plot(t,average_duration,'Linewidth',line_width,'Color',color4,'LineStyle',LineStyle_4); hold on;
plot([pre_t 0],[average_duration(1) average_duration(1)],'Linewidth',line_width,'Color',color4,'LineStyle',LineStyle_4);
scatter(0,average_duration(1),markersize,color4,'o');
plot(0,average_duration(1),'*','Linewidth',line_width,'Color',color_rss_val,'MarkerSize',2);
scatter(0,average_duration(1),markersize,color4,'o','filled'); 
scatter(t(end),average_duration(end),markersize,color4,'>','filled'); 

subplot(3,2,5)
plot(tau,steady.psi_ss(:),'Linewidth',line_width,'Color',color2,'LineStyle',LineStyle_2); hold on;
plot(tau,results.psi_rss(:),'Linewidth',line_width,'LineStyle',LineStyle_4,'Color',color4);
plot(tau,results.psi_n(:,1),'Linewidth',line_width,'Color',color3,'LineStyle',':');
xlabel('Maturity, $\tau$','FontSize',font_size,'interpreter','latex')
title('(e) Bond prices, $\psi(\tau,t)$','FontSize',font_size,'interpreter','latex');
hold on; grid on;
legend({'Default','$\sigma = 0$','$\psi$($\tau$,0), $\sigma = 0$'},...
    'Box','off','FontSize',12,'interpreter','latex');

subplot(3,2,6)
plot(tau,steady.v_ss(:),'Linewidth',line_width,'Color',color2,'LineStyle',LineStyle_2); hold on;
plot(tau,results.v_rss(:),'Linewidth',line_width,'LineStyle',LineStyle_4,'Color',color4);
plot(tau,results.v_n(:,1),'Linewidth',line_width,'Color',color3,'LineStyle',':');
xlabel('Maturity, $\tau$','FontSize',font_size,'interpreter','latex')
title('(f) Domestic valuations, $v(\tau,t)$','FontSize',font_size,'interpreter','latex');
hold on; grid on;
legend({'Default','$\sigma = 0$','$v$($\tau$,0), $\sigma = 0$'},...
    'Box','off','FontSize',12,'interpreter','latex');
