%% Main Plots
% Main Plot Params
pre_t       = -10;
font_size   = 14;
line_width  =2;
base_color  = [0.2 0.2 0.6];

% Main Plots
pre_t       = -10;
pos_t       = 40;

% Parameters
dt = parameters.dt;

%% Plots for the Transitional Dynamics
color_base =[0 0 1]; 
color1   =  [0.4 0.2 0.8];
color2   =  [0.2 0.6 0.8];
color3   =  [0.2 0.2 0.8];
color4   =  [0.8 0.2 0.8];
color_mat=  [color1; color2; color3; color4];
color_base2=[0 0.4 0.1];
color_rss=  [0.6 0.1 0.1];
color_rss_val=  [0.6 0.6 0.1];
markersize   = 40;
LineStyle_2 = '-.';
LineStyle_3 = '--';
LineStyle_4 = ':';


% Ranges for debt levels
year_1=5;
year_2=10;
year_3=15;
%year_4=20;
range_1  = (1:year_1/dt);
range_5  = (year_1/dt+1:year_2/dt);
range_10 = (year_2/dt+1:year_3/dt);
range_20 = (year_3/dt+1:N_tau);

% Time T variables
t         = 0:parameters.dt:parameters.t_max;    % calendar time

% Build Results 
i_1 =  sum(results.iota_n(range_1,:))*dt;
i_5 = sum(results.iota_n(range_5,:))*dt;
i_10 = sum(results.iota_n(range_10,:))*dt;
i_20 = sum(results.iota_n(range_20,:))*dt;
i_1_rss = sum(results.iota_rss(range_1,end))*dt;
i_5_rss = sum(results.iota_rss(range_5,end))*dt;
i_10_rss = sum(results.iota_rss(range_10,end))*dt;
i_20_rss = sum(results.iota_rss(range_20,end))*dt;

% Building Densities
f_1 =  sum(results.f_n(range_1,:))*dt;
f_5 = sum(results.f_n(range_5,:))*dt;
f_10 = sum(results.f_n(range_10,:))*dt;
f_20 = sum(results.f_n(range_20,:))*dt;
f_1_rss = sum(results.f_rss(range_1,end))*dt;
f_5_rss = sum(results.f_rss(range_5,end))*dt;
f_10_rss = sum(results.f_rss(range_10,end))*dt;
f_20_rss = sum(results.f_rss(range_20,end))*dt;

Debt_rss = sum(results.f_rss(:,end))*dt;
Debt     = sum(results.f_n(:,:))*dt;

% Duration
[average_duration, average_term_to_maturity] = Portfolio_duration(results, parameters, paths);


% %% Plot the Yield Curve - RSS
% figure(1)
% subplot(2,2,1)
% plot(tau,results.psi_rss,'Linewidth',line_width,'Color',color_base)
% xlabel('Maturity $\tau$','FontSize',font_size,'interpreter','latex')
% title('Price of debt $\psi_{rss}(\tau)$','FontSize',font_size,'interpreter','latex')
%  
% 
% subplot(2,2,2)
% plot(tau,results.yield_rss ,'Linewidth',line_width,'Color',color_base)
% xlabel('Maturity $\tau$','FontSize',font_size,'interpreter','latex')
% title('Yield curve $\Psi_{rss}(\tau)$','FontSize',font_size,'interpreter','latex')
% hold on; axis tight; grid on; 
% 
% subplot(2,2,3:4)
% plot(t,paths.r_bar_n*100,'Linewidth',3,'Color',color_base); hold on;
% plot([pre_t 0],[parameters.r_bar_ss parameters.r_bar_ss]*100,'Linewidth',3,'Color',color_base); axis tight;
% scatter(0,parameters.r_bar_ss*100,markersize,2,'o','MarkerEdgeColor',color_base);
% scatter(0,paths.r_bar_n(1)*100,markersize,2,'o','filled','MarkerEdgeColor',color_base,'MarkerFaceColor',color_base);
% xlabel('Time $t$ (years)','FontSize',font_size,'interpreter','latex')
% title('Short rate $\bar{r}_t$ (bps)','FontSize',font_size,'interpreter','latex')
% hold on; grid on; % ylim([0 max(parameters.r_bar_n)]);
% axis tight;
% 
% %% Steady State - Comparisons
% figure
% if parameters.mode == 'PF'
%     subplot(221)
%     plot(tau,results.f_n(:,end),'Linewidth',line_width,'Color',base_color); hold on;
%     xlabel('Maturity, $\tau$','FontSize',font_size,'interpreter','latex')
%     title('Outstanding debt, $f(\tau)$','FontSize',font_size,'interpreter','latex')
%     hold on; grid on;
%     subplot(222)
%     plot(tau,results.iota_n(:,end),'Linewidth',line_width,'Color',base_color); hold on;
%     xlabel('Maturity, $\tau$','FontSize',font_size,'interpreter','latex')
%     title('Issuances $\iota(\tau)$','FontSize',font_size,'interpreter','latex')
%     hold on; grid on;
%     subplot(223)
%     plot(tau,results.f_n(:,end)./abs(sum(results.f_n(:,end))*dt),'Linewidth',line_width,'Color',base_color); hold on;
%     xlabel('Maturity, $\tau$','FontSize',font_size,'interpreter','latex')
%     title('Maturity dist., $f(\tau) / \int_0^T f(\tau) d\tau$','FontSize',font_size,'interpreter','latex')
%     hold on; grid on;
%     subplot(224)
%     plot(tau,paths.psi_n(:,end),'Linewidth',line_width,'Color',base_color); hold on;
%     plot(tau,results.v_n(:,end),'Linewidth',line_width,'Color',base_color,'LineStyle',LineStyle_2);
%     xlabel('Maturity, $\tau$','FontSize',font_size,'interpreter','latex')
%     title('$\psi(\tau)$ vs. $v(\tau)$','FontSize',font_size,'interpreter','latex');
%     hold on; grid on;
%     legend({'$\psi(\tau)$','$v(\tau)$'},...
%         'Box','off','FontSize',14,'interpreter','latex');
% else
%     subplot(221)
%     plot(tau,results.f_n(:,end),'Linewidth',line_width,'Color',base_color); hold on;
%     plot(tau,results.f_rss(:),'Linewidth',line_width,'Color',color_rss)
%     xlabel('Maturity, $\tau$','FontSize',font_size,'interpreter','latex')
%     title('Outstanding debt $f(\tau)$','FontSize',font_size,'interpreter','latex')
%     hold on; grid on;
%     subplot(222)
%     plot(tau,results.iota_n(:,end),'Linewidth',line_width,'Color',base_color); hold on;
%     plot(tau,results.iota_rss(:),'Linewidth',line_width,'Color',color_rss);
%     xlabel('Maturity, $\tau$','FontSize',font_size,'interpreter','latex')
%     title('Issuances $\iota(\tau)$','FontSize',font_size,'interpreter','latex')
%     hold on; grid on;
%     subplot(223)
%     plot(tau,results.f_n(:,end)./abs(sum(results.f_n(:,end))),'Linewidth',line_width,'Color',base_color); hold on;
%     plot(tau,results.f_rss(:)./abs(sum(results.f_rss(:))),'Linewidth',line_width,'Color',color_rss)      ;
%     xlabel('Maturity, $\tau$','FontSize',font_size,'interpreter','latex')
%     title('Maturity distribution','FontSize',font_size,'interpreter','latex')
%     hold on; grid on;
%     subplot(224)
%     plot(tau,paths.psi_n(:,end),'Linewidth',line_width,'Color',base_color); hold on;
%     plot(tau,results.psi_rss(:),'Linewidth',line_width,'Color',color_rss); hold on;
%     plot(tau,results.v_n(:,end),'Linewidth',line_width,'Color',base_color,'LineStyle',LineStyle_2);
%     plot(tau,results.v_rss(:),'Linewidth',line_width,'LineStyle',LineStyle_2,'Color',color_rss);
%     plot(tau,results.v_n(:,1),'Linewidth',line_width,'Color',base_color/2,'LineStyle','-.');
%     xlabel('Maturity, $\tau$','FontSize',font_size,'interpreter','latex')
%     title('$\psi(\tau)$ vs. $v(\tau)$','FontSize',font_size,'interpreter','latex');
%     hold on; grid on;
%     legend({'$\psi(\tau,ss)$','$\psi(\tau,rss)$','$v(\tau,ss)$','$v(\tau,rss)$','$v(\tau,0)$'},...
%         'Box','off','FontSize',14,'interpreter','latex');
% end

%% TRANSITIONAL DYNAMICS AFTER The SHOCK
figure
subplot(2,3,1)
plot(t,results.r_n*100,'Linewidth',line_width,'Color',color_base); hold on;
plot(t,paths.r_bar_n*100,'Linewidth',line_width,'Color',color_base2,'LineStyle',LineStyle_2); 
scatter(0,parameters.rho*100,markersize,color_base,'o');
scatter(0,parameters.r_bar_ss *100,markersize,color_base2,'o');
plot([pre_t 0],[parameters.rho parameters.rho]*100,'Linewidth',line_width,'Color',color_base);
scatter(t(end),results.r_n(end)*100,markersize,color_base,'>','filled'); 
plot([pre_t 0],[parameters.r_bar_ss  parameters.r_bar_ss ]*100,'Linewidth',line_width,'Color',color_base2,'LineStyle',LineStyle_2);
scatter(t(end),paths.r_bar_n(end)*100,markersize,color_base,'>','filled'); 
plot(0,parameters.rho*100,'*','Linewidth',line_width,'Color',color_rss_val,'MarkerSize',2);
scatter(0,results.r_n(1)*100,markersize,color_base,'o','filled'); 
plot(0,parameters.r_bar_ss *100,'*','Linewidth',line_width,'Color',color_rss_val,'MarkerSize',2);
scatter(0,paths.r_bar_n(1)*100,markersize,color_base2,'o','filled'); 
title('(a) Discount, $r(t)$, and rates, $\bar{r}(t)$','FontSize',font_size,'interpreter','latex')
xlabel('Time, $t$ (years)','FontSize',font_size,'interpreter','latex')
ylabel('$\%$','FontSize',font_size,'interpreter','latex')
xlim([pre_t pos_t]); grid on;
legend({'$r(t)$','$\bar{r}(t)$'},'interpreter','latex','FontSize',font_size,...
    'Location','Best','Box','off','Orientation','horizontal');

subplot(2,3,2)
plot(t,results.c_n,'Linewidth',line_width,'Color',color_base); hold on;
plot(t,paths.y_n,'Linewidth',line_width,'Color',color_base2,'LineStyle',LineStyle_2);
scatter(0,results.c_n(1),markersize,color_base,'o','filled'); 
scatter(0,paths.y_n(1),markersize,color_base2,'o','filled'); 
scatter(t(end),results.c_n(end),markersize,color_base,'>','filled');
plot([pre_t 0],[results.c_rss results.c_rss],'Linewidth',line_width,'Color',color_base);
scatter(t(end),paths.y_n(end),markersize,color_base2,'>','filled');
plot([pre_t 0],[paths.y_n(end) paths.y_n(end)],'Linewidth',line_width,'Color',color_base2,'LineStyle',LineStyle_2);
scatter(0,results.c_rss,markersize,color_base,'o');
scatter(0,paths.y_n(end),markersize,color_base2,'o');
plot(0,results.c_rss,'*','Linewidth',4,'Color',color_rss_val,'MarkerSize',2);
plot(0,paths.y_n(end),'*','Linewidth',4,'Color',color_rss_val,'MarkerSize',2);
grid on; xlabel('Time, $t$ (years)','FontSize',font_size,'interpreter','latex');
title('(b) Consumption, $c(t)$, and output, $y(t)$','FontSize',font_size,'interpreter','latex');
xlim([pre_t pos_t]); 
legend({'$c(t)$','y(t)'},'interpreter','latex','FontSize',font_size,...
    'Location','Best','Box','off','Orientation','horizontal');

subplot(2,3,3)
plot(t,Debt*100,'Linewidth',line_width,'Color',color_base); hold on;
plot([pre_t 0],[Debt_rss Debt_rss]*100,'Linewidth',line_width,'Color',color_base);
scatter(0,Debt_rss*100,markersize,color_base,'o');
plot(0,Debt_rss*100,'*','Linewidth',line_width,'Color',color_rss_val,'MarkerSize',2);
scatter(0,Debt(1)*100,markersize,color_base,'o','filled'); 
scatter(t(end),Debt(end)*100,markersize,color_base,'>','filled'); 

xlabel('Time, $t$ (years)','FontSize',font_size,'interpreter','latex')
ylabel('$\%$ $y_{ss}$','FontSize',font_size,'interpreter','latex')
title('(c) Total debt, $b(t)$','FontSize',font_size,'interpreter','latex')
hold on; xlim([pre_t pos_t]); grid on;

subplot(2,3,4)
plot(t,average_duration,'Linewidth',line_width,'Color',color_base); hold on;

plot([pre_t 0],[average_duration(1) average_duration(1)],'Linewidth',line_width,'Color',color_base);
scatter(0,average_duration(1),markersize,color_base,'o');
plot(0,average_duration(1),'*','Linewidth',line_width,'Color',color_rss_val,'MarkerSize',2);
scatter(0,average_duration(1),markersize,color_base,'o','filled'); 
scatter(t(end),average_duration(end),markersize,color_base,'>','filled'); 
xlabel('Time, $t$ (years)','FontSize',font_size,'interpreter','latex')
ylabel('Years','FontSize',font_size,'interpreter','latex')
title('(d) Average duration','FontSize',font_size,'interpreter','latex')
hold on; xlim([pre_t pos_t]); grid on;


subplot(2,3,5);
plot(t,(i_1'),'Linewidth',2,'Color',color_mat(1,:)); hold on;
plot(t,(i_5'),'Linewidth',2,'Color',color_mat(2,:),'LineStyle',LineStyle_2);
plot(t,(i_10'),'Linewidth',2,'Color',color_mat(3,:),'LineStyle',LineStyle_3);
 plot(t,(i_20'),'Linewidth',2,'Color',color_mat(4,:),'LineStyle',LineStyle_4);
plot([pre_t 0],[i_1_rss i_1_rss],'Linewidth',2,'Color',color_mat(1,:));
plot([pre_t 0],[i_5_rss i_5_rss],'Linewidth',2,'Color',color_mat(2,:),'LineStyle',LineStyle_2);
plot([pre_t 0],[i_10_rss i_10_rss],'Linewidth',2,'Color',color_mat(3,:),'LineStyle',LineStyle_3);
 plot([pre_t 0],[i_20_rss i_20_rss],'Linewidth',2,'Color',color_mat(4,:),'LineStyle',LineStyle_4);
scatter(0,i_1(1),markersize,color_mat(1,:),'o','filled');
scatter(0,i_5(1),markersize,color_mat(2,:),'o','filled'); 
scatter(0,i_10(1),markersize,color_mat(3,:),'o','filled'); 
scatter(0,i_20(1),markersize,color_mat(4,:),'o','filled');
scatter(t(end),i_1(end),markersize,color_mat(1,:),'>','filled');
scatter(t(end),i_5(end),markersize,color_mat(2,:),'>','filled'); 
scatter(t(end),i_10(end),markersize,color_mat(3,:),'>','filled'); 
 scatter(t(end),i_20(end),markersize,color_mat(4,:),'>','filled');
scatter(0,i_1_rss,markersize,color_mat(1,:),'o');  plot(0,i_1_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
scatter(0,i_5_rss,markersize,color_mat(2,:),'o');  plot(0,i_5_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
scatter(0,i_10_rss,markersize,color_mat(3,:),'o'); plot(0,i_10_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
 scatter(0,i_20_rss,markersize,color_mat(4,:),'o'); plot(0,i_20_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
scatter(0,i_1_rss,markersize,color_mat(1,:),'o');  plot(0,i_1_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
scatter(0,i_5_rss,markersize,color_mat(2,:),'o');  plot(0,i_5_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
scatter(0,i_10_rss,markersize,color_mat(3,:),'o'); plot(0,i_10_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
 scatter(0,i_20_rss,markersize,color_mat(4,:),'o'); plot(0,i_20_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
xlim([pre_t pos_t]);
% legend({'0-5 year','5-10 years','10-15 years','15 $\geq$ years'},'Interpreter','Latex','FontSize',font_size)
xlabel('Time, $t$ (years)','FontSize',font_size,'interpreter','latex')
ylabel('$\%$ $y_{ss}$','FontSize',font_size,'interpreter','latex')
title('(e) Issuances rate, $\iota(t)$','FontSize',font_size,'interpreter','latex')
grid on;

subplot(2,3,6)
plot(t,(f_1'),'Linewidth',2,'Color',color_mat(1,:)); hold on;
plot(t,(f_5'),'Linewidth',2,'Color',color_mat(2,:),'LineStyle',LineStyle_2);
plot(t,(f_10'),'Linewidth',2,'Color',color_mat(3,:),'LineStyle',LineStyle_3);
plot(t,(f_20'),'Linewidth',2,'Color',color_mat(4,:),'LineStyle',LineStyle_4);
plot([pre_t 0],[f_1_rss f_1_rss],'Linewidth',2,'Color',color_mat(1,:));
plot([pre_t 0],[f_5_rss f_5_rss],'Linewidth',2,'Color',color_mat(2,:),'LineStyle',LineStyle_2);
plot([pre_t 0],[f_10_rss f_10_rss],'Linewidth',2,'Color',color_mat(3,:),'LineStyle',LineStyle_3);
plot([pre_t 0],[f_20_rss f_20_rss],'Linewidth',2,'Color',color_mat(4,:),'LineStyle',LineStyle_4);
scatter(0,f_1(1),markersize,color_mat(1,:),'o','filled');
scatter(0,f_5(1),markersize,color_mat(2,:),'o','filled'); 
scatter(0,f_10(1),markersize,color_mat(3,:),'o','filled'); 
scatter(0,f_20(1),markersize,color_mat(4,:),'o','filled');
scatter(t(end),f_1(end),markersize,color_mat(1,:),'>','filled');
scatter(t(end),f_5(end),markersize,color_mat(2,:),'>','filled'); 
scatter(t(end),f_10(end),markersize,color_mat(3,:),'>','filled'); 
scatter(t(end),f_20(end),markersize,color_mat(4,:),'>','filled');
scatter(0,f_1_rss,markersize,color_mat(1,:),'o');  plot(0,f_1_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
scatter(0,f_5_rss,markersize,color_mat(2,:),'o');  plot(0,f_5_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
scatter(0,f_10_rss,markersize,color_mat(3,:),'o'); plot(0,f_10_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
scatter(0,f_20_rss,markersize,color_mat(4,:),'o'); plot(0,f_20_rss,'*','Linewidth',4,'Color',color_rss_val,'MarkerSize',2);
xlim([pre_t pos_t]);
  legend({'0-5 (years)','5-10 ','10-15','15-20'},...
        'Position',[0.332239224848658 -0.00773908315260428 0.373750358484675 0.0524046428167998],'Interpreter','Latex','FontSize',font_size,'Location','SouthEast','Box','off','Orientation','horizontal');

xlabel('Time, $t$ (years)','FontSize',font_size,'interpreter','latex')
ylabel('$\%$ $y_{ss}$','FontSize',font_size,'interpreter','latex')
title('(f) Debt outstanding, $f(t)$','FontSize',font_size,'interpreter','latex')
grid on;
