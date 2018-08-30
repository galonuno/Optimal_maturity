% Code to compute the Risky Steady State of the model with default
%(c) Bigio, Nuno & Passadore model
% This version codes Galo's solution using standard finite difference methods.
% The Code Solves for the Risky Steady State with default
close all; % clear;
plotit=1;

%% Code parameters
% Tolerance
parameters.tol_path   = 1e-3;  % tolerance parameter optimal path
parameters.tol        = 1e-5;  % tolerance parameter
parameters.max_iter   = 5000;  % maximum number of iterations
parameters.relax      = 0.1 ;  % relaxation coefficient
parameters.relax_path = 0.01;  % relaxation coefficient optimal path

% Numerical parameters
parameters.t_max     = 60    ; % maximum number of years (time)
parameters.tau_max   = 20    ; % maximum number of years (maturity)
parameters.dt        = 1/12  ; % monthly steps

% Method to Solve the model
parameters.solve_method='iter'; % iter for iterative, solver for solver

% Number of Shocks
N_r=   2;
N_y=   3;
parameters.N_r=N_r;
parameters.N_y=N_y;

%% Model parameters
% Time Period annual
% Preference Parameters
parameters.gamma   = 2.00     ; % risk aversion
parameters.rho     = 0.0416   ; % discount factor - Model is Quarterly

% Output parameters
parameters.y_ss         = 1.0    ; % output at steady state
parameters.y_0(1:N_y,1)   = [0.99; 1.0; 0.95 ]; % [s*] initial value, output after a shock
parameters.rho_y(1:N_y,1) = [0.25; 0.45; 0.22]; % [s*] persistence output after shock - half life 1/rho months
 
%Interest Rate Shock Parameters
parameters.r_bar_ss           = 0.04                ; % steady-state short rate
parameters.delta              = parameters.r_bar_ss ; % coupon = r_ss so we have bonds at par
parameters.r_bar_0(1:N_r,1)   = [0.06; 0.08]        ; % [s*] initial value [0.042 ; 0.043] 
parameters.rho_r_bar(1:N_r,1) = [0.25; 0.35]        ; % [s*] persistence (in years)

%Risky Steady State Parameters
parameters.phi  = 0.05       ; % Poisson arrival rate

% adjustment cost
parameters.lambda_bar = 7.0240  ; %

% Default parameters - Uniform distribution
parameters.defmodel='logit'     ; % default model
parameters.V_low = -30;
parameters.V_high = -30;
% parameters.V_low = -9;
% parameters.V_high = 1;

% Update to a matrix of joint probabilities - discrete matrix of marginals
% - sum of row and columns must add up to 1
parameters.prob_mat=zeros(N_r,N_y); % [s*]
parameters.prob_mat=[0.1 0.25 0.25; 0.3 0.05 0.05]; % [s*] % joint probability distribution, correlation is built in
if sum(sum(parameters.prob_mat))~=1;
    error('Prob matrix won''t add to one');
end

% Price of Risk
parameters.U         =  @(c) (c.^(1-parameters.gamma)-1)/(1-parameters.gamma);
parameters.U_p_ratio =  @(c_a,c_b) (c_a/c_b).^(-parameters.gamma);
parameters.dU        =  @(c) (c.^(-parameters.gamma));

%% Pre-Allocation
time_preallocate;

%% Paths
% Path for the Short Rate and output
exo_path_construct;

% Autarky Values
autarkyV_construct;
parameters.V_a=kron(ones(N_r,1),V_a)     ;

%% Save Exogenous Paths
paths.r_bar_n = r_bar_n;
paths.y_n     = y_n    ;

%% Step 1: Bond Prices in Steady State - index of r shock in columns
psi_n=zeros(N_tau,N_t,N_r); % [s*]
psi_0_mat=zeros(N_tau,N_r); % [s*]
for rr=1:N_r
    psi_n(:,:,rr)          = solve_HJB_path(r_bar_n(rr,:),r_bar_ss,parameters) ; % [s*]
    psi_0_mat(:,rr)        = psi_n(:,1,rr)                                     ; % [s*]
end
psi_ss         = solve_HJB_ss(r_bar_ss,parameters)                             ; % [s*]
rss.psi_0_mat  = psi_0_mat                                                     ; % [s*]
% results
results.psi_ss  = psi_ss                                                       ; % [s*]

%% Step 2: Value functions in Steady State 
steady=solve_steady(psi_ss,parameters)      ;
% Save Steady States
steady.psi_ss  = psi_ss                     ;
paths.psi_n    = psi_n                      ;

%% Step 2: Risky-Steady State
path_out        = solve_rss_path(rss,steady,paths,parameters); % [s*<-]
rss.psi_rss     = path_out.psi_rss  ;
results.psi_rss = path_out.psi_rss;
paths.psi_rss   = path_out.psi_rss             ;
results.yield_rss = parameters.delta-log(results.psi_rss)./(tau');

%% Step 3: Pack Results
results.f_n         = path_out.f_n          ;
results.f_rss       = path_out.f_out(:,1)   ;
results.iota_n      = path_out.iota_n       ;
results.iota_rss    = path_out.iota_out(:,1);
results.v_n         = path_out.v_n          ; 
results.v_rss       = path_out.v_out(:,1)   ; 
results.r_n         = path_out.r_n          ; 
results.c_rss       = path_out.c_out(:,1)   ; 
results.c_n         = path_out.c_n          ;

%% Main Plots
% Main Plot Params
pre_t       = -10;
font_size   = 20;
line_width  =4;
LineStyle_2 = ':';
base_color  = [0.2 0.2 0.6];

% Main Plots
pre_t       = -10;

% Parameters
dt = parameters.dt;

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
markersize   = 40;

% Ranges for debt levels
year_1=2;
year_2=10;
year_3=15;
% year_4=20;
range_1  = (1:year_1/dt);
range_5  = (year_1/dt+1:year_2/dt);
range_10 = (year_2/dt+1:N_tau);
% range_10 = (year_2*/dt+1:year_3*/dt);
% range_20 = (year_3*/dt+1:N_tau);

% Time T variables
t         = 0:parameters.dt:parameters.t_max;    % calendar time

% Build Results 
i_1 =  sum(results.iota_n(range_1,:))*dt;
i_5 = sum(results.iota_n(range_5,:))*dt;
i_10 = sum(results.iota_n(range_10,:))*dt;
% i_20 = sum(results.iota_n(range_20,:))*results.dt;
i_1_rss = sum(results.iota_rss(range_1,end))*dt;
i_5_rss = sum(results.iota_rss(range_5,end))*dt;
i_10_rss = sum(results.iota_rss(range_10,end))*dt;
% i_20_rss = sum(results.iota_rss(range_20,end))*results.dt;

% Building Densities
f_1 =  sum(results.f_n(range_1,:))*dt;
f_5 = sum(results.f_n(range_5,:))*dt;
f_10 = sum(results.f_n(range_10,:))*dt;
% f_20 = sum(results.f_n(range_20,:))*results.dt;
f_1_rss = sum(results.f_rss(range_1,end))*dt;
f_5_rss = sum(results.f_rss(range_5,end))*dt;
f_10_rss = sum(results.f_rss(range_10,end))*dt;

%% Plot the Yield Curve - RSS
figure(1)
subplot(2,2,1)
plot(tau,results.psi_rss,'Linewidth',line_width,'Color',color_base)
xlabel('Maturity $\tau$','FontSize',font_size,'interpreter','latex')
title('Price of Debt $\psi_{rss}(\tau)$','FontSize',font_size,'interpreter','latex')
hold on; axis tight; grid on; 

subplot(2,2,2)
plot(tau,results.yield_rss ,'Linewidth',line_width,'Color',color_base)
xlabel('Maturity $\tau$','FontSize',font_size,'interpreter','latex')
title('Yield Curve $\Psi_{rss}(\tau)$','FontSize',font_size,'interpreter','latex')
hold on; axis tight; grid on; 

subplot(2,2,3:4)
plot(t,paths.r_bar_n*100,'Linewidth',3,'Color',color_base); hold on;
plot([pre_t 0],[parameters.r_bar_ss parameters.r_bar_ss]*100,'Linewidth',3,'Color',color_base); axis tight;
scatter(0,parameters.r_bar_ss*100,markersize,2,'o','MarkerEdgeColor',color_base);
scatter(0,paths.r_bar_n(1)*100,markersize,2,'o','filled','MarkerEdgeColor',color_base,'MarkerFaceColor',color_base);
xlabel('Calendar Time $t$','FontSize',font_size,'interpreter','latex')
title('Short Rate $\bar{r}_t$ (bps)','FontSize',font_size,'interpreter','latex')
hold on; grid on; % ylim([0 max(parameters.r_bar_n)]);
axis tight;

%% Steady State - Comparisons
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
plot(tau,paths.psi_n(:,end),'Linewidth',line_width,'Color',base_color); hold on;
plot(tau,results.psi_rss(:),'Linewidth',line_width,'Color',color_rss); hold on;
plot(tau,results.v_n(:,end),'Linewidth',line_width,'Color',base_color,'LineStyle',LineStyle_2); 
plot(tau,results.v_rss(:),'Linewidth',line_width,'LineStyle',LineStyle_2,'Color',color_rss); 
plot(tau,results.v_n(:,1),'Linewidth',line_width,'Color',base_color/2,'LineStyle','-.');
xlabel('Maturity, $\tau$','FontSize',font_size,'interpreter','latex')
title('$\psi(\tau)$ vs. $v(\tau)$','FontSize',font_size,'interpreter','latex');
hold on; grid on;
legend({'$\psi(\tau,ss)$','$\psi(\tau,rss)$','$v(\tau,ss)$','$v(\tau,rss)$','$v(\tau,0)$'},...
    'Box','off','FontSize',14,'interpreter','latex');    

%% TRANSITIONAL DYNAMICS AFTER The SHOCK
figure(3)
subplot(2,2,1)
plot(t,results.r_n*100,'Linewidth',line_width,'Color',color_base); hold on;
plot([pre_t 0],[parameters.rho parameters.rho]*100,'Linewidth',3,'Color',color_base);
scatter(0,parameters.rho*100,markersize,color_base,'o');
plot(0,parameters.rho*100,'*','Linewidth',line_width,'Color',color_rss_val,'MarkerSize',2);
scatter(0,results.r_n(1)*100,markersize,color_base,'o','filled'); 
scatter(t(end),results.r_n(end)*100,markersize,color_base,'>','filled'); 
xlabel('Time, $t$','FontSize',font_size,'interpreter','latex')
ylabel('$\%$','FontSize',font_size,'interpreter','latex')
title('Time Discount $r(t)$ (bps)','FontSize',font_size,'interpreter','latex')
hold on; xlim([pre_t t(end)]); grid on;

subplot(2,2,2)
plot(t,results.c_n,'Linewidth',line_width,'Color',color_base); hold on;
plot(t,paths.y_n,'Linewidth',line_width,'Color',color_base2);
% scatter(0,results.c_ss,markersize,color_base,'o');
scatter(0,results.c_n(1),markersize,color_base,'o','filled'); 
scatter(t(end),results.c_n(end),markersize,color_base,'>','filled');
plot([pre_t 0],[results.c_rss results.c_rss],'Linewidth',3,'Color',color_base);
% plot([pre_t 0],[steady.c_ss steady.c_ss],'Linewidth',3,'Color',color_base);
% scatter(0,steady.y_ss,markersize,color_base2,'o');
scatter(0,paths.y_n(1),markersize,color_base2,'o','filled'); 
scatter(t(end),paths.y_n(end),markersize,color_base2,'>','filled');
plot([pre_t 0],[paths.y_n(end) paths.y_n(end)],'Linewidth',3,'Color',color_base2);
scatter(0,results.c_rss,markersize,color_base,'o');
plot(0,results.c_rss,'*','Linewidth',4,'Color',color_rss_val,'MarkerSize',2);
grid on; xlabel('Time, $t$','FontSize',font_size,'interpreter','latex');
title('$c(t)$ vs. $y(t)$','FontSize',font_size,'interpreter','latex');
xlim([pre_t t(end)]); legend({'$c(t)$','y(t)'},'interpreter','latex');

subplot(2,2,3);
plot(t,(i_1'),'Linewidth',2,'Color',color_mat(1,:)); hold on;
plot(t,(i_5'),'Linewidth',2,'Color',color_mat(2,:));
plot(t,(i_10'),'Linewidth',2,'Color',color_mat(3,:));
% plot(t,(i_20'),'Linewidth',2,'Color',color_mat(4,:));
plot([pre_t 0],[i_1_rss i_1_rss],'Linewidth',2,'Color',color_mat(1,:));
plot([pre_t 0],[i_5_rss i_5_rss],'Linewidth',2,'Color',color_mat(2,:));
plot([pre_t 0],[i_10_rss i_10_rss],'Linewidth',2,'Color',color_mat(3,:));
% plot([pre_t 0],[i_20_rss i_20_rss],'Linewidth',2,'Color',color_mat(4,:));
scatter(0,i_1(1),markersize,color_mat(1,:),'o','filled');
scatter(0,i_5(1),markersize,color_mat(2,:),'o','filled'); 
scatter(0,i_10(1),markersize,color_mat(3,:),'o','filled'); 
% scatter(0,i_20(1),markersize,color_mat(4,:),'o','filled');
scatter(t(end),i_1(end),markersize,color_mat(1,:),'>','filled');
scatter(t(end),i_5(end),markersize,color_mat(2,:),'>','filled'); 
scatter(t(end),i_10(end),markersize,color_mat(3,:),'>','filled'); 
% scatter(t(end),i_20(end),markersize,color_mat(4,:),'>','filled');
scatter(0,i_1_rss,markersize,color_mat(1,:),'o');  plot(0,i_1_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
scatter(0,i_5_rss,markersize,color_mat(2,:),'o');  plot(0,i_5_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
scatter(0,i_10_rss,markersize,color_mat(3,:),'o'); plot(0,i_10_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
% scatter(0,i_20_rss,markersize,color_mat(4,:),'o'); plot(0,i_20_rss,'*','Linewidth',4,'Color',color_rss_val,'MarkerSize',8);
scatter(0,i_1_rss,markersize,color_mat(1,:),'o');  plot(0,i_1_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
scatter(0,i_5_rss,markersize,color_mat(2,:),'o');  plot(0,i_5_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
scatter(0,i_10_rss,markersize,color_mat(3,:),'o'); plot(0,i_10_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
% scatter(0,i_20_rss,markersize,color_mat(4,:),'o'); plot(0,i_20_rss,'*','Linewidth',4,'Color',color_rss_val,'MarkerSize',8);
xlim([pre_t t(end)]);
% legend({'0-5 year','5-10 years','10-15 years','15 $\geq$ years'},'Interpreter','Latex','FontSize',font_size)
xlabel('Time, $t$','FontSize',font_size,'interpreter','latex')
ylabel('$\%$ $y_{ss}$','FontSize',font_size,'interpreter','latex')
title('Issuances Rate, $\iota(t)$','FontSize',font_size,'interpreter','latex')
grid on;

subplot(2,2,4)
% f_20_rss = sum(results.f_rss(range_20,end))*results.dt;
plot(t,(f_1'),'Linewidth',2,'Color',color_mat(1,:)); hold on;
plot(t,(f_5'),'Linewidth',2,'Color',color_mat(2,:));
plot(t,(f_10'),'Linewidth',2,'Color',color_mat(3,:));
% plot(t,(f_20'),'Linewidth',2,'Color',color_mat(4,:));
plot([pre_t 0],[f_1_rss f_1_rss],'Linewidth',2,'Color',color_mat(1,:));
plot([pre_t 0],[f_5_rss f_5_rss],'Linewidth',2,'Color',color_mat(2,:));
plot([pre_t 0],[f_10_rss f_10_rss],'Linewidth',2,'Color',color_mat(3,:));
% plot([pre_t 0],[f_20_rss f_20_rss],'Linewidth',2,'Color',color_mat(4,:));
scatter(0,f_1(1),markersize,color_mat(1,:),'o','filled');
scatter(0,f_5(1),markersize,color_mat(2,:),'o','filled'); 
scatter(0,f_10(1),markersize,color_mat(3,:),'o','filled'); 
% scatter(0,f_20(1),markersize,color_mat(4,:),'o','filled');
scatter(t(end),f_1(end),markersize,color_mat(1,:),'>','filled');
scatter(t(end),f_5(end),markersize,color_mat(2,:),'>','filled'); 
scatter(t(end),f_10(end),markersize,color_mat(3,:),'>','filled'); 
% scatter(t(end),f_20(end),markersize,color_mat(4,:),'>','filled');
scatter(0,f_1_rss,markersize,color_mat(1,:),'o');  plot(0,f_1_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
scatter(0,f_5_rss,markersize,color_mat(2,:),'o');  plot(0,f_5_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
scatter(0,f_10_rss,markersize,color_mat(3,:),'o'); plot(0,f_10_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
% scatter(0,f_20_rss,markersize,color_mat(4,:),'o'); plot(0,f_20_rss,'*','Linewidth',4,'Color',color_rss_val,'MarkerSize',8);
xlim([pre_t t(end)]);
  legend({'0-1 (years)','1-5 ','5 to 20'},...
        'Position',[0.332239224848658 -0.00773908315260428 0.373750358484675 0.0524046428167998],'Interpreter','Latex','FontSize',16,'Location','SouthEast','Box','off','Orientation','horizontal');

xlabel('Time, $t$','FontSize',font_size,'interpreter','latex')
ylabel('$\%$ $y_{ss}$','FontSize',font_size,'interpreter','latex')
title('Debt Outstanding, $f(t)$','FontSize',font_size,'interpreter','latex')
grid on;
