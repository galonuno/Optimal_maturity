% Code to compute the Risky Steady State
%(c) Bigio, Nuno & Passadore model
% Solves the frictionless case with only income shocks
close all; clear;
plotit=1;

%% Code parameters
% Tolerance
parameters.tol_path   = 1e-3;  % tolerance parameter optimal path
parameters.tol        = 1e-8;  % tolerance parameter
parameters.max_iter   = 5000;  % maximum number of iterations
parameters.relax      = 0.1 ;  % relaxation coefficient
parameters.relax_path = 0.01;  % relaxation coefficient optimal path

% Numerical parameters
parameters.t_max     = 1050    ; % maximum number of years (time)
parameters.tau_max   = 20    ; % maximum number of years (maturity)
parameters.dt        = 1/12  ; % monthly steps

% Method to Solve the model
parameters.solve_method='iter'; % iter for iterative, solver for solver

% Number of Shocks
N_r=1           ;
N_y=1           ;

%% Model parameters
% Time Period annual
% Preference Parameters
parameters.gamma   = 2.00     ; % risk aversion
parameters.rho     = 0.0545   ; % discount factor - Model is Quarterly

% Output parameters
parameters.y_ss         = 1.0    ; % output at steady state
parameters.y_0(1:N_y)   = 0.9    ; % initial value, output after a shock
parameters.rho_y(1:N_y) = 0.9975 ; % persistence output after shock - half life 1/rho months

%Interest Rate Shock Parameters
parameters.r_bar_ss           = 0.04                ; % steady-state short rate
parameters.delta              = parameters.r_bar_ss ; % coupon = r_ss so we have bonds at par
parameters.r_bar_0(1:N_r,1)   = parameters.r_bar_ss ; % initial value
parameters.rho_r_bar(1:N_r,1) = 1.00                ; % persistence (in months)

%Risky Steady State Parameters
parameters.phi  = 0.125       ; % Poisson arrival rate

% there's only one state
parameters.prob_vec=1         ;

% Price of Risk
if parameters.gamma~=1
    parameters.U=@(c) (c.^(1-parameters.gamma)-1)./(1-parameters.gamma);
else
    parameters.U=@(c) log(c);
end
parameters.U_p_ratio=@(c_a,c_b) (c_a/c_b).^(-parameters.gamma);
parameters.Upinv=@(x) (x).^(-1/parameters.gamma);

%% Pre-Allocation
time_preallocate;

%% Paths
% Path for the Short Rate
exo_path_construct;

%% Save Exogenous Paths
paths.r_bar_n = r_bar_n;
paths.y_n     = y_n    ;

%% Computations
% Step 1: Build NPV
% Unpack Parameters
dt     = parameters.dt          ;  
rho    = parameters.rho         ;  
gamma  = parameters.gamma       ;  
U_p_ratio= parameters.U_p_ratio ;  
phi= parameters.phi             ; 
U=parameters.U                  ;

% Unpack Paths
r_bar_n= paths.r_bar_n;
y_n    = paths.y_n    ;

% Defintions
npv=@(r_vec,z_vec) sum(z_vec.*exp(-cumsum([0 r_vec(1:end-1)].*dt)).*dt)   ; % NPV formula
npv_cg=@(r_vec) sum(ones(1,N_t).*exp(cumsum([0 r_vec(1:end-1)].*dt)).*dt) ; % NPV of constant
cg_d=@(r_vec) (r_vec*(1-gamma)-rho)/gamma                                 ; % Consumption Growth discounted

% Applying formulas
npv_y_n=npv(r_bar_n,y_n); % NPV income
cg_rate_n=cg_d(r_bar_n);  % Consumption Growth Rate - Discount Rate
npv_c=npv_cg(cg_rate_n);  % Chi function (NPV of C_t as function of C_0)

% Step 2: Back Consumption after shock
fc_rss=@(B) y_ss+r_bar_ss*B                                     ;
fc_0=@(B) (B+npv_y_n)/npv_c                                     ;
B_res=@(B) rho-(r_bar_ss+phi*(U_p_ratio(fc_0(B),fc_rss(B))-1))  ;
B_rss = fzero(@(B) B_res(B),0.0)                                ;
c_rss = fc_rss(B_rss);
c_0   = fc_0(B_rss);
c_n   = c_0*exp(cumsum((([0 r_bar_n(1:end-1)]-rho)/gamma).*dt)) ;
D_rss =-B_rss;

% Packaging the Result
results.c_n  = c_n;
results.c_rss= c_rss;

% Value after Shock
V_path=@(c_n) sum(U(c_n).*exp(-cumsum(rho*[0 1+0*c_n(1:end-1)].*dt)).*dt)  ;
V_0=V_path(c_n)                                                            ;
V_rss=(U(c_rss)+phi*V_0)/(rho+phi)                                         ;

%% Solve Full Value Function
% Value at new path
N_b=2000                    ; 
plotiter=1                  ;
U  = parameters.U           ;
N_t=length(t)               ;
Upinv=parameters.Upinv      ;

% Formulas
c_path=@(B_vec,r_vec) fc_0(B_vec).*(ones(1,N_t).*exp(cumsum([0 r_vec(1:end-1)].*dt)));
c_test=c_path(B_rss,((r_bar_n-rho)/gamma));
tol=parameters.tol    ;
gamma=parameters.gamma;

b_vec=linspace(-npv_y_n+0.1*npv_y_n, 0.1*y_ss/r_bar_ss, N_b)';
V    =U(fc_rss(b_vec))/rho;
V_0_b=zeros(N_b,1);
for bb=1:N_b
    V_0_b(bb)=V_path(c_path(b_vec(bb),((r_bar_n-rho)/gamma)));    
end
cond=2*parameters.tol;
sigma = zeros(N_b,1);
iter = 0;
% Initial Guess
dvF = zeros(N_b,1); dvB = zeros(N_b,1);
while cond>tol % value function iteration
    V_in =V;
    
    % Finite difference approximation
    dvF(1:end-1) = (V(2:end) - V(1:end-1))./(b_vec(2:end) - b_vec(1:end-1)); % Forward
    dvF(N_b)       = (y_ss + r_bar_ss*b_vec(end))^(-gamma);
    dvB(2:end)   = (V(2:end) - V(1:end-1))./(b_vec(2:end) - b_vec(1:end-1)); % Backward
    dvB(1)       = (y_ss + r_bar_ss*b_vec(1))^(-gamma);
    
    % Consumption Back and Forward
    cF = Upinv(dvF)           ; % First-order condition
    cB = Upinv(dvB)           ; % First-order condition
    
    % Update Drifts and Volatilities  
    % Assuming we work with risk-free technology at borrowing constraint | % Adopting Risky Technology Away from borrowing Constraint
    muF = r_bar_ss.*b_vec + (y_ss - cF)  ;
    muB = r_bar_ss.*b_vec + (y_ss - cB)  ;
    HcF       = U(cF) + dvF.*muF;      % forward Hamiltonian
    HcB       = U(cB) + dvB.*muB;      % backward Hamiltonian
    c0 = r_bar_ss.*b_vec + y_ss; 

    % Handling with non-convexities
    Ineither = (1-(muF>0)) .* (1-(muB<0));
    Iunique  = (muB<0).*(1-(muF>0)) + (1-(muB<0)).*(muF>0);
    Iboth    = (muB<0).*(muF>0);
    Ib       = Iunique.*(muB<0) + Iboth.*(HcB>=HcF);
    If       = Iunique.*(muF>0) + Iboth.*(HcF>=HcB);
    I0       = Ineither;   
    Ib(N_b) = 1; If(N_b) = 0; I0(N_b) = 0;        
    
    c_vec = cF.*If + cB.*Ib + c0.*I0;
      
    U_vec =U(c_vec)+phi*V_0_b;
    
    % Solve HJB_implicit
    V=HJB_implicit(U_vec,b_vec,muF,muB,sigma,rho+phi,dt,V,Ib,If);
    
    
    % Update condition
    cond=max(abs(V./V_in-1));
    
    % Plot if needed
    if plotiter==1
        if rand>0.95
            plot(b_vec,V); drawnow; hold on;
        end
    end
    iter=iter+1;
end
% Find Steady State
[~,index_min]=min(abs(-c_vec+r_bar_ss*b_vec+y_ss));
b_rss_im=b_vec(index_min);
c_rss_val=0;
disp(['Difference in implicit and direct method: ' num2str(abs(b_rss_im-B_rss)/B_rss)*100 '\% away']);

% Comparison

%% Main Plot
pre_t=-10;
font_size=20;
line_width=4;
LineStyle_2=':';
base_color=[0.2 0.2 0.6];

% Main Plots
pre_t=-round(0.1*parameters.t_max);

% Parameters
dt=parameters.dt;

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
markersize   = 60;

figure
plot(b_vec,V,'LineWidth',3); hold on; 
scatter(B_rss,V_rss,markersize,2,'MarkerFaceColor',color_rss_val);
axis tight; line([B_rss B_rss],[ylim],'LineWidth',3);

% Plot Solutions
figure
plot(b_vec,muF,'LineWidth',3); hold on; plot(b_vec,muB,'LineWidth',3); axis tight; grid on;
plot(b_vec,-c_vec+r_bar_ss*b_vec+y_ss,'LineWidth',3);
axis tight; line([B_rss B_rss],[ylim],'LineWidth',3)
[~,index_min]=min(abs(-c_vec+r_bar_ss*b_vec+y_ss));
b_rss_im=b_vec(index_min);
c_rss_val=0;
disp(['Difference in implicit and direct method: ' num2str(abs(b_rss_im-B_rss)/B_rss)]);

figure
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


% psi_n          = solve_HJB_path(r_bar_n,r_bar_ss,parameters) ;
% psi_0_mat      = psi_n(:,1)                                  ;
% psi_ss         = solve_HJB_ss(r_bar_ss,parameters)           ;
% psi_rss        = solve_HJB_rss(r_bar_ss,psi_0_mat,parameters);
% 
% % results
% results.psi_ss=psi_ss  ;
% results.psi_rss=psi_rss;
% % results.yield_rss = parameters.delta-log(results.psi_rss)./(tau'); % not
% % sure why coupon was here.
% results.yield_rss = 0-log(results.psi_rss)./(tau');
% 
% %% Step 2: Solve Steady State Values
% steady=solve_steady(psi_ss,parameters);
% 
% % Save Steady States
% steady.psi_ss=psi_ss;
% rss.psi_rss=psi_rss  ;
% paths.psi_n=psi_n   ;
% 
% % Saving Variables
% paths.psi_rss= psi_rss             ;
% if plotit==1
%    surf(t,tau,psi_n,'edgecolor','none','FaceAlpha',0.5); axis tight;
%    xlabel('t','Interpreter','latex'); ylabel('$\tau$','Interpreter','latex'); 
%    hold on;
%    plot3(1+tau*0,tau,psi_rss,'LineWidth',2,'Color','k','LineStyle','-');
% end
% 
% %% Step 2: Risky-Steady State
% path_out = solve_rss_path(rss,steady,paths,parameters);
% 
% %% Step 3: Pack Results
% results.f_n         = path_out.f_n          ;
% results.f_rss       = path_out.f_out(:,1)   ;
% results.iota_n      = path_out.iota_n       ;
% results.iota_rss    = path_out.iota_out(:,1);
% results.v_n         = path_out.v_n          ; 
% results.v_rss       = path_out.v_out(:,1)   ; 
% results.r_n         = path_out.r_n          ; 
% results.c_rss       = path_out.c_out(:,1)   ; 
% results.c_n         = path_out.c_n          ;
% 
% %% Main Plots

% 
% % Ranges for debt levels
% year_1=2;
% year_2=10;
% year_3=15;
% % year_4=20;
% range_1  = (1:year_1/dt);
% range_5  = (year_1/dt+1:year_2/dt);
% range_10 = (year_2/dt+1:N_tau);
% % range_10 = (year_2*/dt+1:year_3*/dt);
% % range_20 = (year_3*/dt+1:N_tau);
% 
% % Time T variables
% t         = 0:parameters.dt:parameters.t_max;    % calendar time
% 
% % Build Results 
% i_1 =  sum(results.iota_n(range_1,:))*dt;
% i_5 = sum(results.iota_n(range_5,:))*dt;
% i_10 = sum(results.iota_n(range_10,:))*dt;
% % i_20 = sum(results.iota_n(range_20,:))*results.dt;
% i_1_rss = sum(results.iota_rss(range_1,end))*dt;
% i_5_rss = sum(results.iota_rss(range_5,end))*dt;
% i_10_rss = sum(results.iota_rss(range_10,end))*dt;
% % i_20_rss = sum(results.iota_rss(range_20,end))*results.dt;
% 
% % Building Densities
% f_1 =  sum(results.f_n(range_1,:))*dt;
% f_5 = sum(results.f_n(range_5,:))*dt;
% f_10 = sum(results.f_n(range_10,:))*dt;
% % f_20 = sum(results.f_n(range_20,:))*results.dt;
% f_1_rss = sum(results.f_rss(range_1,end))*dt;
% f_5_rss = sum(results.f_rss(range_5,end))*dt;
% f_10_rss = sum(results.f_rss(range_10,end))*dt;
% 
% %% Plot the Yield Curve - RSS
% figure(1)
% subplot(2,2,1)
% plot(tau,results.psi_rss,'Linewidth',line_width,'Color',color_base)
% xlabel('Maturity $\tau$','FontSize',font_size,'interpreter','latex')
% title('Price of Debt $\psi_{rss}(\tau)$','FontSize',font_size,'interpreter','latex')
% hold on; axis tight; grid on; 
% 
% subplot(2,2,2)
% plot(tau,results.yield_rss ,'Linewidth',line_width,'Color',color_base)
% xlabel('Maturity $\tau$','FontSize',font_size,'interpreter','latex')
% title('Yield Curve $\Psi_{rss}(\tau)$','FontSize',font_size,'interpreter','latex')
% hold on; axis tight; grid on; 
% 
% subplot(2,2,3:4)
% plot(t,paths.r_bar_n*100,'Linewidth',3,'Color',color_base); hold on;
% plot([pre_t 0],[parameters.r_bar_ss parameters.r_bar_ss]*100,'Linewidth',3,'Color',color_base); axis tight;
% scatter(0,parameters.r_bar_ss*100,markersize,2,'o','MarkerEdgeColor',color_base);
% scatter(0,paths.r_bar_n(1)*100,markersize,2,'o','filled','MarkerEdgeColor',color_base,'MarkerFaceColor',color_base);
% xlabel('Calendar Time $t$','FontSize',font_size,'interpreter','latex')
% title('Short Rate $\bar{r}_t$ (bps)','FontSize',font_size,'interpreter','latex')
% hold on; grid on; % ylim([0 max(parameters.r_bar_n)]);
% axis tight;
% 
% %% Steady State - Comparisons
% figure(2)
% subplot(221)
% plot(tau,results.f_n(:,end),'Linewidth',line_width,'Color',base_color); hold on;
% plot(tau,results.f_rss(:),'Linewidth',line_width,'Color',color_rss)
% xlabel('Maturity, $\tau$','FontSize',font_size,'interpreter','latex')
% title('Outstanding Debt $f(\tau)$','FontSize',font_size,'interpreter','latex')
% hold on; grid on;
% subplot(222)
% plot(tau,results.iota_n(:,end),'Linewidth',line_width,'Color',base_color); hold on;
% plot(tau,results.iota_rss(:),'Linewidth',line_width,'Color',color_rss);
% xlabel('Maturity, $\tau$','FontSize',font_size,'interpreter','latex')
% title('Issuances $\iota(\tau)$','FontSize',font_size,'interpreter','latex')
% hold on; grid on;
% subplot(223)
% plot(tau,results.f_n(:,end)./abs(sum(results.f_n(:,end))),'Linewidth',line_width,'Color',base_color); hold on;
% plot(tau,results.f_rss(:)./abs(sum(results.f_rss(:))),'Linewidth',line_width,'Color',color_rss)      ;
% xlabel('Maturity, $\tau$','FontSize',font_size,'interpreter','latex')
% title('Maturity Distribution','FontSize',font_size,'interpreter','latex')
% hold on; grid on;
% subplot(224)
% plot(tau,paths.psi_n(:,end),'Linewidth',line_width,'Color',base_color); hold on;
% plot(tau,results.psi_rss(:),'Linewidth',line_width,'Color',color_rss); hold on;
% plot(tau,results.v_n(:,end),'Linewidth',line_width,'Color',base_color,'LineStyle',LineStyle_2); 
% plot(tau,results.v_rss(:),'Linewidth',line_width,'LineStyle',LineStyle_2,'Color',color_rss); 
% plot(tau,results.v_n(:,1),'Linewidth',line_width,'Color',base_color/2,'LineStyle','-.');
% xlabel('Maturity, $\tau$','FontSize',font_size,'interpreter','latex')
% title('$\psi(\tau)$ vs. $v(\tau)$','FontSize',font_size,'interpreter','latex');
% hold on; grid on;
% legend({'$\psi(\tau,ss)$','$\psi(\tau,rss)$','$v(\tau,ss)$','$v(\tau,rss)$','$v(\tau,0)$'},...
%     'Box','off','FontSize',14,'interpreter','latex');    
% 
% %% TRANSITIONAL DYNAMICS AFTER The SHOCK
% figure(3)
% subplot(2,2,1)
% plot(t,results.r_n*100,'Linewidth',line_width,'Color',color_base); hold on;
% plot([pre_t 0],[parameters.rho parameters.rho]*100,'Linewidth',3,'Color',color_base);
% scatter(0,parameters.rho*100,markersize,color_base,'o');
% plot(0,parameters.rho*100,'*','Linewidth',line_width,'Color',color_rss_val,'MarkerSize',2);
% scatter(0,results.r_n(1)*100,markersize,color_base,'o','filled'); 
% scatter(t(end),results.r_n(end)*100,markersize,color_base,'>','filled'); 
% xlabel('Time, $t$','FontSize',font_size,'interpreter','latex')
% ylabel('$\%$','FontSize',font_size,'interpreter','latex')
% title('Time Discount $r(t)$ (bps)','FontSize',font_size,'interpreter','latex')
% hold on; xlim([pre_t t(end)]); grid on;
% 
% subplot(2,2,2)
% plot(t,results.c_n,'Linewidth',line_width,'Color',color_base); hold on;
% plot(t,paths.y_n,'Linewidth',line_width,'Color',color_base2);
% % scatter(0,results.c_ss,markersize,color_base,'o');
% scatter(0,results.c_n(1),markersize,color_base,'o','filled'); 
% scatter(t(end),results.c_n(end),markersize,color_base,'>','filled');
% plot([pre_t 0],[results.c_rss results.c_rss],'Linewidth',3,'Color',color_base);
% % plot([pre_t 0],[steady.c_ss steady.c_ss],'Linewidth',3,'Color',color_base);
% % scatter(0,steady.y_ss,markersize,color_base2,'o');
% scatter(0,paths.y_n(1),markersize,color_base2,'o','filled'); 
% scatter(t(end),paths.y_n(end),markersize,color_base2,'>','filled');
% plot([pre_t 0],[paths.y_n(end) paths.y_n(end)],'Linewidth',3,'Color',color_base2);
% scatter(0,results.c_rss,markersize,color_base,'o');
% plot(0,results.c_rss,'*','Linewidth',4,'Color',color_rss_val,'MarkerSize',2);
% grid on; xlabel('Time, $t$','FontSize',font_size,'interpreter','latex');
% title('$c(t)$ vs. $y(t)$','FontSize',font_size,'interpreter','latex');
% xlim([pre_t t(end)]); legend({'$c(t)$','y(t)'},'interpreter','latex');
% 
% subplot(2,2,3);
% plot(t,(i_1'),'Linewidth',2,'Color',color_mat(1,:)); hold on;
% plot(t,(i_5'),'Linewidth',2,'Color',color_mat(2,:));
% plot(t,(i_10'),'Linewidth',2,'Color',color_mat(3,:));
% % plot(t,(i_20'),'Linewidth',2,'Color',color_mat(4,:));
% plot([pre_t 0],[i_1_rss i_1_rss],'Linewidth',2,'Color',color_mat(1,:));
% plot([pre_t 0],[i_5_rss i_5_rss],'Linewidth',2,'Color',color_mat(2,:));
% plot([pre_t 0],[i_10_rss i_10_rss],'Linewidth',2,'Color',color_mat(3,:));
% % plot([pre_t 0],[i_20_rss i_20_rss],'Linewidth',2,'Color',color_mat(4,:));
% scatter(0,i_1(1),markersize,color_mat(1,:),'o','filled');
% scatter(0,i_5(1),markersize,color_mat(2,:),'o','filled'); 
% scatter(0,i_10(1),markersize,color_mat(3,:),'o','filled'); 
% % scatter(0,i_20(1),markersize,color_mat(4,:),'o','filled');
% scatter(t(end),i_1(end),markersize,color_mat(1,:),'>','filled');
% scatter(t(end),i_5(end),markersize,color_mat(2,:),'>','filled'); 
% scatter(t(end),i_10(end),markersize,color_mat(3,:),'>','filled'); 
% % scatter(t(end),i_20(end),markersize,color_mat(4,:),'>','filled');
% scatter(0,i_1_rss,markersize,color_mat(1,:),'o');  plot(0,i_1_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
% scatter(0,i_5_rss,markersize,color_mat(2,:),'o');  plot(0,i_5_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
% scatter(0,i_10_rss,markersize,color_mat(3,:),'o'); plot(0,i_10_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
% % scatter(0,i_20_rss,markersize,color_mat(4,:),'o'); plot(0,i_20_rss,'*','Linewidth',4,'Color',color_rss_val,'MarkerSize',8);
% scatter(0,i_1_rss,markersize,color_mat(1,:),'o');  plot(0,i_1_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
% scatter(0,i_5_rss,markersize,color_mat(2,:),'o');  plot(0,i_5_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
% scatter(0,i_10_rss,markersize,color_mat(3,:),'o'); plot(0,i_10_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
% % scatter(0,i_20_rss,markersize,color_mat(4,:),'o'); plot(0,i_20_rss,'*','Linewidth',4,'Color',color_rss_val,'MarkerSize',8);
% xlim([pre_t t(end)]);
% % legend({'0-5 year','5-10 years','10-15 years','15 $\geq$ years'},'Interpreter','Latex','FontSize',font_size)
% xlabel('Time, $t$','FontSize',font_size,'interpreter','latex')
% ylabel('$\%$ $y_{ss}$','FontSize',font_size,'interpreter','latex')
% title('Issuances Rate, $\iota(t)$','FontSize',font_size,'interpreter','latex')
% grid on;
% 
% subplot(2,2,4)
% % f_20_rss = sum(results.f_rss(range_20,end))*results.dt;
% plot(t,(f_1'),'Linewidth',2,'Color',color_mat(1,:)); hold on;
% plot(t,(f_5'),'Linewidth',2,'Color',color_mat(2,:));
% plot(t,(f_10'),'Linewidth',2,'Color',color_mat(3,:));
% % plot(t,(f_20'),'Linewidth',2,'Color',color_mat(4,:));
% plot([pre_t 0],[f_1_rss f_1_rss],'Linewidth',2,'Color',color_mat(1,:));
% plot([pre_t 0],[f_5_rss f_5_rss],'Linewidth',2,'Color',color_mat(2,:));
% plot([pre_t 0],[f_10_rss f_10_rss],'Linewidth',2,'Color',color_mat(3,:));
% % plot([pre_t 0],[f_20_rss f_20_rss],'Linewidth',2,'Color',color_mat(4,:));
% scatter(0,f_1(1),markersize,color_mat(1,:),'o','filled');
% scatter(0,f_5(1),markersize,color_mat(2,:),'o','filled'); 
% scatter(0,f_10(1),markersize,color_mat(3,:),'o','filled'); 
% % scatter(0,f_20(1),markersize,color_mat(4,:),'o','filled');
% scatter(t(end),f_1(end),markersize,color_mat(1,:),'>','filled');
% scatter(t(end),f_5(end),markersize,color_mat(2,:),'>','filled'); 
% scatter(t(end),f_10(end),markersize,color_mat(3,:),'>','filled'); 
% % scatter(t(end),f_20(end),markersize,color_mat(4,:),'>','filled');
% scatter(0,f_1_rss,markersize,color_mat(1,:),'o');  plot(0,f_1_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
% scatter(0,f_5_rss,markersize,color_mat(2,:),'o');  plot(0,f_5_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
% scatter(0,f_10_rss,markersize,color_mat(3,:),'o'); plot(0,f_10_rss,'*','Linewidth',2,'Color',color_rss_val,'MarkerSize',2);
% % scatter(0,f_20_rss,markersize,color_mat(4,:),'o'); plot(0,f_20_rss,'*','Linewidth',4,'Color',color_rss_val,'MarkerSize',8);
% xlim([pre_t t(end)]);
%   legend({'0-1 (years)','1-5 ','5 to 20'},...
%         'Position',[0.332239224848658 -0.00773908315260428 0.373750358484675 0.0524046428167998],'Interpreter','Latex','FontSize',16,'Location','SouthEast','Box','off','Orientation','horizontal');
% 
% xlabel('Time, $t$','FontSize',font_size,'interpreter','latex')
% ylabel('$\%$ $y_{ss}$','FontSize',font_size,'interpreter','latex')
% title('Debt Outstanding, $f(t)$','FontSize',font_size,'interpreter','latex')
% grid on;
