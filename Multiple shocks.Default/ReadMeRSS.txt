RISKY STEADY STATE: Dynamics and Plots

1) main_rss. 

Parameters and plots. Calls the functions risky risky_SS_rss.m,bond_prices, bond_prices_rss, risky_SS.

2) risky_SS_rss.m defines the parameters, call functions bond_prices, bond_prices_rss, risky_SS. 

Este codigo es entero tuyo. 


3) bond_prices (parameters). 

Computes the dynamics of arbitrage free bond prices after a shock. 
The shock is permanent and expected, government by a AR(1) proces for the rates. parameters.r_bar_0 is the value 
after the shock. parameters.rho_r_bar is the speed of mean reversion. parameters.r_bar is the steady state short rate. 

Parameters: t_max, tau_max, dt, delta, r_bar,r_bar_0, rho_r_bar. 


4) bond_prices_rss (parameters). Computes the prices before the shock. This prices are a solution of a ODE with a JUMP. 


Parameters: r_bar,phi (prob shock), delta. 

5) Other codes. There are two codes that are called yield curve, that are just for plots. Not used.

