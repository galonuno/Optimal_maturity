%% Computes value at autarky
penalty  = parameters.penalty ;
y_n_penalty = y_n .* (1-penalty);
V_a    = compute_value_functional([y_n_penalty y_ss*ones(N_y,1)],parameters) ;
