function Phi = phi_cdf(parameters,V,V_a)

% Galo, you are imposing a uniform solution
defmodel=parameters.defmodel;
% defmodel='exo';
switch defmodel
    case 'exo'
        b = parameters.V_low;
        a = parameters.V_high;
        if V<b
            PHI = 0;
        else if V>a
                PHI =1;
            else
                PHI = (V-b)/(a-b); 
            end
        end
    case 'probit'
        
    case 'logit'
        % parameters.muprobit=1;
        muprobit=10;
        Phi=exp(muprobit*V)./(exp(muprobit*V)+exp(muprobit*V_a));
end