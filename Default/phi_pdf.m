function Phi_p = phi_pdf(parameters,V,V_a)
defmodel=parameters.defmodel;
muprobit=parameters.muprobit;

% defmodel='exo';
switch defmodel
    case 'exo'
        b = parameters.V_low;
        a = parameters.V_high;
        if V<b
            PHI_prime = 0;
        else if V>a
                PHI_prime = 0;
            else
                PHI_prime = 1/(a-b);
            end
        end
    case 'probit'
        
    case 'logit'
        %Phi_p=muprobit*(exp(muprobit*V)./(exp(muprobit*V)+exp(muprobit*V_a))-exp(muprobit*V)./(exp(muprobit*V)+exp(muprobit*V_a).^(2)));
        Phi_p=muprobit*(exp(muprobit*V).*exp(muprobit*V_a))./((exp(muprobit*V)+exp(muprobit*V_a)).^2);
end

