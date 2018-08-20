function PHI_prime = phi_pdf(parameters,V)
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