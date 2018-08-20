function PHI = phi_cdf(parameters,V)
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
