% Berechnung der gewöhnlichen Legendre-Polynome 2. Art per Rekurrenz-Formel
function Q = Legendre_Q(n, x)

    if n == 0       
        Q = atanh(x);        
    elseif n == 1        
        Q = x .* atanh(x) - 1;        
    else        
        Qn1 = atanh(x);             % Startwert: n = 0
        Qn = x .* atanh(x) - 1;     % Startwert: n = 1
        for j = 2:n
            Q = ((2*j-1) * x .* Qn - (j-1) * Qn1)./(j);
            Qn1 = Qn;
            Qn = Q;            
        end
    end
end
