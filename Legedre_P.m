% Berechnung der gewöhnlichen Legendre-Polynome 1. Art per Rekurrenz-Formel
function P = Legendre_P(n, x)

    if n == 0       
        P = 0 * x + 1;       
    elseif n == 1        
        P = x;       
    else        
        Pn1 = 1;            % Startwert: n = 0
        Pn = x;             % Startwert: n = 1
        for j = 2:n
            P = ((2*j-1) * x .* Pn - (j-1) * Pn1)./(j);
            Pn1 = Pn;
            Pn = P;            
        end
    end
end
