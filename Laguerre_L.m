% Berechnung der Laguerre-Polynome 1. Art per Rekurrenz-Formel
function L = Laguerre_L(n, x)

    if n == 0
        L = 0 * x + 1;
    elseif n == 1
        L = 1 - x;
    else
       Ln2 = 1;
       Ln1 = 1 - x;
       for j = 2:n
           L = (2*j-1-x)/j .* Ln1 - (j-1)/j * Ln2;
           Ln2 = Ln1;
           Ln1 = L;
       end  
    end
end
