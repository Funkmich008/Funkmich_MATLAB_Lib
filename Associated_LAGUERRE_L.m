% Berechnung der assoziierten Laguerre-Polynome 1. Art per Rekurrenz-Formel
function L = Associated_LAGUERRE_L(n, k, x)
    if n == 0
        L = 0 * x + 1;
    elseif n == 1
        L =  1 + k - x;
    else
       Ln2 =  1;
       Ln1 =  1 + k - x;
       for j = 2:n
           L = (2*j-1+k-x)/j .* Ln1 - (j+k-1)/j * Ln2;
           Ln2 = Ln1;
           Ln1 = L;
       end  
    end
end