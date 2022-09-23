%% Bulirsch Complete Elliptic Integral of the first kind
% Paper: Numerical Calculation of Elliptic Integrals and Elliptic Functions III
% Autor: R. Bulirsch
% Algorithm by Bartky
function cel = BulirschElliptic_cel(k_c, p, a, b)

    mu_old = 1;
    nu_old = abs(k_c);
    a_old = a;
    b_old = b./sqrt(p);
    p_old = sqrt(p);   

    K = 20;
    for k = 0:K
              
       mu_new = nu_old + mu_old;
       nu_new = 2*sqrt(mu_old .* nu_old);
       p_new = p_old + (mu_old .* nu_old)./p_old;
       a_new = a_old + b_old./p_old;
       b_new = 2*(b_old + (mu_old .* nu_old)./p_old .* a_old);
        
       mu_old = mu_new;
       nu_old = nu_new;
       p_old  = p_new;
       a_old  = a_new;
       b_old  = b_new;
       
    end

    cel = pi/2 * (a_new .* mu_new + b_new)./(mu_new .* (mu_new + p_new));
    
end
