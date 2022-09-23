function RJ = CarlsonElliptic_RJ(x_0, y_0, z_0, p_0)
    
    A_0 = (x_0 + y_0 + z_0 + 2*p_0)/5;
    
    x = x_0;
    y = y_0;
    z = z_0;
    p = p_0;
    A = A_0; 
    
    delta = (p-x) .* (p-y) .* (p-z);
    
    K = 25;
    V = 0;
    for k = 0:K
       
        lambda = sqrt(x.*y) ...
               + sqrt(y.*z) ...
               + sqrt(x.*z);
        
        sp = sqrt(p);
        d = (sp + sqrt(x)) .* (sp + sqrt(y)) .* (sp + sqrt(z));
        e = 4^(-3*k)./d.^2 .* delta;   
        V = V + 4^(-k)./d .* Carlson_RC(1, 1+e);
        
        x = 1/4 * (x + lambda);
        y = 1/4 * (y + lambda);
        z = 1/4 * (z + lambda);
        p = 1/4 * (p + lambda);
        A = 1/4 * (A + lambda);
           
    end
    
    X = (A_0 - x_0)./(4^K .* A);
    Y = (A_0 - y_0)./(4^K .* A);
    Z = (A_0 - z_0)./(4^K .* A);
    P = (-X -Y -Z)/2;
    
    E_2 = X.*Y + X.*Z + Y.*Z - 3*P.^2;
    E_3 = X.*Y.*Z + 2*E_2.*P + 4*P.^3;
    E_4 = (2*X.*Y.*Z.* + E_2.*P + 3*P.^3) .* P;
    E_5 = X.*Y.*Z.*P.^2;
    
    RJ = 1./(4^K * A.^(3/2)) .* (1 - 3/14 * E_2 ...
                                   + 1/6 * E_3 ...
                                   + 9/88 * E_2.^2 ...
                                   - 3/22 * E_4 ...
                                   - 9/52 * E_2.*E_3 ...
                                   + 3/26 * E_5) ...
                                   + 6 * V;

end
