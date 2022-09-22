%% Carlson Elliptic Integral RD
% Paper: Numerical computation of real or complex elliptic integrals
% Author: B. C. Carlson
function RD = Carlson_RD(x_0, y_0, z_0)
    
    A_0 = (x_0 + y_0 + 3*z_0)/5;

    x = x_0;
    y = y_0;
    z = z_0;
    A = A_0;
    
    K = 30;
    V = 0;
    for k = 0:K      
        
        lambda = sqrt(x.*y) ...
               + sqrt(x.*z) ...
               + sqrt(y.*z);
        V = V + 4^(-k)./(sqrt(z) .* (z + lambda));  
        x = 1/4 * (x + lambda);
        y = 1/4 * (y + lambda);
        z = 1/4 * (z + lambda);
        A = 1/4 * (A + lambda);

    end
     
    X = (A_0 - x_0)./(4^K * A);
    Y = (A_0 - y_0)./(4^K * A);
    Z = - (X + Y)/3;
    
    E_2 = X.*Y - 6 * Z.^3;
    E_3 = (3 * X.*Y - 8*Z.^2) .* Z;
    E_4 = 3*(X.*Y - Z.^2) .* Z.^2;
    E_5 = X.*Y.*Z.^3;
    
    
    RD = 4^(-K)./sqrt(A).^3 .* (1 - 3/14 * E_2 ...
                                  + 1/6 * E_3 ...
                                  + 9/88 * E_2.^2 ...
                                  - 3/22 * E_4 ...
                                  - 9/52 * E_2 .* E_3 ...
                                  + 3/26 * E_5) ...
                                  + 3 * V;

end