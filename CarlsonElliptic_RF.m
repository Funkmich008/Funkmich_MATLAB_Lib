%% Carlson Elliptic Integral RF
% Paper: Numerical computation of real or complex elliptic integrals
% Author: B. C. Carlson
function RF = CarlsonElliptic_RF(x_0, y_0, z_0)

    A_0 = (x_0+y_0+z_0)/3;
    
    x = x_0;
    y = y_0;
    z = z_0;
    A = A_0;
    
    K = 20;
    for k = 0:K      
        
        lambda = sqrt(x.*y) ...
               + sqrt(y.*z) ...
               + sqrt(x.*z);
           
        x = 1/4 * (x + lambda);
        y = 1/4 * (y + lambda);
        z = 1/4 * (z + lambda);
        A = 1/4 * (A + lambda);
        
    end
     
    X = (A_0 - x_0)./(4^K .* A);
    Y = (A_0 - y_0)./(4^K .* A);
    Z = - X - Y;
    
    E_2 = X .* Y - Z.^2;
    E_3 = X.*Y.*Z;
    
    RF = 1./sqrt(A) .* (1 - 1/10 * E_2 ...
                          + 1/14 * E_3 ...
                          + 1/24 * E_2.^2 ...
                          - 3/44 * E_2 .* E_3);

end
