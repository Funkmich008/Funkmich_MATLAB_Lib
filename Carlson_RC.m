%% Carlson Elliptic Integral RC
% Paper: Numerical computation of real or complex elliptic integrals
% Author: B. C. Carlson
function RC = Carlson_RC(x_0, y_0)
    
    A_0 = (x_0 + 2*y_0)/3;

    x = x_0;
    y = y_0;
    A = A_0;
    
    K = 20;
    for k = 0:K      
        
        lambda = 2*sqrt(x.*y) + y;
           
        x = 1/4 * (x + lambda);
        y = 1/4 * (y + lambda);
        A = 1/4 * (A + lambda);

    end
     
    s = (y_0 - A_0)./(4^K * A);
    
    RC = 1./sqrt(A) .* (1 + 3/10 * s.^2 ...
                          + 1/7 * s.^3 ...
                          + 3/8 * s.^4 ...
                          + 9/22 * s.^5 ...
                          + 159/208 * s.^6 ...
                          + 9/8 * s.^7);

end