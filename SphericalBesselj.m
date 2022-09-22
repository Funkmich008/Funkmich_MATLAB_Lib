% Reccurence Relation for the spherical Bessel functions of the first kind
% for integer orders of n
function j = SphericalBesselj(n, x)   
    if n>= 0
        %% Reccurence for positive n
        if n == 0
            j = sin(x)./x;
        elseif n == 1
            j = sin(x)./x.^2 - cos(x)./x;
        else
            
            j0 = sin(x)./x;
            j1 = (sin(x)./x.^2 - cos(x)./x);
    
            j = 0;
            for i=2:n
                j = (2*(i-1)+1)./x .* j1 - j0;
                j0 = j1;
                j1 = j;
            end
        end
    else
        %% Reccurence for negative n
        j0 = sin(x)./x;
        j1 = (sin(x)./x.^2 - cos(x)./x);

        j = 0;
        for i=0:abs(n)-1
            j = (-2*i+1)./x .* j0 - j1;          
            j1 = j0;
            j0 = j;
        end
    end

    % limit at x = 0
    if n == 0
        j(x==0) = 1;
    else
        j(x==0) = 0;
    end
end
