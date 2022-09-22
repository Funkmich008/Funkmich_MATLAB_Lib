% Reccurence Relation for the spherical Bessel functions of the second kind
% for positive and negative integers of n
function y = SphericalBessely_REC_SINGLE(n, x)
    
    if n>= 0
        %% Reccurence for positive n
        if n == 0
            y = -cos(x)./x;
        elseif n == 1
            y = -cos(x)./x.^2 - sin(x)./x;
        else
            
            y0 = - cos(x)./x;
            y1 = -cos(x)./x.^2 - sin(x)./x;
    
            y = 0;
            for i=2:n
                y = (2*(i-1)+1)./x .* y1 - y0;
                y0 = y1;
                y1 = y;
            end
        end
    else
        %% Reccurence for negative n
        y0 = -cos(x)./x;
        y1 = -cos(x)./x.^2 - sin(x)./x;

        y = 0;
        for i=0:abs(n)-1
            y = (-2*i+1)./x .* y0 - y1;          
            y1 = y0;
            y0 = y;
        end
    end

    % limit at x = 0
    if n == -1
        y(x==0) = 1;
    elseif n < -1
        y(x==0) = 0;
    end

end