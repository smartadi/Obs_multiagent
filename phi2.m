function [y,s] = phi2(x)   
    
    k = 0.1*(norm(x-[0.5;0.5;0.5]) + 0.1);
    % r = 2*(0.5-rand(size(x)));
    r = (0.5-rand(size(x)));
    if norm(r)>= k
        r = k*r/norm(r);
    end
    y = x + r;


end