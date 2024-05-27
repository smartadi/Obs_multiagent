function [y,s] = phi(x)   
    
    f=1;
    s = norm(x-[0.5;0.5;0.5]);
    while f == 1
    % r = 2*(0.5-rand(size(x)));
    r = 2*(0.5-rand(size(x)));

    r = r/norm(r);
    if norm(r) <= s
        f = 0;
        y = x + (norm(x-[0.5;0.5;0.5]) + 0.1)*r;
    end
    end


end