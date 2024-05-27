function [y,k,L] = phi_set_L(x)   
    
    k = 0.1*((x-[0.5;0.5;0.5])'*(x-[0.5;0.5;0.5]) + 0.1);
    % r = 2*(0.5-rand(size(x)));
    r = (0.5-rand(size(x)));
    if r'*r > k
        r = k*r./(r'*r);
    end
    y = x + r;
    L = 0.1*norm(x-[0.5;0.5;0.5]);

end