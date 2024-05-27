function k = phi_set_k_cvx(x)   
    p = [1;0;0];
    k = 0.1*((x-p)'*(x-p) + 0.1);
    % r = 2*(0.5-rand(size(x)));
    r = (0.5-rand(size(x)));
    % if r'*r > k
    %     r = k*r./(r'*r);
    % end
    y = x + r;
    L = 0.1*(x-p)'*(x-p);

end