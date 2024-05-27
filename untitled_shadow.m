
A=diag([-3,-1,-2]);

cvx_begin sdp
variable P(3,3) 
A'*P + P*A <= -10*eye(3)
P >= eye(3)
cvx_end

