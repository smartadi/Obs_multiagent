%% data informativity test
clear all;
close all;
clc;
A = [-0.2414 -0.8649 0.6277;
      0.3192 -0.0301 1.0933;
      0.3129 -0.1649 1.1093];

B = [1 0;
     0 2;
     1 1];

C = [0 0 1];

A=[ 0.850, -0.038, -0.380;
    0.735, 0.815, 1.594;
    -0.664, 0.697, -0.064];

B = [ 1.431, 0.705;
    1.620 -1.129;
    0.913 0.369];

n=3;
m=2;
p=1;

T = 20;

I = [];

x0 = normrnd(0,1,n,1);
U = normrnd(0,1,m,T);

X = x0;
W=[];

e = 2.4;
for i = 1:T
    W = [W,sqrt(e)*rand(n,1)/sqrt(n)];
    X = [X, A*X(:,end) + B*U(:,i) + W(:,end)];
    
end
Y = C*X;


eps = max(vecnorm(W,2,1).^2);
ph1 = T*eps*eye(n,n);
ph12 = 0*zeros(n,T);
ph2 = -eye(T,T);
% figure()
% plot(X(1,:));hold on;
% plot(X(2,:));
% plot(X(3,:));

Xm = X(:,1:end-1);
Ym = C*Xm;
Xp = X(:,2:end);
Um = U;
Zm = [Ym;Um];
Hz = [eye(p);
    zeros(1,1);
    zeros(1,1);
    zeros(1,1)];
Hu = eps*eye(p);
Rm = eye(T);                                                                         

Q11 = T*Hu;
Q12 = zeros(size(Hz,2),T);
Q22 = -eye(T);

Qq = [eye((p+m)*l), Hz*Ym*Rm';
    zeros(n+m,n), -Zm*Rm';
    zeros(m,n), -Um*Rm';
    zeros(n,n), zeros(n,T)];
Qe=[Hz*Q11*Hz^T, H_z*Q12;
    Q12^T*Hz^T, Q22];
% q = [eye(n), Xp;
%     zeros(n,n), -Xm;
%     zeros(m,n), -Um];
% Phi = [ph1  ,ph12;
%        ph12',ph2];

N = Qq*Qe*Qq';
echo off


cvx_begin sdp quiet
    variable P(n,n) symmetric
    variable L(m,n)
    variable b nonnegative
    variable a nonnegative 
    P >= 0.000001*eye(n);
    [P - b*eye(n), -J1*P-J2*L, zeros(n,m),  J1*P+J2*L;
     -P*J1'-L*J2',         -P,        -L', zeros(n,n);
       zeros(m,n),         -L, zeros(m,m),          L;
     P*J1'+L'*J2', zeros(n,n),         L',          P] - a*N >= 0*eye(3*n+m);
cvx_end

K = L*inv(P);

t = max(eig(A+B*K));

I = [I,t];

%%
x0 = normrnd(0,1,n,1);
U = normrnd(0,1,m,T);

H = zeros(n*T,m*T);

X = x0;
W=[];


for i = 1:T
    W = [W,sqrt(e)*rand(n,1)/sqrt(n)];
    X = [X, (A+B*K)*X(:,end) + W(:,end)];
end


figure()
plot(X(1,:));hold on;
plot(X(2,:));
plot(X(3,:));
%%
