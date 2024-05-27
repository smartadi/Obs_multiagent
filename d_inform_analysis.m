%% data informativity test
clear all;
close all;
clc;


A=[ 0.850, -0.038, -0.380;
    0.735, 0.815, 1.594;
    -0.664, 0.697, -0.064]

B = [ 1.431, 0.705;
    1.620 -1.129;
    0.913 0.369];

n=3;
m=2;

T = 20;

I = [];
for iter=1:500

    iter
x0 = normrnd(0,1,n,1);
U = normrnd(0,1,m,T);

X = x0;
W=[];

e = 2.4;
for i = 1:T
    W = [W,sqrt(e)*rand(n,1)/sqrt(n)];
    X = [X, A*X(:,end) + B*U(:,i) + W(:,end)];
end

eps = max(vecnorm(W,2,1).^2);
ph1 = T*eps*eye(n,n);
ph12 = 0*zeros(n,T);
ph2 = -eye(T,T);
% figure()
% plot(X(1,:));hold on;
% plot(X(2,:));
% plot(X(3,:));



Xm = X(:,1:end-1);
Xp = X(:,2:end);
Um = U;
Q = [eye(n), Xp;
    zeros(n,n), -Xm;
    zeros(m,n), -Um;
    zeros(n,n), zeros(n,T)];

q = [eye(n), Xp;
    zeros(n,n), -Xm;
    zeros(m,n), -Um];
Phi = [ph1  ,ph12;
       ph12',ph2];

N = Q*Phi*Q';
NN = q*Phi*q';
echo off
cvx_begin sdp quiet
    variable P(n,n) symmetric
    variable L(m,n)
    variable b nonnegative
    variable a nonnegative 
    P >= 0.000001*eye(n);
    [P - b*eye(n), zeros(n,n), zeros(n,m), zeros(n,n);
     zeros(n,n),         -P,        -L', zeros(n,n);
     zeros(m,n),         -L, zeros(m,m),          L;
     zeros(n,n), zeros(n,n),         L',          P] - a*N >= 0*eye(3*n+m);
cvx_end

K = L*inv(P);

t = max(eig(A+B*K));

I = [I,t];
end
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

sum(I < 1)/5