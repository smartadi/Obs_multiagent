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





x0 = normrnd(0,1,n,1);
U = normrnd(0,1,m,T);

H = zeros(n*T,m*T);

X = x0;
W=[];

e = 0.01;
for i = 1:T
    W = [W,e*rand(n,1)];
    X = [X, A*X(:,end) + B*U(:,i) + W(:,end)];
end

eps = max(vecnorm(W,2,1).^2)*1.05;
ph1 = T*eps*eye(n,n);
ph12 = 0*zeros(n,T);
ph2 = -eye(T,T);
figure()
plot(X(1,:));hold on;
plot(X(2,:));
plot(X(3,:));



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
% cvx_begin sdp
%     variable M(3*n+m,3*n+m) symmetric
%     M - N >= -0*eye(3*n+m)
% cvx_end

%% noise cov
C = W*(eye(T) - ones(T,T)/T)*W'*1/T;

ph22 = -(eye(T) - ones(T,T)/T)*1/T;

Phi2 = [ph1  ,ph12;
       ph12',ph22];

nn = norm(Phi-Phi2)



%% Assumption 1 *holds
AA = [eye(n);W']
BB = AA'*Phi2*AA;
eig(BB)
%% Remark 1 **does not hold

R = [ph1,W;W',-inv(ph22)];
e = eig(R)

r = eig(-inv(ph22) - W'*inv(ph1)*W)

%% equation 8 *holds
clc
AAA = [eye(n);A';B'];
S = eig(AAA'*q*Phi2*q'*AAA)
%% Theorem 13
clc
N22 = eig([Xm;Um]*ph22*[Xm;Um]')
%% Theorem 14
nnn = eig(NN)
%%
N = Q*Phi*Q';
% N = Q*Phi2*Q';

cvx_begin sdp
    variable P(n,n) symmetric
    variable L(m,n)
    variable b nonnegative
    variable a nonnegative 
    P >= 0.000001*eye(n)
    [P - b*eye(n), zeros(n,n), zeros(n,m), zeros(n,n);
     zeros(n,n),         -P,        -L', zeros(n,n);
     zeros(m,n),         -L, zeros(m,m),          L;
     zeros(n,n), zeros(n,n),         L',          P] - a*N >= 0*eye(3*n+m)
cvx_end
% 
% cvx_begin sdp
%     variable sA(2n+m,2n+m) symmetric
%     variable sB(2n+m,n)
%     variable sC(n,n) symmetric
% 
%     sC > 0*eye(n)
%     [P - b*eye(n), zeros(n,n), zeros(n,m), zeros(n,n);
%      zeros(n,n),         -P,        -L', zeros(n,n);
%      zeros(m,n),         -L, zeros(m,m),          L;
%      zeros(n,n), zeros(n,n),         L',          P] - a*N >= 0*eye(3*n+m)
% cvx_end

%%
K = L*inv(P)

eig(A+B*K)
%%
% 
% FS = [P - b*eye(n), zeros(n,n), zeros(n,m), zeros(n,n);
%      zeros(n,n),         -P,        -L', zeros(n,n);
%      zeros(m,n),         -L, zeros(m,m),          L;
%      zeros(n,n), zeros(n,n),         L',          P] - a*N;
% 
% eig(FS)

%%
% N = Q*Phi*Q';
% NN = q*Phi*q';
% eps = 1e-9
% cvx_begin sdp
%     variable M(3*n+m,3*n+m) symmetric
%     variable a nonnegative
%     variable b nonnegative
%     M - a*N >= -eps*eye(3*n+m)
%     M(1:n,n+1:end) == 0
%     M(n+1:2*n,end-n+1:end) == 0
%     M(2*n+1:2*n+m,2*n+1:2*n+m) == 0
%     M(n+1:2*n,n+1:2*n) - M(end-n+1:end,end-n+1:end) ==0
%     M(2*n+1:2*n+m,n+1:2*n) - M(2*n+1:2*n+m,end-n+1:end) == 0
%     M(1:n,1:n) - M(end-n+1:end,end-n+1:end) + b*eye(n) == 0
%     M(end-n+1:end,end-n+1:end) > 0
% cvx_end
% 
% Mm = full(M)
% a
% L = M(2*n+1:2*n+m,n+1:2*n)
% P = M(end-n+1:end,end-n+1:end)
% b
% K = L*inv(P)
% eig(A+B*K)
%%
x0 = normrnd(0,1,n,1);
U = normrnd(0,1,m,T);

H = zeros(n*T,m*T);

X = x0;
W=[];

e = 0.01;
for i = 1:T
    W = [W,e*rand(n,1)];
    X = [X, A*X(:,end) + B*U(:,i) + W(:,end)];
end


figure()
plot(X(1,:));hold on;
plot(X(2,:));
plot(X(3,:));