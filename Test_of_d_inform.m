%% data informativity test
clear all;
close all;
clc;


A=1

B=1;


n=1;
m=1;

T = 3;


eps = .5;
ph1 = 1;
ph12 = 0*zeros(n,T);
ph2 = -eye(T,T);



% 
% x0 = normrnd(0,1,n,1);
% U = normrnd(0,1,m,T);
% 
% H = zeros(n*T,m*T);
% 
% X = x0;
% W=[];
% for i = 1:T
%     W = [W,0.1*rand(n,1)];
%     X = [X, A*X(:,end) + B*U(:,i) + W(:,end)];
% end
% 
% 
% figure()
% plot(X(1,:));hold on;
% plot(X(2,:));
% plot(X(3,:));
% 


Xm = [0 0 1];
Xp = [0 1 0];
Um = [-1/2 1/2 -3/2];
W = [1/2 1/2 1/2]
Q = [eye(n), Xp;
    zeros(n,n), -Xm;
    zeros(m,n), -Um;
    zeros(n,n), zeros(n,T)];

q = [eye(n), Xp;
    zeros(n,n), -Xm;
    zeros(m,n), -Um];
Phi = [ph1  ,ph12;
       ph12',ph2];
% 
% N = Q*Phi*Q';
% NN = q*Phi*q';
% cvx_begin sdp
%     variable M(3*n+m,3*n+m) symmetric
%     M - N >= -0*eye(3*n+m)
% cvx_end
% 
% 
% %%
% 
% 
% %%
% % P
% % L
% % K = L*inv(P)
% % 
% % eig(A)
% % eig(A + B*K)
% %% noise cov
% C = W*(eye(T) - ones(T,T)/T)*W'*1/T;
% 
% ph22 = -(eye(T) - ones(T,T)/T)*1/T;
% 
% Phi2 = [ph1  ,ph12;
%        ph12',ph22];
% 
% nn = norm(Phi-Phi2)
% 
% 
% 
% %% Assumption 1 *holds
% AA = [eye(n);W']
% BB = AA'*Phi2*AA;
% eig(BB)
% %% Remark 1 **does not hold
% 
% R = [ph1,W;W',-inv(ph22)];
% e = eig(R)
% 
% r = eig(-inv(ph22) - W'*inv(ph1)*W)
% 
% %%
% % clc
% % eig(A)
% % eig(A + B*K)
% % QQ = eig(P - (A+B*K)*P*(A+B*K)')
% %% equation 8 *holds
% clc
% AAA = [eye(n);A';B'];
% S = eig(AAA'*q*Phi2*q'*AAA)
% %% Theorem 13
% clc
% N22 = eig([Xm;Um]*ph22*[Xm;Um]')
% %% Theorem 14
% nnn = eig(NN)
%%
N = Q*Phi*Q';
% N = Q*Phi2*Q';

cvx_begin sdp
    variable P(n,n) symmetric
    variable L(m,n)
    variable b nonnegative
    variable a nonnegative
    P > 0.001*eye(n)
    [P - b*eye(n), zeros(n,n), zeros(n,m), zeros(n,n);
     zeros(n,n),         -P,        -L', zeros(n,n);
     zeros(m,n),         -L, zeros(m,m),          L;
     zeros(n,n), zeros(n,n),         L',          P] - a*N >= -0.001*eye(3*n+m)
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
K = L*inv(P);

eig(A+B*K)
%%

FS = [P - b*eye(n), zeros(n,n), zeros(n,m), zeros(n,n);
     zeros(n,n),         -P,        -L', zeros(n,n);
     zeros(m,n),         -L, zeros(m,m),          L;
     zeros(n,n), zeros(n,n),         L',          P] - a*N;

eig(FS)


%%
N = Q*Phi*Q';
NN = q*Phi*q';
cvx_begin sdp
    variable M(3*n+m,3*n+m) symmetric
    variable a nonnegative
    variable b nonnegative
    M - a*N >= 0.001*eye(3*n+m)
    M(3,3) == 0
    M(1,2) == 0
    M(1,3) == 0
    M(1,4) == 0
    M(2,4) == 0
    M(2,2) == -M(4,4)
    M(3,2) == -M(3,4)
    M(1,1) - M(4,4) +b <= 0
cvx_end

full(M)
a
L = M(3,4)
P = M(4,4)
b
K = L/P
