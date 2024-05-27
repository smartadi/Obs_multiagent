clc;
close all;
clear all;


s=[];
Y=[];
K=[];
d = [0;0;0];
for i = 1:10000
    [y,k] =phi_set(d);
    Y=[Y,y];
    K=[K,k];
    s = [s,norm(y-d)];
    i
end

figure()
plot(s);hold on;
plot(K)
max(s)

figure()
plot3(Y(1,:),Y(2,:),Y(3,:),'ro')

%%
n= 6;
m= 3;
l= 3;

N= 1000;
x= zeros(n,N);
u= zeros(m,N);
y= zeros(l,N);


A = [0,0,0,1,0,0;
     0,0,0,0,1,0;
     0,0,0,0,0,1;
     0,0,0,0,0,0;
     0,0,0,0,0,0;
     0,0,0,0,0,0];

% A = [0,0,0,1,0,0;
%      0,0,0,0,1,0;
%      0,0,0,0,0,1;
%      0,0,0,0,0,0;
%      0,0,0,0,0,0;
%      0,0,0,0,0,0] + 0.01*diag([ones(6,1)])

eig(A)

B = [0 0 0;
     0 0 0;
     0 0 0;
     1 0 0;
     0 1 0;
     0 0 1];

C = [1,0,0,0,0,0;
     0,1,0,0,0,0;
     0,0,1,0,0,0];
 
D = [];

dt = 0.01;

sys1 = ss(A,B,C,D);
sys2 = c2d(sys1,dt);

At = sys2.A;
Bt = sys2.B;
Ct = sys2.C;


Qx = 1*eye(3);
Qv = 0.1*eye(3);
Q = [Qx,zeros(3,3);
    zeros(3,3),Qv];
R = 0.1*eye(3);

[Ks,S,e] = dlqr(At,Bt,0.001*Q,0.001*R,[]);

%stabilized system
As = At-Bt*Ks;


Qx = 1*eye(3);
Qv = 0.1*eye(3);
Q = [Qx,zeros(3,3);
    zeros(3,3),Qv];
R = 0.1*eye(3);

[K,S,e] = dlqr(As,Bt,Q,R,[]); 

x0 = [1,1,1,0,0,0]';  
x(:,1) = x0; 
for i=1:N
    u(:,i) = -K*x(:,i);
    x(:,i+1) = As*x(:,i) + Bt*u(:,i);
    y(:,i) = Ct*x(:,i);
end


figure()
plot3(x(1,:),x(2,:),x(3,:));
title('trajectory')
grid on;


%% 

xest = [1,1,1,0,0,0]'; 
% window size 
N = 130;
A = [];
H = [];
for i = 1:N
    A = [A;As^(i)];
    G=[];
    for j= 1:N
        if i >= j
        G = [G,As^(i-j)*Bt];
        else
        G = [G,zeros(n,m)];
        end
    end

    H = [H;G];

end

xobs = xest;

%%
for iter = 1:1
    x0 = xest;
    u = zeros(m,N);
    x = zeros(n,N);
    y = zeros(l,N);
    x(:,1) = x0;
    Lxt=[];
    % Generate optimal trajectory
    for i=1:N
        u(:,i) = -K*x(:,i);
        x(:,i+1) = As*x(:,i) + Bt*u(:,i);
        [yp(:,i),k,lxt] = phi_set_L(Ct*x(:,i));
        y(:,i) = Ct*x(:,i);
    end

    X_bar = x(:);
    X_bar(1:6)=[];

    U_bar = u(:);


    ly_exp=[];
    eps=1;
    % for i=1:N
    %     ly_exp = [ly_exp;eps*(min(real(eig((As-B*K)^i))))];
    % end
        
    % g = [2.^(N:-1:1)];

    cvx_begin
        variable U(m*N)
        variable X(n*(N))
 %       variable Lxt(N) nonnegative
 %       variable Qxt(N) nonnegative
        %maximize((1 - sum(*Lxt)) - sum(abs(Lxt'*X)))
        L=[];
        Kl=[];
        obj=0;
        for i=1:N
            lx = phi_set_L_cvx(Ct*X((i-1)*(n)+1:i*n));
            kx = phi_set_k_cvx(Ct*X((i-1)*(n)+1:i*n));
            L=[L,lx];
            Kl=[Kl,kx];

%            Lxt(i) = l;
%            Qxt(i) = k;
            ly = eps*(min(real(eig((As-Bt*K)^i))));
            obj = obj + ly'*ly*(1 - lx) - 2*kx;
        end       

        x0=zeros(n,1);

        maximize(obj)
        subject to
            X  == X_bar + A*x0 + H*(U);
            % X(1:n) == X_bar(1:n);
            (X-X_bar)'*(X-X_bar) <= 1e2
            X(end-n:end) == X_bar(end-n:end);
            U'*U <= 1e2

    cvx_end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u_obs = U(1:m)
    xest = As*x0 + Bt*u_obs;
    xobs=[xobs,xest];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%
xx = reshape(X,6,[]);

%
close all;
figure()
plot3(xx(1,:),xx(2,:),xx(3,:));hold on;
plot3(xx(1,1),xx(2,1),xx(3,1),'ko');
plot3(xx(1,end),xx(2,end),xx(3,end),'bo');
plot3(1,0,0,'ro')
plot3(0,0,0,'go')
xlabel('x');
ylabel('y');
zlabel('z');

figure()
plot3(x(1,:),x(2,:),x(3,:));hold on;
plot3(x(1,1),x(2,1),x(3,1),'ko');
plot3(x(1,end),x(2,end),x(3,end),'bo');
plot3(1,0,0,'ro')
plot3(0,0,0,'go')
xlabel('x');
ylabel('y');
zlabel('z');