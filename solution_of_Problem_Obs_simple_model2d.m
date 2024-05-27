clc;
close all;
clear all;

%%
n= 4;
m= 2;
l= 2;

N= 1000;
x= zeros(n,N);
u= zeros(m,N);
y= zeros(l,N);


A = [0,0,1,0;
     0,0,0,1;
     0,0,0,0;
     0,0,0,0];

eig(A)

B = [0 0;
     0 0;
     1 0;
     0 1;];

C = [1,0,0,0;
     0,1,0,0];
 
D = [];

dt = 0.01;

sys1 = ss(A,B,C,D);
sys2 = c2d(sys1,dt);

At = sys2.A;
Bt = sys2.B;
Ct = sys2.C;


Qx = 1*eye(2);
Qv = 0.1*eye(2);
Q = [Qx,zeros(2,2);
    zeros(2,2),Qv];
R = 0.1*eye(2);

[Ks,S,e] = dlqr(At,Bt,0.001*Q,0.001*R,[]);

%stabilized system
As = At-Bt*Ks;


Qx = 1*eye(2);
Qv = 0.1*eye(2);
Q = [Qx,zeros(2,2);
    zeros(2,2),Qv];
R = 0.1*eye(2);

[K,S,e] = dlqr(As,Bt,Q,R,[]); 

x0 = [1,1,0,0]';  
x(:,1) = x0; 
for i=1:N
    u(:,i) = -K*x(:,i);
    x(:,i+1) = As*x(:,i) + Bt*u(:,i);
    y(:,i) = Ct*x(:,i);
end


figure()
plot(x(1,:),x(2,:));
title('trajectory')
grid on;


%% 

xest = [1,1,-1,0]'; 
% window size 
N = 150;
N = 170;
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
        y(:,i) = Ct*x(:,i);
    end

    X_bar = x(:);
    X_bar(1:4)=[];

    U_bar = u(:);


    ly_exp=[];
    eps=1;
    

    cvx_begin
        variable U(m*N)
        variable X(n*(N))
        L=[];
        Kl=[];
        obj=0;
        for i=1:N
            lx = phi_set_L_cvx2(Ct*X((i-1)*(n)+1:i*n));
            kx = phi_set_k_cvx2(Ct*X((i-1)*(n)+1:i*n));
            L=[L,lx];
            Kl=[Kl,kx];

            ly = eps*(min(real(eig((As-Bt*K)^i))));
            obj = obj + ly'*ly*(1 - lx) - 2*kx;
        end       

        x0=zeros(n,1);

        maximize(obj)
        subject to
            %X  == X_bar + A*x0 + H*(U);
            X  == X_bar + A*x0 + H*(U);

            % X(1:n) == X_bar(1:n);
            (X-X_bar)'*(X-X_bar) <= 1e2*2e1
            X(end-n:end) == X_bar(end-n:end);
            U'*U <= 1e2*2e1

    cvx_end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % u_obs = U(1:m)
    % xest = As*x0 + Bt*u_obs;
    % xobs=[xobs,xest];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%
xx = reshape(X,n,[]);
uu = reshape(U,m,[]);

%%
close all;
figure()
plot(xx(1,:),xx(2,:),'k','LineWidth',2);hold on;
plot(xx(1,1),xx(2,1),'ko');
plot(xx(1,end),xx(2,end),'bo');
plot(x(1,:),x(2,:),'b','LineWidth',2);hold on;
plot(x(1,1),x(2,1),'ko');
plot(x(1,end),x(2,end),'bo');
plot(1,0,'ro')
plot(0,0,'go')
xlabel('x');
ylabel('y');

dx = xx-x(:,1:end-1);

figure()
plot(dx(1:2,:))

figure()
plot(x(1,:),x(2,:));hold on;
plot(x(1,1),x(2,1),'ko');
plot(x(1,end),x(2,end),'bo');
plot(1,0,'ro')
plot(0,0,'go')
xlabel('x');
ylabel('y');

uuu = uu+u;
figure()
plot(u(1,:),'-b','LineWidth',2);hold on;
plot(u(2,:),'-r','LineWidth',2)
plot(uuu(1,:),'b')
plot(uuu(2,:),'r')
plot(uu(1,:),'bo')
plot(uu(2,:),'ro')