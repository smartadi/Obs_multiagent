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

As = At-Bt*Ks;


Qx = 1*eye(3);
Qv = 0.1*eye(3);
Q = [Qx,zeros(3,3);
    zeros(3,3),Qv];
R = 0.1*eye(3);

[K,S,e] = dlqr(As,Bt,Q,R,[]); 

x0 = [1,1,1,1,0,0]'; 
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

%% Add non linear sectors
close all;
N=1000;
x= zeros(n,N);
u= zeros(m,N);
yp= zeros(l,N);
y= zeros(l,N);

cost=[];
costy=[];
c=0;
cy=0;
costyy=[];
cyy=0;
x(:,1) = x0;

for i=1:N
    
    u(:,i) = +2*K*x(:,i);
    x(:,i+1) = As*x(:,i) + Bt*u(:,i);
    yp(:,i) = phi2(Ct*x(:,i));
    y(:,i) = Ct*x(:,i);
    c = c+x(:,i)'*x(:,i) + u(:,i)'*u(:,i);
    cost = [cost; c];

    cy =  cy + yp(:,i)'*yp(:,i);
    costy = [costy; cy];
    
    cyy =  cyy + y(:,i)'*y(:,i);
    costyy = [costyy; cyy];
end


figure()
plot3(x(1,:),x(2,:),x(3,:))
title('trajectory')


figure()
plot3(yp(1,:),yp(2,:),yp(3,:));hold on;
plot3(y(1,:),y(2,:),y(3,:));
legend('error','normal')
title('output')


figure()
% plot(cost);hold on;
plot(costy);hold on;
plot(costyy)
legend('uncertain','exact')
title('cost')


output_error = vecnorm(y - yp,2,1);

out = abs(costy - costyy);
%
figure()
plot(output_error)
title('output error')

figure()
plot(out)
title('cost difference')

%% Output based controller 
x= zeros(n,N);
uy= zeros(m,N);
yp= zeros(l,N);

Ky = 0.5*eye(3);
Kyd = 1*eye(3); 
x(:,1) = x0; 
ypre = yp(:,1);
for i=1:N
    
    yp(:,i) = Ct*x(:,i);
    err = (yp(:,i) - ypre)/dt;
    uy(:,i) = -Ky*yp(:,i) - Kyd*err;
    x(:,i+1) = At*x(:,i) + Bt*uy(:,i);
    ypre = yp(:,i);
end

figure()
plot3(x(1,:),x(2,:),x(3,:))


figure()
plot(x(1,:));hold on;
plot(x(2,:));
plot(x(3,:));

%% Output based control sector based
close all;
x= zeros(n,N);
uy= zeros(m,N);
yp= zeros(l,N);

Ky = 0.5*eye(3);
Kyd = 1*eye(3);
x(:,1) = x0; 
ypre = yp(:,1);
for i=1:N
    
    yp(:,i) = phi2(Ct*x(:,i));
    err = (yp(:,i) - ypre)/dt;
    uy(:,i) = -Ky*yp(:,i) - Kyd*err;
    x(:,i+1) = At*x(:,i) + Bt*uy(:,i);
    ypre = yp(:,i);
end

figure()
plot3(x(1,:),x(2,:),x(3,:))

figure()
plot(x(1,:));hold on;
plot(x(2,:));
plot(x(3,:));


%% 

xest = [1,1,1,1,0,0]'; 
% window size 
N = 100;
A = [];
H = [];
for i = 1:N
    A = [A;As^(i)];
    G=[];
    for j= 1:N
        if i <= j
        G = [G,As^(i-j)*B];
        else
        G = [G,zeros(n,m)];
        end
    end

    H = [H;G];

end

A = [eye(n);A];
H = [zeros(n,m*N);H];
% while error > 1e-3
for iter = 1:1
    x0 = xest;
    u = zeros(m,N);
    x = zeros(n,N+1);
    y = zeros(l,N);
    x(:,1) = x0;
    Lxt=[];
    % Generate optimal trajectory
    for i=1:N
        u(:,i) = K*x(:,i);
        x(:,i+1) = As*x(:,i) + Bt*u(:,i);
        [yp(:,i),k,lxt] = phi_set_L(Ct*x(:,i));
        y(:,i) = Ct*x(:,i);
    end

    Xt = x(:);
    %U = u(:);
    ly_exp=[];
    for i=1:N
        ly_exp = [ly_exp;eps*(min(eig((As-B*K)^i)))];
    end
        
    g = [2.^(N:-1:1)];
    cvx_begin
        variables U(m*N) Lxt(N) XT(n*N) Xt(n*N) X(N)
        %maximize((1 - sum(*Lxt)) - sum(abs(Lxt'*X)))
        maximize obj

        obj = max(0,g)

        subject to
            X == A*x0 + H*[zeros(m,1),U];
            Xt == X(1:end-n);
            XT == X(n+1:end);
            for i=1:N
                Lxt(i) == phi_set_L*(Xt(i));
                X(i) = norm(XT((i-1)*n+1:i*n)-Xt((i-1)*n+1:i*n));
            end
    cvx_end

end

%%

s=[];
Y=[];
K=[];
D=[];
for i = 1:10000
    d = 0.5-rand(3,1);
    D=[D,d];
    [y,k] =phi_set(d);
    Y=[Y,y];
    K=[K,k];
    s = [s,norm(y-d)];
    i
end


%%
[a1 b1] = min(K);
[a2 b2] = max(K);


c1=D(:,b1)
c2=D(:,b2)


%%
(a1-a2)/norm(c1-c2)

%%
D(:,1)
D(:,4)
K(1)
K(4)
abs(K(1) - K(4))/norm(D(:,1)-D(:,4))


%% Lipschitz constant is invariant
d1=[0.5,0.5,0.5];
[y1,k1] =phi_set(d1)

d2=[0.6,0.6,0.6];
[y2,k2] =phi_set(d2)


diff = (k2-k1)/norm(d2-d1)


d3=[0.7,0.7,0.7];
[y3,k3] =phi_set(d3)


diff = (k3-k1)/norm(d3-d1)




d4=[1,2,3];
[y4,k4] =phi_set(d4)


diff = (k4-k1)/norm(d4-d1)
k4-k1
norm(d4-d1)

%%
DD = vecnorm(D -[0.5;0.5;0.5],2,1);

[DD I] = sort(DD)

KK = K(I);

figure()
plot(DD,KK)

%%
DD=[];
Lkm=[];
for j=1:1000
d = 0.5-rand(3,1);
d = 2*e/norm(e);
DD = [DD,d];
% d = [1;1;1];
Lk=[];
for i= 1:100
e = 0.5-rand(3,1);
e = 0.1*e/norm(e);
dd = d + e;
[y,k] =phi_set(d);
[yy,kk] =phi_set(dd);

Lk = [Lk, norm(kk-k)];
end
Lkm = [Lkm,max(Lk)];
end
DD = vecnorm(DD -[0.5;0.5;0.5],2,1);

[DD I] = sort(DD)

Lkm = Lkm(I);

figure()
plot(DD,Lkm)


figure()
plot(Lk)
%%
clc;
close all;
DD=[];
Lk=[];
for j=1:1000
d = 0.5-rand(3,1);
d = 2*d/norm(d);
DD = [DD,d];
% d = [1;1;1];
[y,k,L] =phi_set_L(d);
Lk = [Lk,L];
end
DD = vecnorm(DD -[0.5;0.5;0.5],2,1);

[DD I] = sort(DD);

Lk = Lk(I);

figure()
plot(DD,Lk)


