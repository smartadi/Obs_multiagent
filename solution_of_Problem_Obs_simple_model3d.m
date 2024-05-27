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

xest = [1,1,1,0,0,0]'; 
% window size 
N = 125;
A = [];
H = [];
for i = 1:N
    A = [A;As^(i)];
    G=[];
    for j= 1:N
        if i >= j
        G = [G,As^(i-j)*B];
        else
        G = [G,zeros(n,m)];
        end
    end

    H = [H;G];

end

A = [eye(n);A];
H = [zeros(n,m*N);H];
xobs = xest;
% while error > 1e-3
%%
for iter = 1:1
    x0 = xest;
    u = zeros(m,N);
    x = zeros(n,N+1);
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
    %U = u(:);
    ly_exp=[];
    eps=1;
    % for i=1:N
    %     ly_exp = [ly_exp;eps*(min(real(eig((As-B*K)^i))))];
    % end
        
    % g = [2.^(N:-1:1)];

    cvx_begin
        variable U(m*N)
        variable X(n*(N+1))
 %       variable Lxt(N) nonnegative
 %       variable Qxt(N) nonnegative
        %maximize((1 - sum(*Lxt)) - sum(abs(Lxt'*X)))
        L=[];
        Kl=[];
        obj=0;
        for i=1:N
            l = phi_set_L_cvx(C*X((i-1)*(n)+1:i*n));
            k = phi_set_k_cvx(C*X((i-1)*(n)+1:i*n));
            L=[L,l];
            Kl=[Kl,k];

%            Lxt(i) = l;
%            Qxt(i) = k;
            ly = eps*(min(real(eig((As-Bt*K)^i))));
            obj = obj + ly'*ly*(1 - l) - 2*k;
        end       

%        obj =
% (ly_exp)'*diag(ones(N,1) - Lxt)*(ly_exp) - sum(Qxt);
        maximize(obj)
        subject to
            X  == X_bar + A*x0 + H*U;
            % X(1:n) == X_bar(1:n);
            X(end-n:end) == X_bar(end-n:end);
            U'*U <= 1e1

    cvx_end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u_obs = U(1:m)
    xest = As*x0 + Bt*u_obs;
    xobs=[xobs,xest];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%%
xx = reshape(X,6,[]);

%%
close all;
figure()
plot3(xx(1,:),xx(2,:),xx(3,:));hold on;
plot3(xx(1,1),xx(2,1),xx(3,1),'ko');
plot3(xx(1,end),xx(2,end),xx(3,end),'bo');
plot3(0.5,0,0.5,'ro')
plot3(0,0,0,'go')


figure()
plot3(x(1,:),x(2,:),x(3,:));hold on;
plot3(x(1,1),x(2,1),x(3,1),'ko');
plot3(x(1,end),x(2,end),x(3,end),'bo');
plot3(0.5,0,0.5,'ro')
plot3(0,0,0,'go')