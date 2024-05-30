clc;
close all;
clear all;
%% Non MPC Solution
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
N = 200;



%
s = 6;
Qx = 1*eye(2);
Qv = 0.1*eye(2);
Q = [Qx,zeros(2,2);
    zeros(2,2),Qv];

QS = kron(eye(s),Q);
RS = kron(eye(s),R);


AS = kron(eye(s),As);
BS = kron(eye(s),Bt);
CS = kron(eye(s),Ct);
[KS,S,e] = dlqr(AS,BS,QS,RS,[]);

% AA = kron(eye(s),A);
% HS = kron(eye(s),H);
%%
AA=[];
HS = [];
for i = 1:N
    AA = [AA;AS^(i)];
    G=[];
    for j= 1:N
        if i >= j
        G = [G,AS^(i-j)*BS];
        else
        G = [G,zeros(n*s,m*s)];
        end
    end

    HS = [HS;G];

end
%%
xest=[]

for i=1:s
    xest = [xest;0.1*(0.5-rand(n,1))];
end

xest =  kron(ones(s,1),[1;1;0;0])+xest;
xobs = xest;

x0 = xest;
u = zeros(s*m,N);
x = zeros(s*n,N);
y = zeros(s*l,N);
x(:,1) = x0;

% Generate optimal trajectory
for i=1:N
    u(:,i) = -KS*x(:,i);
    x(:,i+1) = AS*x(:,i) + BS*u(:,i);
    y(:,i) = CS*x(:,i);
end

X_bar = x(:);
X_bar(1:n*s)=[];

U_bar = u(:);

ly_exp=[];
eps=1;
delta=1e2;   
L=[];
Kl=[];

Lb=[];
Klb=[];
obj=0; 
cvx_begin 
    variable U(m*s*N)
    variable X(n*s*N)
        
               
        
    % Xre = reshape(X,n*s,N)



        for i=1:N
            ds=[];

            % L=[];
            % Kl=[];
            ls = [];
            ks = [];
            ok=0;
            ol=0;
            ly = eps*(min(abs(eig((As-Bt*K)^i))));

            for j=1:s
                 lx = phi_set_L_cvx2(Ct*X((i-1)*s*n + (j-1)*n+1 : (i-1)*s*n +(j)*n));
                 kx = phi_set_k_cvx2(Ct*X((i-1)*s*n + (j-1)*n+1 : (i-1)*s*n +(j)*n));

                 lb = phi_set_L_cvx2(Ct*X_bar((i-1)*s*n + (j-1)*n+1 : (i-1)*s*n +(j)*n));
                 kb = phi_set_k_cvx2(Ct*X_bar((i-1)*s*n + (j-1)*n+1 : (i-1)*s*n +(j)*n));



                % lx = phi_set_L_cvx2(Ct*Xre((j-1)*n+1 : (j)*n,i));
                % kx = phi_set_k_cvx2(Ct*Xre((j-1)*n+1 : (j)*n,i));

            L=[L;lx];
            Kl=[Kl;kx];

            Lb=[Lb;lb];
            Klb=[Klb;kb];

            ls=[ls;lx];
            ks=[ks;kx];

            % ly = eps*(min(abs(eig((As-Bt*K)^i))));

            % ly = eps*(min(abs((eig((As-Bt*K)^i)))));

            
            % ok = ok  ;
            obj = obj + ly'*ly*(1-lx) - 2*kx;
      
            end

        end
        for j = 1:2n
            
        
        end


        dx0=zeros(n*s,1);

        

        maximize(obj)
        subject to
            %X  == X_bar + A*x0 + H*(U);
            X  == X_bar + AA*dx0 + HS*(U);

            X(1:n*s) == X_bar(1:n*s);
            
            (X-X_bar)'*(X-X_bar) <= s^2*1e3

            X(end-s*n:end) == X_bar(end-s*n:end);
            %U'*U <= 1e5

            U'*U <= U_bar'*U_bar + delta*s^2


    cvx_end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % u_obs = U(1:m)
    % xest = As*x0 + Bt*u_obs;
    % xobs=[xobs,xest];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    u_obs = U(1:s*m) + U_bar(1:s*m);
    % xest = X_bar(1:n) + As*dx0 + Bt*u_obs;
    xest = kron(eye(s),As)*x0 + kron(eye(s),Bt)*u_obs;

%%
xx = reshape(X,s*n,[]);
uu = reshape(U,s*m,[]);
xb = reshape(X_bar,s*n,[]);
ub = reshape(U_bar,s*m,[]);

ls = reshape(L,s,[]);
ks = reshape(Kl,s,[]);

%%
close all;
xe = xx-xb
figure()
for i =1:s
plot(xx((i-1)*n+1,:),xx((i-1)*n+2,:));hold on;
plot(xx((i-1)*n+1,1),xx((i-1)*n+2,1),'ko');
plot(xx((i-1)*n+1,end),xx((i-1)*n+2,end),'bo');
plot(xb((i-1)*n+1,:),xb((i-1)*n+2,:));hold on;
plot(xb((i-1)*n+1,1),xb((i-1)*n+2,1),'ko');
plot(xb((i-1)*n+1,end),xb((i-1)*n+2,end),'bo');
end
plot(1,0,'ro')
plot(0,0,'go')
xlabel('x');
ylabel('y');
title("obs")
%

figure()
for i =1:s
plot(xe((i-1)*n+1,:),xe((i-1)*n+2,:));hold on;
plot(xe((i-1)*n+1,1),xe((i-1)*n+2,1),'ko');
plot(xe((i-1)*n+1,end),xe((i-1)*n+2,end),'bo');
end
plot(1,0,'ro')
plot(0,0,'go')
xlabel('x');
ylabel('y');
title("obs")


figure()
for i =1:s
plot(xe((i-1)*n+1,:));hold on;
end
title("obs")

