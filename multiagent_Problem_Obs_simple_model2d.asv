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
%N = 150;
N = 15;
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
I = 200;
xiter = [];
s = 4;
xest=[]

for i=1:s
    xest = [xest;0.5*(0.5-rand(4,1))];
end

xest =  kron(ones(4,1),[1;1;0;0])+xest;


Q = kron(eye(s),Q);
R = kron(eye(s),R);

[Kk,S,e] = dlqr(kron(eye(s),As),kron(eye(s),Bt),Q,R,[]); 


Q = kron(eye(N),Q);
R = kron(eye(N),R);


% for iter = 1:I
for iter = 1:2
    x0 = xest;
    u = zeros(s*m,N);
    x = zeros(s*n,N);
    y = zeros(s*l,N);
    x(:,1) = x0;
    
    
    % Generate optimal trajectory
    for i=1:N
        u(:,i) = -Kk*x(:,i);
        x(:,i+1) = kron(eye(s),As)*x(:,i) + kron(eye(s),Bt)*u(:,i);
        y(:,i) = kron(eye(s),Ct)*x(:,i);
    end


    % Xref = kron(ones(N),[1])

    % cvx_begin quiet
    %     variable Uu(m*N)
    %     variable Xx(n*(N))
    % 
    %     obj= (Xx)'*Q*(Xx) + Uu'*R*Uu;       
    % 
    %     dx0=zeros(n,1);
    % 
    %     minimize(obj)
    %     subject to
    %         %X  == X_bar + A*x0 + H*(U);
    %         Xx  ==  A*x0 + H*(Uu);
    % 
    % 
    % cvx_end
    % 
    % 
    % X_bar = Xx;
    %X_bar(1:4)=[];

    X_bar = x(:);
    X_bar(1:12)=[];

    U_bar = u(:);


    ly_exp=[];
    eps=1;
    delta=1e1
    

    cvx_begin quiet
        variable U(m*s*N)
        variable X(n*s*(N))
        L=[];
        Kl=[];
        obj=0;            
        
        Xre = reshape(X,n*s,N)
        for i=1:N
            ds=[];
            for j=1:s
                lx = phi_set_L_cvx2(Ct*Xre((j-1)*n+1 : (j)*n,i));
                kx = phi_set_k_cvx2(Ct*Xre((j-1)*n+1 : (j)*n,i));

                ds = [ds;kron(eye(s),Ct)*Xre((j-1)*n+1 : (j)*n,i)];
            
            
            %lx = phi_set_L_cvx2(Ct*X((i-1)*s*n +1+ (j-1)*s : i*s*n + (j-1)*s));
            %kx = phi_set_k_cvx2(Ct*X((i-1)*s*n +1+ (j-1)*s : i*s*n + (j-1)*s));
            L=[L,lx];
            Kl=[Kl,kx];



            
            end

            Ds = norm([ 1 -1  0;
                   0  1 -1;
                  -1  0  1]*ds)^2;

            obj = obj+Ds;

            ly = eps*(min(real(eig((As-Bt*K)^i))));

            obj = obj + ly'*ly*(1 - lx) - 2*kx - 2*kx;

        end
        

        dx0=zeros(n,1);

        maximize(obj)
        subject to
            %X  == X_bar + A*x0 + H*(U);
            X  == X_bar + A*dx0 + H*(U);

            % X(1:n) == X_bar(1:n);
            (X-X_bar)'*(X-X_bar) <= 1e2
            X(end-n:end) == X_bar(end-n:end);
            U'*U <= 1e2

            U'*U <= U_bar'*U_bar + delta


    cvx_end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u_obs = U(1:m) + U_bar(1:m);
    % xest = X_bar(1:n) + As*dx0 + Bt*u_obs;
    xest = As*x0 + Bt*u_obs;
    xobs=[xobs,xest];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iter
    xiter = [xiter;reshape(X,n,[])];


end

%
xx = reshape(X,n,[]);
uu = reshape(U,m,[]);

%%
close all;
figure()
plot(xx(1,:),xx(2,:),'k','LineWidth',2);hold on;
plot(x(1,:),x(2,:),'b','LineWidth',2);hold on;
plot(xx(1,1),xx(2,1),'ko');
plot(xx(1,end),xx(2,end),'bo');
plot(x(1,:),x(2,:),'b','LineWidth',2);hold on;
plot(x(1,1),x(2,1),'ko');
plot(x(1,end),x(2,end),'bo');
plot(1,0,'r*')
plot(0,0,'go')
xlabel('x');
ylabel('y');
legend('optimal','improved')
title("observability based trajectory")


%%
%
figure()
plot(xobs(1,:),xobs(2,:));hold on;
plot(xobs(1,1),xobs(2,1),'ko');
plot(xobs(1,end),xobs(2,end),'bo');
plot(1,0,'ro')
plot(0,0,'go')
xlabel('x');
ylabel('y');
title("MPC")

%
figure()
plot(xobs(1,:),xobs(2,:),'r','LineWidth',3);hold on;
for i=1:I
    plot(xiter((i-1)*n+1,:),xiter((i-1)*n+2,:),'--k')
    plot(xiter((i-1)*n+1,1),xiter((i-1)*n+2,1),'bo')
    plot(xiter((i-1)*n+1,end),xiter((i-1)*n+2,end),'go')
end
title('MPC')