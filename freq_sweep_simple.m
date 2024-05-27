clear all;
close all;
clc;


dt=1/35;

%% start with 1hz for 5 seconds
init = 1;

% min and max freq 
fmin = 0.1;
fmax = 15 %% signal amplitude

% amplitude must be chosen such that its not high enough to cause
% saturation and not low enough to encounter dead zones
T = 30;
a = 1;
t = 0:dt:T;
w_min = fmin*2*pi;
w_max = fmax*2*pi;
C1 = 4;
C2 = 0.0187;

C1 = 1;
C2 = 0.1;
K = C2*(exp(C1*t/T)-1);

w = w_min + K*(w_max-w_min);
th = cumsum(w);

y = a* sin(th); 

figure()
plot(y);

figure()
plot(w)


figure()
plot(th)