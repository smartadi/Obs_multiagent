%% Experiment 2 single sequence length = 500, repeats = 250 


clc;
clear all;
close all;
T = 500
L = 100
n = 10
m = 1
trials = 250

% signal duration corresponds to the frame rate
% so dt = 1/35 is fixed

%% Signal Generator
r = 452134


inputs2=[];
for i = 1:trials
    U = u_signal(T,r+i);
    inputs2=[inputs2;U];
end



H = hankel(U);
Hu = double(H(1:n+L,1:T-n-L-1));
rank(Hu)

figure()
plot(U)

%% 
% save inputs2.mat inputs2



