clc;
%% Experiment 1 single sequence length = 1000, repeats 150 

clear all;
close all;
T = 1000
L = 100
n = 10
m = 1
r = 452134;
trials = 150

% signal duration corresponds to the frame rate
% so dt = 1/35 is fixed

%% Signal Generator
% random signal
% U = rand(1,T)>0.8

% spaced signal
% th = 0.4
% 
% avg = 5;
% U=[];
% for i = 1:T
%     rng(r+i)
%     u = double(rand(1)>th);
%     U = [U,u];
%     if i <= avg
%         U(:,end) = 0;
%     elseif U(:,end-1) == 0 && mean(U(:,end-4:end)) > 1/avg
%         U(:,end) = 0;
%     end
% 
% end

U = u_signal(T,r)

H = hankel(U);
Hu = double(H(1:n+L,1:T-n-L-1));
rank(Hu)

figure()
plot(U)

%% repeat sequence 250 times

inputs1 = repmat(U,trials,1);

save inputs1.mat inputs1



