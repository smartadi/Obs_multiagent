clear all;
close all;
clc;
T=30;
f=35;
% chirp = dsp.Chirp(...
%     'Type','Logarithmic',...
%     'SweepDirection', 'Bidirectional', ...
%     'TargetFrequency', 25, ...
%     'InitialFrequency', 0.01,...
%     'TargetTime', T/2, ...
%     'SweepTime', T/2, ...
%     'SamplesPerFrame', T*f, ...
%     'SampleRate', f, ...
%     'InitialPhase',pi/2);

chirp = dsp.Chirp(...
    'Type','logarithmic',...
    'SweepDirection', 'Unidirectional', ...
    'TargetFrequency', 15, ...
    'InitialFrequency', 0.01,...
    'TargetTime', T, ...
    'SweepTime', T, ...
    'SamplesPerFrame', T*f, ...
    'SampleRate', f, ...
    'InitialPhase',pi/2);

y= chirp();
t = 0:1/f:T;

figure()
plot(t(1:end-1),y);
xlabel('time')
ylabel('amp')

% init=35;
% y = [zeros(init,1);y];
% 
% shift = 1.5;
% 
% y = (shift+y)/2;
% 
% figure()
% plot(y);
% xlabel('time')
%% white noise
% 
% rng(5); 
% e = normrnd(0,0.01,T*f+init,1);
% ye = y+e;
% 
% figure()
% plot(ye);

%% experiment repeat

% input3 = repmat(e',200,1);

% save input3.mat input3