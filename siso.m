% MATLAB code for SISO-ANC system
clc;
clear all;
close all;
% loading the Input data
load('foetal_ecg.dat');
x=foetal_ecg;
% time signal;
% given sampling frequency
Fs=500;
timesig=x(:,1);
% abdnomial signals
abdomin1=x(:,2);
abdomin2=x(:,3);
abdomin3=x(:,4);
abdomin4=x(:,5);
abdomin5=x(:,6);
%thoriad signals
thoirad1=x(:,7);
thoirad2=x(:,8);
thoirad3=x(:,9);
figure
subplot(3,1,1);
plot(timesig,abdomin1);
title('abdomin1');
xlabel('time[sec]');
ylabel('Amplitude [mV]');
subplot(3,1 ,2);
plot(timesig,abdomin2);
title('abdomin2');
xlabel('time[sec]');
ylabel('Amplitude [mV]');
subplot(3,1,3);
plot(timesig,abdomin3);
title('abdomin3');
xlabel('time[sec]');
ylabel('Amplitude [mV]');
figure
subplot(2,1,1);
plot(timesig,abdomin4);
title('abdomin4');
xlabel('time[sec]');
ylabel('Amplitude [mV]');
subplot(2,1,2);
plot(timesig,abdomin5);
title('abdomin5');
xlabel('time[sec]');
ylabel('Amplitude [mV]');
figure
subplot(3,1,1);
plot(timesig,thoirad1,'r');
title('thoirad1');
xlabel('time[sec]');
ylabel('Amplitude [mV]');
subplot(3,1,2);
plot(timesig,thoirad2,'r');
title('thoirad2');
xlabel('time[sec]');
ylabel('Amplitude [mV]');
subplot(3,1,3);
plot(timesig,thoirad3,'r');
title('thoirad3');
xlabel('time[sec]');
ylabel('Amplitude [mV]');
%% LMS application for fetus Extraction
% Taking the average of Fetus+ Mothers signal
d=(abdomin1+abdomin2+abdomin3+abdomin4+abdomin5)/5;
% Taking the Average of mothers signal mostly
a=(thoirad3+thoirad1+thoirad2)/3;
% Intialising the Step size
mue= 0.00000002;
% nth order of the filter
nord=12;
X=convm(a,nord);
%Applying LMS algorithm using lms basic function.
[A1,E1,y1] = lms(X,d,mue,nord);
%% Applying NLMS algorithm using nlms basic function
% Defining the beta value
beta=0.009;
nord=12;
X=convm(a,nord);
%Applying nLMS algorithm using lms basic function.
[A2,E2,y2] = nlms(X,d,beta,nord);
%% Applying for LLMS using llms basic function
mu=0.0000002;
gammax=0.001;
nord=12;
X=convm(a,nord);
%Applying LMS algorithm using lms basic function.
[A4,E4,y4] = llms(X,d,mu,gammax,nord);
%% Fetus + mothers;
figure
subplot(2,1,1)
plot(timesig,d(1:2500),'r-');
hold on;
plot(timesig,d(1:2500),'k--');
hold on;
plot(timesig,d(1:2500),'g-.');
hold on;
title('Fetus+mothers ECG for SISO');
xlabel('time[sec]');
ylabel('Amplitude [mV]');
legend('SISO-LMS','SISO-NLMS','SISO-LLMS');
%% Mothers signal
subplot(2,1,2)
plot(timesig,a(1:2500),'r-');
hold on;
plot(timesig,a(1:2500),'k--');
hold on;
plot(timesig,a(1:2500),'g-.');
hold on;
title('Mothers ECG for SISO');
xlabel('time[sec]');
ylabel('Amplitude [mV]');
legend('SISO-LMS','SISO-NLMS','SISO-LLMS');
%% Filter output
figure
subplot(2,1,1)
plot(timesig,y1(1:2500),'r-');
hold on;
plot(timesig,y2(1:2500),'k--');
hold on;
plot(timesig,y4(1:2500),'g-.');
hold on;
title('filter output');
xlabel('time[sec]');
ylabel('Amplitude [mV]');
%axis([0 1 -15 15])
legend('SISO-LMS','SISO-NLMS','SISO-LLMS')
%% plotting final fetus signalsignal.
subplot(2,1,2);
plot(timesig,E1(1:2500),'r-');
hold on;
plot(timesig,E2(1:2500),'k--');
hold on;
plot(timesig,E4(1:2500),'g-.');
hold on;
title('Fetus ECG for SISO');
xlabel('time[sec]');
ylabel('Amplitude [mV]');
%axis([0 1 -15 15])
legend('SISO-LMS','SISO-NLMS','SISO-LLMS');