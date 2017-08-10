%Script to plot MSE vs SNR with CRB

close all, clear all, clc;
N = 2000;
n = 0:N-1;
%clean signal
%example 1
%y = cos(2*0.4*pi.*n + 0.1*pi) + 0.5*cos(2*0.5*pi.*n+0.3*pi) + 0.2*cos(2*0.6*pi.*n);
%example 2
y = cos(2*0.25*pi.*n) + 0.5*cos(2*0.26*pi.*n + 0.25*pi);
theta = [0.25,1,0,0.26,0.5,0.25*pi];
%example 3
%y = cos(0.04.*n) + 0.5*cos(0.05.*n);
%nornalize signal power to 0dB
y_norm = y./max(abs(y));
nsignals = 2;

snr = -20:5:50;
bounds = zeros(3*nsignals,length(snr));
for k = 1:length(snr)
    %signal+noise
    x = awgn(y_norm, snr(k));
    sigma_z = var(y_norm - x);
    %get CRB
    bounds(:,k) = crb(nsignals,N,theta,sigma_z);
end

figure;
%plot CRB
plot(snr, log10(bounds(1:3:end,:)));grid on;
xlabel('SNR in dB');ylabel('log_{10} Mean squared error');
legend('CRB');
title('Frequency of sinusoid = 0.25Hz');


    
    
