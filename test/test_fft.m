%Script to test my fft function

close all,clear all, clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fs = 1;
% L = 20;
% n = 0:L-1;
% y = cos(0.4*pi.*n + 0.1*pi) + 0.5*cos(0.5*pi.*n+0.3*pi) + ...
%     0.2*cos(0.6*pi.*n);
% %y = [1,zeros(1,L-1)];
% 
% Y_m = fftshift(mixed_radix_fft(y', L));
% Y = fftshift(fft(y',L));
% 
% freqs = linspace(-fs/2,fs/2,L)*2*pi;
% figure;
% plot(freqs/pi,abs(Y));hold on;grid on;
% plot(freqs(L/2+1:end)/pi,abs(Y_m));hold off;grid on;
% legend('Matlab fft','My fft');
% xlabel('Frequency in radians normalized by pi');
% ylabel('Magnitude');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = 50:100:5000;
nL = length(L);
t = zeros(4,nL);
err = zeros(1,nL);

for k = 1:nL
    n = 0:L(k)-1;
    y = cos(0.05.*n) + cos(0.04.*n + 0.25*pi);
    tic;
    Y_cor = fft(y',L(k));
    t(1,k) = toc;
    tic;
    Y_est =  mixed_radix_fft(y', L(k));
    t(2,k) = toc;
    tic;
    Y_split = resamp_splitradix_fft(y',2^nextpow2(L(k)),L(k));
    t(3,k) = toc;
    tic;
    dft_mat = exp(-1j*2*pi/L(k)*(0:L(k)-1)'*(0:L(k)-1));
    Y_dft = dft_mat(L(k))*y';
    t(4,k) = toc;
    err(k) = norm(abs(Y_est) - abs(Y_cor));
    
end

figure;
plot(L, t(1,:));hold on;grid on;
plot(L, t(2,:));hold on;grid on;
plot(L, t(3,:));hold on;grid on;
plot(L, t(4,:));hold off;grid on;
xlabel('Length of fft');ylabel('Time in seconds');
legend('Matlab fft','Mixed radix FFT', 'Resampled split radix FFT', 'DFT');
title('FFT speed comparisons');
    
figure;
plot(L, log10(err));grid on;
xlabel('Length of fft');ylabel('Error (log_{10})');
title('Error between my mixed radix fft and Matlab''s fft');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%let's look at correct and resampled fft plots to see what's going on
%This works with fftshift

% freqsycor = linspace(-0.5,0.5,length(Y_cor));
% fs = 2^nextpow2(L(k))/L(k);
% freqsysplit = linspace(-fs/2,fs/2,length(Y_split));
% figure(3);
% subplot(211);plot(freqsycor*2*pi,fftshift(abs(Y_cor)));
% subplot(212);plot(freqsysplit*2*pi, fftshift(abs(Y_split)));
% xlabel('Frequency in rad');ylabel('Amplitude');

