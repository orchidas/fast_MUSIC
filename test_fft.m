%Script to test my fft function

close all, clear all, clc;
% fs = 1;
% L = 20;
% n = 0:L-1;
% y = cos(0.4*pi.*n + 0.1*pi) + 0.5*cos(0.5*pi.*n+0.3*pi) + 0.2*cos(0.6*pi.*n);
% %y = [1,zeros(1,L-1)];
% 
% Y_m = fftshift(mixed_radix_fft(y', L));
% Y = fftshift(fft(y',L));
% 
% freqs = linspace(-fs/2,fs/2,L)*2*pi;
% figure;
% plot(freqs/pi,abs(Y));hold on;grid on;
% plot(freqs/pi,abs(Y_m));hold off;grid on;
% legend('Matlab fft','My fft');
% xlabel('Frequency in radians normalized by pi');
% ylabel('Magnitude');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 50:100:3000;
nL = length(L);
t = zeros(3,nL);
err = zeros(1,nL);

for k = 1:nL
    n = 0:L(k)-1;
    %y = cos(0.4*pi.*n + 0.1*pi) + 0.5*cos(0.5*pi.*n+0.3*pi) + 0.2*cos(0.6*pi.*n);
    y = cos(2*0.25*pi.*n) + cos(2*0.26*pi.*n + 0.25*pi);
    tic;
    Y_cor = fft(y',L(k));
    t(1,k) = toc;
    tic;
    Y_est =  mixed_radix_fft(y', L(k));
    t(2,k) = toc;
    tic;
    dft_mat = exp(-1j*2*pi/L(k)*(0:L(k)-1)'*(0:L(k)-1));
    Y_dft = dft_mat(L(k))*y';
    t(3,k) = toc;
    err(k) = norm(abs(Y_est) - abs(Y_cor));
end

figure;
plot(L, t(1,:));hold on;grid on;
plot(L, t(2,:));hold on;grid on;
plot(L, t(3,:));hold off;grid on;
xlabel('Length of fft');ylabel('Time in seconds');
legend('Matlab fft','My FFT', 'DFT');
title('FFT speed comparisons');
    
figure;
plot(L, log10(err));grid on;
xlabel('Length of fft');ylabel('Error (log_{10})');
title('Error between my fft and Matlab''s fft');
