%Script to test MUSIC algorithm

close all, clear all, clc;
n = 0:1999;
%clean signal
%example 1
%y = cos(0.4*pi.*n + 0.1*pi) + 0.5*cos(0.5*pi.*n+0.3*pi) + 0.2*cos(0.6*pi.*n);
%example 2
y = cos(0.5*pi.*n) + cos(0.52*pi.*n + 0.25*pi);
%nornalize signal power to 0dB
y_norm = y./max(abs(y));
snr = 10;
%signal+noise
x = awgn(y_norm, snr);

% %fast MUSIC
% freqs_fast = fast_music(x,2,500,'fft','fft');
% sort(freqs_fast/pi)
% 
% %MUSIC
% freqs = music(x,2,500,'gram_schmidt','fft');
% sort(freqs/pi)

% %measuring computation speed and accuracy
nbins = 50:100:1500;
nsig = 2;
nmethods = 5;
M = 100;
t = zeros(length(nbins), nmethods);
err = zeros(length(nbins),nmethods);
%sig_freqs = [-0.6, -0.5, -0.4, 0.4, 0.5, 0.6];
sig_freqs = [-0.52,-0.5,0.5,0.52];

%frequencies detected by MUSIC with basic QR
for n = 1:length(nbins)
    tic;
    freqs = music(x, nsig, nbins(n), 'gram_schmidt','fft');
    t(n,1) = toc;
    err(n,1) = norm(sort(freqs/pi) - sig_freqs);
end

%this is needed so all variables used by the function are cleared
clearvars -except x nbins nsig nmethods t err sig_freqs snr M

%frequencies detected by MUSIC with hessenberg QR 
for n = 1:length(nbins)
    tic;
    freqs = music(x, nsig, nbins(n), 'hess','fft');
    t(n,2) = toc;
    err(n,2) = norm(sort(freqs/pi) - sig_freqs);
end

clearvars -except x nbins nsig nmethods t err sig_freqs snr M

for n = 1:length(nbins)
    tic;
    freqs = music(x, nsig, nbins(n), 'tridiag','fft');
    t(n,3) = toc;
    err(n,3) = norm(sort(freqs/pi) - sig_freqs);
end

clearvars -except x nbins nsig nmethods t err sig_freqs snr M 

%frequencies detected by fast_MUSIC with fft
for n = 1:length(nbins)
    tic;
    freqs_fast = fast_music(x, nsig, nbins(n), 'fft', 'fft');
    t(n,4) = toc;
    err(n,4) = norm(sort(freqs_fast/pi) - sig_freqs);
end

clearvars -except x nbins nsig nmethods t err sig_freqs snr M

%frequencies detected by fast_MUSIC with dft
for n = 1:length(nbins)
    tic;
    freqs_fast = fast_music(x, nsig, nbins(n), 'dft', 'fft');
    t(n,5) = toc;
    err(n,5) = norm(sort(freqs_fast/pi) - sig_freqs);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
for k = 1:nmethods
    plot(nbins, t(:,k));hold on;grid on;
end
hold off;
xlabel('Number of bins in search space');
ylabel('Time in seconds');
legend('MUSIC basic QR','MUSIC hess QR','MUSIC implicit QR', ...
    'fast MUSIC fft','fast MUSIC dft');
title(strcat('Order of autocorrelation matrix =', num2str(M)));

figure;
for k = 1:nmethods
    plot(nbins, log10(err(:,k)+eps));hold on;grid on;
end
hold off;
xlabel('Number of bins in search space');
ylabel('Estimation error (log_{10})');
legend('MUSIC basic QR','MUSIC hess QR','MUSIC implicit QR', ...
    'fast MUSIC fft','fast MUSIC dft');
title(strcat('Order of autocorrelation matrix =', num2str(M)));


