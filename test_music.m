%Script to test MUSIC and fast_MUSIC 

close all, clear all, clc;

%number of available data points
N = 8000;
n = 0:N-1;
%clean signal
%example 1
%y = cos(2*0.4*pi.*n + 0.1*pi) + 0.5*cos(2*0.5*pi.*n+0.3*pi) + 0.2*cos(2*0.6*pi.*n);
%example 2
%y = cos(2*0.24*pi.*n) + 0.5*cos(2*0.26*pi.*n + 0.25*pi);
%example 3
y = cos(0.04.*n) + 0.5*cos(0.05.*n);
%normalize signal power to 0dB
y_norm = y./max(abs(y));
snr = 10;
%signal+noise
x = awgn(y_norm, snr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fast MUSIC - split radix does not work - see what's going on
freqs_fast = fast_music(x,2,2000,'default', 'fft');
sort(freqs_fast)

%MUSIC
freqs = music(x,2,2000,'default','fft',200);
sort(freqs)

%%QIFFT
freqs = qifft(x,4096,'black',5,2);
sort(freqs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%measuring computation speed and accuracy
%nbins = 50:50:500;
nbins = 500;
nsig = 2;
nmethods = 6;
M = 50:150:3500;
L = length(M);
t = zeros(L, nmethods);
err = zeros(L,nmethods);
freqs = zeros(L, nmethods, 2*nsig);
%sig_freqs = [-0.6, -0.5, -0.4, 0.4, 0.5, 0.6];
sig_freqs = [-0.26,-0.24,0.24,0.26]*2*pi;
%sig_freqs = [-0.05,-0.04,0.04,0.05];
f = zeros(1,2*nsig);

%frequencies detected by MUSIC with basic QR
for n = 1:L
    tic;
    freqs(n,1,:) = music(x, nsig, nbins, 'gram_schmidt','fft',M(n));
    t(n,1) = toc;
    f(1,:) = freqs(n,1,:);
    err(n,1) = norm(sort(f) - sig_freqs);
end

%this is needed so all variables used by the function are cleared
clearvars -except x nbins nsig nmethods t err sig_freqs snr M freqs f L

%frequencies detected by MUSIC with hessenberg QR 
for n = 1:L
    tic;
    freqs(n,2,:) = music(x, nsig, nbins, 'hess','fft',M(n));
    t(n,2) = toc;
    f(1,:) = freqs(n,2,:);
    err(n,2) = norm(sort(f) - sig_freqs);
end

clearvars -except x nbins nsig nmethods t err sig_freqs snr M freqs f L

for n = 1:L
    tic;
    freqs(n,3,:) = music(x, nsig, nbins, 'implicit','fft',M(n));
    t(n,3) = toc;
    f(1,:) = freqs(n,3,:);
    err(n,3) = norm(sort(f) - sig_freqs);
end

clearvars -except x nbins nsig nmethods t err sig_freqs snr M freqs f L

%frequencies detected by fast_MUSIC with fft
for n = 1:L
    tic;
    freqs(n,4,:) = fast_music(x, nsig, nbins, 'mixed_radix', 'fft',M(n));
    t(n,4) = toc;
    f(1,:) = freqs(n,4,:);
    err(n,4) = norm(sort(f) - sig_freqs);
end

clearvars -except x nbins nsig nmethods t err sig_freqs snr M freqs f L

for n = 1:L
    tic;
    freqs(n,5,:) = fast_music(x, nsig, nbins, 'resample_split_radix','fft',M(n));
    t(n,5) = toc;
    f(1,:) = freqs(n,5,:);
    err(n,5) = norm(sort(f) - sig_freqs);
end

clearvars -except x nbins nsig nmethods t err sig_freqs snr M freqs f L

%frequencies detected by fast_MUSIC with dft
for n = 1:L
    tic;
    freqs(n,6,:) = fast_music(x, nsig, nbins, 'dft', 'fft',M(n));
    t(n,6) = toc;
    f(1,:) = freqs(n,6,:);
    err(n,6) = norm(sort(f) - sig_freqs);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

markers = 'odv*xs';
%plot computation time and MSE plots
figure;
for k = 1:nmethods
    %plot(nbins, t(:,k));hold on;grid on;
    p = plot(M, log(t(:,k)), strcat('-',markers(k)));hold on;grid on;
    p.MarkerSize = 8;
end
hold off;
xlabel('Order of autocorrelation matrix');
%xlabel('Number of bins in search space');
ylabel('Time in seconds (log)');
legend('MUSIC basic QR','MUSIC hess QR','MUSIC implicit QR', ...
    'fast MUSIC mixed radix fft','fast MUSIC resampled split radix fft',...
    'fast MUSIC dft');
%title(strcat('Order of autocorrelation matrix =', num2str(M)));
title(strcat('Number of bins in search space =', num2str(nbins)));

figure;
for k = 1:nmethods
    %plot(nbins, log10(err(:,k)+eps));hold on;grid on;
    p = plot(M, log10(err(:,k)+eps), strcat('-',markers(k)));hold on;grid on;
    p.MarkerSize = 8;
end
hold off;
%xlabel('Number of bins in search space');
xlabel('Order of autocorrelation matrix');
ylabel('Mean squared error in Hz (log_{10})');
legend('MUSIC basic QR','MUSIC hess QR','MUSIC implicit QR', ...
    'fast MUSIC mixed radix fft','fast MUSIC resampled split radix fft',...
    'fast MUSIC dft');
%title(strcat('Order of autocorrelation matrix =', num2str(M)));
title(strcat('Number of bins in search space =', num2str(nbins)));


