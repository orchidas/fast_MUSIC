%Script to test MUSIC algorithm

close all,clear all, clc;
n = 0:1999;
%clean signal
%example 1
y = cos(0.4*pi.*n + 0.1*pi) + 0.5*cos(0.5*pi.*n+0.3*pi) + 0.2*cos(0.6*pi.*n);
%example 2
%y = cos(0.5*pi.*n) + cos(0.52*pi.*n + 0.25*pi);
%nornalize signal power to 0dB
y_norm = y./max(abs(y));
snr = 10;
%signal+noise
x = awgn(y_norm, snr);

%fast MUSIC
freqs_fast = fast_music(x,3,500,'fft','fft');
sort(freqs_fast/pi)

%MUSIC
freqs = music(x,3,500,'fft');
sort(freqs/pi)

%measuring computation speed and accuracy
nbins = 50:100:1500;
nsig = 3;
t1 = zeros(length(nbins),1);
t2 = zeros(length(nbins),1);
err1 = zeros(length(nbins),1);
err2 = zeros(length(nbins),1);
sig_freqs = [-0.6, -0.5, -0.4, 0.4, 0.5, 0.6];

%frequencies detected by fast_MUSIC
for n = 1:length(nbins)
    tic;
    freqs_fast = fast_music(x, nsig, nbins(n), 'fft', 'fft');
    t1(n) = toc;
    err1(n) = norm(sort(freqs_fast/pi) - sig_freqs);
end

%this is needed so all variables used by the function are cleared
clearvars -except x nbins nsig t1 err1 sig_freqs snr freqs_fast

%frequencies detected by MUSIC
for n = 1:length(nbins)
    tic;
    freqs = fast_music(x, nsig, nbins(n), 'fft');
    t2(n) = toc;
    err2(n) = norm(sort(freqs/pi) - sig_freqs);
end

figure;
plot(nbins, t1);hold on;grid on;
plot(nbins, t2 );hold off;grid on;
xlabel('Number of bins in search space');
ylabel('Time in seconds');
legend('fast MUSIC','MUSIC');

figure;
plot(nbins, log10(err1+eps));hold on;grid on;
plot(nbins, log10(err2+eps));hold off;grid on;
xlabel('Number of bins in search space');
ylabel('Estimation error (logarithm)');
legend('fast MUSIC','MUSIC');


% %let's time computation speeds
% 
% for n= 1:length(nbins);
%   f1 = @()music(x,3,500,'fft');
%   t1(n) = timeit(f1);
%  
%   f2 = @()fast_music(x,3,500,'fft', 'fft');
%   t2(n) = timeit(f2);
% end
