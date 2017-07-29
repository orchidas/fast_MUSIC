%Script to test MUSIC algorithm

close all, clc;
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

%frequencies detected by MUSIC
freqs = music(x, 3, 500, 'fft');
sort(freqs/pi)

%frequencies detected by fast_MUSIC
freqs_fast = fast_music(x,3,500,'fft');
sort(freqs_fast/pi)

%let's time computation speeds
% f1 = @()music(x,3,1000);
% t1 = timeit(f1)
% 
% f2 = @()fast_music(x,3,1000);
% t2 = timeit(f2)
