%Script to test MUSIC algorithm

close all,clear all,clc;
n = 0:2000;
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

freqs = music(x, 5, 1000);
%frequencies detected by MUSIC
freqs/pi