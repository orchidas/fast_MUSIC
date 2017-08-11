%Script to plot MSE vs SNR with CRB
close all, clc;

N = 1000;
n = 0:N-1;
%clean signal
%example 1
%y = cos(2*0.4*pi.*n + 0.1*pi) + 0.5*cos(2*0.5*pi.*n+0.3*pi) + 0.2*cos(2*0.6*pi.*n);
%example 2
y = cos(2*0.25*pi.*n) + 0.5*cos(2*0.26*pi.*n + 0.25*pi);
%example 3
%y = cos(0.04.*n) + 0.5*cos(0.05.*n + 0.25*pi);

%theta = [0.4,1,0.1*pi,0.5,0.5,0.3*pi,0.6,0.2,0];
theta = [0.25,1,0,0.26,0.5,0.25*pi];
%theta = [0.04/(2*pi), 1, 0, 0.05/(2*pi), 0.5, 0.25*pi];
%nornalize signal power to 0dB
y_norm = y./max(abs(y));
nsig = 2;

snr = -20:5:50;
nbins = 1000;
bounds = zeros(3*nsig,length(snr));
err_music = zeros(2,length(snr));
err_fmusic = zeros(2,length(snr));
%sig_freqs = [-0.05,-0.04,0.04,0.05];
sig_freqs = [-0.26,-0.25,0.25,0.26]*2*pi;


for k = 1:length(snr)
    %signal+noise
    x = awgn(y_norm, snr(k));
    sigma_z = var(y_norm - x);
    %get CRB
    bounds(:,k) = crb(nsig,N/2,theta,sigma_z);
    freqs_music = sort(music(x, nsig, nbins, 'default','fft'));
    freqs_fmusic = sort(fast_music(x, nsig, nbins, 'default', 'fft'));
    err_music(1,k) = (freqs_music(3) - sig_freqs(3))^2;
    err_fmusic(1,k) = (freqs_fmusic(3) - sig_freqs(3))^2;
    err_music(2,k) = (freqs_music(4) - sig_freqs(4))^2;
    err_fmusic(2,k) = (freqs_fmusic(4) - sig_freqs(4))^2;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot MSE with CRB
figure;
plot(snr, log10(bounds(1,:)+eps));grid on;hold on;
plot(snr, log10(err_music(1,:)+eps));grid on;hold on;
plot(snr, log10(err_fmusic(1,:)+eps));grid on;hold off;
xlabel('SNR in dB');ylabel('log_{10} Mean squared error');
title('Frequency of sinusoid = 0.25Hz');
legend('CRB','MUSIC','fast MUSIC');

figure;
plot(snr, log10(bounds(4,:)+eps));grid on;hold on;
plot(snr, log10(err_music(2,:)+eps));grid on;hold on;
plot(snr, log10(err_fmusic(2,:)+eps));grid on;hold off;
xlabel('SNR in dB');ylabel('log_{10} Mean squared error');
title('Frequency of sinusoid = 0.26Hz');
legend('CRB','MUSIC','fast MUSIC');


    
    
