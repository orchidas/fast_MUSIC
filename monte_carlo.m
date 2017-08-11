%Monte Carlo simulations to get MSE
close all, clear all, clc;

nsims = 100;
%uniformly sampled random phase between [-pi,pi]
phi = -pi + 2*pi*rand(nsims,1);
N = 400;
n = 0:N-1;
nsig = 2;
snr = -20:5:50;
nbins = 500;
bounds = zeros(nsims,3*nsig);
err_music = zeros(nsims,2);
err_fmusic = zeros(nsims,2);
mse_music = zeros(2,length(snr));
mse_fmusic = zeros(2,length(snr));
sig_freqs = [-0.25,-0.26,0.25,0.26]*2*pi;
%sig_freqs = [-0.05,-0.04,0.04,0.05];
crb_bounds = zeros(3*nsig,length(snr));

for k = 1:length(snr)
    sigma_z = 10^(-snr(k)/10);
    for l = 1:nsims
        y = cos(2*0.25*pi.*n) + 0.5*cos(2*0.26*pi.*n + phi(l));
        %y = cos(0.04.*n) + 0.5*cos(0.05.*n + phi(1));
        theta = [0.25,1,0,0.26,0.5,phi(l)];
        %theta = [0.04/(2*pi), 1, 0, 0.05/(2*pi), 0.5, phi(l)];
        y_norm = y./max(abs(y));
        x = awgn(y_norm, snr(k));
        %get CRB
        bounds(l,:) = crb(nsig,N/2,theta,sigma_z);
        freqs_music = sort(music(x, nsig, nbins, 'default','fft',200));
        freqs_fmusic = sort(fast_music(x, nsig, nbins, 'default', 'fft'));
        err_music(l,1) = (freqs_music(3) - sig_freqs(3))/(2*pi);
        err_fmusic(l,1) = (freqs_fmusic(3) - sig_freqs(3))/(2*pi);
        err_music(l,2) = (freqs_music(4) - sig_freqs(4))/(2*pi);
        err_fmusic(l,2) = (freqs_fmusic(4) - sig_freqs(4))/(2*pi);
    end
    mse_music(1,k) = norm(err_music(:,1));
    mse_music(2,k) = norm(err_music(:,2));
    mse_fmusic(1,k) = norm(err_fmusic(:,1));
    mse_fmusic(2,k) = norm(err_fmusic(:,2));
    crb_bounds(:,k) = mean(bounds,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot MSE with CRB
figure;
plot(snr, log10(crb_bounds(1,:)+eps));grid on;hold on;
plot(snr, log10(mse_music(1,:)+eps));grid on;hold on;
plot(snr, log10(mse_fmusic(1,:)+eps));grid on;hold off;
xlabel('SNR in dB');ylabel('log_{10} Mean squared error');
title('Frequency of sinusoid = 0.25Hz');
legend('CRB','MUSIC','fast MUSIC');

figure;
plot(snr, log10(crb_bounds(4,:)+eps));grid on;hold on;
plot(snr, log10(mse_music(2,:)+eps));grid on;hold on;
plot(snr, log10(mse_fmusic(2,:)+eps));grid on;hold off;
xlabel('SNR in dB');ylabel('log_{10} Mean squared error');
title('Frequency of sinusoid = 0.26Hz');
legend('CRB','MUSIC','fast MUSIC');


        