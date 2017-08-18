%Monte Carlo simulations to get MSE
close all, clear all, clc;

nsims = 100;
%uniformly sampled random phase between [-pi,pi]
phi = -pi + 2*pi*rand(nsims,1);
N = 2000;
n = 0:N-1;
nsig = 2;
snr = -20:5:50;
nbins = 2000;
bounds = zeros(nsims,3*nsig);
%crb is in Hz
crb_bounds = zeros(3*nsig,length(snr));
%crb_bounds = zeros(nsig, length(snr));
%the following should be in Hz
err_music = zeros(nsims,nsig);
err_fmusic = zeros(nsims,nsig);
err_qifft = zeros(nsims,nsig);
mse_music = zeros(nsig,length(snr));
mse_fmusic = zeros(nsig,length(snr));
mse_qifft = zeros(nsig,length(snr));
%this is in radians
%sig_freqs = [-0.26,-0.24,0.24,0.26]*2*pi;
sig_freqs = [-0.05,-0.04,0.04,0.05];
%sig_freqs = [-0.26,0.26]*2*pi;
%theta = [0.24,1,0,0.26,0.5,0];
theta = [0.04/(2*pi), 1, 0, 0.05/(2*pi), 0.5, 0];
%theta = [0.26 1 0];
snr_s = zeros(nsig, length(snr));

for k = 1:length(snr)
    
    sigma_z = 10^(-snr(k)/10);
    for l = 1:nsims
        %y = cos(2*0.24*pi.*n) + 0.5*cos(2*0.26*pi.*n + phi(l));
        y = cos(0.04.*n) + 0.5*cos(0.05.*n + phi(1));
        theta(end) = phi(l);
        y_norm = y./max(abs(y));
        x = awgn(y_norm, snr(k));
        bounds(l,:) = crb(nsig,N,theta,sigma_z);
        freqs_music = sort(music(x, nsig, nbins, 'default','fft',200));
        freqs_fmusic = sort(fast_music(x, nsig, nbins, 'default', 'fft'));
        freqs_qifft = sort(qifft(x,4096,'win',5,nsig));
        for m = 1:nsig
            err_music(l,m) = (freqs_music(nsig + m) - sig_freqs(nsig + m))/(2*pi);
            err_fmusic(l,m) = (freqs_fmusic(nsig + m) - sig_freqs(nsig + m))/(2*pi);
            err_qifft(l,m) = (freqs_qifft(nsig + m) - sig_freqs(nsig + m))/(2*pi);
        end
    end
    
    for m = 1:nsig
        %total error in DFS bins
        mse_music(m,k) = N*norm(err_music(:,m))/nsims;
        mse_fmusic(m,k) = N*norm(err_fmusic(:,m))/nsims;
        mse_qifft(m,k) = N*norm(err_qifft(:,m))/nsims;
        %CRB in DFS bins - these bounds work if the spacing between 2
        %frequencies is 1.8 bins or more
        %spectrum SNR
        snr_s(m,k) = N*(theta(3*(m-1)+2)^2)/(sigma_z);
        %crb_bounds(m,k) = 1.5/(pi^2 * snr_s(m,k));
    end
    crb_bounds(:,k) = mean(bounds,1);
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot MSE with CRB

for m = 1:nsig
    figure(m);
    plot(10*log10(snr_s(m,:)), log10(crb_bounds(3*(m-1)+1,:)*(N^2)+eps));grid on;hold on;
    %plot(10*log10(snr_s(m,:)), log10(crb_bounds(m,:)+eps));grid on;hold on;
    plot(10*log10(snr_s(m,:)), log10(mse_music(m,:)+eps));grid on;hold on;
    plot(10*log10(snr_s(m,:)), log10(mse_fmusic(m,:)+eps));grid on;hold on;
    plot(10*log10(snr_s(m,:)), log10(mse_qifft(m,:)+eps));grid on;hold off;
    xlabel('Spectrum SNR in dB');ylabel('Root mean squared error in bins (loq_{10}');
    title(strcat('Frequency of sinusoid = ',num2str(sig_freqs(nsig+m)),'rad'));
    %axis([-20,80,-100,40]);
    legend('CRB','MUSIC','fast MUSIC','QIFFT');
end



        