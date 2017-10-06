%Monte Carlo simulations to get MSE
close all, clear all, clc;

nsims = 100;
%uniformly sampled random phase between [-pi,pi]
phi = -pi + 2*pi*rand(nsims,1);
N = 1000;
n = 0:N-1;
nsig = 2;
snr = -20:10:100;
nbins = 1000;
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
%sig_freqs = [-2.6,-2.4,2.4,2.6];
sig_freqs = [-0.05,-0.04,0.04,0.05];
%sig_freqs = [-0.26,0.26]*2*pi;
%sig_freqs = [-0.05,0.05];

%theta = [2.4/(2*pi),1,0,2.6/(2*pi),0.5,0];
theta = [0.04/(2*pi), 1, 0, 0.05/(2*pi), 0.5, 0];
%theta = [0.26,1,0];
%theta = [0.05/(2*pi) 1 0];
snr_s = zeros(nsig, length(snr));

for k = 1:length(snr)
    
    sigma_z = 10^(-snr(k)/10);
    for l = 1:nsims
        %y = cos(2.4.*n) + 0.5*cos(2.6.*n + phi(l));
        y = cos(0.04.*n) + 0.5*cos(0.05.*n + phi(l));
        %y = cos(2*pi*0.26.*n + phi(l));
        %y = cos(0.05.*n + phi(l));
        theta(end) = phi(l);
        y_norm = y./max(abs(y));
        x = awgn(y_norm, snr(k));
        bounds(l,:) = crb(nsig,N,theta,sigma_z);
        freqs_music = sort(music(x, nsig, nbins, 'default','fft',300));
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
        %mse_music(m,k) = var(err_music(:,m));
        mse_music(m,k) = mean(err_music(:,m).^2);
        %mse_fmusic(m,k) = var(err_fmusic(:,m));
        mse_fmusic(m,k) = mean(err_fmusic(:,m).^2);
        %mse_qifft(m,k) = var(err_qifft(:,m));
        mse_qifft(m,k) = mean(err_qifft(:,m).^2);
        %CRB in DFS bins - these bounds work if the spacing between 2
        %frequencies is 1.8 bins or more
        %spectrum SNR
        %snr_s(m,k) = N*(theta(3*(m-1)+2)^2)/(sigma_z);
        %crb_bounds(m,k) = 1.5/(pi^2 * snr_s(m,k));
    end
    crb_bounds(:,k) = mean(bounds,1);
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot MSE with CRB

% for m = 1:nsig
%     figure(m);
%     %plot(10*log10(snr_s(m,:)), 10*log10(crb_bounds(3*(m-1)+1,:)*(N^2)+eps));grid on;hold on;
%     plot(10*log10(snr_s(m,:)), 10*log10(crb_bounds(m,:)+eps));grid on;hold on;
%     plot(10*log10(snr_s(m,:)), 10*log10(mse_music(m,:)+eps));grid on;hold on;
%     plot(10*log10(snr_s(m,:)), 10*log10(mse_fmusic(m,:)+eps));grid on;hold on;
%     plot(10*log10(snr_s(m,:)), 10*log10(mse_qifft(m,:)+eps));grid on;hold off;
%     xlabel('Spectrum SNR in dB');ylabel('Mean squared error in bins (dB)');
%     title(strcat('Frequency of sinusoid = ',num2str(sig_freqs(nsig+m)),'rad'));
%     %axis([-20,80,-100,40]);
%     legend('CRB','MUSIC','fast MUSIC','QIFFT');
% end

for m = 1:nsig
    figure(m);
    plot(snr, 10*log10((crb_bounds(3*(m-1)+1,:))+eps), '-s','MarkerSize',8);grid on;hold on;
    plot(snr, 10*log10(mse_music(m,:)+eps),'-d','MarkerSize',8);grid on;hold on;
    plot(snr, 10*log10(mse_fmusic(m,:)+eps),'-x','MarkerSize',8);grid on;hold on;
    plot(snr, 10*log10(mse_qifft(m,:)+eps),'-v','MarkerSize',8);grid on;hold off;
    xlabel('SNR in dB');ylabel('Mean squared error (dB)');
    %title(strcat('Frequency of sinusoid = ',num2str(sig_freqs(nsig+m)),'rad'));
    legend('CRB','MUSIC','fast MUSIC','QIFFT');
    %save figure as eps file
    print(strcat('monte_carlo_',num2str(sig_freqs(nsig+m)),'rad'),'-deps');
end



        