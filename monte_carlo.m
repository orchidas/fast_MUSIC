%Monte Carlo simulations to get MSE
close all;

nsims = 1000;
%uniformly sampled random phase between [-pi,pi]
phi = -pi + 2*pi*rand(nsims,1);
N = 2500;
n = 0:N-1;
nsig = 2;
snr = -20:2:60;
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
sig_freqs = [0.004,0.005]*2*pi;

theta = [0.004, 1, 0, 0.005, 0.5, 0];
snr_s = zeros(nsig, length(snr));

for k = 1:length(snr)
    
    sigma_z = 10^(-snr(k)/10);
    for l = 1:nsims
        y = cos(2*pi*sig_freqs(1).*n) + 0.5*cos(2*pi*sig_freqs(2).*n + phi(l));
        theta(end) = phi(l);
        y_norm = y./max(abs(y));
        x = awgn(y_norm, snr(k));
        bounds(l,:) = crb(nsig,N,theta,sigma_z);
        [~,freqs_fmusic,M] = fast_music(x,1, nsig, nbins, 'default', 'fft','');
        freqs_fmusic = sort(freqs_fmusic);
        [peaks,freqs_music] = music(x,1, nsig, nbins, 'default','fft','',200);
        freqs_music = sort(freqs_music);
        [peaks,freqs_qifft] = qifft(x,1,4096,'win',5,nsig);
        freqs_qifft = sort(freqs_qifft);
        for m = 1:nsig
            err_music(l,m) = freqs_music(m) - sig_freqs(m)/(2*pi);
            err_fmusic(l,m) = freqs_fmusic(m) - sig_freqs(m)/(2*pi);
            err_qifft(l,m) = freqs_qifft(m) - sig_freqs(m)/(2*pi);
        end
    end
    
    for m = 1:nsig
        %total error in DFS bins
        mse_music(m,k) = mean(err_music(:,m).^2);
        mse_fmusic(m,k) = mean(err_fmusic(:,m).^2);
        mse_qifft(m,k) = mean(err_qifft(:,m).^2);
    end
    crb_bounds(:,k) = mean(bounds,1);
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot MSE with CRB

for m = 1:nsig
    figure(m);
    plot(snr, 10*log10((crb_bounds(3*(m-1)+1,:))+eps), '-s','MarkerSize',8,...
        'MarkerIndices',1:5:length(snr));grid on;hold on;
    plot(snr, 10*log10(mse_music(m,:)+eps),'-d','MarkerSize',8,...
        'MarkerIndices',1:5:length(snr));grid on;hold on;
    plot(snr, 10*log10(mse_fmusic(m,:)+eps),'-x','MarkerSize',8,...
        'MarkerIndices',1:5:length(snr));grid on;hold on;
    plot(snr, 10*log10(mse_qifft(m,:)+eps),'-v','MarkerSize',8,...
        'MarkerIndices',1:5:length(snr));grid on;hold off;
    xlabel('SNR in dB');ylabel('Mean squared error (dB)');
    title(strcat('Frequency of sinusoid = ',num2str(sig_freqs(m)/(2*pi)),'Hz'));
    legend('CRB','MUSIC','fast MUSIC','QIFFT');
    legend('CRB','MUSIC','fast MUSIC');
    %save figure as eps file
    print(strcat('monte_carlo_',num2str(sig_freqs(m))/(2*pi),'Hz'),'-deps');
end



        
