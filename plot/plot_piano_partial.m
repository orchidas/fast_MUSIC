% this script reads in an audio file containing a piano note,
% bandpass filters it so that it only has the 10th partial, then runs
% FAST MUSIC on it to detect the peaks (there should be 2). It overlays it 
% on top of the FFT plot.

%clear all;
close all;
path = '../piano data/A3/';
file = 'A3.wav';
fs = 44100;
[b,a] = butter(4,[2500 2900]/(fs/2));
ftpower = 14:18;
nsamps = round(2.^ftpower);
nns = length(nsamps);
nbins = 2^19;
tfs = zeros(2*nbins,nns);
f = (0:nbins)/nbins*fs/2;
freqs_foundfm_L = zeros(nns,1);
freqs_foundfm_U = zeros(nns,1);
npeaks = 2;

%read signal
[x,fs] = audioread(strcat(path,file));
%filter signal
y = filtfilt(b,a,x);

for n = 1:nns
    nsamp = nsamps(n);
    ys = y(fs+(1:nsamp));
    tfs(:,n) = fft(ys, 2*nbins);
    [peaksfm, freqs_foundfm] = fast_music(ys',fs,npeaks,500000,'default','fft',...
        strcat(file,'_',num2str(ftpower(n))));
    freqs_foundfm = sort(freqs_foundfm);
    freqs_foundfm_L(n) = freqs_foundfm(1);
    freqs_foundfm_U(n) = freqs_foundfm(2);

end

% figure(1); 
% plot(f, 20*log10(abs(tfs(1:nbins+1,:))/max(abs(tfs(:))))); grid; hold on;
% %plot([freqs_foundfm_L';freqs_foundfm_L'],[-60;0]*ones(1,nns),'-o');hold on;
% %plot([freqs_foundfm_U';freqs_foundfm_U'],[-60;0]*ones(1,nns),'-o');hold off;
% ylim([-60 0]);
% xlim([2660,2680]);

%load peak_plots;
% heights = zeros(nns,npeaks);
for n = 1:nns
    %heights(n,1) = interp1(f',20*log10(abs(tfs(1:nbins+1,n))/max(abs(tfs(:)))), freqs_foundfm_L(n));
    %heights(n,2) = interp1(f',20*log10(abs(tfs(1:nbins+1,n))/max(abs(tfs(:)))), freqs_foundfm_U(n));
    figure(n);
    plot(f',20*log10(abs(tfs(1:nbins+1,n))/max(abs(tfs(:)))));hold on;grid;
    plot([freqs_foundfm_L(n);freqs_foundfm_L(n)],[-60;0],'b-o');hold on;
    plot([freqs_foundfm_U(n);freqs_foundfm_U(n)],[-60;0],'b-o');hold off;
    title(strcat('Data size = ', num2str(nsamps(n))));
    xlabel('Frequency in Hz','fontsize',16);
    ylabel('Magnitude in dB','fontsize',16);
    ylim([-60 0]);
    xlim([2660,2680]);
end

%plot([freqs_foundfm_L';freqs_foundfm_U'], [heights(:,1)';heights(:,2)'],'o');
%plot(freqs_foundfm_U, heights(:,2),'o');hold off;
