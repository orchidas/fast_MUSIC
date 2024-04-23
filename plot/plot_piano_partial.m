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

%load peak_plots;
for n = 1:nns
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

