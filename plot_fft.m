

close all;
%clear all;
path = '../piano data/A3/';
files = {'A3.wav    ','A3_1D.wav ','A3_2D.wav ','A3_3D.wav ','A3_23D.wav','A3_13D.wav','A3_12D.wav'};
fs = 44100;
nfiles = length(files);
[b,a] = butter(4,[2500 2900]/(fs/2));
nsamp = 6*fs;
window = blackman(nsamp);
%nbins = 2^nextpow2(nsamp);
nbins = 2^19;
tf = zeros(2*nbins,nfiles);
npeaks = 2;
all_freqs_found_qifft = zeros(npeaks,nfiles);
all_freqs_found_fmusic = zeros(npeaks,nfiles);
fmusic_size = 2^14;
f = (0:nbins)/nbins*fs/2;


for k = 1:nfiles
 
    [x,fs] = audioread(strcat(path,files{k}));  
    y = filtfilt(b,a,x);
    ys = y(fs+(1:nsamp));
    xs = x(fs+(1:nsamp));
    tf(:,k) = fft(ys.*window,2*nbins);% y(fs+(1:nsamp)).*window],2*nbins);
    figure(1); plot(f, 20*log10(abs(tf(1:nbins+1,k))/max(abs(tf(:,k))))); grid; 
    ylim([-60 0]);
    %xlim([2660,2680]);
    xlim([2400,2900]);
    %xlim([2.66,2.68]);
    hold on;

    %plot spectrogram to see beats
    myspecgram(xs',fs,2047,512,4096);
    ylim([0,5000]);
    
    %do QIFFT
    [peaks,freqs_found] = qifft(ys',fs,nbins,'black',5,npeaks);
    freqs_found_Hz = (freqs_found/pi) * (fs/2);
    all_freqs_found_qifft(k,1:npeaks) = freqs_found_Hz;
    figure(1);
    plot(freqs_found_Hz, peaks, '*'); grid on;hold on;
    
    %do fast Music
    [peaksfm, freqs_foundfm] = fast_music(ys(1:fmusic_size)',fs,npeaks,50000,'default','fft',files{k});
    freqs_foundfm_Hz = freqs_foundfm/pi * (fs/2);
    all_freqs_found_fmusic(k,1:npeaks) = freqs_foundfm_Hz;
    
    figure(1);plot(freqs_foundfm_Hz,zeros(2,1),'o');
    

end

% figure(1); plot((0:nbins)/nbins*fs/2000,20*log10(abs(tf(1:nbins+1,:))/max(abs(tf(:))))); 
% ylim([-80 0]);
% xlim([2.66,2.68]);
% grid on;hold off;
legend(files);