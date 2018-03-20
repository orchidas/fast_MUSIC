%test tuning of piano strings
close all,clc;
path = '../piano data/A3/';
%files = {'A3.wav    ', 'A3_23D.wav','A3_13D.wav','A3_12D.wav'};
files = {'A3.wav    ','A3_1D.wav ','A3_2D.wav ', 'A3_3D.wav ', 'A3_23D.wav','A3_13D.wav', 'A3_12D.wav'};
%files = {'C4.wav    ','C4_1D.wav ','C4_2D.wav ', 'C4_3D.wav ', 'C4_23D.wav','C4_13D.wav', 'C4_12D.wav'};
fftsize = 8192;
fc = 300;
fs = 44100;
start = 1.5*fs;
stop = start + fftsize -1;
%second order butterworth lowpass filter 
[b,a] = butter(6,fc/(fs/2));
%create frequency axis
freqs = linspace(0, fs/2, fftsize/2);
[v,cutoff] = min(abs(freqs - fc));
npeaks = 5;
%all_freqs_found_qifft = [];
%all_freqs_found_music = [];
%all_freqs_found_fmusic = [];

all_freqs_found_qifft = zeros(length(files),npeaks);
all_freqs_found_fmusic = zeros(length(files),npeaks);
all_freqs_found_music = zeros(length(files),npeaks);


for k = 1:length(files)
    
    [x,fs] = audioread(strcat(path,files{k}));   
    %take 4096 samples from middle of the file (ignore attack)
    xs = x(start:stop);
    
    %Do a simple FFT
    x_win = xs .* hann(fftsize);
    %take fft
    X = fftshift(fft(x_win,fftsize));
    X = abs(X)/max(abs(X));
    %preserve positive half of spectrum only
    X = X(fftsize/2+1:end);
    %convert to dB
    Xmag = 20*log10(X);
    %plot FFT of signal
    figure;
    plot(freqs(1:cutoff),Xmag(1:cutoff));hold on;grid on;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %filter X so that it only has a fundamental
    ys = filter(b,a,xs);
    y_win = ys .* blackman(fftsize);
    %take fft
    Y = fftshift(fft(y_win,fftsize));
    Y = abs(Y)/max(abs(Y));
    %preserve positive half of spectrum only
    Y = Y(fftsize/2+1:end);
    %convert to dB
    Ymag = 20*log10(Y);
    %plot FFT of filtered signal 
    plot(freqs(1:cutoff),Ymag(1:cutoff),'r');hold on;grid on;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %apply QIFFT to it
    %should see a total of 3 peaks for coupled strings
    if(k == 1)
        npeaks = 3;
    elseif(k > 1 && k <= 4)
        npeaks = 2;
    %should see just 1 peak for uncoupled strings
    else
        npeaks = 1;
    end

    
    [peaks,freqs_found] = qifft(ys',fs,fftsize,'black',5,npeaks);
    freqs_found_Hz = (freqs_found/pi) * (fs/2);
    all_freqs_found_qifft(k,1:npeaks) = freqs_found_Hz;
    
    %plot QIFFT
    plot(freqs_found_Hz, peaks, 'k*'); grid on;hold off;
    xlabel('Frequency in Hz');ylabel('DFT magnitude in decibels');
    axis([0, fc+100, -80, max(peaks)+10]);    
    legend('Original signal','Filtered signal','QIFFT peaks');
    title(files{k});
    
    %apply MUSIC
    [peaksm, freqs_foundm] = music(ys',fs,npeaks,20000,'default','fft',files{k});
    freqs_foundm_Hz = freqs_foundm/pi * (fs/2);
    all_freqs_found_music(k,1:npeaks) = freqs_foundm_Hz;
    
    %apply fast MUSIC
    [peaksfm, freqs_foundfm] = fast_music(ys',fs,npeaks,20000,'default','fft',files{k});
    freqs_foundfm_Hz = freqs_foundfm/pi * (fs/2);
    all_freqs_found_fmusic(k,1:npeaks) = freqs_foundfm_Hz;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

save(strcat(path,'A3'),'all_freqs_found_qifft','all_freqs_found_music','all_freqs_found_fmusic');