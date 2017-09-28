function [freqs] = qifft(x,N,win,zpf,npeaks)
%Quadratically Interpolated FFT (QIFFT)
%method for estimating sinusoidal parameters from
%peaks in spectral magnitude data

%x - block of signal
%N - length of fft (power of 2)
%win - type of window, rect, hann, ham, black
%zpf - zero padding factor
%npeaks - number of peaks in spectrum

%length of fft
%window length
M = round(N/zpf);
w = zeros(M,1);
freqs = zeros(1,npeaks);

if(strcmp(win,'rect'))
    w = ones(M,1);
elseif(strcmp(win,'ham'))
    w = hamming(M);
elseif(strcmp(win,'hann'))
    w = hann(M);
else
    w = blackman(M);
end

%step 1 - Calculate the amplitude and phase spectrum
%of audio data, by using an appropriately zero-
%padded FFT with an appropriate window of an appropriate length
xwin = x(1:M).*w';
X = fftshift(fft(xwin,N));
Xmag = abs(X)/max(abs(X));
%since signal is real, spectrum is symmetric. Retain half of it.
Xmag = 20*log10(Xmag(N/2+1:end));
fbins = linspace(0,pi,N/2);

for k = 1:npeaks
    %Find the bin number of the maximum peak magnitude
    [val,kmax] = max(Xmag);

    %Quadratically interpolate the log-amplitude of
    %the peak using two neighboring samples
    if (kmax > 1 && kmax < N/2-1)
        [A_hat, del_hat] = parabolic_interpolation(Xmag(kmax-1), Xmag(kmax),...
        Xmag(kmax+1));
    else
        del_hat = 0;
    end

    %Estimate the peak frequency in bins
    freqs(k) = kmax + del_hat;
    freqs(k) = (freqs(k)-1)*(2*pi/N);
    
%     figure;
%     plot(fbins/pi,Xmag);hold on;grid on;
%     plot(freqs(k)/pi,A_hat,'*');hold off;grid on;
%     title('Magnitude spetrum');xlabel('Frequency in radians/pi');
    
    %Subtract the peak from the FFT data for sub-
    %sequent processing.
    start = find(Xmag(1:kmax-1) < -6, 1,'last');
    stop = find(Xmag(kmax+1:end) < -6, 1,'first');
    Xmag(start:kmax+stop) = -100;
    
end

%since the data is real,spectrum is symmetric
freqs = [-freqs, freqs];

    function[peak,pos] = parabolic_interpolation(a,b,c)
        %given 3 points, it returns the result of their parabolic interpolation
        pos = 0.5 * ((a-c)/(a - 2*b + c));
        peak = b - 0.25*(a-c)*pos;
    end

end
