function [peaks,freqs,M] = fast_music(x, fs, nsignals, nbins, method_eig, method_autocorr,file,M)

%Replace eigenvalue decomposition in MUSIC with FFT
%x - signal
%fs - sampling frequency, required for plotting pseudospectrum
%nsignals - number of real sinusoids in signal
%nbins - number of points in search space
%method_eig - calculate eigenvalues with dft or fft
%method_autocorr - method for estimating autocorrelation function, direct
%or fft
%M - autocorrelation matrix order (ideally should be calcuated from ACF
%periodicity, but included just for plotting accuracy vs M).

if nargin == 4
    method_eig = 'default';
elseif nargin == 5
    method_autocorr = 'fft';
end

N = length(x);
%estimate autocorrelation function
R = estimate_autocorrelation_function(x, N, method_autocorr);
%take last half of autocorrelation only
R = R(N+1:end);
%resampling will shift the spectrum
shift = 1;

%M is the number of antenna, or the dimension of the autocorrelation matrix
%in our case.
if nargin == 7
    period = find_periodicity(R,0.05);
    %period = 198;
    %take more periods for better estimation
    M =period*floor(N/period);
    %M = period;
    %if signal is not periodic, or too short to be periodic
    if(M < 1)
        M = N;
    end
end
R = R(1:M);


%since we know for order = M, autocorrelation matrix will be circulant - 
%no need to do explicit eigenvalue decomposition, just multiply 
%autocorrelation function with DFT matrix

if strcmp(method_eig,'dft')
     dftm = exp(-1j*2*pi/M*(0:M-1)'*(0:M-1));
     eigvals = dftm * R';
     
%direct multiplication with DFT matrix is expensive -- 
%reduce computation time by using FFT

elseif strcmp(method_eig, 'default')
    %Matlab's built-in FFT (much faster)
    eigvals = fft(R',M);
elseif strcmp(method_eig,'mixed_radix')
    %My mixed radix FFT
    eigvals = mixed_radix_fft(R', M);
elseif strcmp(method_eig,'resample_split_radix')
    %My resampled split radix fft
    eigvals = resamp_splitradix_fft(R',2^nextpow2(M), M);
    %oversampling factor
    shift = 2^nextpow2(M)/M;
    M = 2^nextpow2(M);
end

[eig_vals_sorted, inds] = sort(abs(eigvals),'descend');

p = 2*nsignals;
noise_eigvals_pos = inds(p+1:M);
sig_eigvals_pos = inds(1:p);

%since the signal is real, our search space can be over positive
%frequencies only
%k = -nbins/2+1:nbins/2;
k = 0:nbins/2-1;
P = zeros(nbins/2,1);
meps = 10^-6;

%alternative pseudospectrum estimate from closed-form solution
for m = 1:length(k)
    
%    %vectorized code
%     curn = noise_eigvals_pos-1;
%     curk = k(m);
%     pos = find(curn*nbins == curk*M);
%     if(~isempty(pos))
%         P(m) =  M;
%         curn = noise_eigvals_pos([1:pos-1 pos+1:end])-1;
%     end
%     P(m) = P(m) +  (sum(abs(sin(pi.*(curk/nbins - curn/M)*M)./...
%            sin(pi.*(curk/nbins - curn/M)))));
%     P(m) = 1./P(m);

     %noise subspace approximation
%      curn = noise_eigvals_pos-1;
%      curk = k(m);
%      P(m) = M./(sum((abs(sin(pi.*(curk/nbins - curn/M)*M)./...
%            sin(pi.*(curk/nbins - curn/M)))).^2));
%      if isnan(P(m))
%         P(m) = 1/(M*(M-p));
%      end  
    
    %signal subspace approximation
    curn = sig_eigvals_pos-1;
    curk = k(m);
    P(m) = 1./(M - (1/M)*sum((sin(pi.*(curk/nbins - curn/M)*M)./...
           sin(pi.*(curk/nbins - curn/M))).^2));
    if isnan(P(m))
        P(m) = 1/meps;
    end
    
end

%shift spectrum if resampled 
P = P * shift;
%frequency estimates
[peaks,freqs] = find_peaks(P,nsignals);
freqs = (freqs-1)/length(P)*fs/2;
 
h = figure;
plot(k/(nbins/2) * (fs/2), P);hold on;grid on;
plot(freqs, peaks, '*');hold off;grid on;
xlim([0,0.01]);
%xlim([2400,2900]);ylim([0,1.1*max(peaks)]);
%xlim([2660,2680]);
ylabel('Pseudospectrum','fontsize',16);
xlabel('Frequency in Hz','fontsize',14);
%title(strcat('Fast MUSIC ', file));
%print(strcat('../Figures/',file,'.eps'),'-deps');
%savefig(h,strcat('../piano data/A3/fmusic-',file,'.fig'));


end

