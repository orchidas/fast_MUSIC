function [peaks,freqs,M] = fast_music(x, fs, nsignals, nbins, varargin)

%% 
% Replace eigenvalue decomposition in MUSIC with FFT
% Inputs:
% x - signal
% fs - sampling frequency, required for plotting pseudospectrum
% nsignals - number of real sinusoids in signal
% nbins - number of points in search space
% Returns:
% peaks : pseudospectrum peaks values
% freqs : frequencies of the peaks in Hz
% M : periodicity of signal in samples
%%

switch nargin
    case 4
        method_eig = 'default';
        method_autocorr = 'fft';
        file = '';
        plot_spec = 0;
    case 5
        method_eig = varargin{1};
        method_autocorr = 'fft';
        file = '';
        plot_spec = 0;
    case 6
        method_eig = varargin{1};
        method_autocorr = varargin{2};
        file = '';
        plot_spec = 0;
    case 7
        method_eig = varargin{1};
        method_autocorr = varargin{2};
        file = varargin{3};
        plot_spec = 0;
    case 8
        method_eig = varargin{1};
        method_autocorr = varargin{2};
        file = varargin{3};
        plot_spec = varargin{4};
    otherwise
        error('Wrong number of inputs');
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
M = find_periodicity(R,0.05);
%if signal is not periodic, or too short to be periodic
if(M < 1)
    M = N;
end

M = M*floor(N/M);
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


%this is important if peaks are not significantly close to each other,
%then the peak next to the maximum is selected as the second largest but it may
%not be so.
p = 2*nsignals;
[sig_eigvals, sig_eigvals_pos] = find_peaks(abs(eigvals),p, 'n');
for i = 1:p
    %clearly there is no second peak in this case
    if(abs(sig_eigvals(i)) < 0.1*max(abs(sig_eigvals)))
        [~, inds] = sort(abs(eigvals),'descend');    
        %noise_eigvals_pos = inds(p+1:M);
        sig_eigvals_pos = inds(1:p);
    break;
    end 
end

%since the signal is real, our search space can be over positive
%frequencies only
k = 0:nbins/2-1;
P = zeros(nbins/2,1);
meps = 10^-6;

%alternative pseudospectrum estimate from closed-form solution
for m = 1:length(k)
        
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
[peaks,freqs] = find_peaks(P,nsignals,'y');
freqs = (freqs-1)/length(P)*fs/2;
 

if plot_spec
    figure;
    plot(k/(nbins/2) * (fs/2), P);hold on;grid on;
    plot(freqs, peaks, '*');hold off;grid on;
    ylabel('Pseudospectrum','fontsize',16);
    xlabel('Frequency in Hz','fontsize',14);
    title(strcat('Fast MUSIC ', file));
end

end

