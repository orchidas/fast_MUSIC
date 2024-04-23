function [peaks,freqs] = music(x, fs, nsignals, nbins, varargin)

%% 
% MUSIC algorithm for sinusoid parameter estimation.
% Inputs:
% x - signal corrupted with white noise
% fs - sampling frequency - required for plotting pseudospectrum
% nsignals - number of real sinusoids in signal
% nbins - number of bins in search space
% method_eig (optional) - method used for eigenvalue decomposition
% method_autocorr (optional) - method used for calculating autocorrelation
%                               function
% file (optional) - name of the signal to be plotted
% plot_spec (optional) - whether to plot the pseudospectrum
% Returns:
% peaks - array of peak values
% freqs - frequencies at which the peaks are found, in Hz
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

%M is the number of antenna, or the dimension of the autocorrelation matrix
%in our case.
period = find_periodicity(R,0.05);
%take more periods for better estimation
M =period*floor(N/period);
%if signal is not periodic, or too short to be periodic
if(M < 1)
    M = N;
end

%get autocorrelation matrix
Rx = toeplitz(R(1:M));

%get eigenvalues

%matlab's built-in function
if strcmp(method_eig,'default')
    [eig_vec, eig_vals] = eig(Rx);

%my eigen decomposition function
else
    %not all methods need same number of iterations to converge
    if strcmp(method_eig,'gram_schmidt')
        niter = 50;
    elseif strcmp(method_eig,'hess')
        niter = 30;
    else
        niter = 30;
    end
    [eig_vec,eig_vals] = eig_decomp(Rx,method_eig,niter);
end

[~, inds] = sort(abs(diag(eig_vals)),'descend');

%twice the number of real sinusoids
p = 2*nsignals;
noise_eigvals_pos = inds(p+1:M);
%eigenvectors spanning noise subspace
noise_eigvec = eig_vec(:,noise_eigvals_pos);
noise_subspace = noise_eigvec*noise_eigvec';

%since the signal is real, our search space can be over positive
%frequencies only
omega = linspace(0,pi,nbins/2+1);
omega = omega(1:end-1);
k = 0:M-1;
P = zeros(length(omega),1);

for n = 1:length(omega)
    a = exp(1i*omega(n).*k');
    %pseudospectrum estimation
    P(n) = 1/(a'*noise_subspace*a);
end
%frequency estimates
[peaks,freqs] = find_peaks(abs(P),nsignals);
freqs = (freqs-1)/length(P) * fs/2;

if plot_spec
    plot(omega/pi * (fs/2), abs(P));hold on;grid on;
    plot(freqs, peaks, '*');hold off;grid on;
    xlim([0,0.01]);
    ylabel('Pseudospectrum');
    xlabel('Frequency in Hz');
    title(strcat('MUSIC ', file));
end



