function [freqs] = music(x, nsignals, nbins, method_eig, method_autocorr)

%MUSIC algorithm for sinusoid parameter estimation
%x - signal corrupted with white noise
%nsignals - number of real sinusoids in signal
%nbins - number of bins in search space
%method eig - algorithm for eigenvalue decomposition
%method_autocorr - method for calculating autocorrelation function, direct
%or fft

if nargin == 3
    method_eig = 'default';
elseif nargin == 4
    method_autocorr = 'fft';
end

N = length(x);
%estimate autocorrelation function
R = estimate_autocorrelation_function(x, N/2, method_autocorr);

%M is the number of antenna, or the dimension of the autocorrelation matrix
%in our case.
M = find_periodicity(R,0.05);
%if signal is not periodic, or too short to be periodic
if(M < 1)
    M = N;
end
%get autocorrelation matrix
%method 2
Rx = toeplitz(R(1:M));

%get eigenvalues
%matlab's built-in function
if strcmp(method_eig,'default')
    [eig_vec, eig_vals] = eig(Rx);
%my eigen decomposition function
else
    %not all methods need same no of iterations to converge
    if strcmp(method_eig,'gram_schmidt')
        niter = 100;
    elseif strcmp(method_eig,'hess')
        niter = 100;
    else
        niter = 100;
    end
    [eig_vec,eig_vals] = eig_decomp(Rx,method_eig,niter);
end

[eig_vals_sorted, inds] = sort(abs(diag(eig_vals)),'descend');

% figure;
% stem(1:M, eig_vals_sorted);hold off;
% title('Sorted eigenvalues');

%nsignals = determine_number_of_sinusoids(eig_vals_sorted, max_signals);
%twice the number of real sinusoids
p = 2*nsignals;
noise_eigvals_pos = inds(p+1:M);
%eigenvectors spanning noise subspace
noise_eigvec = eig_vec(:,noise_eigvals_pos);
noise_subspace = noise_eigvec*noise_eigvec';

omega = linspace(-pi,pi,nbins);
k = 0:M-1;
P = zeros(length(omega),1);

for n = 1:length(omega);
    a = exp(1i*omega(n).*k');
    %pseudospectrum estimation
    P(n) = 1/(a'*noise_subspace*a);
end
%frequency estimates
[peaks,freqs] = find_peaks(abs(P),p);
freqs = -pi + (freqs-1)*(2*pi/length(P));

% figure;
% plot(omega/pi, abs(P));hold on;grid on;
% plot(freqs/pi, peaks, '*');hold off;grid on;
% ylabel('Pseudospectrum');
% xlabel('Frequency in radians normalized by pi');
% title('MUSIC');


end



