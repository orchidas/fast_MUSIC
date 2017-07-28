function [freqs] = music(x, nsignals, nbins)
%MUSIC algorithm for sinusoid parameter estimation
%x - signal corrupted with white noise
%nsignals - number of real sinusoids in signal
%nbins - number of bins in search space

N = length(x);
%estimate autocorrelation function
R = estimate_autocorrelation_function(x, N);

%M is the number of antenna, or the dimension of the autocorrelation matrix
%in our case.
M = find_periodicity(R,0.05);
%if signal is not periodic, or too short to be periodic
if(M < 1)
    M = N;
end
%get autocorrelation matrix
%method 1
%[corr, Rx] = corrmtx(x, M-1);
%method 2
Rx = toeplitz(R(1:M));

%get eigenvalues
[eig_vec, eig_vals] = eig(Rx);
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

figure;
plot(omega/pi, abs(P));hold on;grid on;
plot(freqs/pi, peaks, '*');hold off;grid on;
ylabel('Pseudospectrum');
xlabel('Frequency in radians normalized by pi');
title('MUSIC');


end



