function [freqs] = fast_music(x, nsignals, nbins)
%Replace eigenvalue decomposition in MUSIC with FFT

N = length(x);
%estimate autocorrelation function
R = estimate_autocorrelation_function(x, N);

% figure;
% plot(0:N-1,R);title('Autocorrelation function');
% xlabel('Lags');

%M is the number of antenna, or the dimension of the autocorrelation matrix
%in our case.
M = find_periodicity(R,0.05);
%if signal is not periodic, or too short to be periodic
if(M < 1)
    M = N;
end

% %since we know for order M, autocorrelation matrix will be circulant - 
% %no need to do explicit eigenvalue decomposition, just multiply 
% %autocorrelation function DFT matrix
% dftm = dftmtx(M);
% eigvals = dftm * R';

%direct multiplication with DFT matrix is more expensive than doing
%eigenvalue decomposition if M is large - reduce computation time by using
%FFT
%M = 2^nextpow2(M);
%use fft to reduce computation
eigvals = fft(R(1:M), M);
[eig_vals_sorted, inds] = sort(abs(eigvals),'descend');

% figure;
% stem(1:M, eig_vals_sorted);hold off;
% title('Sorted eigenvalues');

p = 2*nsignals;
noise_eigvals_pos = inds(p+1:M);
%eigenvectors spanning noise subspace
%noise_eigvec = 1/sqrt(M) .* dftm(:,noise_eigvals_pos);
noise_eigvec = exp(2*pi*1i*(0:M-1)'*(noise_eigvals_pos-1)/M);

omega = linspace(-pi,pi,nbins);
omega = omega(1:end-1);
k = 0:M-1;
P = zeros(length(omega),1);

for n = 1:length(omega);
    a = exp(1i*omega(n).*k');
    %pseudospectrum estimation
    P(n) = 1/(a'*(noise_eigvec*noise_eigvec')*a);
end
%frequency estimates
[peaks,freqs] = find_peaks(abs(P),p);

% figure;
% plot(omega/pi, abs(P));hold on;grid on;
% plot(freqs/pi, peaks, '*');hold off;grid on;
% ylabel('Pseudospectrum');
% xlabel('Frequency in radians normalized by pi');
% title('Fast MUSIC');

end

