function [freqs] = fast_music(x, nsignals, nbins, method)
%Replace eigenvalue decomposition in MUSIC with FFT
%x - signal
%nsignals - number of real sinusoids in signal
%nbins - number of points in search space
%method - dft or fft

if(nargin < 4)
    method = 'fft';
end

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

% %since we know for order M, autocorrelation matrix will be circulant - 
% %no need to do explicit eigenvalue decomposition, just multiply 
% %autocorrelation function DFT matrix
if(strcmp(method,'dft'))
     dftm = dftmtx(M);
     eigvals = dftm * R';
%direct multiplication with DFT matrix is expensive -- 
%reduce computation time by using FFT
%M = 2^nextpow2(M);
%use fft to reduce computation
else
    eigvals = fft(R(1:M), M);
end

[eig_vals_sorted, inds] = sort(abs(eigvals),'descend');

% figure;
% stem(1:M, eig_vals_sorted);hold off;
% title('Sorted eigenvalues');

p = 2*nsignals;
noise_eigvals_pos = inds(p+1:M);

%eigenvectors spanning noise subspace
% if(strcmp(method, 'dft'))
%     noise_eigvec = 1/sqrt(M) .* dftm(:,noise_eigvals_pos);
% else
%     noise_eigvec = exp(2*pi*1i*(0:M-1)'*(noise_eigvals_pos-1)/M);
% end
% noise_subspace = noise_eigvec*noise_eigvec';

if mod(nbins,2) == 0
    k = -nbins/2+1:nbins/2;
else
    k = -(nbins-1)/2:(nbins-1)/2;
end
P = zeros(nbins,1);

% omega = linspace(-pi,pi,nbins);
% P = zeros(length(omega),1);
% k = 0:M-1;
% for n = 1:length(omega);
%     a = exp(1i*omega(n).*k');
%     %pseudospectrum estimation
%     %P(n) = 1/(a'*noise_subspace*a);
%     P(n) = 1/sum(a'*noise_subspace);
% end

%alternative pseudospectrum estimate from closed-form solution
for m = 1:length(k);
%     for n = 1:length(noise_eigvals_pos)
%         curn = noise_eigvals_pos(n)-1;
%         %to avoid NaN error, evaluated using L'Hospital's rule
%         if(k(m)*M == curn*nbins)
%             P(m) = P(m) + M;
%         else
%             P(m) = P(m) + 1/sqrt(M)*(abs(sin(pi*(k(m)/nbins - curn/M)*M)/...
%             sin(pi*(k(m)/nbins - curn/M))));
%         end
%     end
%     P(m) = 1/P(m);

    % vectorized code
    P(m) =  1/(sum(abs(sin(pi.*(k(m)/nbins - (noise_eigvals_pos-1)/M)*M)./...
           sin(pi.*(k(m)/nbins - (noise_eigvals_pos-1)/M)))));
    if isnan(P(m))
         P(m) = 1/M;
    end
end

%frequency estimates
[peaks,freqs] = find_peaks(P,p);

figure;
plot(2*k/nbins, P);hold on;grid on;
plot(freqs/pi, peaks, '*');hold off;grid on;
ylabel('Pseudospectrum');
xlabel('Frequency in radians normalized by pi');
title('Fast MUSIC');

end

