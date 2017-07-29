function [freqs] = fast_music(x, nsignals, nbins, method_eig, method_autocorr)
%Replace eigenvalue decomposition in MUSIC with FFT
%x - signal
%nsignals - number of real sinusoids in signal
%nbins - number of points in search space
%method_eig - calculate eigenvalues with dft or fft
%method_autocorr - method for estimating autocorrelation function, direct
%or fft

if(nargin == 3)
    method_eig = 'fft';
elseif nargin == 4
    method_autocorr = 'direct';
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

% %since we know for order M, autocorrelation matrix will be circulant - 
% %no need to do explicit eigenvalue decomposition, just multiply 
% %autocorrelation function DFT matrix
if strcmp(method_eig,'dft')
     dftm = dftmtx(M);
     eigvals = dftm * R';
%direct multiplication with DFT matrix is expensive -- 
%reduce computation time by using FFT
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
%
% omega = linspace(-pi,pi,nbins);
% P = zeros(length(omega),1);
% k = 0:M-1;
% for n = 1:length(omega);
%     a = exp(1i*omega(n).*k');
%     %pseudospectrum estimation
%     P(n) = 1/(a'*noise_subspace*a);
% end


if mod(nbins,2) == 0
    k = -nbins/2+1:nbins/2;
else
    k = -(nbins-1)/2:(nbins-1)/2;
end
P = zeros(nbins,1);


%alternative pseudospectrum estimate from closed-form solution
for m = 1:length(k);
%     for n = 1:length(noise_eigvals_pos)
%         curn = noise_eigvals_pos(n)-1;
%         %this is required to make the function symmetric
%         curk = abs(k(m));
%         %to avoid NaN error, evaluated using L'Hospital's rule
%         if(curk*M == curn*nbins)
%             P(m) = P(m) + M;
%         else
%             P(m) = P(m) + 1/sqrt(M)*(abs(sin(pi*(curk/nbins - curn/M)*M)/...
%             sin(pi*(curk/nbins - curn/M))));
%         end
%     end
%     P(m) = 1/P(m);

    % vectorized code
    curn = noise_eigvals_pos-1;
    curk = abs(k(m));
    inds = find(curn*nbins == curk*M);
    if(~isempty(inds))
        P(m) =  M;
        curn = noise_eigvals_pos([1:inds-1 inds+1:end])-1;
    end
    P(m) = P(m) +  (sum(abs(sin(pi.*(curk/nbins - curn/M)*M)./...
           sin(pi.*(curk/nbins - curn/M)))));
    P(m) = 1./P(m);
    
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

