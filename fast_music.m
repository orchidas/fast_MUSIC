function [freqs] = fast_music(x, nsignals, nbins, method_eig, method_autocorr,M)

%Replace eigenvalue decomposition in MUSIC with FFT
%x - signal
%nsignals - number of real sinusoids in signal
%nbins - number of points in search space
%method_eig - calculate eigenvalues with dft or fft
%method_autocorr - method for estimating autocorrelation function, direct
%or fft
%M - autocorrelation matrix order (ideally should be calcuated from ACF
%periodicity, but included just for plotting accuracy vs M).

if nargin == 3
    method_eig = 'fft';
elseif nargin == 4
    method_autocorr = 'fft';
end

N = length(x);
%estimate autocorrelation function
R = estimate_autocorrelation_function(x, N/2, method_autocorr);

%M is the number of antenna, or the dimension of the autocorrelation matrix
%in our case.
if nargin == 5
    M = find_periodicity(R,0.05);
    %if signal is not periodic, or too short to be periodic
    if(M < 1)
        M = N/2;
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

else
    %My FFT
    %eigvals = mixed_radix_fft(R', M);
    %Matlab's built-in FFT (much faster)
    eigvals = fft(R',M);
end

[eig_vals_sorted, inds] = sort(abs(eigvals),'descend');

% figure;
% stem(1:M, eig_vals_sorted);hold off;
% title('Sorted eigenvalues');

p = 2*nsignals;
noise_eigvals_pos = inds(p+1:M);

k = 0:nbins/2-1;
%k = -nbins/2+1:nbins/2;
P = zeros(nbins/2,1);

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
%

    %vectorized code
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

     curn = noise_eigvals_pos-1;
     curk = k(m);
     P(m) = M./(sum((abs(sin(pi.*(curk/nbins - curn/M)*M)./...
           sin(pi.*(curk/nbins - curn/M)))).^2));
    if isnan(P(m))
        P(m) = 1/(M*(M-p));
    end  
    
end

%frequency estimates
[peaks,freqs] = find_peaks(P,nsignals);
freqs = (freqs-1)*(pi/length(P));

% figure;
% plot(k/nbins, P);hold on;grid on;
% plot(freqs/(2*pi), peaks, '*');hold off;grid on;
% ylabel('Pseudospectrum');
% xlabel('Frequency in Hz');
% title('Fast MUSIC');

%since the signal is real, the spectrum will be symmetric
freqs = [-freqs,freqs];
pos_nan = find(isnan(freqs));
freqs(pos_nan) = 0;

end

