function [R_hat,M] = estimate_autocorrelation_function(data, nlags, method)

%%
% Unbiased estimation of autocorrelation function
% Inputs:
% data - input time domain signal
% nlags - number of lags for which to calculate R_hat
% method - direct/fft, FFT is much faster
% Outputs:
% R_hat - autocorrelation function as a function of lages
% M - number of lags
%%

if(nargin < 3)
    method = 'direct';
end

M = nlags;
R_hat = zeros(1,M);
n = length(data);

if strcmp(method,'direct')
    meu = mean(data);
    for k = 0:M-1
        for t = 1:n-k
            R_hat(k+1) = R_hat(k+1) + (data(t)-meu)*(conj(data(t+k))-meu);
        end
        R_hat(k+1) = R_hat(k+1)/(n-k);
    end
    
elseif strcmp(method,'fft')
    L = M;
    nfft = 2^nextpow2(2*M-1);
    R_hat = ifft( fft(data,nfft) .* conj(fft(data,nfft)) );
    % rearrange and keep values corresponding to lags: -(len-1):+(len-1)
    R_hat = [R_hat(end-M+1:end) , R_hat(1:M)];
    l = [M-1:-1:0, 0:M-1];
    % remove bias (remove triangular weightings)
    R_hat = R_hat./(L-l);
end



end

