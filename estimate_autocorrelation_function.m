function [R_hat,M] = estimate_autocorrelation_function(data, nlags, method)
%Unbiased estimation of autocorrelation function from data
%nlags - number of lags for which to calculate R_hat
%method - direct/fft

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
else
    L = 2*M-1;
    L = 2^nextpow2(L);
    R_hat = ifft(fft(data,L).*fft(fliplr(conj(data)),L));
    R_hat = R_hat./M;
    %since ACF is symmetric, preserve positive lags only
    R_hat = R_hat(L/2+1:end);
    M = L/2;   
    %TO DO - need to get rid of tapering in ACF
    
end

% figure;
% plot(0:M-1,R_hat);title('Autocorrelation function');
% xlabel('Lags');


end

