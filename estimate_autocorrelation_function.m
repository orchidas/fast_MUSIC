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
elseif strcmp(method,'fft')
    L = 2*M-1;
    R_hat = fftshift(ifft(fft(data,L).*fft(fliplr(conj(data)),L)));
    %filter out bias - remove triangular weighting
    l = -(L-1)/2:L/2;
    R_hat = R_hat./(L-abs(l));
    %since ACF is symmetric, preserve positive lags only
    R_hat = R_hat((L+1)/2:end);
    M = (L+1)/2;   
    
end

% figure;
% plot(0:M-1,R_hat);title('Autocorrelation function');
% xlabel('Lags');


end

