function [R_hat] = estimate_autocorrelation_function(data, nlags)
%Unbiased estimation of autocorrelation function from data
%nlags - number of lags for which to calculate R_hat

M = nlags;
R_hat = zeros(1,M);
n = length(data);
meu = mean(data);
for k = 0:M-1
    for t = 1:n-k
        R_hat(k+1) = R_hat(k+1) + (data(t)-meu)*(data(t+k)-meu);
    end
    R_hat(k+1) = R_hat(k+1)/(n-k);
end


end

