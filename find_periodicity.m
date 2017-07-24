function [period] = find_periodicity(x, thresh)
%finds periodicity in autocorrelation function using AMDF

tau_min = 0;
tau_max = length(x)-1;

D = zeros(1,tau_max - tau_min + 1);
L = length(x);

for n = tau_min:tau_max
    for k = n+1:L
      D(n-tau_min+1) = D(n-tau_min+1) + abs(x(k) - x(k-n));
    end
    D(n-tau_min+1) = D(n-tau_min+1)/L;
end

% figure;
% plot(tau_min:tau_max, D);grid on;
% xlabel('Lags in samples');ylabel('AMDF function');

minPos = -1;
 for i = 2:tau_max-tau_min-1
     if(D(i-1) > D(i) && D(i+1) > D(i) && D(i) <= thresh)
        minPos = i-1;
        break;
    end
 end
                 
period = minPos + tau_min;


end

