function [period] = find_periodicity(x, thresh, varargin)

%%
% Finds periodicity in autocorrelation function using AMDF
% x - input time domain signal
% thresh - threshold for peak calculation
%%

switch nargin
    case 2
        plot_func = 0;
    case 3
        plot_func = varargin{1};
    otherwise
        error("Wrong number of inputs");
end
    
tau_min = 0;
tau_max = length(x)-1;
period = -1;

D = zeros(1,tau_max - tau_min + 1);
L = length(x);

for n = tau_min:tau_max
    for k = n+1:L
      D(n-tau_min+1) = D(n-tau_min+1) + abs(x(k) - x(k-n));
    end
    D(n-tau_min+1) = D(n-tau_min+1)/L;
end

if plot_func
    figure;
    plot(tau_min:tau_max, D);grid on;
    xlabel('Lags in samples');ylabel('AMDF function');
end

local_minima = zeros(length(D),2);
k = 1;
minPos = -1;
 for i = 2:tau_max-tau_min-1
     if(D(i-1) > D(i) && D(i+1) > D(i) && D(i) <= thresh)
        local_minima(k,1) = i-1;
        local_minima(k,2) = D(i);
        k = k+1;
    end
 end
 
local_minima = local_minima(1:k-1,:);

if isempty(local_minima)
    return;
end

%again find first minimum among all minima in AMDF
%first local minimum considerably smaller than others
if(local_minima(1,2) < mean(local_minima(:,2)) && ...
    local_minima(1,2) > 10)
    minPos = local_minima(1,1);
else
for i = 2:k-2
    if(local_minima(i-1,2) > local_minima(i,2) && local_minima(i+1,2) > local_minima(i,2))
        minPos = local_minima(i,1);
        if minPos < 10
            continue;
        else
            break;
        end
    end
end
end
  
period = round(minPos);


end

