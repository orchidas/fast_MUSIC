function [num_sine] = determine_number_of_sinusoids(eigvals_sorted, max_signals)

%determines number of sinusoids in signal from a list of sorted eigenvalues
%eigvals_sorted - eigenvalues of autocorrelation matrix sorted in
%descending order
%max_signals - maximum number of possible sinusoids

%we make use of the fact that in case of signal of r harmonics + AWGN, 
%the last M-r eigenvalues will be same (denoting the noise variance factor)

eig_odd = eigvals_sorted(1:2:end);
pos = find(abs(diff(eig_odd)) < 0.05);
d = diff(pos);
count = zeros(length(d),2);
k = 1;
i = 1;
while(i < length(d))
    if(d(i) ~= 1)
        i = i+1;
        continue;
    else
        count(k,1) = i;
        while(d(i) == 1 && i<length(d))
            count(k,2) = count(k,2)+1;
            i = i+1;
        end
        k = k+1;
    end
end

ind = find(count(:,2) > max_signals,1);
num_sine = pos(count(ind,1))-1;
    
end

