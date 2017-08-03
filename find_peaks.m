function [peaks,freqs] = find_peaks(Xwdb, maxPeaks)

%finds peaks in signal/plot

allPeaks = [];
indPos = [];
k = 1;

for i=2:length(Xwdb)-1
    if(Xwdb(i) >= Xwdb(i-1) && Xwdb(i) >= Xwdb(i+1))
        allPeaks(k) = Xwdb(i);
        indPos(k) = i;
        k = k+1;
    end
end

%-- Order from largest to smallest magnitude, keep only maxPeaks of them --%
%sort along first column only
[peaks, pos] = sort(allPeaks,'descend'); 
if(length(peaks) < maxPeaks)
    peaks = zeros(1,maxPeaks);
    freqs = zeros(1,maxPeaks);
    return;
end
peaks = peaks(1:maxPeaks);
inds = indPos(pos);
inds = inds(1:maxPeaks);
freqs = zeros(size(peaks));

%-- Do parabolic interpolation in dB magnitude to find more accurate peak --%
%-- and frequency estimates --%

for i=1:maxPeaks
    %idx=find(Xwdb==peaks(i));
    idx = inds(i);
    %parabolic interpolation
    a=Xwdb(idx-1);
    b=Xwdb(idx);
    c=Xwdb(idx+1);
    p = 0.5*((a-c)/(a+c-2*b));
    peaks(i) = b - (0.25*(a-c)*p);
    freqs(i) = (idx + p); %in bins
end


end

