function [peaks,peak_pos] = find_peaks(Xwdb, maxPeaks, varargin)

%%
% Finds peaks in signal.
% Inputs:
% Xwdb - input signal (can be in dB)
% maxPeaks - maximum number of peaks we are looking for
% do_parabolic_interp (optional) - wheter to do parabolic interpolation (y/n)
% Returns:
% peaks - array of peaks (maxPeaks x 1)
% peak_pos - x values where peaks are found (in samples)
%%
if nargin == 2
    do_parabolic_interp = 'n';
else
    do_parabolic_interp = varargin{1};
end
allPeaks = [];
indPos = [];
k = 1;

%find all local maxima
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
    peak_pos = zeros(1,maxPeaks);
    return;
end
peaks = peaks(1:maxPeaks);
inds = indPos(pos);
inds = inds(1:maxPeaks);
peak_pos = inds;

%-- Do parabolic interpolation in dB magnitude to find more accurate peak --%
%-- and frequency estimates --%

if strcmp(do_parabolic_interp, 'y')
    for i=1:maxPeaks
        %idx=find(Xwdb==peaks(i));
        idx = inds(i);
        %parabolic interpolation
        a=Xwdb(idx-1);
        b=Xwdb(idx);
        c=Xwdb(idx+1);
        p = 0.5*((a-c)/(a+c-2*b));
        peaks(i) = b - (0.25*(a-c)*p);
        peak_pos(i) = (idx + p); %in bins
    end
end



