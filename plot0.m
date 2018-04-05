close all;
path = '../piano data/A3/';
files = {'A3.wav    ','A3_1D.wav ','A3_2D.wav ','A3_3D.wav ','A3_23D.wav','A3_13D.wav','A3_12D.wav'};
fs = 44100;
nfiles = length(files);
[b,a] = butter(4,[2500 2900]/(fs/2));
lengths = [0.5,1,2,4,6];
nl = length(lengths);
nbins = 2^19;
tfs = zeros(2*nbins,nl);
f = (0:nbins)/nbins*fs/2;
%read signal
[x,fs] = audioread(strcat(path,files{1}));  

for n = 1:nl

nsamp = round(lengths(n)*fs);
window = blackman(nsamp);
tfs(:,n) = fft(x(fs+(1:nsamp)).*window, 2*nbins);

end

figure(1); plot(f, 20*log10(abs(tfs(1:nbins+1,:))/max(abs(tfs(:))))+20); grid;
ylim([-60 0]);
xlim([2660,2680]);
%xlim([2400,2900]);