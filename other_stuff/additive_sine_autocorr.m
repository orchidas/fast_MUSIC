%Author - Orchisama Das

%this script generates additive sinusoids artificially and compares their
%theoretical and estimated autocorrelation function and matrices.
%Additionally, the circulant structure of the matrices under certain
%conditions and the corresponding eigenvalues are explored.

close all, clear all, clc;
n = 0:2000;
%clean signal
%example 1
%y = cos(0.4*pi.*n + 0.1*pi) + 0.5*cos(0.5*pi.*n+0.3*pi) + 0.2*cos(0.6*pi.*n);
%example 2
y = cos(0.5*pi.*n) + cos(0.52*pi.*n + 0.25*pi);
%nornalize signal power to 0dB
y_norm = y./max(abs(y));
snr = 10;
%signal+noise
x = awgn(y_norm, snr);
var_w = var(x - y_norm);

figure;
plot(n(1:200),x(1:200));
xlabel('Samples');title('Signal+noise');

%theoretical autocorrelation
%example 1
%R = 0.5*cos(n.*0.4*pi) + 0.125*cos(0.5*pi.*n) + 0.02*cos(0.6*pi.*n);
%example 2
R = 0.5*cos(0.5*pi.*n) + 0.5*cos(0.52*pi.*n);
%normalize autocorrelation function
R = R./(max(abs(y))^2);
%compensating for additive white noise
R(1) = R(1) + var_w;
%find periodicity of autocorrelation function
M = find_periodicity(R,0.001);
R = R(1:M);
%estimate autocorrelation function from data
R_hat = estimate_autocorrelation_function(x,M);
%this will determine if the resulting estimated autocorrelation matrix will
%be circulant
symm = is_symmetric(R_hat)

figure;
plot(1:M,R);hold on;
plot(1:M,R_hat);hold off;grid on;
xlabel('Lags');legend('Theoretical ACF','Estimated ACF');
title('Autocorrelation function');

%estimating power spectrum
win_R = hamming(M)'.*R;
win_Rhat = hamming(M)'.*R_hat;
Nfft = 2^nextpow2(M);
psd_thr = fftshift(fft(win_R,Nfft));
psd_est = fftshift(fft(win_Rhat,Nfft));
freq = linspace(-pi,pi,Nfft);
psd_thr = 20*log10(abs(psd_thr)/max(abs(psd_thr)));
psd_est = 20*log10(abs(psd_est)/max(abs(psd_est)));
[peaks_thr, freqs_thr] = find_peaks(psd_thr,4);
[peaks_est, freqs_est] = find_peaks(psd_est,4);

%plot power spectrums
figure;
plot(freq/pi, psd_thr);hold on;grid on;
plot(freq/pi, psd_est);hold on;grid on;
plot(freqs_thr/pi, peaks_thr,'*b');hold on;
plot(freqs_est/pi, peaks_est,'*r');hold off;
xlabel('Frequency in pi*rad/s');ylabel('Mag in dB');
legend('Theoretical PSD','Estimated PSD');
title('Power spectrum');

%%Comments - we can see that the autocorrelation function is periodic. Let
%%the period be denoted as p. If we choose M = np, where n is any positive
%%integer, the resulting autocorrelation matrix will be circulant. 
%%The minimum order which will yield a circulant matrix is M=p %%

%get theoretical autocorrelation matrix
auto_mat = toeplitz(R);
%try with matlab's built-in function - this is not circulant
[corr_mat,auto_mat_est] = corrmtx(x, M-1);
%try converting estimated autocorrelation function to toeplitz - this is
%not exactly circulant but is very close
%auto_mat_est = toeplitz(R_hat);
% get eigenvalues of theoretical autocorrelation matrix
[eigvec_R, eigval_R] = eig(auto_mat);
%get DFT matrix
dftm = dftmtx(M);
%get eigenvalues of theoretical autocorrelation matrix via dft
eigval_R_dft = dftm * R';
% get eigenvalues of estimated autocorrelation matrix via dft
eigval_R_est_dft =  dftm * R_hat';

%find a corresponding circulant matrix
%toeplitz(col, row)
circ = toeplitz(circshift(flipud(R(:)),1),R);
[eigvec_circ, eigval_circ] = eig(circ);

%check if created circulant matrix is indeed circulant for sanity check
flag0 = is_circulant(circ)
%check if theoretical autocorrelation matrix is circulant
flag1 = is_circulant(auto_mat)
%check if estimated autocorrelation matrix is circulant
flag2 = is_circulant(auto_mat_est)

%plot eigenvalues
figure;
stem(1:M, sort(abs(diag(eigval_R)), 'descend'));hold on;
stem(1:M, sort(abs(diag(eigval_circ)), 'descend'));hold on;
stem(1:M, sort(abs(eigval_R_dft), 'descend'));hold on;
stem(1:M, sort(abs(eigval_R_est_dft), 'descend'));hold off;
legend('eigenvalues of theoretical autocorrelation matrix',...
    'eigenvalues of circulant matrix', 'eigenvalues of theoretical matrix calcuated with DFT', 'eigenvalues of estimated matrix calculated with DFT');

%mean squared error between theoretical eigenvalues
mse_theor = mean(abs((sort(abs(diag(eigval_R)) - sort(abs(eigval_R_dft))))).^2)
%mean squared error between estimated eigenvalue
mse_est = mean(abs((sort(abs(diag(eigval_R)) - sort(abs(eigval_R_est_dft))))).^2)

%find eigenvectors corresponding to top eigenvalues
[eig_sort, inds] = sort(abs(eigval_R_est_dft), 'descend');
pos = determine_number_of_sinusoids(eig_sort, 5);
top_eigvals = eig_sort(1:2*pos);
top_eigvals_pos = inds(1:2*pos);
%eigenvectors are corresponding columns in DFT matrix
eigvec_top = dftm(:, top_eigvals_pos);
%frequency resolution
f_res = 2*pi/M;
%frequency present in signal according to positions of top eigenvalues
f_detected = -pi + (top_eigvals_pos-1)*f_res;
%frequency detected normalized by pi radians
f_detected/pi 

