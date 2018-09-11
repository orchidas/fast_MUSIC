# FAST MUSIC
<h3>A more efficient MUSIC algorithm for sinusoidal parameter estimation.</h3>

<p>This repository contains code to evaluate FAST MUSIC, an efficient extension of the MUSIC algorithm for approximately periodic signals.
It is particularly useful for resolving sinusoids close in frequency. Eigenvalue decomposition in MUSIC is replaced by the FFT by making the
autocorrelation matrix circulant. For more details see :

<i> <b> <a href = "http://dafx2018.web.ua.pt/papers/DAFx2018_paper_43.pdf">FAST MUSIC -- An efficient implementation of the MUSIC algorithm for frequency estimation of approximately periodic signals </a></b>- O. Das, J. Abel, J.O Smith </i> in Proc of 21st International Conference on Digital Audio Effects (DAFx 18), Aveiro, Portugal.

Additionally, it contains scripts to do eigenvalue decomposition (QR factorization, reduction to Hessenberg form), FFTs (Self sorting mixed radix and resampled split radix) for any composite order.
FAST MUSIC is evaluated against traditional MUSIC and QIFFT.
