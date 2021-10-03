# FAST MUSIC
<h3>A more efficient MUSIC algorithm for sinusoidal parameter estimation.</h3>

<p>This repository contains code to evaluate FAST MUSIC, an efficient extension of the MUSIC algorithm for parameter estimation of approximately periodic signals.
It is useful for resolving close-frequency modes of a system. Eigenvalue decomposition in MUSIC is replaced by the FFT by making the
autocorrelation matrix circulant. For more details see

<i> <b> <a href = "http://dafx2018.web.ua.pt/papers/DAFx2018_paper_43.pdf">FAST MUSIC - An efficient implementation of the MUSIC algorithm for frequency estimation of approximately periodic signals </a></b> </i>- Orchisama Das, Jonathan S. Abel, Julius O. Smith in Proc. of International Conference on Digital Audio Effects, DAFx 2018.</b>

Additionally, it contains scripts to do eigenvalue decomposition (QR factorization, reduction to Hessenberg form), FFTs (Self sorting mixed radix and resampled split radix) for any composite order.
FAST MUSIC is evaluated against traditional MUSIC and QIFFT.
