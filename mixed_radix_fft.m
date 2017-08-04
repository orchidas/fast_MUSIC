function [X] = mixed_radix_fft(x, N)

%Function to compute mixed radix FFT for any composite length N
%according to algorithm given in 'Self-sorting mixed radix fast Fourier
%transforms' - C Temperton, published in Journal of Computational Physics
%X - DFT of x
%x - Input vector of length(N), can be real or complex

C = zeros(N,1);
A = x;

%get prime factors of N
ifax = factor(N);
nfax = length(ifax);
trigs = exp(-2*1i*pi.*(0:N-1)/N);

%decimation in frequency
LA = 1;
for n = 1:nfax
    ifac = ifax(n);
    %reverse roles of A and C
    [C,A] = pass(A,C,trigs,ifac,LA,N);
    LA = LA*ifac;
end

    function [A,C]  = pass(A,C,trigs,ifac,LA,N)
        %self-sorting decimation in frequency
        M = N/ifac;
        ia = (0:ifac-1)*M;
        ja = (0:ifac-1)*LA;
       
        i = 1;j = 1;
        jump = (ifac-1)*LA;
        for k = 0:LA:M-LA
            for l = 1:LA
                if(k == 0)
                    omega = eye(ifac);
                else
                    omega = diag(trigs((0:(ifac-1))*k+1));
                end
                %to do - simplify matrix calculation
                C(ja+j) = omega * (dftmtx(ifac)*A(ia+i));
                i = i+1;
                j = j+1;
            end
        j = j+jump;
        end
    end

X = A;
end


