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

%need to have 4 and 6 as factors as well
pos2 = find(ifax == 2);
pos3 = find(ifax == 3);
minl = min(length(pos2), length(pos3));

%find all factors of 6
for n = 1:minl
    ifax(pos2(n)) = ifax(pos2(n))*ifax(pos3(n));
    ifax(pos3(n)) = 0;
end

%find all factors of 4
pos2 = find(ifax == 2);
for n = 1:2:length(pos2)-1
    ifax(pos2(n)) = ifax(pos2(n))*ifax(pos2(n+1));
    ifax(pos2(n+1)) = 0;
end

%final array of factors in ascending order including 4 and 6
ifax = sort(ifax(ifax~=0));    
nfax = length(ifax);
trigs = exp(-2*1i*pi.*(0:N-1)/N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        for k = 0:LA:(M-LA)
            for l = 1:LA
                
               %unsimplified matrix calculation
%                 if(k == 0)
%                     omega = eye(ifac);
%                 else
%                     omega = diag(trigs((0:(ifac-1))*k+1));
%                 end                 
%                 C(ja+j) = omega * (dftmtx(ifac)*A(ia+i));
%                
%               %simplified matrix calculation - fastest
                for m = 1:ifac
                    C(ja(m)+j) = trigs(k*(m-1)+1) * sum(A(ia + i) .* ...
                        exp(-2*pi*1i*(0:ifac-1)*(m-1)/ifac).'); 
%                   C(ja(m)+j) = exp(-2*pi*1i*k*(m-1)/N) * sum(A(ia + i) .* ...
%                         exp(-2*pi*1i*(0:ifac-1)*(m-1)/ifac).'); 
                end

%                 %vectorized version
%                 C(ja + j) =  trigs((0:ifac-1)*k+1).' .* ...
%                     (dftmtx(ifac)*A(ia+i));
                
                i = i+1;
                j = j+1;
            end
        j = j+jump;
        end
    end

X = A;
end


