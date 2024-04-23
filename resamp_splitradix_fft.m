function [X] = resamp_splitradix_fft(x,N,D)

%%
% Resample data so that it is periodic with a power of 2
% and perform splitradix fft on it
% Inputs:
% x - data to be resampled
% N - required fft length (power of 2)
% D - actual fft length
% Returns:
% X - Fourier transform of x of size N
%%
[P,Q] = rat(N/D);
x_resamp = resample(x,P,Q);

%splitradix fft (some of this code is taken from Ilias Konsoulas
%"Sprint race for fast butterflies" Matlab package).
X = zeros(1,N);
X = split_radix(x_resamp,N);

    function [X] = split_radix(x,N)
        
        %base cases for N=1 and N=2
        if(N == 1)
            X(1) = x(1);
        elseif(N == 2)
            X(1) = x(1) + x(2);
            X(2) = x(1) - x(2);
        else
            ind1 = 1:2:N-1;
            ind2 = 2:4:N-2;
            ind3 = 4:4:N;

            g = x(ind1); %N/2 length fft
            h = x(ind2); %N/4 length fft
            i = x(ind3); %N/4 length fft
        
            G = split_radix(g,N/2);
            H = split_radix(h, N/4);
            I = split_radix(i,N/4);
            
            %create twiddle factors
            k = 0:N/4-1;
            W = exp(-1i*2*pi*k/N);
            
            %apply formulae
            X(1:N/4) = G(1:N/4) + W.*H + (W.^3).*I;
            X(N/4+1:N/2) = G(N/4+1:N/2) - 1i*W.*H + 1i*(W.^3).*I;
            X(N/2+1:3*N/4) = G(1:N/4) - W.*H - (W.^3).*I;
            X(3*N/4+1:N) = G(N/4+1:N/2) + 1i*W.*H - 1i*(W.^3).*I;
        end
    end

%%to test
% fft([1,1,1,1,0,0,0,0]')
% 
% ans =
% 
%    4.00000 + 0.00000i
%    1.00000 - 2.41421i
%    0.00000 + 0.00000i
%    1.00000 - 0.41421i
%    0.00000 + 0.00000i
%    1.00000 + 0.41421i
%    0.00000 - 0.00000i
%    1.00000 + 2.41421i

end

