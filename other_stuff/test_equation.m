%test closed form equation for pseudospectrum
close all, clc;

M = 40;
%n = 0:M-5;
N = 50;
k = -N/2+1:N/2;
y = zeros(N,1);
n = 5;

for m = 1:length(k)
    curk = k(m);
%    inds = find(n*N == curk*M);
%     if(~isempty(inds))
%         y(m) =  M;
%         curn = n([1:inds-1 inds+1:end]);
%     else
%         curn = n;
%     end
%     y(m) = (sum(abs(sin(pi.*(curk/N - curn/M)*M)./...
%            sin(pi.*(curk/N - curn/M)))));
    if(n*N == curk*M)
        y(m) = M;
    else
        y(m) = abs(sin(pi.*(curk/N - n/M)*M)./...
           sin(pi.*(curk/N - n/M)));
    end

end

figure;
plot(k,y);grid on;