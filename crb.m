function [bounds] = crb(p,N,theta,sigma_z)

%%
% Function to calculate Cramer Rao bound for a sum of real sinusoids
% in noise for unbiased estimator - See "Modern Spectral Estimation" 
% by Steven Kay, Section 13.4
% p - number of real sinusoids
% N - number of data points available
% theta - a vector consisting of all amplitudes, frequencies and phases,
%         in order - [f1,a1,phi1,...,fp,ap,phip]
% sigma_z = variance of noise corrupting the signal
% bounds - array containing CRB for each parameter in theta
%% 

M = zeros(3*p,3*p);
bounds = zeros(3*p,1);
n = 0:N-1;

%construct the M matrix
for i = 1:3:3*p
    for j = 1:3:3*p
        val = del(i,j,n,theta);
        M(i,j) = sum((n.^2).*cos(val));
        M(i,j+1) = -sum(n.*sin(val));
        M(i,j+2) = sum(n.*cos(val));
        M(i+1,j) = -M(i,j+1);
        M(i+1,j+1) = sum(cos(val));
        M(i+1,j+2) = sum(sin(val));
        M(i+2,j) = M(i,j+2);
        M(i+2,j+1) = -M(i+1,j+2);
        M(i+2,j+2) = M(i+1,j+1);
    end
end

M_inv = inv(M);
for i = 1:3:3*p
    %frequency bound
    bounds(i) = (sigma_z * M_inv(i,i))/(2*(2*pi*theta(i+1))^2);
    %amplitude bound
    bounds(i+1) = (sigma_z * M_inv(i,i))/2;
    %phase bound
    bounds(i+2) = (sigma_z* M_inv(i,i))/(2*(theta(i+1)^2));
end

    function [val] = del(i,j,n,theta)
        val = 2*pi*(theta(i)-theta(j)).*n + (theta(i+2) - theta(j+2));
    end
        
end