function [H,U] = hessenberg(A)

%Function to compute Hessnberg matrix from A using Householder reflections.
%H - hessenberg matrix, U - eigenvectors

n = length(A);
U = eye(n);
u = zeros(n,n);

%if a matrix is symmetric, Hessenberg reduction gives a tridiagonal 
%matrix that can be computed with fewer flops

%if isequal(A,A')
 
 
%non-symmetric matrix
%else
    for k = 1:n-2
        %column vector
        x = A(k+1:n,k);
        alpha = -sign(x(1))*norm(x);
        %householder reflector, P = I - 2*u*u'
        tmp = x - alpha.*[1; zeros(n-k-1,1)];
        u(1:n-k,k) = tmp./norm(tmp);
        %apply P from left
        A(k+1:n,k:n) = A(k+1:n,k:n) - 2*u(1:n-k,k)* ...
        (u(1:n-k,k)'*A(k+1:n,k:n));
        %apply P from right
        A(1:n,k+1:n) = A(1:n,k+1:n) - 2*(A(1:n,k+1:n)* ...
        u(1:n-k,k))*u(1:n-k,k)';
    end
    %to get eigenvectors
    for k = n-2:-1:1
        %U = P*U
        U(k+1:n,k+1:n) = U(k+1:n,k+1:n) - 2*u(1:n-k,k)* ...
        (u(1:n-k,k)'*U(k+1:n,k+1:n));
    end
    H = A;
    %inds = find(abs(H) < 10^-10);
    H(abs(H) < 10^-10) = 0;
    
%end
        


end

