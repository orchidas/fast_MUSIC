function [H,U] = hessenberg(A,sym)

%%
% Function to compute Hessnberg matrix from A using Householder reflections.
% A - input matrix
% sym - is A symmetric?
% H - hessenberg matrix
% U - eigenvectors
%% 

if(nargin == 1)
    sym = 'nsym';
end
n = length(A);
U = eye(n);
u = zeros(n,n);

%if a matrix is symmetric, Hessenberg reduction gives a tridiagonal 
%matrix that can be computed with fewer flops
if strcmp(sym,'sym')
    
    v = zeros(n,n);
    for k = 1:n-2
        x = A(k+1:n,k);
        alpha = -sign(x(1))*norm(x);
        %householder reflector, P = I - 2*u*u'
        tmp = x - alpha.*[1; zeros(n-k-1,1)];
        u(k+1:n,k) = tmp./norm(tmp);
        %P*AP can be computed in single step by making use of v
        %v(:,k) = 2*A*u(:,k) - 2*u(:,k)*(u(:,k)'*A*u(:,k));
        %A = A - u(:,k)*v(:,k)'- v(:,k)*u(:,k)';
        v(:,k) = 2*A(:,k+1:n)*u(k+1:n,k) - 2*u(:,k)*...
        (u(k+1:n,k)'*A(k+1:n,:)*u(:,k));
        A = A-[zeros(k,n);u(k+1:n,k)*v(:,k)'] - ...
            [zeros(n,k),v(:,k)*u(k+1:n,k)'];
    end  
    %to get eigenvectors
    for k = n-2:-1:1
        %U = P*U
        U(k+1:n,k+1:n) = U(k+1:n,k+1:n) - 2*u(k+1:n,k)* ...
        (u(k+1:n,k)'*U(k+1:n,k+1:n));
    end
 
%non-symmetric matrix
else
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
    
end
 H = A;
 H(abs(H) < 10^-8) = 0;
        
end

