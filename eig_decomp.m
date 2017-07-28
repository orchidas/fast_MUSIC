function [V,D] = eig_decomp(X, niter, method)

%function that does eigenvalue decomposition with the QR algorithm
%V - eigenvectors, D = diagonal matrix with eigenvalues
%X - input matrix (must be square)
%niter - number of iterations of QR algorithm
%method - algorithm to be used :
%gram_schmidt - gram schmidt QR O(n^3)
%hessenberg - hessenberg QR O(n^2)

if(nargin == 2)
    method = 'gram_schmidt';
end

A = X;
N = length(A);
if (strcmp(method,'gram_schmidt')) 
    U = eye(N);
    for i = 1:niter
        [Q,R] = gram_schmidt(A);
        A = R*Q;
        U = U*Q;
    end
    D = diag(diag(A));
    V = U;
end
    
    
    


end

