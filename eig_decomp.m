function [V,D] = eig_decomp(X, niter, method)

%function that does eigenvalue decomposition with the QR algorithm
%V - eigenvectors, D = diagonal matrix with eigenvalues
%X - input matrix (must be square)
%niter = number of iterations
%method - algorithm to be used :
%gram_schmidt - gram schmidt QR O(n^3)
%hessenberg - hessenberg QR O(n^2)
%tridiagonal - hessenberg for symmetric matrices

%notes - in general, U returns the Schur vectors, not the eigenvectors
%but for normal matrices (that includes hermitian matrices), the schur
%vectors are the eigenvectors

if(nargin == 2)
    method = 'gram_schmidt';
end

A = X;
n = length(A);

if (strcmp(method,'gram_schmidt')) 
    U = eye(n);
    for i = 1:niter
        [Q,R] = gram_schmidt(A);
        A = R*Q;
        U = U*Q;
    end
    D = diag(diag(A));
    V = U;
    
elseif(strcmp(method,'hess'))
        [H,U] = hessenberg(A);
        c = zeros(n-1,1);
        s = zeros(n-1,1);
        for i = 1:niter
            for k = 1:n-1
                %H = G(k,k+1,v)*H
                [c(k),s(k)] = givens(H(k,k),H(k+1,k));
                rot_mat = [c(k), -s(k);s(k), c(k)];
                H(k:k+1,k:n) = rot_mat * H(k:k+1,k:n);
            end
            for k = 1:n-1
                %apply givens rotation from right
                rot_mat = [c(k), s(k);-s(k), c(k)];
                H(1:k+1,k:k+1) = H(1:k+1,k:k+1)*rot_mat;
                U(:,k:k+1) = U(:,k:k+1)*rot_mat;
            end
        end
        V = U;
        D = diag(diag(H));
end
    
end
    

