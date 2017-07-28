function [Q,R] = gram_schmidt(A)
%QR factorization with gram schmidt orthogonalization
%Q - matrix with orthonormal vectors along the column
%R = upper triangular matrix
N = length(A);
Q = zeros(N,N);
R = zeros(N,N);

%initial values
u = A(:,1);
Q(:,1) = u./norm(u);
R(1,:) = Q(:,1).' * A;

for k = 2:N
    u = A(:,k);
    for j = 1:k-1
        u = u - (u.'*Q(:,j))*Q(:,j);
    end
    Q(:,k) = u./norm(u);
    R(k,k:N) = Q(:,k).'*A(:,k:N);
end

        


end

