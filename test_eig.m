%script to test eigenvalue decomposition algorithm
close all, clc;

% D = diag(1:10);
% rand('seed',36);
% S = rand(10); 
% S = (S-0.5)*2;
% A = S*D/S;

%generate symmetric random matrix
rand('seed',36);
a = rand(5);
A = triu(a) + triu(a,1)';


%reduction to tridiagonal form
% [H1,U1] = hessenberg(A,'sym');
% %reduction to hessenberg form
% [H2,U2] = hessenberg(A);
% %error in computation
% err_H = abs(H2 - H1);
% err_U = abs(U2 - U1);

[Q_e,D_e] = eig_decomp(A,'tridiag',50);
[Q_c,D_c] = eig_decomp(A,'hess',200);

%arrange eigenvalues and eigenvectors in order
[vals_c,inds_c] = sort(diag(D_c), 'descend');
[vals_e,inds_e] = sort(diag(D_e), 'descend');
D_cor = diag(vals_c);
D_est = diag(vals_e);
Q_cor = zeros(size(A));
Q_est = zeros(size(A));

for k = 1:length(A)
    Q_cor(:,k) = Q_c(:,inds_c(k));
    Q_est(:,k) = Q_e(:,inds_e(k));
end

%see error
norm(sort(diag(D_cor)) - sort(diag(D_est)))
err_Q = abs(Q_est-Q_cor);
