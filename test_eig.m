%script to test eigenvalue decomposition algorithm
close all, clear all, clc;

D = diag(1:10);
rand('seed',36);
S = rand(10); 
S = (S-0.5)*2;
A = S*D/S;

%generate symmetric random matrix
% rand('seed',36);
% a = rand(5);
% A = triu(a) + triu(a,1)';

% %my function
% H1 = hessenberg(A);
% %matlab's built-in function
% H2 = hess(A);
% %error in computation
% err_H = abs(H2 - H1);

[Q_e,D_e] = eig_decomp(A,200,'hess');
[Q_c,D_c] = eig(A);

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
