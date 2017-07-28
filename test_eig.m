%script to test eigenvalue decomposition algorithm
close all,clear all, clc;

D = diag(1:10);
rand('seed',36);
S = rand(10); 
S = (S-0.5)*2;
A = S*D*inv(S);
[Q,D_est] = eig_decomp(A,100);
%see error
norm(sort(diag(D)) - sort(diag(D_est)))