function [R] = estimate_covariance_matrix(data, order)
%function that estimates covariance matrix from given data
%data - input data
%order - desired order of covariance matrix

N = length(data);
R = zeros(order,order);
ncols = floor(N/order);
data = data(1:ncols*order);
X = reshape(data,[ncols, order])';
%meu = mean(X,2);

for k = 1:ncols
    R = R +(X(:,k)*X(:,k)');
end
R = R/ncols;


end

