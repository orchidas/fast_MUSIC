function [flag] = is_circulant(X)
%check if given square matrix is circulant
N = length(X);
x = X(1,:);
flag = 1;
th = 10^-2;
for k = 1:N
    for j = 1:N
        if(abs(X(k,j) -x(mod(j-k,N)+1)) > th)
            flag = 0;
            break;
        end
    end
    if (flag == 0)
        break;
    end
end


end

