function [flag] = is_symmetric(x)
%checks if a given signal x is symmetric about its midpoint or not
%we assume that the length of x is always even

L = length(x);
flag = 1;
th = 0.01;
for n = 2:L/2-1
    if(abs(x(n) - x(L-n+2)) > th)
        flag = 0;
        break;
    end
end


end

