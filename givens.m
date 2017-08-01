function [c,s] = givens(x,y)

%givens rotation
c = x/sqrt(abs(x)^2 + abs(y)^2);
s = -y/sqrt(abs(x)^2 + abs(y)^2);

end

