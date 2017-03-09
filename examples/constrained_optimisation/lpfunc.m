function [f,g,H] = lpfunc(x,p)
%
f = sum( (x.^2).^(p/2));
g = p*(x.^2).^((p-1)/2);
H = p*(p-1)*diag( (x.^2).^((p-2)/2));
