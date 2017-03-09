function [dp,e] = DP(g,A,ftrue,sigma,alpha,RH)
% calculate discrepency value
n = size(A,1);
fi = (A'*A + alpha*RH)\(A'*g); % pseudo inverse solution
r = g -A*fi;
dp = r'*r./n - sigma^2;
df = fi - ftrue;
e = df'*df;

