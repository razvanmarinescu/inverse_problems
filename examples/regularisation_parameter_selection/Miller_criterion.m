function [mil,e] = Miller(g,A,ftrue,sigma,alpha,RH)
% calculate discrepency value
n = size(A,1);
fi = (A'*A + alpha*RH)\(A'*g); % pseudo inverse solution
r = g -A*fi;
dp = r'*r./n / sigma^2;
%calculate Prior;
pp = 0.5*fi'*RH*fi/n;
mil = 0.5*dp - pp;
df = fi - ftrue;
e = df'*df;

