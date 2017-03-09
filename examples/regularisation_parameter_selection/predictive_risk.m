function [p,e] = predrisk(g,A,ftrue,sigma,alpha,RH)
% calculate predictive value
n = size(A,1);
fi = (A'*A + alpha*RH)\(A'*g); % pseudo inverse solution
y = A*fi - A*ftrue;
p = y'*y./n;
df = fi - ftrue;
e = df'*df;