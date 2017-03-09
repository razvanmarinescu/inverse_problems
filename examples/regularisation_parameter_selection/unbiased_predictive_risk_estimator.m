function [upre,e] = UPRE(g,A,ftrue,sigma,alpha,RH)
% calculate Unbiased Predictive Risk Estimate 
n = size(A,1);
Adag = inv(A'*A + alpha*RH)*A';
Pdag = A*Adag;
fi = Adag*g; % pseudo inverse solution
r = g -A*fi;
upre = r'*r./n - sigma^2 + 2*sigma^2 *trace(Pdag)/n;
df = fi - ftrue;
e = df'*df;