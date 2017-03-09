function [gcv,e,nf] = GCV(g,A,ftrue,sigma,alpha,RH)
% calculate Generalised Cross Validation value
n = size(A,1);
Adag = inv(A'*A + alpha*RH)*A';
Pdag = A*Adag;
fi = Adag*g; % pseudo inverse solution
r = g -A*fi;
gcv = (r'*r./n) / (trace(eye(n)-Pdag)/n)^2;
df = fi - ftrue;
e = df'*df;
nf = fi'*fi;