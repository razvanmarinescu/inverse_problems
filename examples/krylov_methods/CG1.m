function [x,xi,ri,pi,i] = CG1(x,A,y,tol,n)
% simple CG for symmetric matrix A
% solves y =  A x
m = length(x);
r=A'*y;
p=r;
rsold=r'*r;
r0 = rsold;
xi= zeros(m,n);
pi = zeros(m,n);
ri = zeros(m,n);
for i=1:min(n,size(A,1))
    Ap=A'*A*p;
    alpha=rsold/(p'*Ap);
    x=x+alpha*p;
    xi(:,i) = x;
    ri(:,i) = r;
    pi(:,i) = p;
    r=r-alpha*Ap;
    rsnew=r'*r;
    if sqrt(rsnew/r0)< tol
        break;
    end
    p=r+rsnew/rsold*p;
    rsold=rsnew;
end