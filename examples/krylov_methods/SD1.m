function [x,xi,ri,i] = SD1(x,A,y,tol,n)
% simple SD for symmetric matrix A
% solves y =  A x

nA = length(x);
r=A'*y;
rsold=r'*r;
r0 = rsold;
xi= zeros(nA,n);
ri = zeros(nA,n);
for i=1:n
    Ar = A*r;
    tau = (r'*r)/(Ar'*Ar);
    x = x + tau*r;
    xi(:,i) = x;
    ri(:,i) = r;
    r=r-tau*A'*Ar;
    rsnew=r'*r;
   if sqrt(rsnew/r0)<tol
        break;
    end
end