function [v,step,x] = lm(func,x0)
% script to find minimum of function by Levenberg-Marquardt
% lambda is large to begin with, gets smaller with iteration
niter = 100;
lambda = 1000;
nx = size(x0,1); % number of variables
x = [];%zeros(2,niter);
x(:,1) = x0; % starting point
v = [];%zeros(niter,1);
step = [];
[v(1),g,H] = func(x(:,1));
% p = (H+lambda*eye(nx))\-g;
% xf = x0 + p;
% val = func(xf);
% v(2) = val;
% step(1) = lambda;
% x(:,2) = xf;
step(1) = 0;
tol = 1e-9*abs(v(1));
derr = abs(v(1));
k = 2;
while k < niter && derr > tol
    [val,g,H] = func(x(:,k-1));
    val = v(k-1) + 1;
    while(val - v(k-1) > tol)
        p = (H+lambda*eye(nx))\-g;
        xf = x(:,k-1) + p;
        val = func(xf);
        if(val - v(k-1) < tol)
            v(k) = val;
            step(k) = lambda;
            lambda = lambda/10;
            derr = v(k-1) - val;
        else
            lambda = lambda*10;
        end
    end
    x = [x,xf];
    k = k+1;
end