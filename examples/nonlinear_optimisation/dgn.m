function [v,step,x] = dgn(func,x0,lambda)
% script to find minimum of function by damped Gauss-Newton
% lambda is fixed value to add to Hessian to regularise it
niter = 100;
nx = size(x0,2); % number of variables
x = [];%zeros(2,niter);
x(:,1) = x0; % starting point
v = [];%zeros(niter,1);
step = [];
[v(1),g,H] = func(x(:,1));
p = (H+lambda*eye(nx))\-g;
[v(2),s,xf] = linesearch(@(z) func(z),x(:,1),p);
step(2) = s;
x(:,2) = xf;
tol = 1e-9*abs(v(1));
k = 2;
while k < niter && abs(v(k) - v(k-1)) > tol 
    [v(k),g,H] = func(x(:,k));
    p = (H+lambda*eye(nx))\-g;
    [val,s,xf] = linesearch(@(z) func(z),x(:,k),p);
    v(k+1) = val;
    step(k+1) = s; x = [x,xf]; 
    k = k+1;
end    