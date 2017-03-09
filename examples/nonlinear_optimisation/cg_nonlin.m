function [v,step,x] = cg_nonlin(func,x0)
% script to find minimum of function by nonlinear CG
niter = 100;
x = [];%zeros(2,niter);
x(:,1) = x0; % starting point
v = [];%zeros(niter,1);
step = [];
[v(1),g] = func(x(:,1));
[v(2),s,xf] = linesearch(@(z) func(z),x(:,1),-g);
step(2) = s;
x(:,2) = xf;
tol = 1e-9*abs(v(1));
p = -g;
ngkm = norm(g)^2;
k = 2;
while k < niter && abs(v(k) - v(k-1)) > tol 
    [v(k),g] = func(x(:,k));
    beta = norm(g)^2/ngkm;
    p = -g + beta*p;
    [val,s,xf] = linesearch(@(z) func(z),x(:,k),p);
    v(k+1) = val;
    step(k+1) = s; x = [x,xf]; ngkm = ngkm*beta;
    k = k+1;
end    