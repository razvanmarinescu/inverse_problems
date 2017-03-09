function [v,step,x] = sd2d(func,x0)
% script to find minimum of 2D function by steepest descent
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
k = 2;
while k < niter && abs(v(k) - v(k-1)) > tol 
    [v(k),g] = func(x(:,k));
    [val,s,xf] = linesearch(@(z) func(z),x(:,k),-g);
    v(k+1) = val;
    step(k+1) = s; x = [x,xf];%(:,k) = xf;
    k = k+1;
end    