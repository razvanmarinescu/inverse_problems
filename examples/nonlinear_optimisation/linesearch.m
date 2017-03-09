function [val,step,xf] = linesearch(func,x0,s)
% line search func from x0 in direction s
v0 = func(x0);
tol = 1e-3*abs(v0); % tolerance
lo = 0;
step = 1;
mid = step;
v1 = func(x0 + mid*s);
nup = 1; maxup = 10; % max uphill steps before giving up
while (v1 - v0) > tol && nup < maxup
    step = step/10;      % initial step too large. So shorten it
    mid = step;
    v1 = func(x0 + mid*s);
    nup = nup+1;
end
if(nup == maxup)
    disp('initial line search heads uphill');
    val = v0; step = 0; xf = x0;
else
    hi = 2*step;
    v2 = func(x0 + hi*s);
    while( v2 < v1) % still on decrease
        step = step*2;
        lo = mid;v0 = v1;
        mid = hi;v1 = v2;
        hi = hi+step;v2 = func(x0 + hi*s);        
    end
    % now lo, mid, hi includes the minimum
    % fit a parabola
    A = [lo^2 lo 1; mid^2 mid 1; hi^2 hi 1];
    q = A\[v0;v1;v2];
    step = -0.5*q(2)/q(1);
    xf = x0 + step*s;
    val = func(xf);
end