% Gibbs sampling
function [x] = Gibbsfunc(Pfunc,x1,gam,Niter)
xrange = [-2:0.1:2]*gam;
xlo = min(xrange);
xhi = max(xrange);
nx = length(xrange);
xx = (xhi-xlo)/(nx-1);
n = size(x1,1);
x = zeros(n,Niter+1);
k = 1;
x(:,k) = x1;
while k < Niter
    j = mod(k,n) + 1;  % which dimenstion to search
    ej = zeros(n,1); ej(j) = 1; % unit vector
    xs = x(:,k)*ones(1,41)  + ej*xrange;
    Px = zeros(1,size(xs,2));
    for kx = 1:41
    Px(kx) = Pfunc(xs(:,kx));   % 1D sample of P
    end
    Cx = zeros(size(Px)); % cumulative PDF
    np = size(Px,2);
    for p = 2:np
        Cx(p) = Cx(p-1) + Px(p);
    end
    Cx = Cx/Cx(end);      % normalise
    figure(3); 
    subplot(1,2,1);plot(Px);
    subplot(1,2,2);plot(Cx);
    pause(0.5);
    u = rand; % draw from uniform;
    lo = 1; mid = np/2; hi = np;
    for s = 1:ceil(log2(np))   % binary search with maximum depth
        if(u < Cx(round(mid)))
            hi = mid;
        else
            lo = mid;
        end
        mid = 0.5*(lo + hi);
    end
    dx = mid - round(mid);
    x(:,k+1) = x(:,k) + ej*(xlo + mid*xx);  % only approximate
    k = k+1;
end

