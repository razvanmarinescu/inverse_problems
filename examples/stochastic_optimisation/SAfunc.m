% Simulated Annealing : Metropolis Hastings sampling
function [x] = SAfunc(Pfunc,x1,gam,Niter)

x = zeros(2,Niter+1);
k = 1;
x(:,k) = x1;
while k < Niter
    Px = Pfunc(x(:,k));
    w = gam*randn(size(x(:,k)));
    y = x(:,k) + w;
    Py = Pfunc(y);
    alpha = min(1,Py/Px);
    u = rand; % draw from uniform;
    if u < alpha
        x(:,k+1) = y;
        k = k+1;
    else
        x(:,k+1) = x(:,k);
    end
    if(mod(k,20) == 0)
        gam = 0.98*gam;
    end
end

