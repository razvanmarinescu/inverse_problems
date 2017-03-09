% example of Metropolis Hastings sampling
% Example is from Kaipio and Somersalo 2005, chapter 3.6
x1 = [1.5;-0.8]; % "arbitary" starting point
gam = 0.5;      % randome walk step length
Niter = 100;    % number of iterations

x = Gibbsfunc(@(y) exp(-10*(y(1)^2 - y(2))^2 - (y(2) - 0.25)^4),x1,gam,Niter);

% take mean excluding "burn in". E.g. 2nd half only
burnin = floor(Niter/2);
xs = x(:,burnin:Niter-1);
ns = size(xs,2);
mu = sum(xs,2)/ns;
covar = (xs*xs')/ns - mu*mu';

xlo = -2; dx = 0.05; xhi = 2; 
d = [xlo:dx:xhi];
nd = size(d,2);
im = zeros(nd,nd);
for i = 1:nd
    for j = 1:nd
        xx = [d(i);d(j)];
        im(j,i) = exp(-10*(xx(1)^2 - xx(2))^2 - (xx(2) - 0.25)^4) ;
    end
end
figure(2); clf;
imagesc((im));title(['Gibbs Sampling : \gamma = ',num2str(gam)]);
hold on;
plot((x(1,1:Niter-1)-xlo)/dx+1,(x(2,1:Niter-1)-xlo)/dx+1,'w');
plot((x(1,1)-xlo)/dx+1,(x(2,1)-xlo)/dx+1,'m*');
plot((x(1,Niter-1)-xlo)/dx+1,(x(2,Niter-1)-xlo)/dx+1,'k*');
plot((mu(1)-xlo)/dx+1,(mu(2)-xlo)/dx+1,'ro');
[V,D] = eigs(covar);
w = diag(D);

e1 = [mu - 5*w(1)*V(:,1),mu + 5*w(1)*V(:,1)];
e2 = [mu - 5*w(2)*V(:,2),mu + 5*w(2)*V(:,2)];
plot((e1(1,:)-xlo)/dx+1,(e1(2,:)-xlo)/dx+1,'y');axis equal;
plot((e2(1,:)-xlo)/dx+1,(e2(2,:)-xlo)/dx+1,'y');axis equal;