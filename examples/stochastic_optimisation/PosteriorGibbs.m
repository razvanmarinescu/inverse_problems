% example of Gibbs sampling
% for posterior of Toy problem
eps = 0.202;
A = [1 1+eps; 1 - eps 1];
[U,W,V] = svd(A);
disp(['det A : ',num2str(det(A))]);
disp(['condition number : ',num2str(W(1,1)/W(2,2))]);
sol0 = [1;-0.5];
dat0 = A*sol0;

x1 = [1.5;-0.8]; % "arbitary" starting point
gam = 1.5;      % random walk step length
Niter = 100;    % number of iterations
lambda = 3; %hyperprior;

x = Gibbsfunc(@(y) exp(-norm(dat0 - A*y)^2/2)*MixGaussFunc(y).^(lambda) ,x1,gam,Niter);

% take mean excluding "burn in". E.g. 2nd half only
burnin = floor(Niter/4);
xs = x(:,burnin:Niter-1);
ns = size(xs,2);
mu = sum(xs,2)/ns;
covar = (xs*xs')/ns - mu*mu';

xlo = -4; dx = 0.1; xhi = 4; 
d = [xlo:dx:xhi];
nd = size(d,2);
im = zeros(nd,nd);
for i = 1:nd
    for j = 1:nd
        xx = [d(i);d(j)];
        im(j,i) = exp(-norm(dat0 - A*xx)^2/2)*MixGaussFunc(xx).^(lambda) ;
    end
end
figure(2); clf;
subplot(1,2,1)
imagesc((im));title(['Gibbs Sampling : \gamma = ',num2str(gam)]);
hold on;
plot((x(1,1:Niter-1)-xlo)/dx+1,(x(2,1:Niter-1)-xlo)/dx+1,'w');
plot((x(1,1)-xlo)/dx+1,(x(2,1)-xlo)/dx+1,'m*');
plot((x(1,Niter-1)-xlo)/dx+1,(x(2,Niter-1)-xlo)/dx+1,'k*');
plot((mu(1)-xlo)/dx+1,(mu(2)-xlo)/dx+1,'ro');
[V,D] = eigs(covar);
w = diag(D);

e1 = [mu - 3*w(1)*V(:,1),mu + 3*w(1)*V(:,1)];
e2 = [mu - 3*w(2)*V(:,2),mu + 3*w(2)*V(:,2)];
plot((e1(1,:)-xlo)/dx+1,(e1(2,:)-xlo)/dx+1,'k');axis equal;
plot((e2(1,:)-xlo)/dx+1,(e2(2,:)-xlo)/dx+1,'k');axis equal;
subplot(1,2,2);
imagesc(-log(im));title('negative log posterior');
hold on;

plot((mu(1)-xlo)/dx+1,(mu(2)-xlo)/dx+1,'ro');
plot((e1(1,:)-xlo)/dx+1,(e1(2,:)-xlo)/dx+1,'y');axis equal;
plot((e2(1,:)-xlo)/dx+1,(e2(2,:)-xlo)/dx+1,'y');axis equal;
