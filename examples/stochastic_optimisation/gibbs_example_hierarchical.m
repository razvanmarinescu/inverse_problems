% gibbs sampling for hierarchical model
% Gaussian + smoothness + lam 

clear, close all

%%  ill-posed problem: first-order derivative
n = 100;
m = n;
% A = randn(n,m);
% A = A/norm(A);
A = ones(n);
A = tril(A);
A = A/norm(A);
tt= linspace(0,1,n)';
xe= tt.*(tt<.5) + (1-tt).*(tt>=0.5);
be = A*xe;
sig = 5e-2*max(abs(be));
b  = be + sig*randn(m,1);
figure(10), plot(1:n,xe,'r')

%% smoothness prior
L = spdiags([-1*ones(n,1) ones(n,1)], 0:1,n-1,n);
W = L'*L;

%% augmented Tikhonov 
maxit=20;
lam = 1e-5;
al0 = 1; bt0 = 1e-2;
format short e
for i=1:maxit
    x=(lam*sig^2*W+A'*A)\(A'*b);
    lam = (m/2+al0-1)/(x'*W*x/2+bt0);
    err = norm(x-xe)/norm(xe);
    figure(2),plot(1:m,xe,'k',1:m,x,'r-'),title('joint MAP','fontsize',16)
    pause(0.1)
end


%% Gibbs sampling with hierarchical model
N = 1e4;
Ensx = zeros(m,N);
Enslam = zeros(1,N);
x   = zeros(m,1);      % initial guess for x
lam = 1e-1;
al0 = 1; bt0 = 1e-4;     % parameter for Gamma distribution
Enslam(1) = lam;
tau = 1/sig^2;
for i = 2:N
    % update x componentwise, given x
    for j = 1:m
        a0 = tau*sum(A(:,j).^2) + lam*W(j,j);
        mus = (b-A*x + x(j)*A(:,j));
        mup = sum(W(:,j).*x)+W(j,:)*x-2*x(j)*W(j,j);
        b0 = 2*tau*sum(A(:,j).*mus)-lam*mup;
        mu = b0/(2*a0);
        sigx = sqrt(1/a0);
        x(j) = randn(1)*sigx+mu;
    end
    Ensx(:,i) = x;
    lam = gamrnd(m/2+al0,1/(1/2*x'*W*x+bt0))
    Enslam(i) = lam;
    if rem(i,1000) == 1, i, end
end

%% post processing
x_mean = mean(Ensx(:,1000:N),2);
figure(2), plot(1:n,xe,1:n,x_mean,'linewidth',2), axis square
legend('exact','mean')
x_var  = var(Ensx(:,1000:N),0,2);
figure(3), errorbar(1:n,x_mean,3*sqrt(x_var),'linewidth',2)