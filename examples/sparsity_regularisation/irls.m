% iteratively reweighted least-squares method
n = 1000;
m = round(n/2);
A = randn(m,n);
xe= zeros(n,1);
ind = rand(n,1);
ind = find(ind >.9);
xe(ind) = 10*randn(size(ind));
A  = A/norm(A);
be = A*xe;
b  = be + max(abs(be))*1e-5*randn(m,1);
al = 1e-2;

%% 
%% ill-posed problem: first-order derivative
clear, close all
n = 100;
m = n;
N = 500;
A = ones(n);
A = tril(A);
xe= zeros(n,1);
ind = rand(n,1);
ind = find(ind >.9);
xe(ind) = 10*randn(size(ind));
A  = A/norm(A);
be = A*xe;
b  = be + max(abs(be))*1e-5*randn(m,1);
al = 1e-4;

%%
% iteratively reweighed least-squares
x = zeros(n,1);
res = [];
fcnl= [];
for k = 1:200
    W = diag(2./(abs(x)+1e-10));
    x = (A'*A+al*W)\(A'*b);
    figure(1), plot(1:n,x,'r.',1:n,xe,'k.'), axis square, set(gca,'fontsize',20)
    res(k) = norm(A*x-b);
    fcnl(k)= norm(A*x-b)^2/2+al*norm(x,1);
end

%%
figure(2), plot(1:200,res,'k',1:200,fcnl,'r','linewidth',2), 
legend('residual','functional value')
axis square
set(gca,'fontsize',20)

% fast iterative soft thresholding