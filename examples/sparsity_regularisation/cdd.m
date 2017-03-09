% coordinate descent algorithm, naive implementation

%%
n = 500;
m = n/2;
A = randn(m,n);
xe= zeros(n,1);
ind = rand(n,1);
ind = find(ind >.9);
xe(ind) = 10*randn(size(ind));
A  = A/norm(A);
be = A*xe;
b  = be + max(abs(be))*1e-5*randn(m,1);
al = 1e-2;

%% ill-posed problem: first-order derivative
clear, close all
n = 100;
m = n;
N = 100;
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

% iterative soft thresholding, Gauss-Seidel version
x = zeros(n,1);
r = b-A*x;
res = [];
d = zeros(n,1);
for i = 1:n
    d(i) = norm(A(:,i))^2;
end
for k = 1:N
    for i = 1:n
        t = x(i) + A(:,i)'*r/d(i);
        t = sign(t)*max(abs(t)-al/d(i),0);
        r = r + A(:,i)*(x(i)-t);
        x(i) = t;
    end
    figure(1), plot(1:n,x,'r.',1:n,xe,'k.'), axis square, set(gca,'fontsize',20)
    res(k) = norm(A*x-b);
end

%%
figure(2), plot(1:N,res,'k','linewidth',2)
axis square
set(gca,'fontsize',20)

% %%
% % iterative soft thresholding, Gauss-Seidel version
% x = zeros(n,1);
% r = b-A*x;
% res = [];
% d = zeros(n,1);
% for i = 1:n
%     d(i) = norm(A(:,i))^2;
% end
% for k = 1:300
%     r = b-A*x;
%     for i = 1:n
%         t = x(i) + A(:,i)'*r/d(i);
%         t = sign(t)*max(abs(t)-al/d(i),0);
%         % r = r + A(:,i)*(x(i)-t);
%         x(i) = t;
%     end
%     figure(1), plot(1:n,x,'r.',1:n,xe,'k.'), axis square, set(gca,'fontsize',20)
%     res(k) = norm(A*x-b);
% end
% 
% figure(2), plot(1:300,res,'k','linewidth',2)
% axis square
% set(gca,'fontsize',20)