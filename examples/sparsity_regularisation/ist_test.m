clear, close all

%% simple example: random matrix
n = 500;
m = round(n/2);
N = 300;
A = randn(m,n);
xe= zeros(n,1);
ind = rand(n,1);
ind = find(ind >.9);
xe(ind) = 10*randn(size(ind));
A  = A/norm(A);
be = A*xe;
b  = be + max(abs(be))*1e-5*randn(m,1);
al = 5e-2;

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

% iterative soft thresholding 
x = zeros(n,1);
res = [];
fcnl_ist= [];
for k = 1:N
    x = (x-A'*(A*x-b));
    x = sign(x).*max(abs(x)-al,0);
    figure(1), plot(1:n,x,'r.',1:n,xe,'k.'), axis square, set(gca,'fontsize',20)
    res(k) = norm(A*x-b);
    fcnl_ist(k)= norm(A*x-b)^2/2+al*norm(x,1);
end
x_ist = x;
figure(2), plot(1:N,res,'k',1:N,fcnl_ist,'r','linewidth',2)
axis square
set(gca,'fontsize',20)
x_ls  = (A'*A+al*eye(n))\(A'*b);
figure(1), hold on, plot(1:n,x_ls,'b','linewidth',2)

%% 
% iterative soft thresholding with Cauchy's rule
x = zeros(n,1);
res = [];
fcnl_cauchy= [];
figure(1)
for k = 1:N
    rr = A'*(A*x-b);
    tau = norm(rr)^2/norm(A*rr)^2
    x = (x-tau*A'*(A*x-b));
    x = sign(x).*max(abs(x)-tau*al,0);
    figure(3), plot(1:n,x,'r.',1:n,xe,'k.'), axis square, set(gca,'fontsize',20), title('Cauchy')
    res(k) = norm(A*x-b);
    fcnl_cauchy(k)= norm(A*x-b)^2/2 + al*norm(x,1);
end
figure(4), plot(1:N,res,'k',1:N,fcnl_cauchy,'r','linewidth',2)
axis square
set(gca,'fontsize',20)


%%
% iterative soft thresholding with Barzilai-Borwein rule
x0 = zeros(n,1);
g0 = A'*(A*x0-b);
res = [];
fcnl_bb = [];
figure(1), 
for k = 1:N
    if k == 1
        tau = 1;
        x = x0-tau*g0;
        x = sign(x).*max(abs(x)-tau*al,0);
        x1= x;
        g1= A'*(A*x1-b);
    else
        dx = x1-x0;
        dg = g1-g0;
        % tau = norm(dx)^2/dg'*dx
        tau = dg'*dx/norm(dg)^2;
        tau = max(min(tau,1000),0);  % restrict to the region
        % update the previous step
        x0 = x1;
        g0 = g1;
        x = x1-tau*g1;
        x = sign(x).*max(abs(x)-tau*al,0);
        x1= x;
        g1= A'*(A*x1-b);
    end
    figure(5), plot(1:n,x1,'r.',1:n,xe,'k.'), axis square, set(gca,'fontsize',20)
    title('Barzilai-Borwein')
    res(k) = norm(A*x1-b);
    fcnl_bb(k)= norm(A*x1-b)^2/2+al*norm(x1,1);
end
figure(6), plot(1:N,res,'k',1:N,fcnl_bb,'r','linewidth',2)
axis square
set(gca,'fontsize',20)

%% FIST 
x0 = zeros(n,1);
y  = x0;
t1 = 1;
fcnl_fista = [];
for k = 1:N
    % update x
    x1 = y - A'*(A*y-b);
    x1 = sign(x1).*max(abs(x1)-al,0);
    % update t
    t2 = (1+sqrt(1+4*t1^2))/2;
    % update y
    y = x1 + (t1-1)/t2*(x1-x0);
    x0 = x1;
    t1 = t2;
    figure(7), plot(1:n,x1,'r.',1:n,xe,'k.'), axis square, set(gca,'fontsize',20)
    title('FISTA')
    fcnl_fista(k) = norm(A*x1-b)^2/2+al*norm(x1,1);
end
x = zeros(n,1);

figure(8), semilogy(1:N,abs(fcnl_fista-fcnl_fista(end)),'r',1:N,abs(fcnl_ist-fcnl_fista(end)),'k')

%% figure (9)
figure(9), semilogy(1:N,fcnl_ist,'r',1:N,fcnl_cauchy,'g',1:N,fcnl_bb,'b',1:N,fcnl_fista,'k','linewidth',2)