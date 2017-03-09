addpath('../Lin1D');

n = 64;
sigma = 0.1; % std.dev. of added white noise
f = zeros(n,1);
f(12) = -1; f(20:28) = 1.5; f(32:36) = 2; % a
[y,K] = linblur(f,0.04); 
yn = y + sigma*randn(n,1); % add white noisrbitrary function

figure (1);clf; hold on;
plot(y);
plot(yn,'g');
plot(f,'r');
%%
%xcgn = CGNE1(zeros(size(f)),K,yn,1e-9,floor(n/10) );
alpha = 1e-2;
tol = 1e-3;
niter = floor(n*4/10);
[xcg,xicg,ricg,picg,itcg] = CG1(zeros(size(f)),K'*K+alpha*eye(64),K'*yn,tol,niter);
disp(['CG iterations ',num2str(itcg)]);
[xsd,xisd,risd,itsd] = SD1(zeros(size(f)),K'*K+alpha*eye(64),K'*yn,tol,niter);
disp(['SD iterations ',num2str(itsd)]);


semilogy(sum(ricg.^2,1),'--k');
figure(3); clf; hold on;
title(['reconstruction with \alpha=',num2str(alpha)]);
plot(f,'r');
%plot(xcgn,'k');
plot(xcg,'--k');
plot(xsd,'b');
%% look at residual of normal equations : || A'(y - Ax_*)||
figure(4); clf;
subplot(2,2,1);
semilogy(sum(risd(:,1:itsd).^2,1),'b');
hold on;
semilogy(sum(ricg(:,1:itsd).^2,1),'--k');
title('residual norm ||A^T(y-Ax_*)||');

%% also look at data space norm
yicg = yn*ones(1,itcg)- K*xicg(:,1:itcg);
yisd = yn*ones(1,itsd)-K*xisd(:,1:itsd);

figure(4); 
subplot(2,2,2);
semilogy(sum(yisd.^2,1),'b');
hold on;
semilogy(sum(yicg.^2,1),'--k');
title('residual norm ||(y-Ax_*)||');

%%

figure(5); clf; hold on;
plot(f,'r');
for k = 1:itsd
    plot(xisd(:,k),'b');
    pause(0.5);
end

for k = 1:itcg
    plot(xicg(:,k),'--k');
    pause(0.5);
end

%% ------------ use Matlab versions! -------------
tic;
[xpcg,FLAGpcg,RELRESpcg,ITERpcg,RESVECpcg] = pcg(K'*K+alpha*eye(64),K'*yn,tol,n);
toc;
disp(['PCG iterations ',num2str(ITERpcg)]);

 figure(3);
 plot(xpcg,'g');
 
 figure(4);
 subplot(2,2,3);
 semilogy(RESVECpcg,'g');
 hold on;
 title('residual from Matlab functions');
  %% --- lastly : LSQR

tic;
[xlsqr,FLAGlsqr,RELRESlsqr,ITERlsqr,RESVEClsqr] = lsqr([K; sqrt(alpha)*eye(n)],[yn; zeros(n,1)],tol,n);
toc;
disp(['LSQR iterations ',num2str(ITERlsqr)]);
figure(3);
plot(xlsqr,'m');
 
figure(4);
subplot(2,2,3);
semilogy(RESVEClsqr,'m');
 
