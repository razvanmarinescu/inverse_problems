function ip1(eps)
%eps = 0.202;
A = [1 1+eps; 1 - eps 1];
[U,W,V] = svd(A);
disp(['det A : ',num2str(det(A))]);
disp(['condition number : ',num2str(W(1,1)/W(2,2))]);
nsample = 2000;
sol = zeros(2,nsample);
dat = zeros(2,nsample);
dat0 = [10;20];
sol0 = inv(A)*dat0;
for i = 1:nsample
    dat(:,i) = dat0 + randn(2,1);
    sol(:,i) = inv(A) * dat(:,i);
end

solmean = sum(sol')/nsample;
solcvar = sol*sol' /nsample - solmean' * solmean;
[UP,SP] = eig(solcvar);
sd = 1./diag(sqrt(SP));

figure;hold on; plot(dat0(1),dat0(2),'g+');
subplot(1,2,1);plot(dat(1,:),dat(2,:),'r+'); hold on; plot(dat0(1),dat0(2),'g+');axis equal; 
subplot(1,2,2);plot(sol(1,:),sol(2,:),'r+'); hold on; plot(sol0(1),sol0(2),'g+'); axis equal; 
a1 = V(:,1)*[-1,1]*5/W(1,1);
a2 = V(:,2)*[-1,1]*5/W(2,2);
plot(sol0(1)+a1(1,:),sol0(2)+a1(2,:));axis equal; 
plot(sol0(1)+a2(1,:),sol0(2)+a2(2,:));axis equal; 

disp(['singular vectors of A ',num2str(W(1,1)),' ',num2str(W(2,2))]);
disp(['sd of covariance matrices ',num2str(sd(1)),' ',num2str(sd(2))]);