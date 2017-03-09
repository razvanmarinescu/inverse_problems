% Do regularised inversion using Truncated SVD
%
n = 64;
sigma = 0.1; % std.dev. of added white noise
f = zeros(n,1);
f(12) = -1; f(20:28) = 1.5; f(32:36) = 2; % arbitrary function
[y,K] = linblur(f,0.04);
yn = y + sigma*randn(n,1); % add white noise, 

figure (1);
hold on;
plot(y);
plot(yn,'g');
plot(f,'r');

[U,W,V] = svd(K);
disp(['det K : ',num2str(det(K))]);
disp(['condition number : ',num2str(W(1,1)/W(n,n))]);
figure(2); 
subplot(2,3,1);plot(diag(W));
subplot(2,3,2);plot(U(:,1));
subplot(2,3,3);plot(U(:,2));
subplot(2,3,4);plot(U(:,4));
subplot(2,3,5);plot(U(:,8));
subplot(2,3,6);plot(U(:,16));

figure(3);
plot(f,'r');
hold on;
fi = zeros(size(f));
for k=1:floor(2*n/3)
    fi = fi + U(:,k)'*yn *V(:,k)/W(k,k);
    figure(3);
    plot(fi,'k');
    pause(1);
    figure(1);
    plot(K*fi,'k');
    pause(3);
end
figure(4);
fp = K\yn;
plot(fp,'m');
