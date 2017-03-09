% Compare zero and first-order Tikhonov
%
%
n = 64;
sigma = 0.1; % std.dev. of added white noise
f = zeros(n,1);
f(12) = -1; f(20:28) = 1.5; f(32:36) = 2; % arbitrary function
[y,K] = linblur(f,0.04);
yn = y + sigma*randn(n,1); % add white noise, 
RH0 = eye(n);
D1=lindf(f);
D2 = D1'*D1; % 1D Laplacian;
RH1 = D2;
df = D1*f;
T = 0.5*max(abs(df));
kd = exp(- (df/T).^2); % This is an edge indicator function;
RHA1 = D1'*diag(kd)*D1;

alphaTK0 = 0.0166;  % Got from regln choice using UPRE.
alphaTK1 = 0.06348; % Got from regln choice using UPRE.
alphaATK1 = 0.1; % Got from regln choice using UPRE.

figure (1);
hold on;
plot(y);
plot(yn,'g');
plot(f,'r');

figure(3);
hold on;
plot(f,'r');
plot(df,'--r');
plot(kd,'g');
f0 = (K'*K + alphaTK0*RH0)\(K'*yn);
f1 = (K'*K + alphaTK1*RH1)\(K'*yn);
fa1 = (K'*K + alphaATK1*RHA1)\(K'*yn);
figure(3);
plot(f0,'k');
plot(f1,'m');
plot(fa1,'c');
figure(1);
plot(K*f0,'k');
plot(K*f1,'m');
plot(K*fa1,'c');