% Do regularised inversion using Richardson-Lucy
%
%
addpath('../');
n = 64;
sigma = 0.1; % std.dev. of added white noise
f = zeros(n,1);
f(12) = 1; f(20:28) = 1.5; f(32:36) = 2; % arbitrary function. Must be non-negative
[y,K] = linblur(f,0.04);
yp = 1e10*imnoise(1e-10*y,'poisson');

figure (1);clf
hold on;
plot(y);
plot(yp,'g');
plot(f,'r');


figure(3);clf;
plot(f,'r');
hold on;
fi = ones(size(f)); % initial guess
ks = K'*ones(size(y));
for k=1:10    % the basic MLEM method
    fi = (fi./ks) .* (K'*(yp./(K*fi))); 
    figure(3);title('reconstruction');
    plot(fi,'k');
    pause(1);
    figure(1);
    plot(K*fi,'k');title('data');
    pause(0.1);
end
frl = fi;
% compare to Least-Square 
flsq= (K'*K + 0.02*eye(n))\(K'*yp);

% and weighted Least Square : 
eps = 1;
ws = diag(1./(sqrt(yp)+eps));
fwlsq= (K'*ws.^2*K + 0.02*eye(n))\(K'*ws.^2*yp);

figure(4);clf;
hold on;
plot(frl,'k');
plot(f,'r');
plot(flsq);
plot(fwlsq,'m');

% regularised RL : iterative EM-MAP
D1 = lindf(f); 
D2 = D1'*D1;
D2(n,n) = 2;
fi = ones(size(f)); % initial guess
tau = 0.05;
k = 1; rr = 2; er = 1;
%while abs(rr - 1) > 1e-3 && k < 20
while er > 1e-4 && k < 100
    disp(['iteration ',num2str(k),' relative change ', num2str(rr),' absolute change ', num2str(er)]);
    finext = (fi./(ks+tau*D2*fi)) .* (K'*(yp./(K*fi)));
    rr = norm(finext./fi,1)/n;
    er = norm(finext - fi)/n;
    flast = fi;
    fi = finext;
    figure(3);
    plot(fi,'g');
    pause(1);
    figure(1);
    plot(K*fi,'--g');
    pause(0.1);
    k = k+1;
end
frlr = fi;
figure(4);clf;
hold on;
plot(frl,'k');
plot(f,'r');
plot(flsq);
plot(fwlsq,'m');
plot(frlr,'g');


