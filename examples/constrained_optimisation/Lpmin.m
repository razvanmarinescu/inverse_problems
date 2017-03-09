% minimise LP norm
clear all;
A = [1;2]; b = 5; % constraint
x0 = (pinv(A)*b)';   % pseudo inverse is a feasible point;
%pp = 2; x0 = [0;2.5];
%pp = 1.001; x0 = [5/3;5/3];
%pp = 0.9999;
%pp = 0.0005;
pp= 0.05;
%pp = 20;
[f,g,H] = lpfunc(x0,pp);
x(:,1) = x0;
v(1) = f;
%tau = 0.1;
tau = 1;
KeepGoing = true;
k = 1;
tol = 1e-6;
lambda = 10;
while KeepGoing && k < 30
    % use LM method to make matrix spd
    AA = [H+lambda*eye(2)  A;A'  0];
    z = AA\[-g;0];
    h = z(1:2);
    mu = z(3);
    y = x(:,k) + tau*h;
    [f,g,H] = lpfunc(y,pp);
    if(f <= v(k))
%        tau = max([1,tau*2]);
        lambda = lambda/2;
        v(k+1) = f;
        x(:,k+1) = y;
        if( (v(k)-v(k+1))/v(k) < tol)
            KeepGoing = false;
        end
        k = k+1;
    else
        %tau = tau/2;
        lambda = lambda*2;
    end    
end
    
%--------------------------
figure(4);clf;
subplot(2,2,1); plot(v);title('function value');
%subplot(2,2,2); plot(s);title('steps');
subplot(2,2,3); plot(x(1,:), x(2,:));
title('iterative solutions');
hold on; plot(x(1,:), x(2,:),'r*');

xlo = 0; dx = 0.05; xhi = 3; 
d = [xlo:dx:xhi];
nd = size(d,2);
im = zeros(nd,nd);
for i = 1:nd
    for j = 1:nd
        xx = [d(i);d(j)];
        im(j,i) = lpfunc(xx,pp) ;
    end
end
subplot(2,2,4);
imagesc((im));title('Lp norm Constraint Optimisation');
hold on;
for j = 1:size(A,2) % plot the constraints
   xj = pinv(A(:,j))*b(j);
   p(:,1) = xj' -2*[A(2,j);-A(1,j)];
   p(:,2) = xj' + 2*[A(2,j);-A(1,j)];
   plot((p(1,:)-xlo)/dx+1,(p(2,:)-xlo)/dx+1,'w');
end
% plot the solutions
plot((x(1,:)-xlo)/dx+1,(x(2,:)-xlo)/dx+1,'k')
plot((x(1,:)-xlo)/dx+1,(x(2,:)-xlo)/dx+1,'m*')
%plot((xsol(1)-xlo)/dx+1,(xsol(2)-xlo)/dx+1,'w*');
