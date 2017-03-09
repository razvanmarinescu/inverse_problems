% script to find minimum of rosenbrock by steepest descent
x0 = [3,2];
[v,s,x] = sd2d(@(z) rosenbrock(z),x0);
figure(1);clf;
subplot(2,2,1); plot(v);title('function value');
subplot(2,2,2); plot(s);title('steps');
subplot(2,2,3); plot(x(1,:), x(2,:));
title('iterative solutions');
hold on; plot(x(1,:), x(2,:),'r*');

xlo = -1; dx = 0.05; xhi = 3; 
d = [xlo:dx:xhi];
nd = size(d,2);
im = zeros(nd,nd);
for i = 1:nd
    for j = 1:nd
        im(i,j) = rosenbrock([d(i);d(j)]);
    end
end
subplot(2,2,4);
imagesc((im'));
hold on;
plot((x(1,:)-xlo)/dx+1,(x(2,:)-xlo)/dx+1,'k')
plot((x(1,:)-xlo)/dx+1,(x(2,:)-xlo)/dx+1,'m+')
plot((1-xlo)/dx+1,(1-xlo)/dx+1,'w*');
