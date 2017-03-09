% script to find minimum of rosenbrock by Levenberg-Marquardt
x0 = [3;2];
lambda = 0; % fixed regularisation term
[v,s,x] = lm(@(z) rosenbrockFN(z),x0);
figure(4);clf;
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
        im(i,j) = rosenbrockFN([d(i);d(j)]);
    end
end
subplot(2,2,4);
imagesc((im'));title('Levenberg-Marquardt');
hold on;
plot((x(1,:)-xlo)/dx+1,(x(2,:)-xlo)/dx+1,'k')
plot((x(1,:)-xlo)/dx+1,(x(2,:)-xlo)/dx+1,'m+')
plot((1-xlo)/dx+1,(1-xlo)/dx+1,'w*');
