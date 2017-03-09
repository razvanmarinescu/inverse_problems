% script to find minimum of Scewchuk by Active Sets
%x0 = [sqrt(8);0];  % top right corner
x0 = [2;0];         % top edge
%x0 = [2;0.5];       % an interior point
%x0 = [0;0];         % top left corner
%x0 = [0; sqrt(8)];  % bottom corner
lambda = 0; % fixed regularisation term


P = [11 -9; -9 11];%eye(2);
xsol = [0.5;1.5];%[0.3;1.3];      % specify location of minimum INTERIOR
%xsol = [2.3;1.3];      % specify location of minimum EXTERIOR
%xsol = [1.5;2.5];%[1.3;2.3];      % specify location of minimum EXTERIOR
q = -P*xsol; 
xsol = P\(-q);
A = [-1 0 sqrt(0.5); 0 -1 sqrt(0.5)];
b = [0 ;0 ;3/sqrt(2)];
[v,s,x] = ActiveSet(P,q,A,b,x0);
xAS = x(:,end);
disp(['Starting Point ']);num2str(x0)

disp(['Final from Active Sets ']);num2str(xAS)
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
        xx = [d(i);d(j)];
        im(j,i) = (0.5*xx'*P*xx + q'*xx) ;
    end
end
subplot(2,2,4);
imagesc((im));title('Active Set Constraint Optimisation');
hold on;
for k = 1:size(A,2) % plot the constraints
   h1 = [xlo;xlo]'*A(:,k) - b(k);
   h2 = [xlo;xhi]'*A(:,k) - b(k);
   h3 = [xhi;xlo]'*A(:,k) - b(k);
   h4 = [xhi;xhi]'*A(:,k) - b(k);
   p(:,1) = [xlo;xlo] - h1*A(:,k);
   p(:,2) = [xlo;xhi] - h2*A(:,k);
   p(:,3) = [xhi;xlo] - h3*A(:,k);
   p(:,4) = [xhi;xhi] - h4*A(:,k);
   plot((p(1,:)-xlo)/dx+1,(p(2,:)-xlo)/dx+1,'w');
end
% plot the solutions
plot((x(1,:)-xlo)/dx+1,(x(2,:)-xlo)/dx+1,'k')
plot((x(1,:)-xlo)/dx+1,(x(2,:)-xlo)/dx+1,'m*')
plot((xsol(1)-xlo)/dx+1,(xsol(2)-xlo)/dx+1,'w*');

%im2 = d'*d;
%imagesc(im2);

% fmincon solution
xfm = fmincon(@(x) (0.5*x'*P'*x + q'*x),x0,A',b);
disp(['Final fmincon ']);num2str(xfm)
