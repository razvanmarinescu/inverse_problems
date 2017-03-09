clear all

% number of samples
N = 1000;

% target density
sigma = 0.1;
delta = 1;
Pi = @(x) exp(-0.5/sigma^2*(sqrt(sum(x.^2))-1)^2 - 0.5/delta^2*(x(2)-1)^2);

disp('starting MCMC...')

% the initial sampling point is at the origin
x = [0,0]';         % initial point
n = 2;
gamma = 0.1;       % random walker step size
Ens = zeros(n,N);   % ensemble
Ens(1:n,1) = x;

% target density value of the initial point
Pi0 = Pi(x);
JJ = 2;
acc= 0;
% Metropolis-Hastings loop
for JJ = 2:N
    % random walker proposal
    p = x + gamma*randn(n,1);
    % Compute the target density
    Pi1 = Pi(p);
    % compute the acceptance ratio  --- simplified
    al = log(Pi1)-log(Pi0);
    if (al > log(rand))
        x = p;
        Pi0 = Pi1;
        Ens(1:n,end) = x;
        acc = acc + 1;
    end
    Ens(:,JJ) = x;
end

% postprocessing
mu = mean(Ens(:,N/2:N),2);
figure(1) 
axes('fontsize',12);
hold on
plot([-1.5 1.5],[mu(2), mu(2)],'k','linewidth',2)
plot([mu(1) mu(1)],[-1.5, 1.5],'k','linewidth',2)

N = 200;
[x,y]=meshgrid(linspace(-1.5,1.5,N));
target = zeros(N,N);
for i = 1:N
  for j = 1:N
    target(i,j) = Pi([x(i,j),y(i,j)]);
  end
end
contour(x,y,target);

N = size(Ens,2);
t = linspace(0,2*pi,20);
r = 0.02;
for i = 1:N
  % r = Ens(n+1,i)/scale;
  x = r*cos(t) + Ens(1,i);
  y = r*sin(t) + Ens(2,i);
  fill(x,y,'r')
end
axis equal
axis([-1.5 1.5 -1.5 1.5])

% trace plot for the first component
figure(2), axes('fontsize',16); plot(Ens(1,:))