% This is a reference solution for the exercise 1 of GV08 Optimization
%
% Author: Felix Lucka, f.lucka@ucl.ac.uk

% always invoke the three essential c's to set the stage
clear all
clc
close all

%% 1. Minimum p-Norm Solutions to Underdetermined Problems

% we define A and b such that Ax = b describes the linear problem 
% x_1 + 2 x_2 = 5: 
A = [1,2];
b = 5;

% we define a function handle for Phi
Phi = @(x,p) sum(abs(x).^p);
% we define a vector of values for p
pArray = 1:0.5:5;

% we loop over the different values for p
for iP=1:length(pArray)
    % we call fmincon to solve the problem
    % min Phi(x,p) subject to Ax = b
    % for a fixed p
    % using "@(x) Phi(x,pArray(iP))" as an argument is the same as 
    % defining PhiP = @(x) Phi(x,pArray(iP)) 
    % for a given, fix value for p
    xNorm(iP,:) = fmincon(@(x) Phi(x,pArray(iP)),[0;0],[],[],A,b);
end

%% now we produce a plot like Figure 1 

% the same axis limits 
xLim = [-0.5,3];
yLim = xLim;
% we use this to plot the set of all solutions to Ax = b, (given by x_2 =
% (5 - x_1)/2)
yEquation = (5 - xLim)/2;


% Create figure
figure1 = figure;
% create axis
axes1 = axes('Parent',figure1);
hold(axes1,'on');
% plot a blue line for the set of all solutions
plot(xLim,yEquation)
% plot the minimum p-norm solutions as red dots
plot(xNorm(:,1),xNorm(:,2),'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'Marker','*','LineStyle','none');
% set the axis to the defined limits
xlim(axes1,xLim); ylim(axes1,yLim);
% put a box around the axis
box(axes1,'on');
% set the scaling of the axis to the same values
set(axes1,'DataAspectRatio',[1 1 1],'XGrid','on','YGrid','on');
hold(axes1,'off');

% export the figure as eps, pdf and png using saveas.m
saveas(figure1,'LeastSquaresSolutions.eps','epsc')
saveas(figure1,'LeastSquaresSolutions.pdf','pdf')
saveas(figure1,'LeastSquaresSolutions.png','png')
% export the figure as png using print.m with a higher resolution
print('LeastSquaresSolutionsHiRes.png','-dpng','-r600')

%% compare to other least squares solutions

xPI1        = A' * ((A * A') \ b)
xPI2        = pinv(A) * b
xBackslash  = A\b % gives the sparest least-squares solution (p = 1)


%% 2. Singular Value Decomposition

clear all;close all;clc

%% set up the matrix A by a loop over x and y

n = 100;
y = linspace(0,1,n);
x = linspace(0,1,2*n);
sigma = 0.2;
A = zeros(n,2*n);

for i=1:n
    for j=1:2*n
        A(i,j) = 1/(sqrt(2*pi) * sigma) * exp(-(y(i)-x(j))^2/(2*sigma^2));
    end
end


%% export A as a image 

% using imagesc
figure1 = figure();
imagesc(A);
saveas(figure1,'Aimage1.png','png')

% using imwrite 
Aimg = ceil(A/max(A(:))*256);
colorMap = parula(256);
imwrite(Aimg,colorMap,'Aimage2.png')

%% compute SVD

[U,W,V] = svd(A);
whos('U','W','V')

% check whether A is close to U*W*V'
A_SVD = U*W*V';
norm(A(:)-A_SVD(:),'inf')

% construct the pseudoinverse of W
Wpinv = spdiags(1./diag(W),0,2*n,n);
% plot its structure by using the spy.m function
figure(); spy(Wpinv);

% build the pseudoinverse of A relying on the SVD and by the pinv.m
% function and compare them
APinv_SVD = V*Wpinv*U';
APinv     = pinv(A);
norm(APinv(:)-APinv_SVD(:),'inf') % should be 0 for n = 10 but explode for n > 20


%% plot the first and last 9 singular vectors and the decay of singular values

figure2 = figure();
subplot(1,2,1)
plot(V(:,1:9))
subplot(1,2,2)
plot(V(:,end-9:end))

figure();
plot(log(diag(W)))