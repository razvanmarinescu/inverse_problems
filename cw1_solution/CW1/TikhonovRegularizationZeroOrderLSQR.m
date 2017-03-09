function falpha = TikhonovRegularizationZeroOrderLSQR(Aker,g,alpha)
% TIKHONOVREGULARIZATIONL2 computes the solution of 
% || A f - g ||_2^2 + alpha || f ||_2^2
%
%  falpha = TikhonovRegularizationZeroOrder(Aker,g,alpha)
%
%  INPUT:
%   Aker - blurring kernel
%   g - an image to deblurr
%   alpha - the regularization parameter
%
%  OUTPUTS:
%   falpha - the solution
%
% Author: Felix Lucka

[Nx, Ny] = size(g);


% prepare handle for solving the augmented equation with PCG
AaugHandle = @(f,transposeFlag) Aaug(f,alpha,Aker,transposeFlag,[Nx,Ny]);
rightHandSide = [g(:);zeros(Nx*Ny,1)];


% predefine tolerance and maximal number of iterations
tol = 10^-8;
maxIter = max(Nx,Ny);

[falpha,flag] = lsqr(AaugHandle,rightHandSide,tol,maxIter);

falpha = reshape(falpha,Nx,Ny);
end