function falpha = AnisotropicFiltering(Aker,g,fEdges,T,alpha)
% ANISOTROPICFILTERING computes 
% (A^T A + alpha D^T kappa D)^-1 A^T g
% where D is the gradient operator and kappa is an anisotropic diffusion
% filter given as 
% kappa = exp(- abs(fEdges) / T)
%
%  falpha = AnisotropicFiltering(Aker,g,fEdges,alpha)
%
%  INPUT:
%   Aker - blurring kernel
%   g - an image to deblurr
%   fEdges - an approximation of the image to recover
%   T - a threshold based on the maximum expected edge values in the image.
%   alpha - the regularization parameter
%
%  OUTPUTS:
%   falpha - the solution
%
% Author: Felix Lucka

[Nx, Ny] = size(g);

% prepare function handles
A  = @(x)  imfilter(x,Aker,'circular');
Dx = @(x)  [diff(x,1,2),zeros(size(x,1),1)];
Dy = @(x)  [diff(x,1,1);zeros(1,size(x,2))];
DxT = @(x) [-x(:,1),-diff(x(:,1:end-1),1,2),x(:,end-1)];
DyT = @(x) [-x(1,:);-diff(x(1:end-1,:),1,1);x(end-1,:)];
im2vec = @(image) reshape(image,[],1);
vec2im = @(vector) reshape(vector,Nx,Ny);

kappa = exp(- sqrt(Dx(fEdges).^2+Dy(fEdges).^2)/T);


% reshape f to an image format
ATAalphaGradKappa = @(f) im2vec(A(A(vec2im(f))) + alpha * (DxT(kappa .* Dx(vec2im(f))) + DyT(kappa.* Dy(vec2im(f)))));
rightHandSide = reshape(imfilter(g,Aker,'circular'),[],1);

% predefine tolerance and maximal number of iterations
tol = 10^-8;
maxIter =  max(Nx,Ny);

[falpha,flag] = pcg(ATAalphaGradKappa,rightHandSide,tol,maxIter);

falpha = reshape(falpha,Nx,Ny);
end