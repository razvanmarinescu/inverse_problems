% This is a reference solution for the course work 1 of GV08 Optimization
%
% Author: Felix Lucka, f.lucka@ucl.ac.uk

% always invoke the three essential c's to set the stage
clear all
clc
close all

%% 1. Read an image into Matlab, perform Gaussian convolution, and add noise

% read in the camera man image
f = imread('Cameraman256.png');
% convert it to double precision
f = double(f);
% rescale it to [0,1]
f = f/max(f(:));


% get the size of the image in pixels
[Nx,Ny] = size(f);

% define function handles that convert an image into a vector and vice versa
% look up "Anonymous Functions" in the documentation if your are not
% familiar with this way of defining functions on the fly.
im2vec = @(image) reshape(image,[],1);
vec2im = @(vector) reshape(vector,Nx,Ny); % will only work with the Nx, Ny assigned above!


% construct Gaussian blurring kernel with standart deviation 1.5 (in units of
% pixel length) and a size of 5 standart deviations (must have an odd size,
% otherwise, the transpose is not correct)
sigmaA = 1.5;
sizeKernel = ceil(sigmaA * [5,5]);
sizeKernel = sizeKernel + mod(sizeKernel+1,2);
Aker = fspecial('gaussian',sizeKernel,sigmaA);
% define function handle for A applied to an image (will also return an image)
% we assume periodic boundary conditions
A = @(image) imfilter(image,Aker,'circular');

% in the following, we will use that A^T = A. While this is theoretically true, 
% not all discrete implementations of A fulfill this requirement. Therefore, one
% should check that for two random u, v, we have that <Au,v> = <u,A^Tv> (up
% to numerical inaccuracies
u = randn(size(f));
v = randn(size(f));
Au_dot_v = sum(im2vec(A(u).*v));
u_dot_Av = sum(im2vec(u.*A(v)));
Au_dot_v - u_dot_Av


% generate data
gClean = A(f);
% add Gaussian noise with standart deviation of 0.1
sigmaNoise = 0.1;
gNoisy = gClean + sigmaNoise * randn(size(gClean));
gNoisyVec = im2vec(gNoisy);

% plot the setup
figure(1);
subplot(2,2,1); colormap(gray); imagesc(f); title('original image');
subplot(2,2,2); colormap(gray); imagesc(Aker); title('blurring kernel');
subplot(2,2,3); colormap(gray); imagesc(gClean); title('data, clean');
subplot(2,2,4); colormap(gray); imagesc(gNoisy); title('data, noisy');
drawnow();

%% 2. Deconvolve using linear least squares

% regularization parameter alpha
alpha = 1;
% tolerance of the iterative solver
tol = 10^-8;
% maximal number of iterations
maxIter = 500;

%%% 2.1 Deconvolve using normal equations

% define fucntion handle for (A^T A + alpha * I) * f
ATAalphaId = @(f) ATA(f,alpha,Aker,[Nx,Ny]);
rightHandSide = im2vec(A(gNoisy));

% call pcg and measure the time by tic ... toc
tic
fPCG = pcg(ATAalphaId,rightHandSide,tol,maxIter);
toc
fPCG = vec2im(fPCG);
relativeResiduum = norm(im2vec(A(fPCG))-gNoisyVec)/norm(gNoisyVec)
Energy = norm(im2vec(A(fPCG))-gNoisyVec)^2 + alpha * norm(im2vec(fPCG))^2

% plot the solution
figure('Name','PCG');
subplot(1,3,1); colormap(gray); imagesc(f); title('original image');
subplot(1,3,2); colormap(gray); imagesc(gNoisy); title('data noisy');
subplot(1,3,3); colormap(gray); imagesc(fPCG); title('reconstruction, clean');
drawnow();


% now with gmres
tic
fGMRES = gmres(ATAalphaId,rightHandSide,[],tol,maxIter);
toc
fGMRES = vec2im(fGMRES);
relativeResiduum = norm(im2vec(A(fGMRES))-gNoisyVec)/norm(gNoisyVec)
Energy = norm(im2vec(A(fGMRES))-gNoisyVec)^2 + alpha * norm(im2vec(fGMRES))^2

% plot the solution
figure('Name','GMRES');
subplot(1,3,1); colormap(gray); imagesc(f); title('original image');
subplot(1,3,2); colormap(gray); imagesc(gNoisy); title('data noisy');
subplot(1,3,3); colormap(gray); imagesc(fGMRES); title('reconstruction, clean');
drawnow();


% now with lsqr
AaugHandle = @(f,transposeFlag) Aaug(f,alpha,Aker,transposeFlag,[Nx,Ny]);
rightHandSide = [gNoisyVec;zeros(Nx*Ny,1)];

tic
fLSQR = lsqr(AaugHandle,rightHandSide,sqrt(tol),maxIter);
toc
fLSQR = vec2im(fLSQR);
relativeResiduum = norm(im2vec(A(fLSQR))-gNoisyVec)/norm(gNoisyVec)
Energy = norm(im2vec(A(fLSQR))-gNoisyVec)^2 + alpha * norm(im2vec(fLSQR))^2

% plot the solution
figure('Name','LSQR')
subplot(1,3,1); colormap(gray); imagesc(f); title('original image');
subplot(1,3,2); colormap(gray); imagesc(gNoisy); title('data noisy');
subplot(1,3,3); colormap(gray); imagesc(fLSQR); title('reconstruction, clean');
drawnow();

%% 3.1 Choose a regularization parameter by the discrepancy principle

% we wrapped up the PCG way to solve the Tikhonov regularization problem
% in a regular function and use it to define a function handle
TikhonovHandle = @(alpha) TikhonovRegularizationZeroOrder(Aker,gNoisy,alpha);
% now we define a function handle for the discrepacy function
residuum = @(alpha) norm(im2vec(A(TikhonovHandle(alpha))-gNoisy))^2;
discrepancy = @(alpha) 1/(Nx*Ny) * residuum(alpha) - sigmaNoise^2;

% we first plot the discrepancy for a range of alphas
alpha_vec = 10.^(-2:0.1:0);
discrepancy_vec = zeros(size(alpha_vec));

for i=1:length(alpha_vec)
    i
    discrepancy_vec(i) = discrepancy(alpha_vec(i));
end

figure('Name','Discrepancy Plot')
semilogx(alpha_vec,discrepancy_vec);


% we will use fzero to find the zero of the discrepancy function
fzeroOptions = [];
fzeroOptions.Display = 'iter';
fzeroOptions.TolX = 10^-6;
alphaDiscrepancy = fzero(discrepancy,10.^[-2,0],fzeroOptions)

% plot the solution
figure('Name','alpha discrepancy')
subplot(1,3,1); colormap(gray); imagesc(f); title('original image');
subplot(1,3,2); colormap(gray); imagesc(gNoisy); title('data noisy');
subplot(1,3,3); colormap(gray); imagesc(TikhonovHandle(alphaDiscrepancy)); title('reconstruction, clean');
drawnow();

%% 3.2 Choose a regularization parameter by the Miller criterion

millerAux = @(f) 1/(sigmaNoise^2) * norm(im2vec(A(f)-gNoisy))^2 - alpha * norm(im2vec(f))^2;
miller    = @(alpha) millerAux(TikhonovHandle(alpha));

% we first plot the Miller criterion and the L-curve for a range of alphas
alpha_vec = 10.^(-4:0.25:0);
miller_vec = zeros(size(alpha_vec));

for i=1:length(alpha_vec)
    i
    miller_vec(i) = miller(alpha_vec(i));
end

figure('Name','Miller Plot')
semilogx(alpha_vec,miller_vec);

% we will use fzero to find the zero of the Miller function
alphaMiller = fzero(miller,10.^[-4,0],fzeroOptions)

% plot the solution
figure('Name','alpha miller')
subplot(1,3,1); colormap(gray); imagesc(f); title('original image');
subplot(1,3,2); colormap(gray); imagesc(gNoisy); title('data noisy');
subplot(1,3,3); colormap(gray); imagesc(TikhonovHandle(alphaMiller)); title('reconstruction, clean');
drawnow();

%% 4. use first order regularization; for simplicity, we only repeat the last step

% we wrapped up the computation in a regular function and use it to define a function handle
TikhonovHandle = @(alpha) TikhonovRegularizationFirstOrder(Aker,gNoisy,alpha);
% now we define a function handle for the discrepacy function
discrepancy = @(alpha) 1/(Nx*Ny) * norm(im2vec(A(TikhonovHandle(alpha))-gNoisy))^2 - sigmaNoise^2;

alphaDiscrepancy = fzero(discrepancy,10.^[-1,2],fzeroOptions);
fFirstTikh = TikhonovHandle(alphaDiscrepancy);

% plot the solution
figure('Name','First Order Tikhonov')
subplot(1,3,1); colormap(gray); imagesc(f); title('original image');
subplot(1,3,2); colormap(gray); imagesc(gNoisy); title('data noisy');
subplot(1,3,3); colormap(gray); imagesc(fFirstTikh); title('reconstruction, clean');
drawnow();

%% 5. Construct an anisotropic derivative filter

% define discrete derivatives in x and y direction
Dx = @(x)  [diff(x,1,2),zeros(size(x,1),1)];
Dy = @(x)  [diff(x,1,1);zeros(1,size(x,2))];
% compute the norm of the gradient
gradientNorm = sqrt(Dx(gNoisy).^2+Dy(gNoisy).^2);
Tf = 1;
T = Tf * max(gradientNorm(:));

kappa = exp(- gradientNorm/T);
filterMatrix = spdiags(kappa(:),0,Nx*Ny,Nx*Ny);

%% 6. Hierarchical deblurring
close all

alpha = 5;

% start from the noisy data
fAniso = gNoisy;
nIter = 5;
discrepancy_vec = zeros(nIter,1);

figure('Name','iterative anisotropic deblurring')

for i = 1:nIter
    gradientNorm = sqrt(Dx(fAniso).^2+Dy(fAniso).^2);
    T = 1/alpha * max(gradientNorm(:));
    fAniso = AnisotropicFiltering(Aker,gNoisy,fAniso,T,alpha);
    discrepancy_vec(i) = 1/(Nx*Ny) * norm(im2vec(fAniso-gNoisy))^2 - sigmaNoise^2;
    
    % plot the current iterate and its gradient image
    subplot(2,nIter,i); colormap(gray); imagesc(fAniso); title(['f. iter: ' int2str(i)]);
    subplot(2,nIter,i+nIter); colormap(gray); imagesc(gradientNorm); title(['gradnorm(f), iter: ' int2str(i)]);
    drawnow();
end

figure('Name','iterative anisotropic deblurring, discrepancy')
plot(1:nIter,discrepancy_vec)
