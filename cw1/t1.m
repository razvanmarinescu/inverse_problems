imOrig = imread('Cameraman256.png', 'png');
imOrig = double(imOrig);
subplot(2,2,1);
imagesc(imOrig); 
colormap(gray);
imOrig1d = reshape(imOrig, [numel(imOrig), 1 ]);

sigma = 0.1; % high noise makes the problem also harder 
blurSize = 3; % high blur makes the problem ill-conditioned

%alpha = 0.00;
alpha = 0.00001; % high regularisation makes the edges vetry smooth


subplot(2,2,2);
g = addBlurNoise(imOrig, blurSize, sigma);
g1D = reshape(g, [numel(g), 1 ]);

imagesc(g); 

% solve using normal equations (A^T * A + alpha * I)*fa = A^T * g
% deconvolve using preconditioned conjugate gradients

aTg = reshape(addBlur(g, blurSize), [numel(g), 1 ]);
%solve the system (A^TA + alpha*I)f = A^Tg
tol = 1e-6;
maxit = 1000;
ataFunc = @(x) ATA(x, alpha, blurSize);
[fa,FLAG,RELRES,ITER,RESVEC] = pcg(ataFunc, aTg,... 
tol, maxit, [], [], aTg);
%[fa,FLAG,RELRES,ITER,RESVEC] = pcg(@(x) ATA(x, alpha, blurSize), imgNewBlurred1d);

reconstPCG = reshape(fa, size(g));
subplot(2,2,3);
colormap(gray);
imagesc(reconstPCG); 

subplot(2,2,4);
colormap(gray);
imagesc(imOrig - reconstPCG); 

% solve the system using least squares: 
imgBlurredZeros = [g1D; zeros(size(g1D))];
funcHandle = @(x, transpFlag) AsqrtAI(x, alpha, blurSize, transpFlag);
t = funcHandle(imgBlurredZeros, 'transp');
tp = funcHandle(aTg, 'notransp');
[faLsqr, FLAG_lsqr,RELRES_lsqr,ITER_lsqr,RESVEC_lsqr,LSVEC_lsqr] = lsqr(funcHandle, imgBlurredZeros, [], 1000);

reconstLSQR = reshape(faLsqr, size(g));
subplot(2,2,4);
colormap(gray);
imagesc(reconstLSQR); 


imB = imread('boat.tiff', 'tiff');

