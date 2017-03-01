imOrig = imread('Cameraman256.png', 'png');
imOrig = double(imOrig);
subplot(2,2,1);
imagesc(imOrig); 
colormap(gray);

sigma = 0; % high noise makes the problem also harder 
blurSize = 3; % high blur makes the problem ill-conditioned
alpha = 0; % high regularisation makes the edges vetry smooth


subplot(2,2,2);
imgBlur = addBlurNoise(imOrig, blurSize, sigma);

imagesc(imgBlur); 


% deconvolve using preconditioned conjugate gradients

imgNewBlurred1d = reshape(addBlur(imgBlur, blurSize), [numel(imgBlur), 1 ]);
%solve the system (A^TA + alpha*I)f = A^Tg
tol = 1e-6;
maxit = 10;
[fa,FLAG,RELRES,ITER,RESVEC] = pcg(@(x) ATA(x, alpha, blurSize), imgNewBlurred1d,... 
tol, maxit, [], [], imgNewBlurred1d);

reconstPCG = reshape(fa, size(imgBlur));
subplot(2,2,3);
colormap(gray);
imagesc(imOrig - reconstPCG); 

% solve the system using least squares: 
imgBlurredZeros = [imgNewBlurred1d; zeros(size(imgNewBlurred1d))];
funcHandle = @(x, transpFlag) AsqrtAI(x, alpha, blurSize, transpFlag);
t = funcHandle(imgBlurredZeros, 'transp')
tp = funcHandle(imgNewBlurred1d, 'notransp')
faLsqr = lsqr(funcHandle, imgBlurredZeros);

reconstLSQR = reshape(faLsqr, size(imgBlur));
subplot(2,2,4);
colormap(gray);
imagesc(reconstLSQR); 


imB = imread('boat.tiff', 'tiff');

