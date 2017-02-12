imC = imread('Cameraman256.png', 'png');
imC = double(imC);
subplot(2,2,1);
imagesc(imC); 
colormap(gray);

sigma = 10;
blurSize = 1;
subplot(2,2,2);
imgNew = addBlurNoise(imC, blurSize, sigma);

imagesc(imgNew); 

alpha = 1;
% deconvolve using preconditioned conjugate gradients

imgNewBlurred1d = reshape(addBlur(imgNew, blurSize), [numel(imgNew), 1 ]);
%solve the system (A^TA + alpha*I)f = A^Tg
fa = pcg(@(x) ATA(x, alpha, blurSize), imgNewBlurred1d);

reconstPCG = reshape(fa, size(imgNew));
subplot(2,2,3);
colormap(gray);
imagesc(reconstPCG); 

% solve the system using least squares: 
imgBlurredZeros = [imgNewBlurred1d; zeros(size(imgNewBlurred1d))];
faLsqr = lsqr(@(x) AsqrtAI(x, alpha, blurSize), imgBlurredZeros);

reconstLSQR = reshape(faLsqr, size(imgNew));
subplot(2,2,4);
colormap(gray);
imagesc(reconstLSQR); 


imB = imread('boat.tiff', 'tiff');

