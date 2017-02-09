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

%TODO write A(f) is the addBlurNoise function
fa = pcg(@(x) ATA(x, alpha, blurSize), reshape(addBlur(imgNew, blurSize), [numel(imgNew), 1 ]));

reconstPCG = reshape(fa, size(imgNew));
subplot(2,2,3);
colormap(gray);
imagesc(reconstPCG); 

imB = imread('boat.tiff', 'tiff');

