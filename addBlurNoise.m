function imgNew = addBlurNoise(img, blurSize, sigma)

imgNew = addBlur(img, blurSize);

% add noise
imgNew = imgNew + sigma* randn(size(img));

end