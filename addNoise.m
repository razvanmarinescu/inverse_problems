function imgNew = addNoise(img, sigma)

% need to add noise too

imgNew = img + sigma* randn(size(img));

end