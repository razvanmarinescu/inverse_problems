% This example adds some levels of Poisson and Guassian noise to an image

im = double(phantom(256));
im = im - 2*min(min(im)); % just to ensure positivity.
sqim = sqrt(im);

figure(5);
subplot(2,2,1);
imagesc(im);colorbar;title('original');
subplot(2,2,2);
scale = 1e12;
yn12  = scale*imnoise(im/scale,'poisson');
imagesc(yn12);colorbar;title('Poisson noise 1e12');
subplot(2,2,3);
scale = 1e10;
yn10  = scale*imnoise(im/scale,'poisson');
imagesc(yn10);colorbar;title('Poisson noise 1e10');
subplot(2,2,4);
scale = 1e8;
yn8  = scale*imnoise(im/scale,'poisson');
imagesc(yn8);colorbar;title('Poisson noise 1e8');
%subplot(2,2,3);
%imagesc(yn-im);colorbar;title('Poisson noise');

% compare to Gaussian noise
figure(6);

subplot(2,2,1);
imagesc(im);colorbar;title('original');
subplot(2,2,2);
yg12 = im + sqrt(im).*randn(size(im));
imagesc(yg12);colorbar;title('Gaussian noise 1e12');
subplot(2,2,3);
yg10 = im + 0.1*sqrt(im).*randn(size(im));
imagesc(yg10);colorbar;title('Gaussian noise 1e10');
subplot(2,2,4);
yg8 = im + 0.01*sqrt(im).*randn(size(im));
imagesc(yg8);colorbar;title('Gaussian noise 1e8');

