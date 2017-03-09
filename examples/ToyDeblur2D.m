%%
Nx = 256;
sig = 6;
nlevel = 1;
itol = 1e-16;
MultiplicativeNoise = 1;

im = phantom(Nx);
H = fspecial('gaussian',Nx,sig);
figure(1);clf; hold on;
subplot(2,5,1); imagesc(im);title('true image');
subplot(2,5,2); imagesc(H);title('blurring kernel');

fim = fft2(im); fH = fft2(H);
subplot(2,5,6); imagesc(fftshift(log(abs(fim)))); 
subplot(2,5,7); imagesc(fftshift(log(abs(fH)))); 


%%
fimH = fim.*fH;
imH = ifftshift(ifft2(fimH));
figure(1);
subplot(2,5,3); imagesc(abs(imH));title('blurred image');
subplot(2,5,8); imagesc(fftshift(log(abs(fimH)))); 


%%

fimHn = fft2(imH);
if MultiplicativeNoise
fimHn = fimH.*(1+nlevel*randn(size(im)));
else
fimHn = fimH +nlevel*max(abs(fimH(:)))*randn(size(im));
end
imHn = ifftshift(ifft2(fimHn));
figure(1);
subplot(2,5,3); imagesc(abs(imHn));title('blurred image + noise');
subplot(2,5,8); imagesc(fftshift(log(abs(fimHn)))); 



%%

fimR = (fimH./(fH+itol));
imR = ifft2(fimR);
figure(1);
subplot(2,5,4); imagesc(abs(imR));title('deblurred from noiseless data')
subplot(2,5,9); imagesc(fftshift(log(abs(fimR))));


%%
fimRn = (fimHn./(fH+itol));
imRn = ifft2(fimRn);
figure(1);
subplot(2,5,5); imagesc(abs(imRn));title('deblurred from noisy data')
subplot(2,5,10); imagesc(fftshift(log(abs(fimRn))));

