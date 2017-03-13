function out = ATA(f, alpha, blurSize)
% assumes f is a matrix represented in 1D

dim = sqrt(size(f,1));
f2d = reshape(f, [dim, dim]);

y = addBlur(f2d, blurSize);
z = addBlur(y, blurSize) + alpha*f2d;

% figure(2)
% subplot(2,1,1);
% imagesc(f2d); 
% colormap(gray);
% 
% subplot(2,1,2);
% imagesc(z); 
% colormap(gray);

out = reshape(z, [numel(z), 1]);

end