function out = ATA(f, alpha, blurSize)
% assumes f is a matrix represented in 1D

dim = sqrt(size(f));
f2d = reshape(f, [dim, dim]);

y = addBlur(f2d, blurSize);
z = addBlur(y, blurSize) + alpha*f2d;

out = reshape(z, [numel(z), 1]);

end