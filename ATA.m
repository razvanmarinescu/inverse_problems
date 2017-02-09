function out = ATA(f, alpha, blurSize)
% assumes f is a matrix represented in 1D

dim = sqrt(size(f));
f1d = reshape(f, [dim, dim]);

y = addBlur(f1d, blurSize);
z = addBlur(y, blurSize) + alpha*f1d;

out = reshape(z, [numel(z), 1]);

end