function z = ATA(f,alpha,Aker,sizeImage)
% ATA computes (A^T A + alpha Id) * f for a blurring A
%  z = ATA(f,alpha,Aker,[Nx,Ny])
%
%  INPUT:
%   f - a Nx x Ny image as a vector
%   alpha - the regularization parameter
%   Aker - blurring kernel
%   sizeImage - dimensions of the image
%
%  OUTPUTS:
%   z - the result as a vector
%
% Author: Felix Lucka

% reshape f to an image format
f = reshape(f,sizeImage(1),sizeImage(2));
% y = A(f)
y = imfilter(f,Aker,'circular');
% z = A^T y + alpha * f
z = imfilter(y,Aker,'circular') + alpha * f;
% linearize and return z
z = z(:);

end