function z = ATAgrad(f,alpha,Aker,sizeImage)
% ATA computes (A^T A + alpha DxT Dx + alpha DyT Dy) * f for a blurring A
% and the spatial difference operators Dx and Dy
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

A  = @(x)  imfilter(x,Aker,'circular');
Dx = @(x)  [diff(x,1,2),zeros(size(x,1),1)];
Dy = @(x)  [diff(x,1,1);zeros(1,size(x,2))];
DxT = @(x) [-x(:,1),-diff(x(:,1:end-1),1,2),x(:,end-1)];
DyT = @(x) [-x(1,:);-diff(x(1:end-1,:),1,1);x(end-1,:)];


% reshape f to an image format
f = reshape(f,sizeImage(1),sizeImage(2));

z = A(A(f)) + alpha * (DxT(Dx(f)) + DyT(Dy(f)));
z = z(:);
end