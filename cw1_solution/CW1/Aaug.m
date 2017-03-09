function z = Aaug(f,alpha,Aker,transposeFlag,sizeImage)
% Aaug computes G * f = [A;sqrt(alpha) * Id] * f for a blurring A
%  z = Aaug(f,alpha,Aker,'transp',sizeImage)
%  z = Aaug(f,alpha,Aker,'notransp',sizeImage)
%
%  INPUT:
%   f - a Nx x Ny image as a vector
%   alpha - the regularization parameter
%   Aker - blurring kernel
%   transposeFlag - a string that determines whether G or G' should be applied
%                   'transp' : G'*f is computed
%                   'notransp' : G*f is computed
%   sizeImage - dimensions of the image
%
%  OUTPUTS:
%   z - the result as a vector
%
% Author: Felix Lucka

Nx = sizeImage(1);
Ny = sizeImage(2);


switch transposeFlag
    case 'notransp'
        % reshape f to an image format
        f = reshape(f,Nx,Ny);
        % y = A(f)
        y = imfilter(f,Aker,'circular');
        % linearize and attach sqrt(alpha) f
        z = [y(:);sqrt(alpha)*f(:)];
    case 'transp'
        % splitt f into the two images and reshape them
        f1 = reshape(f(1:Nx*Ny),Nx,Ny);
        f2 = reshape(f((Nx*Ny+1):end),Nx,Ny);
        % y = AT(f)
        y = imfilter(f1,Aker,'circular');
        % linearize and add sqrt(alpha) f
        z = y(:) + sqrt(alpha)*f2(:);
    otherwise
        error('input transposeFlag has to be ''transp'' or ''notransp''')
end

end