function [Coeff,Book] = haardec2(X)
% HAARDEC2 decompose X using a Haar wavelet decomposition.
% the image is assumed to match the area [-1,1]^2 and the coefficients are
% normalized in such a way.
%
% [Coeff,Book] = haardec2(X)
%
%  INPUT:
%   X - a square image, size must be a power of 2
%
%  OUTPUTS:
%   Coeff - the coefficients as a vector 
%   Book  - a bookkeeping of the sizes, see wavedec2
%
% Author: Felix Lucka


[m n] =  size(X);
N = log2(n);
xy_convention = false;

if(m ~= n)
    error('Image must be square')
elseif(mod(N,1) ~= 0)
    error('Imagesize must be multiple of 2')
end

fac = 1/sqrt(2);

Book = [n, n];
Coeff = [];


for i = 1:N
    % low pass in x dir
    L = fac*(X(:,1:2:end)+X(:,2:2:end));
    % high pass in x dir
    H = fac*(X(:,1:2:end)-X(:,2:2:end));
    % low low
    X = fac*(L(1:2:end,:)+L(2:2:end,:));
    % high low
    HL = fac*(H(1:2:end,:)+H(2:2:end,:));
    
    if(xy_convention)
        % low high
        LH = fac*(L(1:2:end,:)-L(2:2:end,:));
        % high high
        HH = fac*(H(1:2:end,:)-H(2:2:end,:));
    else
        % low high
        LH = fac*(L(2:2:end,:)-L(1:2:end,:));
        % high high
        HH = fac*(H(2:2:end,:)-H(1:2:end,:));
    end
    Coeff = [LH(:);HL(:);HH(:);Coeff];
    Book = [size(X);Book];
end

Coeff = [X(:);Coeff];
Book = [size(X);Book];

end
