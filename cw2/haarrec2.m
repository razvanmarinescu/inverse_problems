function X = haarrec2(Coeff)
% HAARREC2 reconstructs an image from the Haar wavelet decomposition given by
% Coeff
%
% X = haarrec2(Coeff)
%
%  INPUT:
%   Coeff - the Haar wavelet coefficients 
%
%  OUTPUTS:
%    X - the reconstructed image
%
% Author: Felix Lucka


xy_convention = false;
n = sqrt(length(Coeff));
N = log2(n);


if(mod(n,1) ~= 0 || mod(N,1) ~= 0)
    error('Coeff must be of length n * n, with n = 2^N, N integer.')
end

fac = 1/sqrt(2);

X = Coeff(1);
Coeff(1) = [];

for i = 0:N-1
    length_details = 2^(2*i);
    LH = Coeff(1:length_details);
    HL = Coeff(length_details+1:2*length_details);
    HH = Coeff(2*length_details+1:3*length_details);
    Coeff(1:3*length_details) = [];
    
    % Reconstruct L and H
    if(xy_convention)
        H_er = fac*(HL-HH);
        H_ur = fac*(HL+HH);
        L_er = fac*(X(:)-LH);
        L_ur = fac*(X(:)+LH);
    else
        H_er = fac*(HL+HH);
        H_ur = fac*(HL-HH);
        L_er = fac*(X(:)+LH);
        L_ur = fac*(X(:)-LH);
    end
    H = [H_ur';H_er'];
    H = reshape(H,2^(i+1),2^i);
    L = [L_ur';L_er'];
    L = reshape(L,2^(i+1),2^i);
    X_uc = fac*(L+H);
    X_ec = fac*(L-H);
    X = zeros(2^(i+1));
    X(:,1:2:end) = X_uc;
    X(:,2:2:end) = X_ec;
end

end
