function out = AsqrtAI(f, alpha, blurSize, transpFlag)
% function that returns (A, sqrt(alpha)I)^T f, where A is addBlur function
% f is a flattened 1d image, followed by the same number of zeros 
if strcmp(transpFlag,'transp')
  % A' = NxM, f=Mx1, M/2 
  dim = sqrt(size(f,2)/2);

  f2d = reshape(f(1:size(f,1)/2), [dim, dim]);

  y = addBlur(f2d, blurSize);
  y1d = reshape(y, [numel(y), 1]);

  out = [y1d; sqrt(alpha)*y1d];
elseif strcmp(transpFlag,'notransp')
  %A = MxN, f=Nx1, M/2
  dim = sqrt(size(f,1)/2);

  f2d = reshape(f(1:size(f,1)/2), [dim, dim]);

  y = addBlur(f2d, blurSize);
  y1d = reshape(y, [numel(y), 1]);

  out = [y1d; sqrt(alpha)*y1d];
end

end