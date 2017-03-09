function out = AsqrtAI(f, alpha, blurSize, transpFlag)
% function that returns (A, sqrt(alpha)I)^T f, where A is addBlur function
% f is a flattened 1d image, followed by the same number of zeros 
if strcmp(transpFlag,'transp')
  % A' = NxM, f=Mx1, out=Nx1, N=M/2 
  N = sqrt(size(f,1)/2); % checked
  
  f1Half = f(1:size(f,1)/2);
  f2Half = f((size(f,1)/2+1):end);
  f1Half2d = reshape(f1Half, [N, N]); 

  f1Half2dBlur = addBlur(f1Half2d, blurSize);
  f1Half2dBlur1d = reshape(f1Half2dBlur, [numel(f1Half2dBlur), 1]);

  out = f1Half2dBlur1d + sqrt(alpha)*f2Half;
  assert(size(out,1) == N*N)
  
elseif strcmp(transpFlag,'notransp')
  %A = MxN, f=Nx1, out=Mx1 N=M/2
  N = sqrt(size(f,1));
  assert(floor(N) == N)
  
  f2d = reshape(f, [N, N]);

  f2dBlur = addBlur(f2d, blurSize);
  f2dBlur1d = reshape(f2dBlur, [numel(f2dBlur), 1]);

  out = [f2dBlur1d; sqrt(alpha)*f];
  
  assert(size(out,1) == 2*N*N)
end

end