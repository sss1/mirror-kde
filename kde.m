% Inputs:
%   X - n-by-d matrix of n d-dimensional data points
%   K - smoothing kernel (1d function handle)
%   h - bandwidth
%
% Outputs:
%   p - kernel density estimate (function handle)
%   

function p = kde(X, K, h)

  p = @(x) point_mirror_kde(x, X, K, h);

end

function p = point_mirror_kde(x, X, K, h)

  [n, d] = size(X);

  dists = bsxfun(@minus, X, x);

  p = 0;
  parfor i=1:n % for each data point X(i,:)
    p = p + prod(arrayfun(@(dist) K(dist/h), dists(i,:)));
  end

  p = p/(n*h^d);

end
