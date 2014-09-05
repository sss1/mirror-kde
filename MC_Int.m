% Uses n-sample Monte Carlo to integrate f over [0,1]^d, given lower and upper
% bounds on f
%
% Inputs:
%   f - function from [0,1]^d to R to be integrated (function handle)
%   n - number of samples to be use for Monte Carlo integration
%   d - dimension of domain of f
%   lower - lower bound of f
%   upper - upper bound of f
%
% Outputs:
%   I - estimated integral of f over [0,1]^d

% TODO: add option to estimate variance
function I = MC_Int(f, n, d, lower, upper)

  % generate sample points
  xs = unifrnd(0, 1, n, d);
  ys = unifrnd(lower, upper, n, 1);

  above = 0;
  for i=1:n
    above = above + (f(xs(i,:)) > ys(i));
  end

  I = (upper - lower)*above/n + lower;

end
