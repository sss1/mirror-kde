% experiment parameters
alpha = 2;
n_trials = 24;
ns = [200 300 400 500 700 1000 1500 2000 2500 3000 5000 8000];
h = 0.1;
kernel = @(x) max((1 - x^2)*3/4, 0); % Epanechnikov kernel
% kernel = @(x) normpdf(x); % Gaussian kernel

% allocate space for results
Is = zeros(n_trials,length(ns));
ts = zeros(n_trials,length(ns));
    
for ni=1:length(ns)
  n = ns(ni);
  parfor trial=1:n_trials

    % generate data
    Zs = unifrnd(0, 1, n, 1);
    Ys = normpdf(Zs, 1);
    Xs = normpdf(Zs, 1);
    
    % Is(trial,ni) = CMI_est(kernel, h, alpha, Xs, Ys, Zs);
    Is(trial,ni) = CMI_est(kernel, h, alpha, Xs, Ys);
  end
  n
end

save('uresults.mat','Is');

% 
% % eliminate some extra computations in certain cases
% if alpha == 0
%   sliceInt = @(z) MC_Int(@(xy:w) (p_XZYZ([xy z]))/(p_Z(z)^2), n_MC, 2, 0, max(2, 2^alpha));
% elseif alpha == 2
%   sliceInt = @(z) MC_Int(@(xy) ((p_XYZ([x y z])^2)/p_XZYZ(x, y, z), n_MC, 2, 0, max(2, 2^alpha));
% else
%   sliceInt = @(z) MC_Int(@(xy) (p_XYZ([x y z])^alpha)*(p_XZYZ(x, y, z)^(1 - alpha))*(p_Z(z)^(alpha - 2)), n_MC, 2, 0, max(2, 2^alpha));
% end

% Monte Carlo formulation
% I = MC_Int(@(z) p_Z(z)*log(sliceInt(z)), n_MC, 1, -10, 10)/(1 - alpha);
% samps = unifrnd(0, 1, 1000, 3);
% tic
% for i=1:500
%   tmp = p_XYZ(samps(i,:));
% end
% toc

% integral/integral2 formulation
% sliceInt = @(z) integral2(@(x, y) (p_XYZ([x y z])^alpha)*((p_XZ(x,z)*p_YZ(y,z))^(1 - alpha))(*p_Z(z)^(alpha - 2)), 0, 1, 0, 1);
% I = integral(@(z) p_Z(z)*log(sliceInt(z)), 0, 1)/(1 - alpha);
