% experiment parameters
alpha = 2;
n_trials = 12;
ns = [500 700 1000 1500 2000 2500 3000 5000 8000]; % 10000];
ns = 24*(floor(ns/24)); % round sample sizes to simplify data splitting
hs = [0.1 0.01];
% kernel = @(x) max((1 - x^2)*3/4, 0); % Epanechnikov kernel
kernel = @(x) normpdf(x); % Gaussian kernel
dim = 3;

% allocate space for results
CMIs = zeros(n_trials,length(ns),length(hs));
MIs = zeros(n_trials,length(ns),length(hs));

for ni=1:length(ns) % for each sample size
  n = ns(ni);

  for hi = 1:length(hs) % for each bandwidth
    h = hs(hi);

    parfor trial=1:n_trials % for each trial
  
      % generate data
      Zs = normrnd(0, 1, n, dim);
      Ys = normrnd(Zs, 1);
      Xs = normrnd(Zs, 1);
      Z2s = normrnd(0, 1, n, dim);

      if trial <= 1
        tic
      end
      % compute estimates of CMI and MI
      CMIs(trial,ni,hi) = CMI_est(kernel, h, alpha, Xs, Ys, Zs);
      MIs(trial,ni,hi) = CMI_est(kernel, h, alpha, Xs, Ys, Z2s);
      if trial <= 1
        toc
      end
    end
    [ni hi] % report completed (sample size, bandwidth) pair

    % save results and parameters
  end
  save('multivar_3d2.mat','CMIs','MIs','ns','hs','alpha','kernel','n_trials');
end
