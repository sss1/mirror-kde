% Estimates the Modified Renyi-alpha Conditional Mutual Information of X and Y
% given Z from respective joint samples Xs, Ys, and Zs of equal length, using
% kernel kernel and bandwidth h
% Presently only works with univariate X, Y, and Z; this will be fixed

function I = CMI_est(kernel, h, alpha, Xs, Ys, Zs)

    data = [Xs Ys Zs];
    n = size(data, 1);

    % split data into 4 pieces (1 for each density to estimate)
    n4 = n/4;
    dat1 = [Xs(1:n4) Ys(1:n4) Zs(1:n4)];
    dat2 = [Xs((n4 + 1):(2*n4)) Zs((n4 + 1):(2*n4))];
    dat3 = [Ys((2*n4 + 1):(3*n4)) Zs((2*n4 + 1):(3*n4))];
    dat4 = Zs((3*n4 + 1):(4*n4));
    
    % define KDEs of appropriate PDFs
    % p_XYZ = mirror_kde(X, kernel, h);
    p_XYZ = kde(dat1, kernel, h);
    p_XZ = kde(dat2, kernel, h);
    p_YZ = kde(dat3, kernel, h);
    p_XZYZ = @(xyz) p_XZ([xyz(1) xyz(3)])*p_YZ([xyz(2) xyz(3)]);
    p_Z = kde(dat4, kernel, h);

    I = 0;

    % Importance sampling with the log outside
    for i = 1:n
      I = I + (p_XZYZ(data(i,:))/(p_XYZ(data(i,:))*p_Z(Zs(i))))^(1 - alpha);
    end

    I = log(I/n)/(1-alpha);
end
