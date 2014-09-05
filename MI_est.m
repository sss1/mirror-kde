% Estimates the Renyi-alpha Mutual Information of X and Y from respective
% joint samples Xs and Ys of equal length, using kernel kernel and bandwidth h
% Presently only works with univariate X and Y; this will be fixed

function I = MI_est(kernel, h, alpha, Xs, Ys)

    data = [Xs Ys];
    n = size(data, 1);

    % split data into 3 pieces (1 for each density to estimate)
    n3 = n/3;
    dat1 = [Xs(1:n3)];
    dat2 = [Ys((n3 + 1):(2*n3))];
    dat3 = [Xs((2*n3 + 1):(3*n3)) Ys((2*n3 + 1):(3*n3))];
    
    % define KDEs of appropriate PDFs
    % p_XYZ = mirror_kde(X, kernel, h);
    p_X = kde(dat1, kernel, h);
    p_Y = kde(dat2, kernel, h);
    p_XY = kde(dat3, kernel, h);

    I = 0;

    % Importance sampling with the log outside
    for i = 1:n
      I = I + (p_X(Xs(i))*p_Y(Ys(i))/(p_XY(data(i,:))))^(1 - alpha);
    end

    I = log(I/n)/(1-alpha);
end
