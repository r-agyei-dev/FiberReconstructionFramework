function Xnorm = normc_custom(X)
    % Normalize each column of X to unit length
    Xnorm = X ./ vecnorm(X);
end