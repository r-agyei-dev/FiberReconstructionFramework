%% normr - normalize rows
function Xnorm = normr_custom(X)
    % Normalize each row of X to unit length
    Xnorm = X ./ vecnorm(X, 2, 2);
end