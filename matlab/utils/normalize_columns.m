function Xn = normalize_columns(X)
%NORMALIZE_COLUMNS Columnwise normalization.
nrm = sqrt(sum(X.^2,1)) + eps;
Xn = X ./ nrm;
end
