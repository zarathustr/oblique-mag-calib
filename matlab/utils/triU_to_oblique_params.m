function [sigma,u2,u3] = triU_to_oblique_params(T)
%TRIU_TO_OBLIQUE_PARAMS Convert upper triangular T into (sigma,u2,u3).
% sigma are log-scales. u2 and u3 are unit vectors on S^1 and S^2.

t11 = T(1,1);
sigma1 = log(abs(t11) + eps);

col2 = [T(1,2); T(2,2)];
s2 = norm(col2) + eps;
u2 = col2 / s2;
if u2(2) < 0, u2 = -u2; end
sigma2 = log(s2);

col3 = [T(1,3); T(2,3); T(3,3)];
s3 = norm(col3) + eps;
u3 = col3 / s3;
if u3(3) < 0, u3 = -u3; end
sigma3 = log(s3);

sigma = [sigma1; sigma2; sigma3];
end
