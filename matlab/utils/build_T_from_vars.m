function T = build_T_from_vars(sigma, u2, u3)
%BUILD_T_FROM_VARS Build upper triangular T from sphere directions and log-scales.
% sigma: [sigma1;sigma2;sigma3] real
% u2: 2x1, ||u2||=1
% u3: 3x1, ||u3||=1
%
% This parameterization is smooth. Diagonal sign conventions can be enforced
% afterward by a diagonal reflection gauge: T <- T*Q, M <- Q*M, where Q is
% orthogonal and diagonal with entries in {-1, +1}.

u2 = u2(:); u3 = u3(:);

t11 = exp(sigma(1));
col2 = exp(sigma(2)) * u2;
col3 = exp(sigma(3)) * u3;

T = [t11, col2(1), col3(1);
     0  , col2(2), col3(2);
     0  , 0      , col3(3)];
end
