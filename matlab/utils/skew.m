function S = skew(w)
%SKEW 3x3 skew-symmetric matrix for cross products.
S = [   0   -w(3)  w(2);
      w(3)   0    -w(1);
     -w(2)  w(1)   0  ];
end
