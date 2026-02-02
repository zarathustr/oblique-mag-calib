function [T,h] = solve_ls_T_h_given_M(Y, M)
%SOLVE_LS_T_H_GIVEN_M Solve min ||Y - h*1^T - T*M||_F with T upper triangular.
% Unknown vector is x = [h1 h2 h3 t11 t12 t13 t22 t23 t33]^T

N = size(Y,2);
A = zeros(3*N, 9);
b = zeros(3*N, 1);

for k = 1:N
    mk = M(:,k);
    yk = Y(:,k);
    r = 3*(k-1);
    b(r+1:r+3) = yk;

    % y1 = h1 + t11*m1 + t12*m2 + t13*m3
    A(r+1,1) = 1;
    A(r+1,4) = mk(1); A(r+1,5) = mk(2); A(r+1,6) = mk(3);

    % y2 = h2 + t22*m2 + t23*m3
    A(r+2,2) = 1;
    A(r+2,7) = mk(2); A(r+2,8) = mk(3);

    % y3 = h3 + t33*m3
    A(r+3,3) = 1;
    A(r+3,9) = mk(3);
end

x = A \ b;
h = x(1:3);
t11 = x(4); t12 = x(5); t13 = x(6); t22 = x(7); t23 = x(8); t33 = x(9);

T = [t11, t12, t13;
     0,   t22, t23;
     0,   0,   t33];
end
