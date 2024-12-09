function x_dot = rhs_fun(~, x, f, P)
%RHS_FUN Summary of this function goes here
%   Detailed explanation goes here

x_dot = [-P.R/P.M         0;
                1        0] * x + ...
        [f; 0];
x_dot(1) = x_dot(1) -P.K/P.M * abs(x(2))^(1/4)*sign(x(2));

end

