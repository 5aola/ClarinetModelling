function U = calc_total_current(Uf, y_dot, params)
%CALC_TOTAL_CURRENT Summary of this function goes here
%   Detailed explanation goes here
U = Uf - params.Sr * y_dot;

end

