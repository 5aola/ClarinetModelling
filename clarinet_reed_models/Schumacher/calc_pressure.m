function p_n = calc_pressure(p_old, U, Z0, r, dt, n)
%CALC_PRESSURE Summary of this function goes here
%   Detailed explanation goes here


calc_n = min(size(r, 2), n);
p_old = p_old(end-calc_n+1:end);
U = U(end-calc_n+1:end);


conv = fliplr(r(1:calc_n)) .* (p_old + Z0 * U)';
p_n = Z0 * U(end) + sum(conv);

end

