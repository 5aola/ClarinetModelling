function uf_new = calc_air_flow_for_newmark(uf, p_in, p_mouth, params, y, dt)
%CALC_AIR_FLOW_FOR_NEWMARK Summary of this function goes here
%   Detailed explanation goes here
ksi = params.H + y;

rw = params.W / ksi;
M_e = (params.rho / (2 * pi * params.W) * (sqrt(rw) + 2 * sqrt(rw) * log(2 * rw)));

uf_dot = p_mouth - p_in - sqrt(abs(uf))^3 * sign(uf) / (sqrt(params.A)^3 * ksi^2);
uf_dot = uf_dot / M_e;
uf_new = uf_dot*dt + uf;
end

