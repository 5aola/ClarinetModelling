function Uf_dot = calc_air_flow(~, U_f, p, params, P_in, ksi)
%CALC_AIR_FLOW Summary of this function goes here
%   Detailed explanation goes here
Uf_dot = P_in - p - sqrt(abs(U_f))^3 * sign(U_f) / (sqrt(params.A)^3 * ksi^2);

rw = params.W / ksi;
% effective mass through the slit:
M_e = (params.rho / (2 * pi * params.W) * (sqrt(rw) + 2 * sqrt(rw) * log(2 * rw))); 

Uf_dot = Uf_dot / M_e;
end

