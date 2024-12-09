function uf_new = calc_bernoulli_flow(uf, p_in, p_mouth, params, y, dt)
%CALC_AIR_FLOW_FOR_NEWMARK Summary of this function goes here
%   Detailed explanation goes here

uf_new = sign(p_mouth - p_in) * params.W  * heaviside(params.H + y)*(params.H + y)*sqrt((2*abs(p_mouth - p_in))/params.rho);
end

