function x_dot = calc_harmonic_osc(~, x, p, U, params, P_in)
%CALC_HARMONIC_OSC Summary of this function goes here
%   Detailed explanation goes here

f = (p - P_in) / params.mu - 0.5 * params.rho / params.W * ...
    U^2/((x(2) + params.H)*tan(params.theta)*params.mu*params.Sr);

x_dot = [-params.gr    -(params.w0^2);
         1              0] * x + ...
        [f; 0];
end

