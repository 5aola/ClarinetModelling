function y_dot = calc_harmonic_osc_with_collision(~, y, p_in, params, p_mouth)
%CALC_HARMONIC_OSC_WITH_COLLISION Summary of this function goes here
%   Detailed explanation goes here

f = params.Sr * (p_in - p_mouth) - ...
    params.k_tg * heaviside(y(2) - params.y_tg)*(y(2) - params.y_tg) ^ params.alpha_tg + ...
    params.k_lay * heaviside(params.y_lay - y(2))*(params.y_lay - y(2)) ^ params.alpha_lay;

ydot_coeff = params.m * params.gamma + ...
    params.k_tg * heaviside(y(2) - params.y_tg)*(y(2) - params.y_tg) ^ params.alpha_tg * params.r_tg + ...
    params.k_lay * heaviside(params.y_lay - y(2))*(params.y_lay - y(2)) ^ params.alpha_lay * params.r_lay;

y_coeff = params.k;

y_dot = [-(ydot_coeff)/params.m    -(y_coeff)/params.m;
         1                                           0] * y + ...
        [f/params.m; 0];

end

