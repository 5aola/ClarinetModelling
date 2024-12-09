function [y_new, y_dot_new, y_dotdot_new] = calc_harmonic_osc_newmark(y, dy, ddy, p, U, params, p_mouth, dt_old)
%CALC_HARMONIC_OSC_NEWMARK Summary of this function goes here
%   Detailed explanation goes here

beta = 0.25;
gamma = 0.5;

nt = 2;
dt = dt_old / nt;

for i = 1 : nt

    % 1. Prediction
    tild_y = y + dt*dy + (1 - 2*beta)*dt^2/2*ddy;
    tild_dy = dy + (1 - gamma)*dt*ddy;
    
    % 2. Solution
    f = (p - p_mouth) / params.mu - 0.5 * params.rho / params.W * ...
    U^2/((y + params.H)*tan(params.theta)*params.mu*params.Sr);
    
    ddy = (-params.w0^2*tild_y - params.gr*tild_dy + f) / (gamma*params.gr*dt + params.w0^2*2*beta*dt^2/2);
    
    % 3. Correction
    y = tild_y + 2*beta*dt^2/2 * ddy;
    dy = tild_dy + gamma*dt * ddy;
end

y_new = y;
y_dot_new = dy;
y_dotdot_new = ddy;

end

