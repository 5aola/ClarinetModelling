function [y_new, y_dot_new, y_dotdot_new, nt] = calc_harmonic_osc_with_collision_newmark(y, dy, ddy, p_in, params, p_mouth, dt_old, dynamic_dt)
%CALC_HARMONIC_OSC_WITH_COLLISION_NEWMARK Summary of this function goes here
%   Detailed explanation goes here

beta = 0.25;
gamma = 0.5;
dt = dt_old;

y_approx = [1, 1];       %[new, old]
dy_approx = [1, 1];
ddy_approx = [1, 1];

y_saved = y;
dy_saved = dy;
ddy_saved = ddy;


nt = 1;
dt = dt/nt;
have_to_iterate = 1;
y_error = 1;

while have_to_iterate == 1
    for i = 1 : nt
        % 1. Prediction
        tild_y = y + dt*dy + (1 - 2*beta)*dt^2/2*ddy;
        tild_dy = dy + (1 - gamma)*dt*ddy;
        
        % 2. Solution
        f = params.Sr * (p_in - p_mouth) - ...
        params.k_tg * heaviside(y - params.y_tg)*(y - params.y_tg) ^ params.alpha_tg + ...
        params.k_lay * heaviside(params.y_lay - y)*(params.y_lay - y) ^ params.alpha_lay;
        
        ydot_coeff = params.m * params.gamma + ...
        params.k_tg * heaviside(y - params.y_tg)*(y - params.y_tg) ^ params.alpha_tg * params.r_tg + ...
        params.k_lay * heaviside(params.y_lay - y)*(params.y_lay - y) ^ params.alpha_lay * params.r_lay;
        
        y_coeff = params.k;
        
        ddy = (- y_coeff*tild_y - ydot_coeff*tild_dy + f) / (params.m + gamma*ydot_coeff*dt + y_coeff*2*beta*dt^2/2);
        
        % 3. Correction
        y = tild_y + 2*beta*dt^2/2 * ddy;
        dy = tild_dy + gamma*dt * ddy;
    end

    if ((y > params.y_tg) || (y < params.y_lay)) && dynamic_dt
        y_approx(2) = y_approx(1);
        dy_approx(2) = dy_approx(1);
        ddy_approx(2) = ddy_approx(1);

        y_approx(1) = y;
        dy_approx(1) = dy;
        ddy_approx(1) = ddy;

        y_error = abs(1 - y_approx(1)/y_approx(2));

        if y_error > 0.0001
            dt = dt / 2;
            nt = nt * 2;

            y = y_saved;
            dy = dy_saved;
            ddy = ddy_saved;    
        else
            have_to_iterate = 0;
        end
    else
        have_to_iterate = 0;
    end
end
y_new = y;
y_dot_new = dy;
y_dotdot_new = ddy;
end

