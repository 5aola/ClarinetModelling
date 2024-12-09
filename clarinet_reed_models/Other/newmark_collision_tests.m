clear;

p_mouth = 5000;
p_in = 0;
params.Sr = 1.46e-4;
% for collisions in portato mode
params.m = 8.62e-6;             % [kg]          effective mass
params.gamma = 8190;            % [1/s]         effective damping 
params.y_lay = -0.385e-3;       % [m]           contact position between reed-mouthpiece -0.385e-3
params.y_tg = 0.2e-3;           % [m]{max=0.63}contact position between tongue-reed
params.k = 1308;                % [N/m]         effective stiffness
params.k_lay = 6.16e6;          % [N/m^2]       effective stiffness
params.k_tg = 1.06e4;           % [N/m^alpha_tg]effective stiffness
params.r_lay = 1.7;             % [s/m]         damping coefficient of mouthpiece
params.r_tg = 1.5;              % [s/m]         damping coefficient of tongue
params.alpha_tg = 1.2;
params.alpha_lay = 2;


y0 = 1.5e-3;                  % Initial displacement [m]
dy0 = 0;                     % Initial velocity [m/s]

dt = 1e-6;
T = .02;
t = (0 : dt : T).';
nt = length(t);

y = y0;
dy = dy0;
ddy = 0;

beta = 0.25;
gamma = 0.5;


y_new = zeros(nt,1);

for i = 1 : nt   
    % 1. Prediction
    tild_y = y + dt*dy + (1 - 2*beta)*dt^2/2*ddy;
    tild_dy = dy + (1 - gamma)*dt*ddy;
    
    % 2. Solution
    f = params.Sr * (p_in - p_mouth) - ...
    params.k_tg * heaviside(y - params.y_tg) ^ params.alpha_tg + ...
    params.k_lay * heaviside(params.y_lay - y) ^ params.alpha_lay;
    
    ydot_coeff = params.m * params.gamma + ...
    0 * params.k_tg * heaviside(y - params.y_tg) ^ params.alpha_tg * params.r_tg + ...
    0 * params.k_lay * heaviside(params.y_lay - y) ^ params.alpha_lay * params.r_lay;
    
    y_coeff = params.k;
    
    ddy = (- y_coeff*tild_y - ydot_coeff*tild_dy + f) / (params.m + gamma*ydot_coeff*dt + y_coeff*2*beta*dt^2/2);
    
    % 3. Correction
    y = tild_y + 2*beta*dt^2/2 * ddy;
    dy = tild_dy + gamma*dt * ddy;
    y_new(i) = y;

end

figure;
plot(y_new);