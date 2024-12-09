clear;

% Lip parameters
k = 1278.8;                 % Spring stiffness [N/m]
m = 1.78e-4;                % Mass [kg]
r = 0;                   % Damper [Ns/m]
A = 1e-4;                   % Lip area [m^2]

w0 = sqrt(k/m);             % Undamped eigenfrequency [rad/s]
xi = r / (2*m*w0);          % Damping factor [-]
w0cs = w0 * sqrt(1 - xi^2); % Damped eigenfrequency [rad/s]
tau = 1 / (xi*w0);

% easy upper
y_uppr = 0.8e-3;
k_uppr = 1.06e4;
r_uppr = 1.5;
alpha_uppr = 2;

% harder lower
y_lowr = -0.5e-3;
k_lowr = 6.16e6;     % 6
r_lowr = 1.7;
alpha_lowr = 1.2;


y0 = 1.5e-3;        % Initial displacement [m]
dy0 = 0;            % Initial velocity [m/s]
f0 = 0;             % Initial force [N]



dt = 1e-6/2;
T = .05;
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
    f =  f0-k_uppr * heaviside(y - y_uppr)*(y - y_uppr) ^ alpha_uppr + ...
    k_lowr * heaviside(y_lowr - y)*(y_lowr - y) ^ alpha_lowr;
    
    velo_c = r + ...
    k_uppr * heaviside(y - y_uppr) * (y - y_uppr) ^ alpha_uppr * r_uppr + ...
    k_lowr * heaviside(y_lowr - y) * (y_lowr - y) ^ alpha_lowr * r_lowr;


    ddy = (-k*tild_y - velo_c * tild_dy + f) / (m + gamma*velo_c*dt + k*2*beta*dt^2/2);
    
    % 3. Correction
    y = tild_y + 2*beta*dt^2/2 * ddy;
    dy = tild_dy + gamma*dt * ddy;
    y_new(i) = y;
end
%%

figure;
plot(1e3*t, 1e3*y_new, 'LineStyle', '-', 'LineWidth', 1.2);
set(gca, 'XGrid', 'On', 'YGrid', 'On', 'Box', 'Off');
xlabel('Idő [ms]');
ylabel('Ajak kitérése [mm]');
legend({'Newmark'});