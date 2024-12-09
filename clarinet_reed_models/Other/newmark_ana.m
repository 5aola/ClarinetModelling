clear;

% Lip parameters
k = 1278.8;                 % Spring stiffness [N/m]
m = 1.78e-4;                % Mass [kg]
r = 9e-2;                   % Damper [Ns/m]
A = 1e-4;                   % Lip area [m^2]

w0 = sqrt(k/m);             % Undamped eigenfrequency [rad/s]
xi = r / (2*m*w0);          % Damping factor [-]
w0cs = w0 * sqrt(1 - xi^2); % Damped eigenfrequency [rad/s]
tau = 1 / (xi*w0);

y0 = 1.5e-3;                  % Initial displacement [m]
dy0 = 0;                     % Initial velocity [m/s]

dt = 1e-5;
T = .02;
t = (0 : dt : T).';
nt = length(t);

y = y0;
dy = dy0;
ddy = 0;

beta = 0.25;
gamma = 0.5;

y_ana = exp(-t/tau).*(y0.*cos(w0cs*t) + (dy0 + y0 / tau)/w0cs*sin(w0cs*t));
y_new = zeros(nt,1);

for i = 1 : nt
   
    % 1. Prediction
    tild_y = y + dt*dy + (1 - 2*beta)*dt^2/2*ddy;
    tild_dy = dy + (1 - gamma)*dt*ddy;
    
    % 2. Solution
    ddy = (-k*tild_y - r*tild_dy) / (m + gamma*r*dt + k*2*beta*dt^2/2); %+ testrehato ero
    
    % 3. Correction
    y = tild_y + 2*beta*dt^2/2 * ddy;
    dy = tild_dy + gamma*dt * ddy;
    y_new(i) = y;
end
err = norm(y_ana-y_new)/norm(y_ana);
fprintf('Relative log10 error: %g\n', log10(err));
%%
fig = figure;
formatfig(fig, [9 6], [1.5 1.5 0.5 0.5]);
hold on;
plot(1e3*t, 1e3*y_new, 'LineStyle', '-', 'LineWidth', 1.2);
plot(1e3*t, 1e3*y_ana, 'LineStyle', '-.', 'LineWidth', 1.2);
set(gca, 'XGrid', 'On', 'YGrid', 'On', 'Box', 'Off');
xlabel('Idõ [ms]');
ylabel('Ajak kitérése [mm]');
legend({'Newmark-séma', 'Analitikus megoldás'});
print(fig, '-dpng', '-r600', 'newmark_analitikus.png');