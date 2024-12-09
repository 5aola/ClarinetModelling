function [p, u, y, f, fs, p_mouth_in_time] = calc_schumacher(note, mouth_pressure, length, opts)
%CALC_SCHUMACHER Summary of this function goes here
%   Detailed explanation goes here

% Check input
narginchk(3, 4);

% Check model
if ~isstring(note)
    error('calc_schumacher:invalid_argument', ...
        'Desired NOTE must be a string');
end

% Check frequency
if ~isnumeric(mouth_pressure)
    error('calc_schumacher:invalid_argument', ...
        'Desired static_mouth_pressure must be a number');
end

% Check frequency
if ~isnumeric(length)
    error('calc_schumacher:invalid_argument', ...
        'Desired length must be a number');
end

% Process calculation options
if nargin == 4
     if ~isstring(opts)
        error('calc_schumacher:invalid_argument', ...
            'If you want to include the reed closing state calculations too then it must be a string');
    end
end

% Parameters
params.mu = 0.0231;             % [kg/m^2]  mass/area
params.Sr = 98.6e-6;            % [m^2]     area               98.6e-6      1.46e-4
params.w0 = 23250;              % [1/s]     resonance freq
params.gr = 3000;               % [1/s]     damping
params.A = 0.0797;              % [mks]     current amplitude parameter   
params.H = 0.4e-3;              % [m]       tip opening
params.W = 0.12;                % [m]       width of the reed (opening)
params.theta = 30 * pi/180;     % [rad]     reed tipping in radian
T = 27.5;
[params.rho, ~] = gas_properties(T, 'air');

[Z_in, Z0, f, fs] = get_clarinet_impedance(note);

params.f_lim = [4000 6000]; % 4-6
[r_t, dt] = calc_reflection(Z_in, Z0, params.f_lim, f);

if (sum(r_t) > -0.9 || sum(r_t) < -1.1)
    error('calc_schumacher:error_in_reflection', ...
                'Total reflection should be (-1) but it is %4.2f ', sum(r_t));
end

fprintf("\nStarting calculations\n");

t_length = length;   % [s] length of calculation
n_sample = floor(t_length / dt);

% initialazing parameters
p_n = zeros(n_sample,1);            % guage pressure
U = zeros(n_sample,1);              % acoustic current 
Uf = 0;
y = zeros(n_sample+1,1);            % reed displacement starting position
y_dot = 0;
y_dotdot = 0;

p_mouth_in_time = zeros(n_sample,1);

for n = 2:n_sample
    if size(mouth_pressure,2) == 3
        tau = mouth_pressure(3)/5;
        p_mouth = (mouth_pressure(2) - mouth_pressure(1)) * (1 - exp( - n*dt/tau)) + mouth_pressure(1);   
    else
        p_mouth = mouth_pressure; % [Pa] static pressure in the player mouth 
    end
    p_mouth_in_time(n) = p_mouth;

    % calculating pressure eq5
    p_n(n) = calc_pressure(p_n(1:n-1), U(1:n-1), Z0, r_t, dt, n-1);

    if nargin == 4
        [y(n+1), y_dot, y_dotdot] = calc_harmonic_osc_newmark(y(n), y_dot, y_dotdot, p_n(n), U(n-1), params, p_mouth, dt);
    else
        % calculating harm osc eq6
        [~, x_new] = ode45(@(t, x)calc_harmonic_osc(t, x, p_n(n), U(n-1), params, p_mouth), [0 dt], [y_dot; y(n)]);
        y_dot = x_new(end,1);
        y(n+1) = x_new(end,2);
    end

    % calculating airflow eq7
    Uf = calc_air_flow_for_newmark(Uf, p_n(n), p_mouth, params, y(n+1), dt);
    %ksi = params.H + y(n+1);
    %[~, Uf_new] = ode45(@(t, U_f)calc_air_flow(t, U_f, p_n(n), params, p_mouth, ksi), [0 dt], Uf);
    %Uf = Uf_new(end);

    % calculating total current eq9
    U(n) = calc_total_current(Uf, y_dot, params);

    if (isnan(y(n+1)) == true)
        fprintf("\nThe reed has exceed the maximum deviation. The pressure is undefinable.\n");
        fprintf("This happened at: %.4f seconds at n = %d \n",dt*n, n);
        break
    end

    if (mod(n,floor(0.1/dt)) == 0)
        fprintf("Calculated time: %4.2f seconds, y = %e \n",dt*n, y(n+1));

    end
end

fprintf("Finished calculations.\n");

p = p_n;
u = U;

end

