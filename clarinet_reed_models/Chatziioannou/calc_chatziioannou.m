function [p, u, y, f, fs, colldetect] = calc_chatziioannou(note, mouth_pressure, length, opts)
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

     if strcmp(opts, "dynamic_dt")
         dynamic_dt = 1;
     end
 else
     dynamic_dt = 0;
 end

% Parameters
params.mu = 0.0231;             % [kg/m^2]  mass/area
params.Sr = 1.46e-4;            % [m^2]     area               98.6e-6
params.w0 = 23250;              % [1/s]     resonance freq
params.gr = 3000;               % [1/s]     damping
params.A = 0.0797;              % [mks]     current amplitude parameter   
params.H = 0.5e-3;              % [m]       tip opening           0.4
params.W = 0.12;            % [m]       width of the reed (opening) INTERNET DATA -> 13.05e-3
params.theta = 30 * pi/180;     % [rad]     reed tipping in radian

% for collisions in portato mode
params.m = 8.62e-6;             % [kg]          effective mass
params.gamma = 8190;            % [1/s]         effective damping 
params.y_lay = -0.3e-3;        % [m]           contact position between reed-mouthpiece -0.385e-3
params.y_tg = 0.1e-3;          % [m]{max=0.63}contact position between tongue-reed
params.y_eq = -0.8e-3;
params.k = 1308;                % [N/m]         effective stiffness
params.k_lay = 6.16e6;          % [N/m^2]       effective stiffness 6
params.k_tg = 1.06e6;           % [N/m^alpha_tg]effective stiffness
params.r_lay = 1.7;             % [s/m]         damping coefficient of mouthpiece
params.r_tg = 1.5;              % [s/m]         damping coefficient of tongue
params.alpha_lay = 2;
params.alpha_tg = 1.3;


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

%dt = 1/(2*fs);
t_length = length;   % [s] length of calculation
n_sample = floor(t_length / dt);

% initialazing parameters
p_n = zeros(n_sample,1);            % guage pressure
U = zeros(n_sample,1);              % acoustic current 
Uf = 0;
y = zeros(n_sample+1,1);  
%y(1:2) = -0.2e-3;                   % reed displacement starting position
y_dot = 0;
y_dotdot = 0;

colldetect = zeros(n_sample+1,1); 

for n = 2:n_sample
    if size(mouth_pressure,2) == 3
        tau = mouth_pressure(3)/5;
        p_mouth = (mouth_pressure(2) - mouth_pressure(1)) * (1 - exp( - n*dt/tau)) + mouth_pressure(1);
    else
        p_mouth = mouth_pressure; % [Pa] static pressure in the player mouth 
    end

    % calculating pressure eq5
    p_n(n) = calc_pressure(p_n(1:n-1), U(1:n-1), Z0, r_t, dt, n-1);
 
    

    [y(n+1), y_dot, y_dotdot, colldetect(n)] = calc_harmonic_osc_with_collision_newmark(y(n), y_dot, y_dotdot, p_n(n), params, p_mouth, dt, dynamic_dt);


    % calculating airflow eq7
    
    %Uf = calc_bernoulli_flow(Uf, p_n(n), p_mouth, params, y(n+1), dt);
    Uf = calc_air_flow_for_newmark(Uf, p_n(n), p_mouth, params, y(n+1), dt);

    % calculating total current eq9
    U(n) = Uf - params.Sr * y_dot;          % ---

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

