clear

%% Parameters
params.mu = 0.0231;    % [kg/m^2]  mass/area
params.Sr = 1.46e-4;   % [m^2]     area
params.w0 = 23250;     % [1/s]     resonance freq
params.gr = 3000;      % [1/s]     damping
params.A = 0.0797;     % [mks]     current amplitude parameter   
params.H = 0.4e-3;     % [m]       tip opening
params.W = 13.05e-3;   % [m]       width of the reed (opening) !INTERNET DATA!
params.theta = 30 * pi/180;     % [rad]     reed tipping in radian

params.P0 = 1670;      % [Pa]      static pressure in the player mouth 

T = 27.5;
[params.rho, ~] = gas_properties(T, 'air');

fs = 40000;

%k = params.w0^2 * params.Sr*params.mu;

%% calculating fs

%frequency scale !!!!! EXPORTED THE IMPEDANCES FOR THIS SCALE,  NOT IN USE
%n_freq = 2^16;
%f_min = 0;
%f_max = fs/2;
%f = linspace(f_min, f_max, n_freq);

%% calculating bore input impedance
note = "A  4";
%[Z_in, Z0, params.w0_bore] = calc_simple_impedance_example(f);

[Z_in, Z0, params.w0_bore, f] = get_clarinet_impedance(note);
%figure 
%plot(f,20*log10(abs(Z_in/Z0)));

%% calculating reflection
params.f_lim = [4000 6000];
[r_t, dt] = calc_reflection(Z_in, Z0, params.f_lim, f);

%% reflection test

reflection_arrival_time = 6*1e-3;

n_reflection_arrival_sample = floor(reflection_arrival_time / dt);
reflection_function_test = zeros(1,n_reflection_arrival_sample*3);
reflection_function_test(1,n_reflection_arrival_sample) = -0.99;

%r_t = reflection_function_test;


%% plotting reflection
time = (0:(size(r_t, 2)-1))*dt;
%figure;
%plot(time, r_t);
%xlabel('Time [s]');
%ylabel('Amp');

%% calculating Schumacher's reed model

t_length = 0.5;   % [s] length of iteration
n_sample = floor(t_length / dt);


% initialazing parameters
p_n = zeros(n_sample,1);            % guage pressure
U = zeros(n_sample,1);              % acoustic current 
Uf = 0;
y = zeros(n_sample+1,1);            % reed displacement starting position
y_dot = 0;


for n = 2:n_sample

    % calculating pressure eq5
    p_n_new = calc_pressure(p_n(1:n-1), U(1:n-1), Z0, r_t, dt, n-1);
    p_n(n) = p_n_new;
    
    % calculating harm osc eq6
    [~, x_new] = ode45(@(t, x)calc_harmonic_osc(t, x, p_n(n), U(n-1), params), [0 dt], [y_dot; y(n)]);
    y_dot = x_new(end,1);
    y(n+1) = x_new(end,2);
    
    % calculating airflow eq7
    ksi = params.H + y(n+1);
    [~, Uf_new] = ode45(@(t, U_f)calc_air_flow(t, U_f, p_n(n), params, ksi), [0 dt], Uf);
    if (ksi < 0)
        Uf = 0;
    else
        Uf = Uf_new(end);
    end
    
    
    % calculating total current eq9
    U_new = calc_total_current(Uf, y_dot, params);
    U(n) = U_new;

    if (mod(n,1000) == 0)
        disp(y(n+1));
    end
end

%% plotting

figure;
plot(p_n);
soundsc(p_n,fs);


