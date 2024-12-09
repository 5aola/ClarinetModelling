function [Z_p, Z0, w_01] = calc_simple_impedance_example(f)
%CALC_SIMPLE_IMPEDANCE_EXAMPLE Summary of this function goes here
%   Detailed explanation goes here

% Temperature where measurements were made
T = 27.5;

bore_length = 0.5;       % [m]
bore_radius = 10 * 1e-3;

% Build model
model = cell(1, 1);
% Add bore: length, start radius, end radius, temp
model{1} = {'dn', 'bore', bore_length, ...
    [bore_radius bore_radius], T};

% Add the radiation impedance at the open end
model{2} = {'dn', 'radiation', bore_radius, T, 'piston'};


% Calculate plane wave impedance
[rho, c] = gas_properties(T, 'air');
Z0 = (rho * c) / (bore_radius^2 * pi);

% Calculating (0,1) mode frequency
w_01 = 1.8412 * c / bore_radius; % (0,1) mode frequency
w_01 = w_01/(2*pi);

% Calculating pipe impedance
Z_p = bore_input_impedance(model, f);
% Calculate input impedance of the recorder
%Z_in = 20*log10(abs(Z_p/Z0))

