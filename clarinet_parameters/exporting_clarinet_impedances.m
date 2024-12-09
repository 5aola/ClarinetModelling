% exporting_clarinet_impedances.m
% This script demonstrates the assembly of the input impedance of a
% Clarinet Bb Yamaha CX1.
clear;

%% Getting the geometry and configurations
% Temperature where measurements were made
T = 27.5;

% sample rate and frequency scale
fs = 96000;
n_freq = 2^16;
f_min = 0;
f_max = fs/2;
f = linspace(f_min, f_max, n_freq);

% Clarinet dimensions
[bores_holes, bores, holes] = get_clarinet_geometry();
N = size(bores_holes,1);

[configs, notes] = get_clarinet_configurations();
n_config = size(configs, 1);

%% Build model
model = cell(N, 1);

for i = 1 : N   
    if bores_holes(i,1) == 0
        % Add bore: length, start radius, end radius, temp
        model{i} = {'dn', 'bore', bores_holes(i,4), [bores_holes(i,2) bores_holes(i,3)], T};
    else
        % Add hole: type, bore radius, temp, hole radius, wall width, 
        model{i} = {'dn', 'hole', '', bores_holes(i,4), T, bores_holes(i,2), bores_holes(i,3)};
    end
end

% Add the radiation impedance at the open end
model{N + 1} = {'dn', 'radiation', bores(end,2), T, 'flanged'};

%% Calculate plane wave impedance

[rho, c] = gas_properties(T, 'air');
Z0 = (rho * c) / (bores(1,1)^2 * pi);

%% Go through each config to export them
Z_all = zeros(n_config,size(f,2));

for i_config = 1 : n_config
    hole_state = configs{i_config, 2};
    hole_cntr = 1;
    for i = 1 : N   
        if bores_holes(i,1) == 1
            switch hole_state(hole_cntr)
                case 'X'
                    model{i}{3} = 'closed';
                case 'O'
                    model{i}{3} = 'open';
            end
           hole_cntr = hole_cntr + 1;
        end
    end
    % Calculate pipe impedance
    Z_p = bore_input_impedance(model, f);
    Z_all(i_config,:) = Z_p;
    
    disp(i_config)
end

%% save the parameters

save("C:\Users\kosty\Documents\MATLAB\clarinet_project\clarinet_parameters\clarinet_input_inpedances.mat", "Z0", "f", "Z_all", "fs");

%% testing the export

%load("C:\Users\kosty\Documents\MATLAB\clarinet_project\clarinet_parameters\clarinet_input_inpedances.mat", "Z0", "f", "Z_all", "fs"); 

%figure;
%plot(f, real(Z_all(i_config-1,:)));

