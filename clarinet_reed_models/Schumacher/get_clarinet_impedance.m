function [Z_p, Z0, f, fs] = get_clarinet_impedance(note)
%CALC_BORE_IMPEDANCE_EXAMPLE Summary of this function goes here
%   Detailed explanation goes here

% Getting configurations
[configs, ~] = get_clarinet_configurations();

load("C:\Users\kosty\Documents\MATLAB\clarinet_project\clarinet_parameters\clarinet_input_inpedances_48k.mat", "Z0", "f", "Z_all", "fs"); 

x = strcmp(configs, note);
Z_p = Z_all(x,:);
    
end

