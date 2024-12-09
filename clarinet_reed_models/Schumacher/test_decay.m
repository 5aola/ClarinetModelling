function [max_pressures, pressure_tests] = test_decay(note, pressure_interval, length)
%TEST_DECAY Summary of this function goes here
%   Detailed explanation goes here

load("C:\Users\kosty\Documents\MATLAB\clarinet_project\clarinet_parameters\clarinet_input_inpedances.mat", "fs"); 
[configs, notes] = get_clarinet_configurations();
x = strcmp(configs, note);
note_freq = notes(x);
samples = floor(fs/note_freq)*4;

pressure_tests = linspace(pressure_interval(1),pressure_interval(2),pressure_interval(3));
max_pressures = zeros(size(pressure_tests,2),2);

for i = 1:size(pressure_tests,2)

    disp(pressure_tests(i));
    [p, ~, ~, ~, ~] = calc_schumacher(note, pressure_tests(i), length, "newmark");
    
    p_mid = p(floor(end/2)-samples:floor(end/2));
    max_pressures(i,1) = max(p_mid) - min(p_mid);

    p_end = p(end-samples:end);
    max_pressures(i,2) = max(p_end) - min(p_end);   

end

end