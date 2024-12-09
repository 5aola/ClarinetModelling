clear

note = "B_4";
mouth_pressure = 1900; % [0 2000 0.1] 0->2000 Pa 0.1 s to get to 5tau
length = 2;
[p, u, y_dynamic, ~, fs, colldetect] = calc_chatziioannou(note, mouth_pressure, length, "dynamic_dt");
%%
figure;
plot((1:size(p,1))*1/fs,p);
audiowrite("Chatziioannou.wav",p/2000,fs)
%%
[p, u, y_dynamic, ~, fs, colldetect] = calc_chatziioannou(note, mouth_pressure, length, "dynamic_dt");
[p, u, y, f, fs, ~] = calc_chatziioannou(note, mouth_pressure, length);
%%
y_div = y_dynamic./y;
figure;
%subplot(3, 1, 1);
plot((1:size(y,1))*1/fs,y_dynamic);
hold on;
plot((1:size(y,1))*1/fs,y);
legend({'Dynamic','Static'},'Location','northeast')
title('Reed Displacement')
xlabel('Time [s]')
ylabel('y [m]')
%subplot(3, 1, 2);
%plot(y_div-1)


%subplot(3, 1, 3);

%plot(colldetect);

soundsc(p,fs);

%%
[configs, notes] = get_clarinet_configurations();
length = 0.5;
pressure_interval = [1500 3000 11];
note_interval = 34:34;


pressure_limits = zeros(size(configs,1),2);

for i = note_interval
    [max_pressures, pressure_tests] = test_decay_chat(convertCharsToStrings(configs{i}), pressure_interval, length);

    select_pressures = max_pressures(:,1)./max_pressures(:,2) < 1.1;
    valid_pressures = pressure_tests(select_pressures);
    if isempty(valid_pressures)
        pressure_limits(i,1) = 0;
        pressure_limits(i,2) = 0;
    else
        pressure_limits(i,1) = min(valid_pressures);
        pressure_limits(i,2) = max(valid_pressures);
    end

    disp(i)
end

%%

figure;
bar(pressure_tests, max_pressures);
title('Input pressures tests @ note C# 6')
xlabel('Input Pressure [Pa]')
ylabel('Peak output pressures @ 2 points [Pa]')
X = categorical({configs{:, 1}});
X = reordercats(X,{configs{:, 1}});

figure;
bar(X, pressure_limits);
title('Input pressures tests')
ylabel('min-max input presurres [Pa]')


%% impedance vs. FFT
[configs, notes] = get_clarinet_configurations();
%mouth_pressures = sum(pressure_limits,2)/2;

n_config = size(configs,1);
deviation = zeros(n_config,1);

starter = 1;
for i_config = starter : (starter + n_config) %n_config

    [p, ~, ~, f, fs] = calc_chatziioannou(convertCharsToStrings(configs{i_config}), 1800, 0.5);
    [Z_p, Z0, ~, ~] = get_clarinet_impedance(convertCharsToStrings(configs{i_config}));

    %if i_config == 34 || i_config == 39 || i_config == 40
    %deviation(i_config) = 0;

    %else
    p_trimed = p(end-2^14-1:end);
    window = hann(size(p_trimed,1));
    windowed_p = p_trimed .* window;
    p_f = fft(windowed_p);  
    f_fft = (0:(size(p_f,1)-1))/size(p_f,1)*fs;
    Z = 20*log10(abs(Z_p/Z0));
    [localpeaks_p, localpeak_locs_p] = findpeaks(abs(p_f));
    [localpeaks_z, localpeak_locs_z] = findpeaks(Z);
        
    freq_of_p = localpeak_locs_p(max(localpeaks_p)==localpeaks_p);
    deviation(i_config) = f_fft(freq_of_p(1)) / f(localpeak_locs_z(1));
    %end

    %dev_in_cents = 1200*log2(1+min_dev);
    
    fprintf('%d \n', i_config);  
end

%% Plotting the deviation

figure;

X = categorical({configs{:, 1}});
X = reordercats(X,{configs{:, 1}});

dev_in_cents = 1200*log2(deviation);

bar(X,dev_in_cents);
title('Alapharmónikusok eltérése a furat bemeneti impedanciájához képest')
%set(gca,'YScale','log')
ylabel('Deviation [cents]');

%% plot

[p, u, y, f, fs] = calc_chatziioannou("A#_4", 1700, 0.5);
[Z_p, Z0, ~, ~] = get_clarinet_impedance("A#_4");

p_trimed = p(end-2^14-1:end);
window = hann(size(p_trimed,1));
windowed_p = p_trimed .* window;
p_f = fft(windowed_p);  
f_fft = (0:(size(p_f,1)-1))/size(p_f,1)*fs;

figure;
plot(windowed_p);


fig = figure('Units', 'normalized', 'Position', [0 0 1 1]);
subplot(2, 1, 1);
loglog(f_fft, abs(p_f));
xlim([10 10000]);

subplot(2, 1, 2);
Z = 20*log10(abs(Z_p/Z0));
semilogx(f,Z);
xlim([10 10000]);
