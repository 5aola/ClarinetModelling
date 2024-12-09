function [r_t, dt] = calc_reflection(Z_in, Z0, f_lim, f)
%CALC_REFLECTION Summary of this function goes here
%   Detailed explanation goes here
Z_in(1) = 0; 

window = zeros(size(f));
window(f < f_lim(1)) = 1;
window(f > f_lim(2)) = 0;

window(f > f_lim(1) & f < f_lim(2)) = 0.5*(1+cos(pi*(f(f > f_lim(1) & f<f_lim(2))-f_lim(1)) / (f_lim(2)-f_lim(1))));

% preparing for reflection function
Z_in_clear = Z_in .* window + Z0*(1-window);
Z_in_clear_ready_for_ifft = cat(2,Z_in_clear,flip(conj(Z_in_clear(2:(end-1)))));

%figure;
%plot(abs(Z_in_clear_ready_for_ifft));
%figure;
%plot(angle(Z_in_clear_ready_for_ifft)*(180/pi));

% ifft
r_t = ifft((Z_in_clear_ready_for_ifft-Z0)./(Z_in_clear_ready_for_ifft+Z0), 'symmetric');
r_t = r_t(1:floor(end/2));

fs = 2*f(size(r_t, 2));
dt = 1/fs;

end

