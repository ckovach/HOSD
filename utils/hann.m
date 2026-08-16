function out = hann(N)

% Hann window. Added to remove dependence on the signal processing toolbox.

x = (0:N-1)'/(N-1);

out = 0.5 - 0.5*cos(2*pi*x);