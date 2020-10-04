% Function to generate sine signal
% Written by Marwin B. Alejo 2020-20221 for EE274_ProgEx01 19/09/2020
function [x, t] = sin_uni(fs, f0, T)
t = 0:1/fs:T;
x = sin(2*pi*f0*t);
end

