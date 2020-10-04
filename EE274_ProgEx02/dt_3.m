% Function for Exercise 2.C.System3 with impulse response and convolution
% method
% Written by Marwin B. Alejo 2020-20221 for EE274_ProgEx02 28/09/2020

function y = dt_3(x)
    a = [1 2 2];
    b = [1.5 0.5];
    h = impz(b,a); %impulse response
    y = conv(h,x);
    figure();
end
