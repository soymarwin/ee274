% Function for Exercise 2.C.System3 with impulse response and convolution
% method
% Written by Marwin B. Alejo 2020-20221 for EE274_ProgEx02 28/09/2020

function y = dt_4(x,L)
    a = [1 -0.5 -0.5];
    b = [1];
    h = impz(b,a); %impulse response
    y = conv(h,x);
end