% Function for Exercise 2.C.System3 with impulse response and convolution
% method
% Written by Marwin B. Alejo 2020-20221 for EE274_ProgEx02 28/09/2020

function y = dt_4(x,L)
    y= zeros(1,length(x));
    b = [1];  % input coefficient n
    a = [1 a_l -0.5 -0.5] % output coefficient n, n-L, n-L-1
    h = impz(b,a); %impulse response
    y = conv(h,x);
end