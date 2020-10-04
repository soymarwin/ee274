% Function to solve section E. Quantization
% Written by Marwin B. Alejo 2020-20221 for EE274_ProgEx01 19/09/2020
function y = adc_uni(x, R, B)
        level = [0:R/(2^B):R-R/(2^B)];
        temp = [-Inf,(level(2:end)-R/(2^(B+1))),Inf];
        y = zeros(1,length(x));
    for i = 1:length(level)
        y = y + (x >= temp(i)).*(x < temp(i+1)).*level(i);
end
