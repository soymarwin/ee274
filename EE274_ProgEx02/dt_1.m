% Function for Exercise 2.C.System1 with recursive method
% Written by Marwin B. Alejo 2020-20221 for EE274_ProgEx02 28/09/2020

function y = dt_1(x)
    % x must be defined in the main function
    y = zeros(1,length(x)); %initialize output signal y
    for n = 1:length(x) %5 time indices
        if n<2
            y(n) = x(n);
        else
            y(n) = 0.5*x(n) + 0.5*x(n-1);
        end
    end
end