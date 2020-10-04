% Function for Exercise 2.C.System2 with recursive method
% Written by Marwin B. Alejo 2020-20221 for EE274_ProgEx02 28/09/2020

function y = dt_2(x)
    % x must be defined in the main function
    y = zeros(1,length(x)); %initialize output signal y
    for n = 1:length(x) %5 time indices
         if n==1
            y(n) = x(n);
         elseif n==2
            y(n) = x(n) - 2*y(n-1);
         else
            y(n) = x(n) - 2*y(n-1)- 2*y(n-2);
         end
    end
end