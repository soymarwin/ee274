% Function that does  signal folding
% Written by Marwin B. Alejo 2020-20221 for EE274_ProgEx01 19/09/2020
function [y,n]=sigfold(x,n)
    y=fliplr(x); n=-fliplr(n);
end