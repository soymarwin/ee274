% Function to generate unit-impulse sequence
% Written by Marwin B. Alejo 2020-20221 for EE274_ProgEx01 19/09/2020
function [x,n]=impseq(n0,n1,n2)
    n = [n1:n2];
    x = [(n-n0) == 0];
end

