% Function to generate unit-step sequence
% Written by Marwin B. Alejo 2020-20221 for EE274_ProgEx01 19/09/2020
function [x,n]=stepseq(n0,n1,n2)
    n = [n1:n2];
    x = [(n0-n) < 0]; %change to less than to satisfy u(n-1)condition.
end