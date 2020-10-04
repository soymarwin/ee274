% Function that adds two discrete time sequences
% Written by Marwin B. Alejo 2020-20221 for EE274_ProgEx01 19/09/2020
function [y,n]=sigadd(x1,n1,x2,n2)
    n = min(min(n1),min(n2)):max(max(n1),max(n2)); % duration of y(n)
    y1 = zeros(1,length(n)); y2 = y1; % initialization
    y1(find((n>=min(n1))&(n<=max(n1))==1))=x1; % x1 with duration of y
    y2(find((n>=min(n2))&(n<=max(n2))==1))=x2; % x2 with duration of y
    y = y1+y2; % sequence addition
end