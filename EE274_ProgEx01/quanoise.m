% Function to generate quantization noise and sqnr for section F. Audio File Formats
% Written by Marwin B. Alejo 2020-20221 for EE274_ProgEx01 19/09/2020
function [y,q,sqnr]=quanoise(y_af,bit)
    y=double(uencode(y_af,bit));
    y=y/max(y);
    y=2*(y-0.5);
    q=y_af-y; %quantization noise
    sqnr=10*log10(norm(y_af)/norm(q));
end