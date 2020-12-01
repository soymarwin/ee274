function hd=ideal_lp(wc,M)
%hd: ideal LPF impulse response between 0 and M-1
%wc: cut-off frequencies in radians
%M: length of the filter
    alpha=(M-1)/2;
    n=[0:M-1];
    m=n-alpha;
    fc=wc/pi;
    hd=fc*sinc(fc*m);
end
