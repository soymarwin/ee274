% Written and performed by Marwin B. Alejo 2020-20221
% Return Impulse Response, Hanning Window, and Magnitude Response plots
% by simply providing the cut-off frequency in pi*rad and transition width.
function hannw(wc,tr_width)
    M=ceil(6.2*pi/tr_width)+1;
    n=[0:M-1];
    alpha=(M-1)/2;
    m=n-alpha;
    fc=wc/pi;
    hd=fc*sinc(fc*m);
    w_hann=(hann(M))';
    b=hd.*w_hann;
    [H,w] = freqz(b,[1],1000,'whole') ;
    H = (H(1:1:501));
    w = (w(1:1:501));
    mag = abs(H);
    db = 20*log10((mag+eps)/max(mag));
    %wvtool(b); % for sidelobe measurement
    figure();
    subplot(3,1,1);stem(n,hd);title('Impulse Response');xlabel('n');ylabel('hd(n)');
    subplot(3,1,2);stem(n,w_hann);title('Hanning Window');xlabel('n');ylabel('w(n)');
    subplot(3,1,3);plot(w/pi,db);title('Magnitude Response in dB');grid;xlabel('frequency in \pi units');ylabel('decibels');
end