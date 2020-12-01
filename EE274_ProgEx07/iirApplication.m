function iirApplication()
    sig_x=ecg(400).';
    sig_y=sgolayfilt(sig_x,0,5);
    [sig_m,sig_n]=size(sig_y);
    sig_fs=200;
    sig_ts=timescope('SampleRate',sig_fs,...
        'TimeSpanSource','Property',...
        'TimeSpan',0.6,...
        'ShowGrid',true,...
        'NumInputPorts',2,...
        'LayoutDimensions',[2 1]);
    sig_ts.ActiveDisplay=1;
    sig_ts.YLimits=[-0.5 1];
    sig_ts.Title = 'Noisy Signal';
    sig_ts.ActiveDisplay=2;
    sig_ts.YLimits=[-0.5 1];
    sig_ts.Title='Filtered Signal';
    % lowpass
    [n,fc]=cheb1ord(40/200,80/200,-20*log10(0.9),-20*log10(0.01));
    [b,a]=cheby1(n,-20*log10(0.9),fc);
    LP=dsp.IIRFilter('Numerator',b,'Denominator',a);
    freqz(LP);
    % highpass
    [b,a]=cheby1(n,-20*log10(0.1),fc, 'high');
    HP=dsp.IIRFilter('Numerator',b,'Denominator',a);
    freqz(HP)
    fprintf('The order of the design IIR filter is %d.\n',n);
    tic;
    while toc < 30
        sig_x=.1*randn(sig_m,sig_n);
        highFreqNoise=HP(sig_x);
        noisySignal=sig_y+highFreqNoise;
        filteredSignal=LP(noisySignal);
        sig_ts(noisySignal,filteredSignal);
    end
    release(sig_ts);
end