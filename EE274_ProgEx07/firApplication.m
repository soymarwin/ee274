function firApplication()
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
    lp_fp=40; lp_fs=80; lp_ap=0.05; lp_sp=0.0009;
    lp_F=[0 lp_fp lp_fs sig_fs/2]/(sig_fs/2);
    lp_A=[1 1 0 0];
    lp_D=[lp_ap lp_sp];
    lp_b=firgr('minorder',lp_F,lp_A,lp_D);
    LP=dsp.FIRFilter('Numerator',lp_b);
    lpord=filtord(lp_b,1);

    % highpass
    hp_fs=40; hp_fp=80; hp_sp=0.05; hp_ap=0.0009;
    hp_F=[0 hp_fs hp_fp sig_fs/2]/(sig_fs/2);
    hp_A=[0 0 1 1];
    hp_D=[hp_sp hp_ap];
    hp_b=firgr('minord',hp_F,hp_A,hp_D);
    HP=dsp.FIRFilter('Numerator',hp_b);
    hpord=filtord(lp_b,1);

    fprintf('The order of the design lowpass and highpass FIR filters are %d and %d.\n',lpord,hpord);

    tic;
    while toc < 50
        sig_x=.1*randn(sig_m,sig_n);
        highFreqNoise=HP(sig_x);
        noisySignal=sig_y+highFreqNoise;
        filteredSignal=LP(noisySignal);
        sig_ts(noisySignal,filteredSignal);
    end

    release(sig_ts);
end
