%% Marwin B. Alejo   2020-20221   EE274_ProgEx06
% Also accessible through
% <http://www.github.com/soymarwin/ee274/EE274_ProgEx06> for history
% tracking.

%% A. Pole-Zero Placement
% This section had been accomplished using filter designer.
% See EE274_ProgEx06.pdf for the complete documentation of the 

%% B. Properties of Windowing Functions

rectw(0.5*pi,0.06*pi); % rectangular window M=31
rectw(0.5*pi,0.03*pi); % rectangular window M=61
rectw(0.5*pi,0.015*pi); % rectangular window M=121
bartw(0.5*pi,0.203*pi); % bartlett window M=31
bartw(0.5*pi,0.102*pi); % bartlett window M=61
bartw(0.5*pi,0.051*pi); % bartlett window M=121
hannw(0.5*pi,0.207*pi); % hanning window M=31
hannw(0.5*pi,0.103*pi); % hanning window M=61
hannw(0.5*pi,0.052*pi); % hanning window M=121
hammw(0.5*pi,0.22*pi); % hamming window M=31
hammw(0.5*pi,0.11*pi); % hamming window M=61
hammw(0.5*pi,0.055*pi); % hamming window M=121
blkw(0.5*pi,0.3667*pi); % blackman window M=31
blkw(0.5*pi,0.183*pi); % blackman window M=61
blkw(0.5*pi,0.0917*pi); % blackman window M=121
kaisw(0.5*pi,0.06*pi,1); % kaiser window M=31 beta 1
kaisw(0.5*pi,0.03*pi,1); % kaiser window M=61 beta 1
kaisw(0.5*pi,0.015*pi,1); % kaiser window M=121 beta 1
kaisw(0.5*pi,0.06*pi,4); % kaiser window M=31 beta 4
kaisw(0.5*pi,0.03*pi,4); % kaiser window M=61 beta 4
kaisw(0.5*pi,0.015*pi,4); % kaiser window M=121 beta 4
kaisw(0.5*pi,0.06*pi,9); % kaiser window M=31 beta 9
kaisw(0.5*pi,0.03*pi,9); % kaiser window M=61 beta 9
kaisw(0.5*pi,0.015*pi,9); % kaiser window M=121 beta 9

%% C. Design a digital FIR lowpass filter
C_wp=0.2*pi; C_ws=0.3*pi; C_tr_width=C_ws-C_wp;
C_M=ceil(6.6*pi/C_tr_width)+1;
C_n=[0:C_M-1];
C_wc=(C_ws+C_wp)/2; %ideal cutoff frequency
C_hd=ideal_lp(C_wc,C_M);
C_w_hamming=(hamming(C_M))';
C_h=C_hd.*C_w_hamming;
[C_db,C_mag,C_pha,C_grd,C_w]=freqz_m(C_h,[1]);
figure();
subplot(3,1,1);stem(C_n,C_h); title('Impulse Response');xlabel('n');ylabel('hd(n)');
subplot(3,1,2);stem(C_n,C_w_hamming);title('Hamming Window');xlabel('n');ylabel('w(n)');
subplot(3,1,3);plot(C_w/pi,C_db);title('Magnitude Response in dB');grid;ylabel('decibel');xlabel('frequency in \pi units');

%% D. Frequency Sampling for grads
D_f = [0 0.125 0.1875 0.25 0.3125 0.375 0.4375 0.5 0.5 0.5625 0.625 0.6875 0.75 0.8125 0.875 0.9375 1];
D_m = [1 -1 1 -1 1 -1 1 -1 0.5 0 0 0 0 0 0 0 0];
D_b1 = fir2(32,D_f,D_m);

D_f = [0 0.5 0.5 1];
D_m = [1 0.5 0 0];
D_b1 = fir2(32,D_f,D_m);
freqz(D_b1,1,[]);