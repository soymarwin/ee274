%% Marwin B. Alejo   2020-20221   EE274_ProgEx05
% Also accessible through
% <http://www.github.com/soymarwin/ee274/EE274_ProgEx05> for history
% tracking.
%%
% _*Disclaimer:* For variable names, I use the section number followed by 
% the actual var name in the original code. This is to ensure that
% variables are not affected by other codes in latter sections on run
% with the publish command. Moreover, all block diagrams required for
% selected sections are in appendix section in lieu to MATLAB's limitation.*_

%% 1. FIR Filter Structure
%% 1.a. Cascade Structure
% tf2zp() require input vectors of the same length hence, the original code
% must be modified as shown below. _Uncomment when needed only_ .
a1_num = input('Numerator coef vector = ');
a1_den = input('Denominator coef vector = ');
[a1_b, a1_a] = eqtflength(a1_num, a1_den);
[a1_z, a1_p, a1_k] = tf2zp(a1_b, a1_a);
a1_sos = zp2sos(a1_z, a1_p, a1_k);

%% 1.b. $H_{1}(z)=2+10z^{-1}+23z^{-2}+34z^{-3}+31z^{-4}+16z^{-5}+4z^{-6}$
b1_num = [2 10 23 34 31 16 4];
b1_den = 1;
[b1_b, b1_a] = eqtflength(b1_num, b1_den);
[b1_z, b1_p, b1_k] = tf2zp(b1_b, b1_a);
b1_sos=zp2sos(b1_z, b1_p, b1_k);
[b1h,b1w]=freqz(b1_num,b1_den);
plot(b1w/pi,20*log10(abs(b1h)));title('Frequency Response of H_{1}(z)');
xlabel('Normalized Frquency (\times\pi rad/sample)'); ylabel('Magnitude(dB)');
fvtool(b1_sos);
b1_cscd=dfilt.dffir(b1_sos);
realizemdl(b1_cscd); % see b1_cscd_mdl.slx for complete cascade generalization diagram of H1(z)
%%
% *$H_{1}z$ is not a linear-phase transfer function due to its coefficients
% not having the required symmetry. Moreover, there is no noticeable
% differences between the magnitude spectrum of $H_{1}z$ as generated using
% freqz() and fvtool() except in the magnitude value near 1pi rad/sample.*

%% 1.c. $H_{2}(z)=6+31z^{-1}+74z^{-2}+102z^{-3}+74z^{-4}+31z^{-5}+6z^{-6}$
c1_num = [6 31 74 102 74 31 6];
c1_den = 1;
[c1_b, c1_a] = eqtflength(c1_num, c1_den);
[c1_z, c1_p, c1_k] = tf2zp(c1_b, c1_a);
c1_sos = zp2sos(c1_z, c1_p, c1_k);
[c1h,c1w]=freqz(c1_num,c1_den);
plot(c1w/pi,20*log10(abs(c1h)));title('Frequency Response of H_{2}(z)');
xlabel('Normalized Frquency (\times\pi rad/sample)'); ylabel('Magnitude(dB)');
fvtool(c1_sos);
c1_cscd=dfilt.dffir(c1_sos);
realizemdl(c1_cscd); % see c1_cscd_mdl.slx for complete cascade generalization diagram of H2(z)

%%
% * $H_{2}z$ is a linear-phase transfer function (type-1) with an odd length
% and even symmetry. Moreover, there is no noticeable differences between 
% the magnitude spectrum of $H_{2}z$ as generated using freqz() and 
% fvtool().* 

%% 2. IIR Filter Structure
% $H(z) = 2(\frac{1+0z^{-1}+z^{-2}}{1+0.8z^{-1}+0.64z^{-2}})(\frac{2-z^{-1}}{1-0.7z^{-1}})(\frac{1+2z^{-1}+z^{-2}}{1+0.81z^{-2}})$

% From TF
fb0_2=2; fb1_2=[1 2 1]; fb2_2=[0 1 -1 2 1];
fa1_2=[1 1 1]; fa2_2=[0.8 0.64 -0.75 0 0.81];
fb_2=fb0_2*conv(fb1_2,fb2_2); fa_2=conv(fa1_2,fa2_2);
[fh_2,fw_2]=freqz(fb_2, fa_2,'whole',101); figure(), plot(fw_2/pi,abs(fh_2));
title('Frequency Response of H(z) with 0 \leq n \leq 100'); 
xlabel('Normalized Frquency (\times\pi rad/sample)'); ylabel('Magnitude(rad)');

%% 2.a. Draw DF1 and DF2 realizations of H(z)
a2_df1=dfilt.df1(fb_2,fa_2);
a2_df2=dfilt.df2(fb_2,fa_2);
realizemdl(a2_df1); % see a2_df1_mdl.slx for complete DF1 realization diagram of H(z)
realizemdl(a2_df2); % see a2_df2_mdl.slx for complete DF2 realization diagram of H(z)

%% 2.b. Draw cascade containing sos DF2 realizations of H(z)
[fb2_2, fa2_2] = eqtflength([1 0 1], [1 0.8 0.64]);
b2_df2=dfilt.df2sos(fb2_2, fa2_2);
[fb2_21, fa2_21] = eqtflength([2 -1], 1 -0.75);
b2_df21=dfilt.df2sos(fb2_21, fa2_21);
[fb2_22, fa2_22] = eqtflength([1 2 1], [1 0.81]);
b2_df22=dfilt.df2sos(fb2_22, fa2_22);
b2_cscd=dfilt.cascade(b2_df2,b2_df21,b2_df22);
realizemdl(b2_cscd); % see b2_df2sos_mdl.slx for complete DF2sos realization diagram

%% 2.c. Draw parallel form containing sos DF2 realization of H(z)
c2_cscd=dfilt.parallel(b2_df2,b2_df21,b2_df22);
realizemdl(c2_cscd); % see c2_df2sos_mdl.slx for complete DF2sos realization diagram

%% 3. Effects of coefficient quantization

% Program P9_1 
% Coefficient Quantization Effects on Direct Form 
% Realization of an IIR Transfer Function 
[b391,a391] = ellip(6,0.05,60,0.4); 
[g391,w391]=Gain(b391,a391); 
bq391 = a2dT(b391,5); aq391 = a2dT(a391,5); 
[gq391,w391] = Gain(bq391, aq391); 
figure(); plot(w391/pi,g391,'b', w391/pi,gq391,'r--'); 
axis([0 1 -80 1]);grid 
xlabel('\omega /\pi');ylabel('Gain, dB'); 
legend('original', 'quantized'); 
figure(); zplane(b391,a391); 
hold on;
pzplot(bq391,aq391); 
hold off;
title('Original pole-zero locations: x, o; New pole-zero locations: +, *') 

% Program P9_2 
% Coefficient Quantization Effects on Cascade 
% Realization of an IIR Transfer Function 
[z923,p923,k923] = ellip(6,0.05,60,0.4); 
[b923,a923] = zp2tf(z923,p923,k923); 
[g923] = Gain(b923,a923); 
sos923 = zp2sos(z923,p923,k923); 
sosq923 = a2dR(sos923,6); 
R1923 = sosq923(1,:);R2923 = sosq923(2,:);R3923 = sosq923(3,:); 
b1923 = conv(R1923(1:3),R2923(1:3));bq923 = conv(R3923(1:3),b1923); 
a1923 = conv(R1923(4:6),R2923(4:6));aq923 = conv(R3923(4:6),a1923); 
[gq923,w923] = Gain(bq923, aq923); 
figure(); plot(w923/pi,g923,'b-',w923/pi,gq923,'r--'); 
axis([0 1 -80 20]);grid 
xlabel('\omega /pi');ylabel('Gain, dB'); 
title('original - solid line; quantized - dashed line'); 
figure(); zplane(b923,a923); 
hold on 
pzplot(bq923,aq923); 
hold off 
title('Original pole-zero locations: x, o; New pole-zero locations: +, *') 

% Modified Program P9_2 
% Coefficient Quantization Effects on Cascade 
% Realization of an IIR Transfer Function 
[z923_1,p923_1,k923_1] = ellip(8,0.1,70,0.55); 
[b923_1,a923_1] = zp2tf(z923_1,p923_1,k923_1); 
[g923_1] = Gain(b923_1,a923_1); 
sos923_1 = zp2sos(z923_1,p923_1,k923_1); 
sosq923_1 = a2dR(sos923_1,5); 
R1923_1 = sosq923_1(1,:); 
R2923_1 = sosq923_1(2,:); 
R3923_1 = sosq923_1(3,:); 
b1923_1 = conv(R1923_1(1:3),R2923_1(1:3));
bq923_1 = conv(R3923_1(1:3),b1923_1); 
a1923_1 = conv(R1923_1(4:6),R2923_1(4:6));
aq923_1 = conv(R3923_1(4:6),a1923_1); 
[gq923_1,w923_1] = Gain(bq923_1, aq923_1); 
figure(); plot(w923_1/pi,g923_1,'b-',w923_1/pi,gq923_1,'r--'); 
axis([0 1 -80 20]); grid; xlabel('\omega /pi');ylabel('Gain, dB'); 
title('original - solid line; quantized - dashed line'); 
figure(); zplane(b923_1,a923_1); hold on; pzplot(bq923_1,aq923_1); 
title('Original pole-zero locations: x, o; New pole-zero locations: +, *') 

%%
% *The frequency response of the unquantized filter (blue) as well as the 
% quantized filter (red). By inspection, the coefficients after 5-bit 
% quantization has adversely affected the frequency response. Frequency 
% response of the quantized filter turned-out differently from that of the 
% original designed filter.*
%%
% *The zero-pole plot of the quantized (red) and unquantized (blue) filter. 
% No poles are present outside the unit circle hence, both the unquantized 
% and quantized filter is stable. Since the filter is quantized on the 
% 8th-order of a cascade of second-order sections, quantized filter’s 
% sensitivity possibly reduced hence, becomes stable.* 