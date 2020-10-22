%% Marwin B. Alejo   2020-20221   EE274_ProgEx03
% Also accessible through
% <http://www.github.com/soymarwin/ee274/EE274_ProgEx03>; for history
% tracking.
%% A.1-2. The Bilateral Z-Transform
%% Sequence (a) $x(n) = (\frac{4}{3})^n u(1-n)$
% *Manual Solution*
%%
% $x(n) = (\frac{4}{3})^n u(-n+1)$
%%
% $X(z) = \sum_{n=-\infty}^{\infty} x(n)z^{-n}$
%%
% $X(z) = \sum_{n=-\infty}^{\infty} (\frac{4}{3})^n u(-n+1)z^{-n}$
%%
% $Let \ k=-n+1 \ and \ n=1-k$
%%
% $X(z) = \sum_{n=-\infty}^{\infty} (\frac{4}{3})^{1-k} u(k)z^{k-1}$
%%
% $X(z) = \sum_{n=0}^{\infty} (\frac{4}{3}) \cdot ((\frac{4}{3})^{-1})^{k} \cdot ((1/z)^{-1})^{k} \cdot z^{-1}$
%%
% $X(z) = (\frac{4z^{-1}}{3}) \ \sum_{n=0}^{\infty} (\frac{3}{4z^{-1}})^{k}$
%%
% $X(z) = {(\frac{4z^{-1}}{3}) \cdot (\frac{1}{1 - \frac{3}{4z^{-1}}}), \ 0 \ < \ \mid {z} \mid \ < \ \frac{4}{3}}$
%%
% $or \ X(z) ={\frac{16z^{-2}}{-9+12z^{-1}}, \ 0 \ < \ \mid {z} \mid \ < \ \frac{4}{3}}$
%%
% $or \ X(z) ={\frac{-16z^{-2}}{9-12z^{-1}}, \ 0 \ < \ \mid {z} \mid \ < \ \frac{4}{3}}$

%%
% *z-plane for 1.(a)*
A1_a_a=[-9, 12, 0];
A1_a_b=[0, 0, -16];
zplane(A1_a_b,A1_a_a);

%%
% *Verification of z-transform v. original sequence with first 8-coef.*
[delta,n]= impseq(0,0,7);
A_a_Xz=filter(A1_a_b,A1_a_a,delta) %A_a_Xz is z-transform sequence
A_a_Xn=[(4/3).^n].*stepseq(1,0,7) 
%A_a_Xn is the original sequence, see stepseq.m

%%
% *Therefore, based on coef values generated from X(z) and x(n),
% the z-transform for sequence(a) is correct.*

%% Sequence (b) $x(n) = 2^{- \mid {n} \mid} + (\frac{1}{3})^{\mid {n} \mid}$
%%
% $X(z) = \sum_{n=0}^{\infty} 2^{-n}z^{-n} + \sum_{n=0}^{\infty} (\frac{1}{3})^{n}z^{-n}$
%%
% $X(z) = \sum_{n=0}^{\infty} (\frac{z^{-1}}{2})^{n} + \sum_{n=0}^{\infty} (\frac{z^{-1}}{3})^{n}$
%%
% $X(z) = \frac{1}{1-\frac{z^{-1}}{2}}+\frac{1}{1-\frac{z^{-1}}{3}}$
%%
% $X(z) = \frac{2}{2-z^{-1}}+\frac{3}{3-z^{-1}}$
%%
% $X(z) = \frac{12-5z^{-1}}{(2-z^{-1})(3-z^{-1})},\ \mid z \mid \ > \ \frac{1}{3} \ \cap \ \mid z \mid \ > \ \frac{1}{2}$
%%
% $or X(z) = \frac{12-5z^{-1}}{6-5z^{-1}+z^{-2}},\ \mid z \mid \ > \ \frac{1}{3} \ \cap \ \mid z \mid \ > \ \frac{1}{2}$
%%
% *z-plane for 1.(b)*
A1_b_a=[6 -5 1];
A1_b_b=[12 -5 0];
zplane(A1_b_b,A1_b_a);

%%
% *Verification of z-transform v. original sequence with first 8-coef.*
[delta,n]= impseq(0,0,7);
A_b_Xz=filter(A1_b_b,A1_b_a,delta) %A_b_Xz is z-transform sequence
A_b_Xn=((2).^(-abs(n)))+((1/3).^(abs(n))) %A_b_Xn is the original sequence

%%
% *Therefore, based on coef values generated from X(z) and x(n),
% the z-transform for sequence(b) is correct.* 

%% A.3. $x(n)=(\frac{1}{3})^{n}u(n-2)+(0.9)^{n-3}u(n)$
%%
% $X(z)={\frac{3z^{-2}}{27-9z^{-1}}}+{\frac{1.3717}{1-0.9z^{-1}}}$
%%
% $X(z)={\frac{{37.0359}-{12.3453z^{-1}}+{3z^{-2}}-{2.7z^{-3}}}{{27}-{33.3z^{-1}}+{8.1z^{-2}}}} \ \mid z \mid \ > \ \frac{1}{3} \ \cap \ \mid z \mid \ > \ {0.9}$

%%
% *z-plane for A.3*
A3_b=[37.0359, -12.3453, 3, -2.7];
A3_a=[27, -33.3, 8.1];
zplane(A3_b,A3_a);

%%
% *Verification of z-transform v. original sequence with first 20-coef.*
[delta,n]= impseq(0,0,19);
A3_Xz=filter(A3_b,A3_a,delta) %A3_Xz is z-transform sequence
A3_Xn=(((1/3).^n).*(stepseq0(2,0,19))+(((0.9).^(n-3)).*(stepseq0(0,0,19)))) 
%A3_Xn is the original sequence, see stepseq0.m

%%
% *Therefore, based on coef values generated from X(z) and x(n),
% the z-transform for sequence in (A.3.) is correct.* 
%% B.4. Inverse Z-Transform
%% Sequence(c) $X(z)={\frac{1-z^{-1}-4z^{-2}+4z^{-3}}{1-\frac{11}{4}z^{-1}+\frac{13}{8}z^{-2}-\frac{1}{4}z^{-3}}}$
B4_b=[1, -1, -4, 4];
B4_a=[1, (-11/4), (13/8), (-1/4)];
[B4_R, B4_p, B4_C]=residuez(B4_b,B4_a);
%%
% $X(z)=\frac{0z}{z-2} - \frac{10z}{z-0.5} + \frac{27z}{z-0.25} - {16}$
%%
% $X(n)=u(-n)-(2^{-2n}(5 \times 2^{n+1}-27)(1-u(-n)))$
%%
% *Verification of z-transform v. ans sequence with first 8-coef.*
%%
% _Disclaimer: First element is a garbage value. Thus, array(2:9)_ 
[delta,n]= impseq(0,0,8);
B4_Xz=filter(B4_b,B4_a,delta); %B4_Xz is z-transform sequence
%B4_Xn is inv. ztrans sequence
B4_Xn=-heaviside(-n)-((2.^(-2*n)).*(5.*(2.^(n+1))-27).*(1-heaviside(-n)));
B4_Xz(2:8)% First 8 coef of B4_Xz - Z-transf 
B4_Xn(2:8)% First 8 coef of B4_Xn - Inv. Z-transf

%% C.5. Signal Generation
%%
% _Generate the periodic even symmetric square pulse x(n) from [0, 1].
% Thepriod pulse is 1 second and a pulse with of 250ms with sampling freq.
% of 8KHz. Plot one period of x(n) and verify if you have the correct
% waveform._ 

C5_t=0:1/8e3:1; % time | x-axis
% C5_x is our x(n)
C5_x=square(cos(4*pi*C5_t)); % since 1pw = 250ms; 1period = 2pw; 4pw in 1s. 
C5_x=(abs(C5_x)+C5_x)/2; % eliminating all -1 with 0.
figure();
plot(C5_t, C5_x);
axis([0 1 -0.2 1.2]); title("1 second plot of the signal: 2 pulses are present"); % 2 periods w/ 250ms pw each 1,0.
xlabel("Time in second (s)");ylabel("Amplitude");
C5_samp_prd = (length(C5_x))/2;
figure();
plot(C5_t(1:C5_samp_prd), C5_x(1:C5_samp_prd));
axis([0 1 -0.2 1.2]); title("A period of the signal");
xlabel("Time in second (s)");ylabel("Amplitude");

%%
% *There are two periods in 1s of specified pulse conf.*

%% 5.a. How many samples in one period?
%% 
sample_periodic=(length(C5_x)-1)/2 % samples in one period
sample_1second=length(C5_x)-1 % samples in 1s
%%
% *There are 4000 samples in one period while 8000 samples for the entire 1s.*
%% 5.b. How many samples with a value of 1?
% -1 was added to negate matlab's indexing rule that starts with 1.
value1_periodic=(sum(C5_x(:)==1)-1)/2 % in one period
value1_1second=(sum(C5_x(:)==1)-1) % in 1s
%%
% *There are 2000 samples with a value of 1 in one period and 4000 samples in 1s.*
%% 5.c. How many zeros?
value0_periodic=sum(C5_x(:)==0)/2 % in one period
value0_1second=sum(C5_x(:)==0) % in 1s
%%
% *There are 2000 samples with a value of 0 in one period and 4000 samples in 1s.*

%% C.6.Fourier Series Analysis Equation
%%
% _Using the analysis equation of the Fourier series, write a program
% that will compute the Fourier series coefficients of the periodic square
% pulse signal. Plot the magnitude and phase of the first 10 Fourier coef._

% plot of 100000 harmonic
C6_t = (0:1/8000:1);
C6_NHarmonics=100000; C6_Ncycles=2; C6_Nsamples=8000;
C6_y(1:C6_Nsamples)=0.5;C6_j=1:C6_Nsamples;
for C6_k=1:C6_NHarmonics
 C6_x(C6_j)=(2*sin(0.5*pi*C6_k)/(pi*C6_k))*cos(C6_k*2*pi*C6_Ncycles*C6_j/C6_Nsamples);
 C6_y=C6_y+C6_x;
end
figure(), plot(C6_t(1:8000),C6_y);axis([0 1 -0.2 1.2]);
title("Illustration of the Fourier synthesis with 100000 harmonics");
xlabel("Time in second (s)"); ylabel("Amplitude");
% plot of 10 harmonic
C6_t = (0:1/8000:1);
C6_NHarmonics=10; C6_Ncycles=2; C6_Nsamples=8000;
C6_y(1:C6_Nsamples)=0.5;C6_j=1:C6_Nsamples;
for C6_k=1:C6_NHarmonics
 C6_x(C6_j)=(2*sin(0.5*pi*C6_k)/(pi*C6_k))*cos(C6_k*2*pi*C6_Ncycles*C6_j/C6_Nsamples);
 C6_y=C6_y+C6_x;
end
figure(), plot(C6_t(1:8000),C6_y);axis([0 1 -0.2 1.2]);
title("Illustration of the Fourier synthesis with 10 harmonics");
xlabel("Time in second (s)"); ylabel("Amplitude");

% compute the magnitude and phase of 10 harmonics signal
C6y_fft = fft(C6_y);
C6_mag = abs(C6y_fft);
C6_ang = angle(C6y_fft);

figure(), stem(C6_mag(1:10)); title("Magnitude Spectrum of first 10 Coef"); xlabel("coef"); ylabel("Magnitude Spectrum");
figure(), stem(C6_ang(1:10)); title("Phase Spectrum of first 10 Coef"); xlabel("coef"); ylabel("Phase Spectrum");

%% 6.a. What is the fundamental frequency of the square pulse?
%%
% *The fundamental frequency is defined by f=1/0.5s as such one complete
% period is 0.5s or 2Hz and having 250ms of on and off.* 

%% 6.b. Enumerate the Magnitude and Phase of first 10 coef.
disp("Magnitude");
C6_mag(1:10)

disp("Phase");
C6_ang(1:10)

%% C.7. Fourier Series Synthesis Equation
%%
% _Using the synthesis equatin for the Foyrier series, synthesiez the
% original square pulse using the first 10 Fourier coefficients. Generate a
% plot of the original square pulse and the synthesized square pulse._ 

figure(), plot(C6_t(1:8000),C6_y, 'color', 'r'); 
hold on; plot(C6_t(1:8000),C5_x(1:8000), 'color', 'b'); 
hold off; title("Original square signal v. 10-harmonics square signal"); 
xlabel("Amplitude"); ylabel("Time in second (s)");

%% 7.a. What is the average MSE of original square pulse vs synthesized pulse?
% *MSE is 1.0100%*
C7_mse_10harm = immse(C5_x(1:8000),C6_y)

%% 7.b. If you use 20 Fourier coef, what will be the MSE?
% *MSE will be 0.51%*

C6_t = (0:1/8000:1);
C6_NHarmonics=20; C6_Ncycles=2; C6_Nsamples=8000;
C6_y(1:C6_Nsamples)=0.5;C6_j=1:C6_Nsamples;
for C6_k=1:C6_NHarmonics
 C6_x(C6_j)=(2*sin(0.5*pi*C6_k)/(pi*C6_k))*cos(C6_k*2*pi*C6_Ncycles*C6_j/C6_Nsamples);
 C6_y=C6_y+C6_x;
end

C7_mse_20harm = immse(C5_x(1:8000),C6_y)

figure(), plot(C6_t(1:8000),C6_y, 'color', 'r'); 
hold on; plot(C6_t(1:8000),C5_x(1:8000), 'color', 'b'); 
hold off; title("Original square signal v. 20-harmonics square signal"); 
xlabel("Amplitude"); ylabel("Time in second (s)");

%% 7.c. What is the effect on the fundamental freq if I increase the pulse width to 300ms? Explain.
% *Fundamental frequency will be shortened to 1.6667Hz from 2.000Hz. One
% period is composed of two pulse width in this activity hence, changing
% the pulse width from 250ms to 300ms yields to an increased period/time.
% Moreover, the original 4 pulses seen in 1 second tends to expand making a
% portion of the original cycle invisble under the 1sec timeframe.*
%% 7.d. What is the effect on the Fourier coef if I change the pulse width?
% *Depending on the rate of change of values of pulse width but when it is
% altered, it inversely alter the value of Fourier coef.*
%% 7.e. What is the effect on the Fourier coef if I change the period?
% *Depending on the rate of change of values of period, being it as
% twice as the pulse width in a periodic signal, but when it is
% altered, it inversely alter the value of Fourier coef.*