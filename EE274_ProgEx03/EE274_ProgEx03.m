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

%% B.4. $X(z)={\frac{1-z^{-1}-4z^{-2}+4z^{-3}}{1-\frac{11}{4}z^{-1}+\frac{13}{8}z^{-2}-\frac{1}{4}z^{-3}}}$
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
