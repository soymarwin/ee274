%% Marwin B. Alejo   2020-20221   EE274_ProgEx03
% Also accessible through
% <http://www.github.com/soymarwin/ee274/EE274_ProgEx03>.

%% A. The Bilateral Z-Transform

%% $(a) \ \ x(n) = (\frac{4}{3})^n u(1-n)$
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
A_a_Xn=[(4/3).^n].*stepseq(1,0,7) %A_A_Xn is the original sequence

%%
% *Therefore, based on coef values generated from X(z) and x(n),
% the z-transform for sequence(a) is correct.*

%% $(b) \ \ x(n) = 2^{- \mid {n} \mid} + (\frac{1}{3})^{\mid {n} \mid}$
%%
% $X(z) = \sum_{n=0}^{\infty} 2^{-n}z^{-n} + \sum_{n=0}^{\infty} (\frac{1}{3})^{n}z^{-n}$
%%
% $X(z) = \sum_{n=0}^{\infty} (\frac{z^{-1}}{2})^{n} + \sum_{n=0}^{\infty} (\frac{z^{-1}}{3})^{n}$
%%
% $X(z) = \frac{1}{1-\frac{z^{-1}}{2}}+\frac{1}{1-\frac{z^{-1}}{3}}$
%%
% $X(z) = \frac{2}{2-z^{-1}}+\frac{3}{3-z^{-1}}$
%%
% $X(z) = \frac{12-5z^{-1}}{(2-z^{-1})(3-z^{-1})},\ \frac{1}{3}\ < \ \mid z \mid \ < \ \frac{1}{2}$
%%
% $or X(z) = \frac{12-5z^{-1}}{6-5z^{-1}+z^{-2}},\ \frac{1}{3}\ < \ \mid z \mid \ < \ \frac{1}{2}$

% z-plane for 1.(b)
A1_b_a=[6 -5 1];
A1_b_b=[12 -5 0];
zplane(A1_b_b,A1_b_a);
