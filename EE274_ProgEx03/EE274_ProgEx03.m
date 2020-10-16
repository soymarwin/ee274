%% Marwin B. Alejo   2020-20221   EE274_ProgEx03
% Also accessible through <http://www.github.com/soymarwin/ee274/EE274_ProgEx03>

%% A. The Bilateral Z-Transform

%% $(a) \ \ x(n) = (\frac{4}{3})^n u(1-n)$
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
% $X(z) = \sum_{n=0}^{\infty} (\frac{4}{3}) \cdot (\frac{4}{3})^{-k} \cdot z^k \times z^{-1}$
%%
% $X(z) = (\frac{4z^{-1}}{3}) \ \sum_{n=0}^{\infty} (\frac{3z}{4})^{k}$
%%
% $X(z) = {(\frac{4z^{-1}}{3}) \cdot (\frac{1}{1 - \frac{3z}{4}}), \ 0 \ < \ \mid {z} \mid \ < \ \frac{4}{3}}$
%%
% $or \ X(z) ={\frac{4z^{-1}}{3 - \frac{9z}{4}}, \ 0 \ < \ \mid {z} \mid \ < \ \frac{4}{3}}$
%%
% $or \ X(z) ={\frac{16z^{-1}}{12 - 9z}, \ 0 \ < \ \mid {z} \mid \ < \ \frac{4}{3}}$
%%
% $or \ X(z) ={\frac{16}{12z - 9z^{2}}, \ 0 \ < \ \mid {z} \mid \ < \ \frac{4}{3}}$
%% $(b) \ \ x(n) = 2^{- \mid {n} \mid} + (\frac{1}{3})^{\mid {n} \mid}$
%%
% $X(z) = \sum_{n=0}^{\infty} 2^{-n}z^{-n} + \sum_{n=0}^{\infty} (\frac{1}{3})^{n}z^{-n}$
%%
% $X(z) = \sum_{n=0}^{\infty} (\frac{z^{-1}}{2})^{n} + \sum_{n=0}^{\infty} (\frac{z^{-1}}{3})^{n}$
%%
% $X(z) = \frac{1}{1-\frac{z^{-1}}{2}}+\frac{1}{1-\frac{z^{-1}}{3}}$
%%
% $X(z) = \frac{2z}{2z-1}+\frac{3z}{3z-1}$
%%
% $X(z) = \frac{z(12z-5)}{(2z-1)(3z-1)},\ \frac{1}{3}\ < \ \mid z \mid \ < \ \frac{1}{2}$
