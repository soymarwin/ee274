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
% $X(z) = \sum_{n=0}^{\infty} (\frac{4}{3}) \times (\frac{4}{3})^{-k} \times z^k \times z^{-1}$
%%
% $X(z) = (\frac{4}{3z}) \ \sum_{n=0}^{\infty} (\frac{3z}{4})^{k}$
%%
% $X(z) = (\frac{4}{3z}) \times (\frac{1}{1 - \frac{3z}{4}}), \ 0 < \ \mid {z} \mid \ < \ \frac{4}{3}$

%% $(b) \ \ x(n) = 2^{- \mid {n} \mid} + (\frac{1}{3})^{\mid {n} \mid}$