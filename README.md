# Polynomial Interpolation of Data Points, Piecewise Linear Interpolation and Cubic Splines

Given data points $(x_0,f(x_0)),\dots,(x_n,f(x_n))$, presumably the output of some unknown, but sought, function $f$, we seek to connect them either by a single polynomial $p(x)$ passing through all of them, or by splicing separate polynomials each interpolating an adjacent pair of points.  The purpose is to *approximating* the $y$-values of $f$ by those of $p$, in order that we may *approximately evaluate* $f(x)$ at other $x$ values besides the ones given, $x_0,\dots,x_n$.  

There is a theoretical justification for this optimism in the form of the **Weierstrass approximation theorem**, which says that polynomials are dense among the continuous functions on a compact interval $[a,b]$ in the supremum norm, $\lVert f\rVert_\infty=\sup_{x\in [a,b]}|f(x)|$.

That is, <span style="text-decoration:overline">$\mathbb{R}[a,b])$</span>$=C([a,b])$[^1] in the technical sense that $\forall \varepsilon>0,\ \exists p(x)\in \mathbb{R}[a,b],\ \lVert f-p\rVert_\infty<\varepsilon$.[^2]  

## Single Polynomial Interpolation of $n+1$ Data Points

Given data points $(x_0,f(x_0)),\dots,(x_n,f(x_n))$, we first seek to connect them by a single polynomial $p(x)$.  

[^1]: $\mathbb{R}[a,b]$ denotes the polynomials $\mathbb{R}[x]$ restricted, as functions, to the interval $[a,b]$.
[^2]:  Theorem 7.26 in Rudin's *Principles of Mathematical Analysis*, or Theorem 8.135 in my *Lectures on Real Analysis*. 