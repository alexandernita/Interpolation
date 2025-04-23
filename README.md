# Interpolation of Data Points

Given data points $(x_0,f(x_0)),\dots,(x_n,f(x_n))$, presumably the output of some unknown, but sought, function $f$, we seek to connect them either by a single polynomial $p(x)$ passing through all of them, or by splicing separate polynomials each interpolating an adjacent pair of points.  The purpose is to *approximating* the $y$-values of $f$ by those of $p$, in order that we may *approximately evaluate* $f(x)$ at other $x$ values besides the ones given, $x_0,\dots,x_n$.  

There is a theoretical justification for this optimism in the form of the **Weierstrass approximation theorem**[^2], which says that polynomials are dense among the continuous functions on a compact interval $[a,b]$ in the supremum norm, $\lVert f\rVert_\infty=\sup_{x\in [a,b]}|f(x)|$.  That is, $\text{Closure}(\mathbb{R}[a,b])=C([a,b])$[^1] in the technical sense that 

$\forall f\in C([a,b]),\ \forall \varepsilon>0,\ \exists p\in \mathbb{R}[a,b],\ \lVert f-p\rVert_\infty<\varepsilon$  

This of course requires the unknown $f$ to at least belong to $C([a,b])$.  If we are using a single polynomial *on the dataset*, then we need $p$ to pass *through* those $n+1$ points, and this puts a hefty restriction on $f$, which then needs to belong to $C^{n+1}([a,b])$.  This is one reason to consider *piecewise* interpolants, which make no demands on $f$ other than continuity.  Instead, the derivation of algorithms for finding piecewise interpolants is harder.  

## Single Polynomial Interpolation of $n+1$ Data Points

Given data points $(x_0,f(x_0)),\dots,(x_n,f(x_n))$, we first seek to connect, or **interpolate** them by a **single polynomial** $p(x)$.  

[^1]: R[a,b] denotes the polynomials R[x] restricted, as functions, to the interval [a,b].
[^2]:  Theorem 7.26 in Rudin's *Principles of Mathematical Analysis*, or Theorem 8.135 in my *Lectures on Real Analysis*. 