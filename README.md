# Interpolation of Data Points

In a world dominated by information, we often find ourselves possessed of collections of **data points**

\begin{array}{l|l}
x           &y\\
\hline
x_0        &y_0\\
\ \vdots   &\ \vdots\\
x_n        &y_n
\end{array}

the second column presumably the outputs $y_i=f(x_i)$ of some unknown, but sought, function $f(x)$ describing a process of interest.  Since we don't know $f(x)$, we try to ***guess*** its other $y$-values by means of some known simple function, like a line or a polynomial.  This ***connecting*** of the dots by a guessed function of convenience is called ***interpolation***.  There are ***two basic methods***:

1. Connecting ***all*** $n+1$ points by a ***single polynomial*** $p(x)$ passing through each of them.
2. Connecting ***each pair*** of data points by ***different polynomials*** (of the same degree) each interpolating one adjacent pair of points.  

The goal, of course, is to ***approximate*** the $y$-values of $f$ by those of our guess, $p$, in order that we may ***approximately evaluate*** $f(x)$ at other $x$ values besides the ones given, $x_0,\dots,x_n$.  

There is a theoretical justification for this optimism in the form of the **Weierstrass approximation theorem**[^2], which says that polynomials are dense among the continuous functions on a compact interval $[a,b]$ in the supremum norm, $\lVert f\rVert_\infty=\sup_{x\in [a,b]}|f(x)|$.  That is, $\text{Closure}(\mathbb{R}[a,b])=C([a,b])$[^1] in the technical sense that 

$$
\forall f\in C([a,b]),\ \forall \varepsilon>0,\ \exists p\in \mathbb{R}[a,b],\ \lVert f-p\rVert_\infty<\varepsilon
$$  

This theorem, however, gives no indication how to *find* such a polynomial for a given $f$, which of course must now belong to $C([a,b])$.  The Lagrange and Hermite polynomials are found *algorithmically*, and exemplify the *practical* side of the mathematical problem.  Such a polynomial needs $p$ to pass through the $n+1$ nodes, and potentially (as with Hermite polynomials) have $p'(x_i)$ agree with $f'(x_i)$ at each node.  This puts extra restrictions on candidate $f$s, which now need to belong to $C^{n+1}([a,b])$ to work.  

This is one reason to consider *piecewise* interpolants, which make no demands on $f$ other than continuity.  Instead, we pay a heavier price in deriving the algorithms for the interpolants.  