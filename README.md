# Interpolation of Data Points

Given data points $(x_0,f(x_0)),\dots,(x_n,f(x_n))$, presumably the output of some unknown, but sought, function $f$, we want to connect them either by a single polynomial $p(x)$ passing through all of them, or by separate polynomials each interpolating an adjacent pair of points.  The goal is to *approximate* the $y$-values of $f$ by those of $p$ *at other points $x$ than the nodes $x_i$*.

There is a theoretical justification for this optimism in the form of the **Weierstrass approximation theorem**[^2], which says that polynomials are dense among the continuous functions on a compact interval $[a,b]$ in the supremum norm, $\lVert f\rVert_\infty=\sup_{x\in [a,b]}|f(x)|$.  That is, $\text{Closure}(\mathbb{R}[a,b])=C([a,b])$[^1] in the technical sense that 

$$
\forall f\in C([a,b]),\ \forall \varepsilon>0,\ \exists p\in \mathbb{R}[a,b],\ \lVert f-p\rVert_\infty<\varepsilon
$$  

This theorem, however, gives no indication how to *find* such a polynomial for a given $f$, which of course must now belong to $C([a,b])$.  The Lagrange and Hermite polynomials are found *algorithmically*, and exemplify the *practical* side of the mathematical problem.  Such a polynomial needs $p$ to pass through the $n+1$ nodes, and potentially (as with Hermite polynomials) have $p'(x_i)$ agree with $f'(x_i)$ at each node.  This puts extra restrictions on candidate $f$s, which now need to belong to $C^{n+1}([a,b])$ to work.  

This is one reason to consider *piecewise* interpolants, which make no demands on $f$ other than continuity.  Instead, we pay a heavier price in deriving the algorithms for the interpolants.  

# Contents 

1. Compute $y$-values of polynomials $p(x)$ in standard form $\sum_{i=0}^n a_ix^i$.  
    ```
    Input = coefficients $a_n,\dots,a_1,a_0$ and $x$-value $x_0$
    Output = $p(x_0)$
    ```
    We do this in two ways, **directly** using ```for``` loops, and via **Horner's method** which uses a ```for``` loop based on the synthetic division algorithm.
2. **Single polynomial interpolation** of $n+1$ data points $(x_0,y_0),\dots,(x_n,y_n)$, presumably on the graph of some unknown but sought function $f(x)$.   
    a. **Lagrange interpolation** is the most common and straightforward.  We illustrtate **three ways** to compute the $y$-value of a Lagrange polynomial:  
        i. **direct** 
        ii. **Neville's method** using divided differences
        iii. **Newton's method** using divided differences
    b. **Hermite interpolation** tries to improve on Lagrange by requiring also the *derivatives* $p'(x_i)$ to agree with $f'(x_i)$ at the nodes.  Its computation is more involved, but the resulting interpolant is more accurate.
3. **Piecewise polynomial interpolation** uses ***different polynomials*** $p_i(x)$ on different subintervals $[x_i,x_{i+1}]$, of the same degree.
    a. **Piecewise linear** interpolation uses different lines between nodes
    b. **piecewise cubic** (aka **cubic splines**) interpolantion uses cubics between nodes