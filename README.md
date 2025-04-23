# Interpolation of Data Points

Given data points $(x_0,f(x_0)),\dots,(x_n,f(x_n))$, presumably the output of some unknown, but sought, function $f$, we seek to connect them either by a single polynomial $p(x)$ passing through all of them, or by splicing separate polynomials each interpolating an adjacent pair of points.  The purpose is to *approximating* the $y$-values of $f$ by those of $p$, in order that we may *approximately evaluate* $f(x)$ at other $x$ values besides the ones given, $x_0,\dots,x_n$.  

There is a theoretical justification for this optimism in the form of the **Weierstrass approximation theorem**[^2], which says that polynomials are dense among the continuous functions on a compact interval $[a,b]$ in the supremum norm, $\lVert f\rVert_\infty=\sup_{x\in [a,b]}|f(x)|$.  That is, $\text{Closure}(\mathbb{R}[a,b])=C([a,b])$[^1] in the technical sense that 

$$
\forall f\in C([a,b]),\ \forall \varepsilon>0,\ \exists p\in \mathbb{R}[a,b],\ \lVert f-p\rVert_\infty<\varepsilon
$$  

This of course requires the unknown $f$ to at least belong to $C([a,b])$.  If we are using a single polynomial *on the dataset*, then we need $p$ to pass *through* those $n+1$ points in $[a,b]$, and this puts a hefty restriction on $f$, which then needs to belong to $C^{n+1}([a,b])$.  This is one reason to consider *piecewise* interpolants, which make no demands on $f$ other than continuity.  Instead, we pay a heavier price to derive the algorithms for the interpolants.  

## Single Polynomial Interpolation of $n+1$ Data Points

### Lagrange Polynomial Interpolation 

Given data points $(x_0,f(x_0)),\dots,(x_n,f(x_n))$, we first seek to connect, or **interpolate** them by a **single polynomial** $p(x)$.  The $n$th **Lagrange polynomial** 

$$
p(x)=f(x_0)L_0(x)+\cdots +f(x_n)L_n(x),
\quad \text{where}\quad L_i(x)=\displaystyle \frac{\prod_{j\neq i}(x-x_j)}{\prod_{j\neq i}(x_i-x_j)}
$$

is the *unique degree* $n$ polynomial passing through the $n+1$ points, since 

$$
L_i(x_j)=
\begin{cases}
0,& \text{if }i\neq j\\ 
1,& \text{if }i=j\\
\end{cases}
\quad\implies \quad
f(x_j)L_i(x_j)=
\begin{cases}
0,& \text{if }i\neq j\\ 
f(x_j),&\text{if }i=j
\end{cases}
$$  

For example, if we only have two data points ($n=1$), say $(x_0,f(x_0))$, $(x_1,f(x_1))$, the 1st Lagrange polynomial is the *unique line through the two points*, 

$$
p(x)    = f(x_0)L_0(x)+f(x_1)L_1(x)
        = f(x_0)\frac{x-x_1}{x_0-x_1}+f(x_1)\frac{x-x_0}{x_1-x_0}
$$

while if we have three data points ($n=2$), say $(x_0,f(x_0))$, $(x_1,f(x_1))$, $(x_2,f(x_2))$, then the 2nd Lagrange polynomial is the *unique parabola passing through the three points*,

$$
\begin{aligned}
p(x)    &= f(x_0)L_0(x)+f(x_1)L_1(x)+f(x_2)L_2(x)\\
        &= f(x_0)\frac{(x-x_1)(x-x_2)}{(x_0-x_1)(x_0-x_2)}+f(x_1)\frac{(x-x_0)(x-x_2)}{(x_1-x_0)(x_1-x_2)}+f(x_2)\frac{(x-x_0)(x-x_1)}{(x_2-x_0)(x_2-x_1)}
\end{aligned}
$$

Uniqueness of $p$ follows from the fact that $(L_0, L_1,\dots, L_n)$ forms a basis for the $(n+1)$-dimensional vector space $\mathbb{R}_n[x]$ of polynomials of degree $\leq n$. 

### Error Bound for Lagrange Interpolation

The **generalized Rolle theorem** says

* If $f\in C^n([a,b])$ and if $f(x_i)=0$ for all $a\leq x_0< x_1<\cdots< x_n\leq b$, then $\exists c\in (x_0,x_n)\subseteq [a,b]$ for which $f^{(n)}(c)=0$.  

Using this we can show---using the function $g(t)=f(t)-p(t)-(f(x)-p(x))\prod_{i=0}^n\frac{(t-x_i)}{(x-x_i)})$---that $\exists c\in (a,b)$ for which 

$$
f(x)-p(x)=\frac{f^{(n+1)(c)}}{(n+1)!}\prod_{i=0}^n(x-x_i)
$$

Then the RHS can be maximized over $[a,b]$.  For $n\geq 4$ the Lagrange polynomial $p(x)$ gives a pretty good approximation to $f(x)$ for *any* $x\in [a,b]$, and clearly more accurately with greater $n$ (which makes $\delta x=\max_i |x-x_i|$ smaller).  

### Computating $y$-Values of Lagrange Polynomials 

Evaluating a polynomial $p(x)=\sum_{k=0}^n a_kx^k$ with known coefficients $a_k$ at a point $x=x_0$ can be done directly with a single for loop:  

```
A=[a_n,...,a_1,a_0]

p = 0
for i in range(0,n+1):
    p = p+A[n-i]*x0**i
```

**Horner's method** uses polynomial division.  Let $q(x)=\sum_{k=0}^{n-1}b_{k+1}x^k$ be the quotient polynomial, and $b_0$ the remainder term, obtained from dividing $p(x)$ by $(x-x_0)$,

$$
p(x)=(x-x_0)q(x)+b_0
$$

Then clearly $p(x_0)=0q(x_0)+b_0=b_0$.  When working by hand, synthetic division is typically used to compute $b_1,\dots, b_n$,

$$
\begin{array}{c|rrr}
        &1&2&1\\
     -1 &&-1&-1\\
        \hline
        &1&1&0
\end{array}
$$

and

| command  | description  |
|--------  |------------  |
| alpha    | so and so    |
| beta     | so and so    | 

[^1]: R[a,b] denotes the polynomials R[x] restricted, as functions, to the interval [a,b].
[^2]:  Theorem 7.26 in Rudin's *Principles of Mathematical Analysis*, or Theorem 8.135 in my *Lectures on Real Analysis*. 