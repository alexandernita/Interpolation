# Interpolation of Data Points

Given data points $(x_0,f(x_0)),\dots,(x_n,f(x_n))$, presumably the output of some unknown, but sought, function $f$, we want to connect them either by a single polynomial $p(x)$ passing through all of them, or by separate polynomials each interpolating an adjacent pair of points.  The goal is to *approximate* the $y$-values of $f$ by those of $p$, in order that we may *approximately evaluate* $f(x)$ at other $x$ values besides the ones given, $x_0,\dots,x_n$.  

There is a theoretical justification for this optimism in the form of the **Weierstrass approximation theorem**[^2], which says that polynomials are dense among the continuous functions on a compact interval $[a,b]$ in the supremum norm, $\lVert f\rVert_\infty=\sup_{x\in [a,b]}|f(x)|$.  That is, $\text{Closure}(\mathbb{R}[a,b])=C([a,b])$[^1] in the technical sense that 

$$
\forall f\in C([a,b]),\ \forall \varepsilon>0,\ \exists p\in \mathbb{R}[a,b],\ \lVert f-p\rVert_\infty<\varepsilon
$$  

This of course requires the unknown $f$ to at least belong to $C([a,b])$.  If we opt for a single polynomial for the dataset, then we need $p$ to pass *through* those $n+1$ points in $[a,b]$, and this puts a hefty restriction on $f$, which then needs to belong to $C^{n+1}([a,b])$ to work.  This is one reason to consider *piecewise* interpolants, which make no demands on $f$ other than continuity.  Instead, we pay a heavier price in deriving the algorithms for the interpolants.  

## Single Polynomial Interpolation of $n+1$ Data Points

### Lagrange Polynomial Interpolation 

Given data points $(x_0,f(x_0)),\dots,(x_n,f(x_n))$, we first seek to connect, or **interpolate** them by a **single polynomial** $p(x)$.  The $n\text{th}$ **Lagrange polynomial** 

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

Then the RHS can be maximized over $[a,b]$.  For $n\geq 4$ the Lagrange polynomial $p(x)$ gives a pretty good approximation to $f(x)$ for *any* $x\in [a,b]$, and clearly more accurately with greater $n$ (which makes $\Delta x=\max_i |x_i-x_{i-1}|$ smaller).  

### Computating $y$-Values of Lagrange Polynomials 

#### Generic Polynomial Evaluation

Evaluating a polynomial $p(x)=\sum_{k=0}^n a_kx^k$ with known coefficients $a_k$ at a point $x=x_0$ can be done directly with a single for loop:  

```
A = [a_n,...,a_1,a_0]

p = 0
for i in range(0,n+1):
    p = p+A[n-i]*x0**i
```

#### Horner's Method 

**Horner's method** uses polynomial division.  Let $q(x)=\sum_{k=0}^{n-1}b_{k+1}x^k$ be the **quotient polynomial**, and $b_0$ the remainder term, obtained from dividing $p(x)$ by $(x-x_0)$,

$$
p(x)=(x-x_0)q(x)+b_0
$$

The most basic fact about this is $p(x_0)=0q(x_0)+b_0=b_0$, so *the remainder term $b_0$ is our desired $y$-value*.  When working by hand, **synthetic division** is typically used to compute $b_1,\dots, b_n$,

$$
\begin{aligned}
&\begin{array}{cccc}
        &a_n&\qquad a_{n-1}\qquad &\ \cdots &\quad a_0\\
        x_0\ &\downarrow&b_n&\ \cdots&\quad \ast
\end{array}\\
&\text{-------------------------------------}\\
&\begin{array}{cccc}
        \quad &a_n&a_{n-1}+b_nx_0&\cdots&\ \ast\\
        &=b_n&=b_{n-1}&\cdots&\ =b_0
\end{array}
\end{aligned}
$$
but the underlying idea, which is what we need to code, anyway, is to define the $b_k$ recursively
$$
\begin{aligned}
b_n&\stackrel{\text{def}}{=}a_n\\
b_{n-1}&\stackrel{\text{def}}{=}a_{n-1}+b_nx_0\\
&\ \vdots\\
b_0&\stackrel{\text{def}}{=}a_0+b_1x_0
\end{aligned}
$$

This can be done with a for loop

```
A = [a_n,...,a_1,a_0]
B = [A[0]]
l = len(A)

for i in range(0,l-1):
    b = B[i]
    B.append(A[i+1]+b*x0)
```

Our desired $y$-value $p(x_0)$ is $b_0=B[l]$.  

#### Direct Evaluation of Lagrange Polynomials

We can also simply **compute the Lagrange polynomial's $y$-value $p(x)$ directly** using nested for loops,

```
A = [a_n,...,a_1,a_0]
x = [x0,x1,...,xn] 
xp = input x-value where to evaluate
yp = 0

for i in range(n):
    p = 1

    for j in range(n):
        if i != j:
            p = p * (xp - x[j])/(x[i] - x[j])
    
    yp = yp + p * y[i]    

```

#### Using Newton's Divided Differences to Find Lagrange Coefficients

Once we have our Lagrange polynomial $p(x)=\sum_{i=0}^n f(x_i)L_i(x)$, the above are easy methods ready to hand to compute its $y$-value at any $x$ in $[a,b]$.  However, computing the coefficients $a_i$ of $p(x)$ from its definition is a tedious matter involving much FOILing.  

If we instead write $p(x)$ in the form $p(x)=\sum_{i=0}^n a_i \prod_{j=0}^{i}(x-x_j)$, then there is a method to compute *these* coefficients $a_i$ 

$$
\begin{aligned}
&f[x_i]=f(x_i)  &\text{zeroth divided difference}\\
&f[x_i,x_{i+1}]=\frac{f[x_{i+1}]-f[x_i]}{x_{i+1}-x_i} &\text{first divided difference}\\
&f[x_i,\dots,x_{i+k}]=\frac{f[x_{i+1},\dots,x_{i+k}]-f[x_i,\dots,x_{i+k-1}]}{x_{i+k}-x_i}&\text{$k$th divided difference}
\end{aligned}
$$

In fact, since $f(x_i)=p(x_i)$ for all $i$, we can recursively solve for the $a_i$,

$$
\begin{aligned}
a_0&=f[x_0]\\
a_1&=f[x_0,x_1]\\
a_k&=f[x_0,\dots,x_k]
\end{aligned}
$$

#### Neville's Method

**Neville's method** uses different $k\text{th}$ Lagrange polynomials, $k=0,\dots,n$, interpolating any $k$ of the data points $(x_0,f(x_0)),\dots,(x_n,f(x_n))$, to **recursively compute** $p(x)$, the $n\text{th}$ **Lagrange polynomial**, at a given $x\in [a,b]$, from all $k\text{th}$ order Lagrange polynomials as $k$ goes from $1$ to $n-1$.  Define the degree $n-1$ polynomial

$$
p_{\hat{i}}(x)=\sum_{{j=0, j\neq i}}^nf(x_j)\prod_{k=0, k\neq j}^n\frac{x-x_j}{x_k-x_j}
$$

Then, for any $i\neq j$ it is straightforward to show

$$
p(x)=
\frac{(x-x_j)p_{\hat{j}}(x)-(x-x_i)p_{\hat{i}}(x)}{x_i-x_j}
$$

This works more generally for any $1\leq k\leq n$, and so allows for a recursion algorithm:  for any multi-index $I_k=`\{i_1<\cdots<i_k `\}`\|$, let $J_k=`\{j_1<\ccdots<j_{n-k}`\::`\:j_\ell\notin I_k`\}$

$$
p_{I_k}
$$


[^1]: R[a,b] denotes the polynomials R[x] restricted, as functions, to the interval [a,b].
[^2]:  Theorem 7.26 in Rudin's *Principles of Mathematical Analysis*, or Theorem 8.135 in my *Lectures on Real Analysis*. 