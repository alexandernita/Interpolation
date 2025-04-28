# 3.5 - Example 2

import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return np.exp(x)

# The inputs and computation of the h_i and alpha_i
n = 3

x = [i for i in range(n+1)]
a = [f(i) for i in range(n+1)]
h = [x[i+1]-x[i] for i in range(n)]
alpha = [3*(a[i+1]-a[i])/h[i]-3*(a[i]-a[i-1])/h[i-1] for i in range(1,n)]

# Solving the tridiagonal matrix equation using Algorithm 6.7
    # Step 1: define l's, m's and z's
l = [1]
m = [0]
z = [0]

for i in range(1,n):
    l.append(2*(x[i+1]-x[i-1])-h[i-1]*m[i-1])
    m.append(h[i]/l[i])
    z.append((alpha[i-1]-h[i-1]*z[i-1])/l[i])

l.append(1)
z.append(0)

    # Step 2: solve for the c's, b's and d's
C = [0]
c = []
b = []
d = []
B = []
D = []

for j in range(n+1):
    C.append(z[n-1-j] - m[n-1-j]*C[j])
    B.append((a[n-j]-a[n-1-j])/h[n-1-j]-h[n-1-j]*(C[j]+2*C[j+1])/3)
    D.append((C[j]-C[j+1])/(3*h[n-1-j]))

for k in range(n):
    c.append(C[n-k])
    b.append(B[n-1-k])
    d.append(D[n-1-k])

a = a[0:3]
for j in range(n):
    print("For S_%d"%j,"(x), the constants are:    a[%d]"%j,"= %0.5f"%a[j],\
          "   b[%d]"%j," = %0.5f"%b[j],\
          "   c[%d]"%j,"= %0.5f"%c[j],\
          "   d[%d]"%j,"= %0.5f"%d[j])

def g(x):
    xi = 0.0
    for i in range(n):
        if xi <= x and x < xi + 1.0:
            return a[i]+b[i]*(x-xi)+c[i]*(x-xi)**2+d[i]*(x-xi)**3

        elif x == 3.0:
            return f(3)
        xi = xi + 1.0
        

x2 = np.linspace(0,3)
Z = []
for j in range(len(x2)):
    rr = x2[j]
    Z.append(g(rr))

# Plot f(x) and g(x)
W = np.linspace(0,3,4)
K, BX = plt.subplots()
B1 = BX.plot(x2,f(x2),color="k",linewidth=1)
B2 = BX.plot(x2,Z,color="r",linewidth=1)
BX.legend(("f(x) = exp(x)","g(x) = the cubic spline"))
BX.set_title("Plot of f(x) against the cubic spline g(x).")
plt.axvline(color="k",linewidth=1)
plt.axhline(color="k",linewidth=1)
B3 = BX.plot(W,f(W),'o',color="b")
plt.show()
        
