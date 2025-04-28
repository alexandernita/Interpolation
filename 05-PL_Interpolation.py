
# Piecewise-Linear Interpolation Example

import numpy as np
import matplotlib.pyplot as plt

# The integer-valued nodes between -1 and 6 
x = np.linspace(1,8)

# The continuous function we'll be approximating
def f(x):
    return x*(x-2)*np.exp(3-x)

# The PL function interpolating f on these nodes
    # Save the y-values of f on the nodes in a list Y
    # to be used in the definition of g(x) following
Y = []
for i in range(1,9):
    Y.append(f(i))
    
def g(x):
    for i in range(1,9):
        if i<= x and x < i+1:
            return (-1)*Y[i-1]*(x-float(i+1))+Y[i]*(x-float(i))
        elif x==8:
            return f(8)

# Save the y-values of g(x) in a list,
# to be used in the plot of g(x) below     
Z = []
for i in range(len(x)):
   Z.append(g(x[i]))

# Plot f(x) and g(x)
W = np.linspace(1,8,8)
B, BX = plt.subplots()
#B1 = BX.plot(x,f(x),color="k",linewidth=1)
B2 = BX.plot(x,Z,color="r",linewidth=1)
BX.legend(("f(x) = x*(x-2)*exp(3-x)","g(x) = the PL function"))
BX.set_title("Plot of f(x) against PL g(x).")
plt.axvline(color="k",linewidth=1)
plt.axhline(color="k",linewidth=1)
B3 = BX.plot(W,f(W),'o',color="b")
plt.show()
