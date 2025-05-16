# Lagrange Interpolation

# Importing NumPy Library
import numpy as np

# Reading number of unknowns
n = int(input('\n\n Enter number of data points, n = '))

# Input the data points

X = []
Y = []

for i in range(n):
    X.append(float(input('\t x['+str(i)+'] = ')))
    Y.append(float(input('\t y['+str(i)+'] = ')))

print("\n\t Data Points Entered")
print("\t x\t y")
print("\t ---------------")
for i in range(n):
    print("\t {:<0.2f}\t {:<0.2f}".format(X[i],Y[i]))

# Reading interpolation point
xp = float(input('\n\t Enter interpolation point: x_0 = '))

# Set interpolated value initially to zero
yp = 0

# Implementing Lagrange Interpolation
for i in range(n):
    
    p = 1
    
    for j in range(n):
        if i != j:
            p = p * (xp - X[j])/(X[i] - X[j])
    
    yp = yp + p * Y[i]    

# Displaying output
print('\n\t The interpolated value is p({0:0.2f}) = {1:.7f}.\n'.format(xp, yp))
