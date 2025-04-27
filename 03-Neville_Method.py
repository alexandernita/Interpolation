# Lagrange Interpolation

import numpy as np

# Reading number of unknowns
# n = int(input('Enter number of data points: '))

x = [-2.0, -1.0, 0.0, 1.0, 2.0]

n = len(x)
Q = [[] for i in range(n)]

# our function and its y-values
def f(x):
    return 3**x

y = [f(x[i]) for i in range(len(x))]

# store y-value and declare x-value at which to evaluate p(x) 
Q[0] = y

# the x-value at which we want to evaluate p(x)
xp = 0.5 

# print title and data points
print("\n\tExercise 3.2.3 (a), Neville's Method:")
print("\n\tGiven %d"%n,"data points on the graph of f,")
print("\n\t x\t f(x) = 3^x")
print("\t------------------------")
for i in range(len(x)):
    print("\t",format(x[i], ".1f"),"\t",format(y[i], ".7f"))
print("\n\tthe degree -1+%d"%n,"Lagrange interpolating polynomial")
print("\tapproximates f(%0.1f)"%xp,"by the value in the bottom right.\n")

# implement Neville's Method
for i in range(1,n):
    for j in range(i,n):
        xi = x[j]
        xij = x[j-i]
        Qij = Q[i]
        R = Q[i-1]
        s = ((xp-xij)*R[j-i+1]-(xp-xi)*R[j-i])/(xi-xij)
        Qij.append(s)
    Q[i]=Qij

# display output
# add zeros to the beginning of all the Q[i] vectors
for i in range(n):
    T = ["" for k in range(i)]
    Q[i]=T+Q[i]

# transpose Q, so that it sits right
S = [[] for i in range(n)]

for i in range(n):
    for j in range(n):
        S[i].append(Q[j][i])

# print out S=transpose(Q)
for i in range(n):
    print("\t",format(S[i][0], ".7f") if type(S[i][0])==float else "",\
          "\t",format(S[i][1], ".7f") if type(S[i][1])==float else "",\
          "\t",format(S[i][2], ".7f") if type(S[i][2])==float else "",\
          "\t",format(S[i][3], ".7f") if type(S[i][3])==float else "",\
          "\t",format(S[i][4], ".7f") if type(S[i][4])==float else "")

# error analysis 
pxp = S[n-1][n-1]
fp = f(xp)
Ep = abs(pxp - fp)
print("\n\tThe approximate and actual values, as well as the error, are:")
print("\n\tp(%0.1f)"%xp,"\t\tf(%0.1f)"%xp,"\t\tE(%0.1f)"%xp)
print("\t--------------------------------------")
print("\t%0.7f"%pxp,"\t%0.7f"%fp,"\t%0.7f"%Ep)
