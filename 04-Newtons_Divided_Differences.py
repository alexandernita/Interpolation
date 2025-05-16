# Coefficients of the nth Lagrange polynomial in the form 
# p(x)=\sum_{i=0}^n a_i \prod_{j=0}^{i}(x-x_j)

# The input data, points on the graph of f
x = [1.0, 1.3, 1.6, 1.9, 2.2]
y = [0.7651977, 0.6200860, 0.4554022, 0.2818186, 0.1103623]
n = len(x)

Q = [[] for i in range(n)]
Q[0] = y

# The interpolation point: 

xp = 1.5

# Implementing Newton's Divided Differences Method
for i in range(1,n):
    for j in range(i,n):
        xi = x[j]
        xij = x[j-i]
        Qij = Q[i]
        R = Q[i-1]
        s = (R[j-i+1]-R[j-i])/(xi-xij)
        Qij.append(s)
    Q[i]=Qij

# Displaying output

# Preamble

print("\n\tExample using Newton's Divided Difference Method:")
print("\n\tGiven %d"%n,"data points on the graph of the unknown f,")
print("\n\t x\t f(x)")
print("\t-------------------")
for i in range(len(x)):
    print("\t",format(x[i], ".1f"),"\t",format(y[i], ".7f"))
print("\n\tthe degree %d"%n,"Lagrange interpolating polynomial in the ")
print("\tform p(x)=sum prod f[x0,..,xk](x-x0)...(x-xk-1) has")
print("\tcoefficients f[x0,..,xk] given along the diagonal of the table:")

# Transpose the matrix Q and print out it's rows
for i in range(n):
    T = ["" for k in range(i)]
    Q[i]=T+Q[i]

S = [[] for i in range(n)]

for i in range(n):
    for j in range(n):
        S[i].append(Q[j][i])

print("\n\tf[xi]\t\tf[xi,xi+1]\tf[xi..xi+2]\tf[xi..xi+3]\tf[xi..xi+4]")
print("\t---------------------------------------------------------------------------")
for i in range(n):
    print("\t",format(S[i][0], ".7f") if type(S[i][0])==float else "",\
          "\t",format(S[i][1], ".7f") if type(S[i][1])==float else "",\
          "\t",format(S[i][2], ".7f") if type(S[i][2])==float else "",\
          "\t",format(S[i][3], ".7f") if type(S[i][3])==float else "",\
          "\t",format(S[i][4], ".7f") if type(S[i][4])==float else "")

# Evaluating at xp using the special form 
yp = S[0][0]

for i in range(1,n):
    p = 1
    for j in range(i):
        p = p*(xp-x[j])
    yp = yp + S[i][i]*p

print("\n\tUsing these, p(%0.1f)"%xp,"approximates f(%0.1f)"%xp,"by:")    
print("\t-------------------")
print("\t p(%0.1f)"%xp," = %0.7f"%yp)
