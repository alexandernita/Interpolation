3# Horner's Method:  if e.g. p(x) = 2x^4 - 3x^2 + 3x - 4, what is p(-2)? Except we'll make it interactive.


n = int(input("\n What is the degree of the polynomial p(x)?\n\n"))
prompt1="\n Enter the coefficients of p(x), in decreasing order starting with a_n, a_n-1,...\n\n"
print(prompt1)

# Make the lists with coefficient names a_k and b_k
D = []
for i in range(0,n+1):
    D.append("a"+str(n-i)+" = ")
    
C = []
for i in range(0,n+1):
    C.append("b"+str(n-i)+" = ") 

# Use the D-list to prompt inputs for the numerical coefficients
A = []
l = len(D)

for i in range(0,l):
    d = D[i]
    A.append(float(input(d)))

# Prompt to enter the specific x-value x0                
x0 = float(input("\n Enter the x-value at which you would like to evaluate p(x).  \n x0 = \n\n"))

# Use the list A and x0 to plug into Horner's method, and generate the b_k's
B = [A[0]]

for i in range(0,l-1):
    b = B[i]
    B.append(A[i+1]+b*x0)

# Verify p(x0) by direct calculation and print out q(x) coefficients
p = 0
for i in range(0,n+1):
    p = p+A[n-i]*x0**i
print("")
print("The qotient polynomial q(x) has coefficients")
print("")
for x,y in zip(C,B):
    print(x,y)
print("")
print("And we remember, of course, that b_0 = p(x_0).")
print("")
print("Let's verify by direct computation:  p(x_0) = ",p)