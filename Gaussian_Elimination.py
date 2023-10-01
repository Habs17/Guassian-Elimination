import numpy as np
import copy
#Gaussian Elimination
def gauss_naive(A,b):
    cp = A
    for k in range(n-1):
        for i in range(k+1, n):
            amut = A[i][k]/A[k][k]
            for j in range(k, n):
                cp[i][j] = A[i][j] - amut*A[k][j]
            b[i] = b[i] - amut*b[k]
            A[i] = cp[i]
    #begin backward substitution
    x = [0 for i in range(n)]
    x[n-1] =  b[n-1]/A[n-1][n-1]
    for i in range(n-2,-1,-1):
        sum = 0
        for k in range(i+1, n):
            sum = sum + A[i][k]*x[k]
        x[i] = b[i] - sum
        x[i] = x[i] / A[i][i]
    return x

def gauss_pivoted(A,b):
    for k in range(n-1):
        #partial pivoting:
        max = A[k][k]
        greatest = k
        for h in range(k,n):
            if abs(A[h][k]) > max:
                greatest = h
                max = abs(A[h][k])
        temp = copy.deepcopy(A[k])
        tempb = copy.deepcopy(b[k])
        A[k] = A[greatest] #row with the greatest absolute val goes on top
        A[greatest] = temp #previous top row moves to the row of maximum element
        b[k] = b[greatest]
        b[greatest] = tempb
        #end of partial pivoting
        cp = copy.deepcopy(A)
        for i in range(k+1, n):
            amut = A[i][k]/A[k][k]
            for j in range(k, n):
                cp[i][j] = A[i][j] - amut*A[k][j]
            b[i] = b[i] - amut*b[k]
            A[i] = copy.deepcopy(cp[i])
    #begin backward substitution
    x = [0 for i in range(n)]
    x[n-1] =  b[n-1]/A[n-1][n-1]
    for i in range(n-2,-1,-1):
        sum = 0
        for k in range(i+1, n):
            sum = sum + A[i][k]*x[k]
        x[i] = b[i] - sum
        x[i] = x[i] / A[i][i]
    return x

def residual_vector(A,x,b):
    r = [0 for i in range(n)]
    for i in range(n):
        asum = 0
        for j in range(n):
            asum += A[i][j]*x[j]
        r[i] = asum - b[j]
    return r

a = [[.729,.81,.9],[1.0,1.0,1.0],[1.331,1.21,1.1]]
b = [.6867,.8338,1.0]
n = len(a)

acp = copy.deepcopy(a)
bcp = copy.deepcopy(b)

part1 = gauss_naive(a,b)
print("Part1: ")
for i in part1:
    print(i, end = " ")
print("\n")
#end of part 1

Rand_Matrix = np.random.randn(n,n)
b2 = [0 for i in range(n)]

Rand_Copy = copy.deepcopy(Rand_Matrix)

Residual_Copy = copy.deepcopy(Rand_Matrix)

for j in range(n):
    b2[j] = j+1

b2cp = copy.deepcopy(b2)
residual_b2cp = copy.deepcopy(b2cp)
print("Part2: ")
part2 = gauss_naive(Rand_Matrix,b2)
for i in part2:
    print(i, end = " ")

print("\n\nPart3: ")
part3 = gauss_pivoted(acp,bcp)
for i in part3:
    print(i, end = " ")

print("\n\nPart4: ")
part4 = gauss_pivoted(Rand_Copy,b2cp)
for i in part4:
    print(i, end = " ")

residual2 = residual_vector(Residual_Copy,part2,residual_b2cp)
print("\n\nResidual Vector of part 2: ")
for i in residual2:
    print(i, end = " ")

residual4 = residual_vector(Residual_Copy,part4,residual_b2cp)
print("\n\nResidual Vector of part 4: ")
for i in residual4:
    print(i, end = " ")

Residual2_Norm = np.linalg.norm(residual2)
Residual4_Norm = np.linalg.norm(residual4)

print("\n\nResiduals:")
print(f"{Residual2_Norm:.20f}")
print(f"{Residual4_Norm:.20f}")
