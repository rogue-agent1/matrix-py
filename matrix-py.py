#!/usr/bin/env python3
"""Matrix operations: multiply, transpose, determinant, inverse, LU decomposition."""
import sys,math

def zeros(r,c):return[[0]*c for _ in range(r)]
def identity(n):return[[1 if i==j else 0 for j in range(n)]for i in range(n)]
def transpose(A):return[[A[j][i]for j in range(len(A))]for i in range(len(A[0]))]
def matmul(A,B):
    r,k,c=len(A),len(A[0]),len(B[0])
    return[[sum(A[i][m]*B[m][j]for m in range(k))for j in range(c)]for i in range(r)]
def det(A):
    n=len(A)
    if n==1:return A[0][0]
    if n==2:return A[0][0]*A[1][1]-A[0][1]*A[1][0]
    d=0
    for j in range(n):
        minor=[[A[i][k]for k in range(n)if k!=j]for i in range(1,n)]
        d+=(-1)**j*A[0][j]*det(minor)
    return d
def inverse(A):
    n=len(A);aug=[A[i]+identity(n)[i]for i in range(n)]
    for i in range(n):
        mx=max(range(i,n),key=lambda r:abs(aug[r][i]))
        aug[i],aug[mx]=aug[mx],aug[i]
        piv=aug[i][i]
        if abs(piv)<1e-12:raise ValueError("Singular")
        aug[i]=[x/piv for x in aug[i]]
        for j in range(n):
            if j!=i:
                f=aug[j][i];aug[j]=[aug[j][k]-f*aug[i][k]for k in range(2*n)]
    return[row[n:]for row in aug]
def lu(A):
    n=len(A);L=identity(n);U=[row[:]for row in A]
    for i in range(n):
        for j in range(i+1,n):
            if abs(U[i][i])<1e-12:continue
            f=U[j][i]/U[i][i];L[j][i]=f
            U[j]=[U[j][k]-f*U[i][k]for k in range(n)]
    return L,U

def main():
    if len(sys.argv)>1 and sys.argv[1]=="--test":
        A=[[1,2],[3,4]];B=[[5,6],[7,8]]
        assert matmul(A,B)==[[19,22],[43,50]]
        assert transpose(A)==[[1,3],[2,4]]
        assert det(A)==-2
        inv=inverse(A)
        prod=matmul(A,inv)
        assert all(abs(prod[i][j]-(1 if i==j else 0))<1e-10 for i in range(2)for j in range(2))
        L,U=lu([[2,1],[4,3]])
        assert abs(L[1][0]-2)<1e-10
        assert det([[1,2,3],[4,5,6],[7,8,9]])==0  # singular
        print("All tests passed!")
    else:
        A=[[1,2],[3,4]];print(f"det={det(A)}")
if __name__=="__main__":main()
