import math
class Matrix:
    def __init__(s,data): s.data=[list(r) for r in data];s.rows=len(data);s.cols=len(data[0])
    def __repr__(s):
        return '\n'.join(' '.join(f'{v:8.3f}' for v in r) for r in s.data)
    def __mul__(s,o):
        r=[[sum(s.data[i][k]*o.data[k][j] for k in range(s.cols)) for j in range(o.cols)] for i in range(s.rows)]
        return Matrix(r)
    def __add__(s,o): return Matrix([[s.data[i][j]+o.data[i][j] for j in range(s.cols)] for i in range(s.rows)])
    def transpose(s): return Matrix([[s.data[j][i] for j in range(s.rows)] for i in range(s.cols)])
    def det(s):
        if s.rows==1: return s.data[0][0]
        if s.rows==2: return s.data[0][0]*s.data[1][1]-s.data[0][1]*s.data[1][0]
        d=0
        for j in range(s.cols):
            minor=Matrix([[s.data[i][k] for k in range(s.cols) if k!=j] for i in range(1,s.rows)])
            d+=(-1)**j*s.data[0][j]*minor.det()
        return d
    def inverse(s):
        n=s.rows;aug=[s.data[i]+[1 if i==j else 0 for j in range(n)] for i in range(n)]
        for i in range(n):
            mx=max(range(i,n),key=lambda r:abs(aug[r][i]));aug[i],aug[mx]=aug[mx],aug[i]
            piv=aug[i][i];aug[i]=[v/piv for v in aug[i]]
            for j in range(n):
                if j!=i: f=aug[j][i];aug[j]=[aug[j][k]-f*aug[i][k] for k in range(2*n)]
        return Matrix([r[n:] for r in aug])
    @staticmethod
    def identity(n): return Matrix([[1 if i==j else 0 for j in range(n)] for i in range(n)])
def demo():
    m=Matrix([[1,2,3],[4,5,6],[7,8,10]])
    print(f"Det: {m.det()}");inv=m.inverse();r=m*inv
    print(f"M*M^-1 diagonal: {[round(r.data[i][i],6) for i in range(3)]}")
if __name__=="__main__": demo()
