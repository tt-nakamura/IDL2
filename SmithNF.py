from ZZlib import XGCD
from numpy import array, eye, nonzero, diag

def SmithNF(A):
    """ Smith Normal Form
    A = array of integers in shape(n,n)
    return d,U where
    d = diagonal components of SmithNF D of A
        in decreasing order, in shape(n,)
    U = numpy array of integers in shape(n,n)
        such that A = VDU, det(U)=1, det(V)=1
    reference: H. Cohen
      "A Course in Computational Algebraic Number Theory"
       Algorithm 2.4.14
    """
    A = array(A, dtype='object')
    U = eye(A.shape[0], dtype='object')
    k = A.shape[0]-1
    while k>0:
        for i in range(k-1,-1,-1):
            if A[i,k] == 0: continue
            d,x,y = XGCD(A[i,k], A[k,k])
            u,v = A[k,k]//d, A[i,k]//d
            A[k,k],A[i,k] = d,0
            A[i,:k],A[k,:k] = u*A[i,:k] - v*A[k,:k],\
                              x*A[i,:k] + y*A[k,:k]
        d = 0
        for j in range(k-1,-1,-1):
            if A[k,j] == 0: continue
            d,x,y = XGCD(A[k,j], A[k,k])
            u,v = A[k,k]//d, A[k,j]//d
            A[k,k],A[k,j] = d,0
            A[:k,j],A[:k,k] = u*A[:k,j] - v*A[:k,k],\
                              x*A[:k,j] + y*A[:k,k]
            U[j],U[k] = y*U[j] - x*U[k],\
                        v*U[j] + u*U[k]
        if d: continue
        i = nonzero(A[:k,:k] % A[k,k])[0]
        if len(i): A[k,:k] = A[i[0],:k]
        else: k -= 1

    return diag(A),U
