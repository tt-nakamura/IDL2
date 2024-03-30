from numpy import array
from ZZlib import XGCD

def HermitNF(A):
    """ Hermite Normal Form
    A = array of integers in shape(m,n)
    return lower triangular matrix in shape(n,n)
    referece: H. Cohen
      "A Course in Computational Algebraic Number Theory"
       Algorithm 2.4.5
    """
    A = array(A, dtype='object')
    if A.shape[0] < A.shape[1]:
        raise RuntimeError("wrong shape in HermitNF")
    k = A.shape[0]-1
    for l in range(A.shape[1]-1, -1, -1):
        for i in range(k-1,-1,-1):
            if A[i,l] == 0: continue
            d,u,v = XGCD(A[i,l], A[k,l])
            s,t = A[k,l]//d, A[i,l]//d
            A[k,l],A[i,l] = d,0
            A[i,:l],A[k,:l] = s*A[i,:l] - t*A[k,:l],\
                              u*A[i,:l] + v*A[k,:l]
        if A[k,l] == 0: continue
        if A[k,l] < 0: A[k,:l+1] = -A[k,:l+1]
        for i in range(k+1, A.shape[0]):
            q, A[i,l] = divmod(A[i,l], A[k,l])
            A[i,:l] -= q*A[k,:l]
        k -= 1

    return A[-A.shape[1]:]
