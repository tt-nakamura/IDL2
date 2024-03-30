from IDL2 import *
from sympy import factorint

def factor(A):
    """ factorize A into prime ideals
    return dictionary of (prime, exponent) pair
    such that product of p**e == A
    A can be int or ZZ2
    reference: T. Takagi
      "Lectures on Elementary Number Theory"
       sections 43,44 (in Japanese)
    """
    if not isinstance(A, IDL2):
        A = principal(A)
    f = {}
    if A.isZero() or A.isUnit(): return f
    B = A.primitive()
    g = factorint(B.norm())
    h = factorint(A.content())

    for p in h:
        k = kron(p)
        if k<0: f[p] = h[p] # inert
        elif k>0: # split
            P = prime(p)
            Q = P.conj()
            f[P],f[Q] = h[p],h[p]
            if p in g:
                if P.divide(B): f[P] += g.pop(p)
                else:           f[Q] += g.pop(p)
        else: # ramify
            P = prime(p)
            f[P] = 2*h[p]
            if p in g: f[P] += g.pop(p)

    for p in g:
        P = prime(p)
        if not P.divide(B): P = P.conj()
        f[P] = g[p]

    return f

def IDL2FromNorm(n):
    """ ideal of given norm
    n = norm (integer)
    return list of all ideals that have norm n
    reference: T. Takagi
      "Lectures on Elementary Number Theory"
       section 50 (in Japanese)
    """
    f = factorint(abs(n))
    J = [IDL2.unit()]
    for p in f:
        l = len(J)
        k = kron(p)
        if k<0 and f[p]&1: return []
        if k<=0:
            q = p**(f[p]>>1)
            if f[p]&1: q *= prime(p)
            for i in range(l): J[i] *= q
        else:
            m = ((f[p] + 1)>>1)<<1
            A = prime(p)
            B = A.conj()
            if f[p] > 1:
                A2 = A*A
                B2 = B*B
                q = p
            C = [A,B] if f[p]&1 else [A2,B2]
            for j in range(2,m,2):
                C.append(C[-2] * A2)
                C.append(C[-2] * B2)
            for j in range(m-4, -1, -2):
                C[j] *= q
                C[j+1] *= q
                q *= p
            if f[p]&1 == 0: C.append(q)
            for D in C[1:]:
                J += [J[i]*D for i in range(l)]
            for i in range(l): J[i] *= C[0]
    return J
