from SmithNF import SmithNF
from numpy import zeros, flatnonzero, prod

def GroupGenerator(L):
    """ generating system of commutative group
    L = candidates of generators (1d array)
    return G,m where
      G = dictionary of (generator, order) pair
      m = order of the group (integer)
    group object must suport following funcitons:
      unit() = unit element of the group
      __mul__(a,b) = group operation a*b
      __pow__(a,e) = exponentiation a**e where e may be e<0
      __eq__(a,b): equality testing a==b
    reference: J. Buchmann and U. Vollmer
      "Binary Quadratic Forms" Algorithm 9.1
    """
    if not L: return [],1
    S = [type(L[0]).unit()]
    F,U,T,m = [],[],[[]],0
    while True:
        for i in S[m:]:
            while i in L: L.remove(i)
        m,l = len(S), len(U)
        if not L: break
        f,H = L[0],[]
        F.append(f)
        while f not in S:
            H.append(f)
            f *= L[0]
        i = S.index(f)
        V = [-T[i][j] for j in range(l)]
        V.append(len(H) + 1)
        U.append(V)
        for i,f in enumerate(H):
            for j in range(m):
                S.append(S[j] * f)
                T.append(T[j] + [i+1])
        for i in range(m):
            T[i].append(0)

    B = zeros((l,l), dtype='object')
    for i in range(l): B[i,:i+1] = U[i]
    D,U = SmithNF(B)
    i = flatnonzero(D>1)
    f = prod(F**U[i], axis=1)
    G = dict(zip(f, D[i]))
    return G,m
