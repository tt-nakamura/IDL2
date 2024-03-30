from IDL2 import *
from IDL2Factoring import IDL2FromNorm
from fractions import gcd
from ZZlib import SqrRoot,XGCD,conductor
from sympy import divisors

def primitive(F):
    """ remove gcd from binary quadratic form
    F = (a,b,c) (tuple of integers)
    return k,G where
      k = gcd(a,b,c) if a>0 else -gcd(a,b,c)
      G = (a/k, b/k, c/k) so that a/k > 0
    """
    a,b,c = F
    k = abs(gcd(gcd(a,b), c))
    if a<0: k=-k
    return k, (a//k, b//k, c//k)

def linsolv2(A,b):
    """ solve Ax=b
    A = 2x2 integer matrix
    b = 2x1 integer vector
    return n,d such that x=n/d
    n = 2x1 integer vector
    d = minimum positive integer scalar
    """
    d = A[0][0]*A[1][1] - A[0][1]*A[1][0]
    if d==0:
        raise RuntimeError("singular matrix in linsolv2")
    s = A[1][1]*b[0] - A[0][1]*b[1]
    t = A[0][0]*b[1] - A[1][0]*b[0]
    g = abs(gcd(gcd(s,t), d))
    if d<0: g=-g
    return (s//g, t//g), d//g
    
def SolveBQE(F,n,N=0):
    """ solve binary quadratic equation in integers
          ax^2 + bxy + cy^2 = n
    F = coefficients (a,b,c) (tuple of integers)
    n = RHS of equation (integer)
    N = maximum number of solutions to search for.
    if N<=0, only fundamental solusions are searched.
    return list of solutions (x,y) in integers.
    if no solution exists, empty list is returned
    reference: T. Takagi
      "Lectures on Elementary Number Theory
       sections 49,52 (in Japanese)
    """
    xy = []
    k,F = primitive(F)
    n,k = divmod(n,k)
    if k: return xy

    a,b,c = F
    D = b**2 - 4*a*c
    if D<0 and n<0: return xy # positive definite
    if D==0: # uninteresting case
        if n<0: return xy
        s = a*n
        r = SqrRoot(s)
        if r*r != s: return xy
        b>>=1
        d,s,t = XGCD(a,b)
        r,k = divmod(r,d)
        if k: return xy
        x,y = s*r, t*r
        xy.append((x,y))
        if N>1: a,b = a//d, b//d
        for i in range(1,N):
            j = (i+1)>>1
            if i&1: j=-j
            xy.append((x - j*b, y + j*a))
        return xy

    f,d = conductor(D)
    try: _p = IDL2.init(d, True)
    except:# D is square (uninteresting case)
        s = SqrRoot(D)
        A = [[a, (b+s)>>1],
             [a, (b-s)>>1]]
        d = divisors(a*n)
        for s in zip(d, d[::-1]):
            v,t = linsolv2(A,s)
            if t==1: xy.append(v)
        return xy
    # quadratic irrationality (interesting case)
    e,s = FundUnit()
    m = 1
    if f>1:
        e1 = e if s>0 else e*e # norm(e1)==1
        if d>0:
            q = e1
            while q.y % f: q,m = q*e1,m+1
        elif d==-4: m = 2
        elif d==-3: m = 3

    r = (b-f if d&1 else b)>>1
    A = principal(a) + ZZ2(r,f)
    for J in IDL2FromNorm(n):
        q = (A*J).generator()
        if q is None: continue
        if (q.norm()>0)^(n>0):
            if s>0: continue
            else: q *= e
        for _ in range(m):
            y,l = divmod(q.y, f)
            if l==0:
                x,l = divmod(q.x - r*y, a)
                if l==0: xy.append((x,y))
                if len(xy)==N: return xy
            if m>1: q *= e1

    if len(xy)==0 or N<1 or D<-4: return xy
    # search for associate solutions
    IDL2.init(D)
    e,s = FundUnit()
    if s<0: e *= e
    r = (b-1 if D&1 else b)>>1
    if   D==-4: N = min(N, 2*len(xy))
    elif D==-3: N = min(N, 3*len(xy))
    p,p1,e1,i = 1,1,e.conj(),0
    xy0 = xy.copy()
    while True:
        if i&1: p *= e
        else: p1 *= e1
        u = p if i&1 else p1
        for x,y in xy0:
            q = ZZ2(a*x + r*y, y)*u
            xy.append(((q.x - r*q.y)//a, q.y))
            if len(xy)==N: return xy
        i += 1
    return xy

def AssocSol(xy0, F, k):
    """ associate solutions of ax^2 + bxy + cy^2 = n
    Assume b^2 - 4ac > 0.
    xy0 = output of SolveBQE
    F = coefficients (a,b,c) (tuple of integers)
    k = exponent of unit e^k to be multiplied to
        quadratic integer that corresponds to sol0.
    return list of associate solutions (x,y).
    """
    xy = []
    _,F = primitive(F)
    a,b,c = F
    D = b**2 - 4*a*c
    if D<=0: return xy

    try: t = IDL2.init(D, True)
    except: return xy # D is square
    e,s = FundUnit()
    if s<0: e *= e
    if k<0: e = e.conj()
    p = e**abs(k)
    r = (b-1 if D&1 else b)>>1
    for x,y in xy0:
        q = ZZ2(a*x + r*y, y)*p
        xy.append(((q.x - r*q.y)//a, q.y))
    return xy
