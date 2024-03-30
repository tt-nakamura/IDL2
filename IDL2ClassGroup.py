from IDL2 import *
from ZZ2 import *
from ZZlib import power, SqrRoot, IsFundDisc
from sympy import primerange, divisors
from GroupGenerator import GroupGenerator
from fractions import gcd

class ICG2:
    """ Ideal Class Group of Quadratic Fields
    reference: T. Takagi
      "Lectures on Elementary Number Theory"
       section 45 (in Japanese)
    """
    amax = 0 # Minkowski bound for a
    class Push:
        def __init__(self):
            """ save current values of parameters """
            self.p = IDL2.Push()
            self.amax = ICG2.amax

        def __del__(self):
            """ restore parameters when this object is destructed """
            ICG2.amax = self.amax

    @staticmethod
    def init(D, push=False):
        """ set discriminant D
        if push is True, save old values of parameters to
          Push object, and return the object
        """
        if push: p = ICG2.Push()
        if not IsFundDisc(D):
            raise RuntimeError('D is not fundamental')
        IDL2.init(D)
        if D>0: ICG2.amax = IDL2.S >> 1
        else: ICG2.amax = SqrRoot(-3*D)//3
        if push: return p

    def __init__(A, ideal):
        A.i = ideal.reduce()

    def __hash__(A):
        return A.i.__hash__()

    def __repr__(A):
        return A.i.__repr__()

    def __mul__(A,B): # multiply and reduce
        return ICG2((A.i * B.i))

    def __eq__(A,B): # test equivalence
        return A.i.isEquiv(B.i)

    def __pow__(A,e): # A**e, e may be e<0
        if e<0: A,e = A.inv(), -e
        return power(A,e)

    @staticmethod
    def unit(): return ICG2(IDL2.unit()) # principal

    def inv(A): return ICG2(A.i.conj()) # conjugate

def ImClassNum():
    """ class number of imaginary quadratic fields
    reference: H. Cohen
      "A Course in Computational Algebraic Number Theory"
       Algorithm 5.3.5
    """
    h = 0
    for r in range((ICG2.amax>>1) + 1):
        b = (r<<1) + ZZ2.Dm4
        d = divisors(ZZ2(r,1).norm())
        n = (len(d) + 1)>>1
        for a,c in zip(d[:n], d[::-1][:n]):
            if a==b or a==c or b==0: h += 1
            elif a>b: h += 2
    return h

def ReClassNum():
    """ class number of real quadratic fields
    reference: J. Buchmann and U. Vollmer
      "Binary Quadratic Forms" section 6.17
    """
    L,h = [],0
    for r in range(IDL2.W1 + 1):
        s = IDL2.W1 - r
        d = divisors(ZZ2(r,1).norm())
        n = (len(d) + 1)>>1
        for a,c in zip(d[:n], d[::-1][:n]):
            if a>s:
                L.append(IDL2(a,r%a,1))
                if a!=c: L.append(IDL2(c,r%c,1))
    while L:
        A = L[0]
        while A in L:
            L.remove(A)
            A = A.cfrac()
        h += 1
    return h

def ClassNum(D=None):
    """ class number of quadratic fields
    if D is None, ZZ2.D is used for D
    else, ZZ2.D is set to D but
      old ZZ2.D is restored on exit.    
    """
    if D is not None:
        p = ICG2.init(D, True)
    if ZZ2.D < 0: return ImClassNum()
    else: return ReClassNum()

def generator(D=None, min=True):
    """ generating system of class group
    return G,h where
      G = dictionary of (generator, order) pair
      h = class number (order of group)
    if D is None, ZZ2.D is used for D
    else, ZZ2.D is set to D but
      old ZZ2.D is restored on exit.
    if min is True, group representatives are
      chosen from ideals of minimum value of a
    """
    if D is not None:
        t = ICG2.init(D, True)
    P = [ICG2(prime(p)) for p in
         primerange(0, ICG2.amax + 1) if kron(p) >= 0]
    G,h = GroupGenerator(P) # utilize general routine
    if not min: return G,h
    H = {}
    for A in G: # find minimum generators
        B = ICG2.unit()
        C = A
        for j in range(1,G[A]):
            B *= A
            if gcd(j, G[A]) > 1: continue
            if B.i.a < C.i.a: C = B
            if ZZ2.D < 0: continue
            # find a minimum in cfrac cycle for D>0
            E = B.i.cfrac()
            while E != B.i:
                if E.a < C.i.a: C = ICG2(E)
                E = E.cfrac()
        H[C] = G[A]
    return H,h
