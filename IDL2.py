from ZZ2 import *
from ZZlib import power, SqrRoot
from numpy import array
from HermitNF import HermitNF
from sympy import sqrt_mod, legendre_symbol

class IDL2:
    """ Integral Ideals in Quadratic Field of discriminant D
    reference: T. Takagi
      "Lectures on Elementary Number Theory"
       section 42 (in Japanese)
    """
    S,W,W1 = 0,0,0 # these are used only for D>0
    # S = floor(sqrt(D)),
    # W = floor(w), W1 = floor(-conj(w))
    # where w = ((D mod 4) + sqrt(D))/2
    class Push:
        def __init__(self):
            """ save current values of parameters """
            self.p = ZZ2.Push()
            self.S = IDL2.S
            self.W = IDL2.W
            self.W1= IDL2.W1
        def __del__(self):
            """ restore parameters when this object is destructed """
            IDL2.S = self.S
            IDL2.W = self.W
            IDL2.W1= self.W1

    @staticmethod
    def init(D, push=False):
        """ set discriminant D
        if push is True, save old values of parameters to
          Push object, and return the object
        """
        if push: p = IDL2.Push()
        if D>0:
            S = SqrRoot(D)
            if S*S == D:
                raise RuntimeError("D must not be square")
        ZZ2.init(D)
        if D>0:
            IDL2.S = S
            IDL2.W = (S + ZZ2.Dm4)>>1
            IDL2.W1= (S - ZZ2.Dm4)>>1
        if push: return p

    def __init__(A, a=0, x=0, y=0):
        """ set components of ideal A = [a, x+yw]
        assume a>0, y>0, 0<=x<a; y divides both a and x;
               a divides norm(x+yw) """
        A.a = int(a)
        A.b = ZZ2(x,y)

    def __hash__(A):
        return A.a ^ A.b.__hash__()

    def __repr__(A):
        return '(' + str(A.a) + ' ' + str(A.b) + ')'

    def __add__(A,B): # A+B
        if not isinstance(B, IDL2):
            B = principal(B)
        H = HermitNF([[A.a, 0], [A.b.x, A.b.y],
                      [B.a, 0], [B.b.x, B.b.y]])
        return IDL2(H[0,0], H[1,0], H[1,1])

    def __mul__(A,B): # A*B
        if isinstance(B, int):
            b = abs(B)
            return IDL2(b*A.a, b*A.b.x, b*A.b.y)
        if isinstance(B, ZZ2):
            s = A.b * B
            H = [[A.a * B.x, A.a * B.y], [s.x, s.y]]
        else:
            s = A.b * B.b
            if A==B:
                H = [[A.a * A.a, 0],
                     [A.a * A.b.x, A.a * A.b.y],
                     [s.x, s.y]]
            else:
                H = [[A.a * B.a, 0],
                     [A.a * B.b.x, A.a * B.b.y],
                     [A.b.x * B.a, A.b.y * B.a],
                     [s.x, s.y]]
        H = HermitNF(H)
        return IDL2(H[0,0], H[1,0], H[1,1])

    def __pow__(A,e): # A**e, assuming e>=0
        return power(A,e)

    def __eq__(A,B): # A==B, assuming 0<=b.x<a
        if not isinstance(B, IDL2):
            B = principal(B)
        return A.a == B.a and A.b == B.b

    def __radd__(A,a):
        return A.__add__(a)

    def __rmul__(A,a):
        return A.__mul__(a)

    @staticmethod # unit ideal
    def unit(): return IDL2(1,0,1)

    def divide(A,B): return A + B == A # test if A divides B
    def norm(A): return A.a * A.b.y
    def isUnit(A): return A.a == 1
    def isZero(A): return A.a == 0
    def content(A): return A.b.y

    def conj(A): # conjugate of A
        b = -A.b.conj()
        if b.x < 0: b.x += A.a
        return IDL2(A.a, b.x, b.y)

    def primitive(A): # primitive part of A
        if A.b.y == 1: return A
        return IDL2(A.a//A.b.y, A.b.x//A.b.y, 1)

    def normalize(A):
        """ assume A = [a, b+w] is primitive.
        let b += a*q for some q, such that
        if D>0 and a<sqrt(D), -a < conj(b+w) < 0,
        else, -a < 2*b + (D mod 4) <= a;
        compute c = |norm(b+w)/a|, return [c, b+w];
        used only internally by cfrac
        reference: J. Buchmann and U. Vollmer,
          "Binary Quadratic Forms" sections 5.2, 6.1
        """
        if ZZ2.D > 0 and A.a <= A.S:
            b = A.W1 - (A.W1 - A.b.x)%A.a
        elif (A.b.x<<1) + ZZ2.Dm4 > A.a:
            b = A.b.x - A.a
        else: b = A.b.x
        c = abs(ZZ2(b,1).norm()//A.a)
        return IDL2(c,b,1)

    def cfrac(A, red=False):
        """ one step of reduction algorithm
        (a.k.a. continued fraction expansion)
        assume A is primitive.
        using b,c of normalize(A),
        return C = [c, k+w] where
          k+w == -conj(b+w) (mod c) and 0<=k<c
        if red is True, test if A is reduced and
          if A is reduced, return (A, True)
          else, return (A.cfrac(), False)
        if A.infra == [n,d] == [ZZ2, int],
          output C.infra = [n*(b+w), d*c]
        reference: J. Buchmann and U. Vollmer,
          "Binary Quadratic Forms" section 6.4
        """
        B = A.normalize()
        if red and IsReduced(A.a, B.b.x, B.a):
            return A, True
        C = B.conj()
        C.b.x %= C.a
        if hasattr(A, 'infra'):
            C.infra = A.infra * [B.b, B.a]
        if red: return C, False
        else: return C

    def reduce(A):
        """ apply cfrac until reduced """
        A,red = A.primitive(), False
        while not red: A,red = A.cfrac(True)
        return A

    def isEquiv(A,B):
        """ test if A and B are equivalent """
        A = A.reduce()
        B = B.reduce()
        if ZZ2.D < 0: return A==B
        C = A
        while A!=B:
            A = A.cfrac()
            if A==C: return False
        return True

    def isPrincipal(A):
        """ test if A is a principal ideal """
        A = A.reduce()
        if ZZ2.D < 0: return A.isUnit()
        B = A
        while not A.isUnit():
            A = A.cfrac()
            if A==B: return False
        return True

    def generator(A):
        """ test if A is a principal ideal and
        if so, return a such that A = (a) and a.y>0
        else, return None
        reference: T. Takagi
          "Lectures on Elementary Number Theory"
           section 51 (in Japanese)
        """
        B,red = A.primitive(), False
        B.infra = array([ZZ2(1), 1])
        while not red: B,red = B.cfrac(True)
        if ZZ2.D > 0:
            C = B
            while not B.isUnit():
                B = B.cfrac()
                if B==C: return None
        elif not B.isUnit(): return None
        b,d = B.infra
        if b.y < 0: b=-b
        return b/d*A.content()

def IsReduced(a,b,c):
    """ test if (a,b,c) is reduced, where
    (a,b,c) is binary quadratic form;
    assume (a,b,c) are normalized.
    reference: J. Buchmann and U. Vollmer,
      "Binary Quadratic Forms" sections 5.3, 6.2
    """
    if ZZ2.D < 0:
        return a<c or a==c and b>=0
    else: return a-b <= IDL2.W

def principal(a):
    """ principal ideal (a) """
    if not isinstance(a, ZZ2):
        return IDL2(a,0,a)
    b = a.times_w()
    H = HermitNF([[a.x, a.y], [b.x, b.y]])
    return IDL2(H[0,0], H[1,0], H[1,1])

def kron(p):
    """ Kronecker symbol for prime number p """
    if p>2: return legendre_symbol(ZZ2.D, p)
    elif ZZ2.Dm4 == 0: return 0
    elif ZZ2.D&7 == 1: return 1
    else: return -1

def prime(p):
    """ prime ideal above prime number p
    assume kron(p) >= 0
    reference: T. Takagi
      "Lectures on Elementary Number Theory"
       section 44 (in Japanese)
    """
    if p==2: return IDL2(2, (ZZ2.D&7)>>2, 1)
    s = sqrt_mod(ZZ2.D, p)
    if (s&1)^(ZZ2.Dm4): s = p-s
    return IDL2(p, s>>1, 1)

def FundUnit(D=None):
    """ Fundamental Unit e of quadratic field
    (a.k.a Pell's equation)
    return e, norm(e)
    if D is None, ZZ2.D is used for D
    else, ZZ2.D is set to D but
      old ZZ2.D is restored on exit.
    if D<-4, e = 1, else if D<0, e = w
    reference: T. Takagi
      "Lectures on Elementary Number Theory"
       sections 47,48 (in Japanese)
    """
    if D is not None:
        p = IDL2.init(D, True)
    if ZZ2.D < -4: return ZZ2(1), 1
    if ZZ2.D < 0: return ZZ2(0,1), 1
    A,s = IDL2.unit(), 1
    A.infra = array([ZZ2(1), 1])
    while True:
        A,s = A.cfrac(),-s
        if A.isUnit(): break
    b,d = A.infra
    return b/d, s
