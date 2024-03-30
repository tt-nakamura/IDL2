from sympy.ntheory.factor_ import core

def SqrRoot(n):
    """ floor(sqrt(n)) for huge integer n """
    x,y = n+1,n
    while x>y: x,y = y, (y + n//y)//2
    return x

def XGCD(a,b):
    """ extended Euclidean algorithm
    compute greatest common divisor d of a,b
    return d,s,t such that s*a + t*b == d
    """
    x,y = abs(a),abs(b)
    s,t,u,v = 1,0,0,1
    while y:
        q,r = divmod(x,y)
        x,y = y,r
        s,u = u,s-q*u
        t,v = v,t-q*v
    if a<0: s=-s
    if b<0: t=-t
    return x,s,t

def power(a,e): # a**e
    """ assume e>=0 """
    if e==0: return a.unit()
    n = (1<<(int(e).bit_length()-1))>>1
    b = a
    while n:
        b *= b
        if e&n: b *= a
        n>>=1
    return b

def conductor(D):
    """
    D = discriminant, D == 0 or 1 (mod 4)
    return f,d such that f**2 divide D,
      d = D/f**2 == 0 or 1 (mod 4) is fundamental disc.
    if d==1 (mod 4), d is square-free
    if d==0 (mod 4), d/4 is square-free and
                     d/4 == 2 or 3 (mod 4)
    """
    d = core(abs(D))
    if D<0: d=-d
    f = SqrRoot(D//d)
    if d&3 <= 1: return f,d
    return f>>1, d<<2

def IsFundDisc(D):
    """ test if conductor f == 1 """
    if D==0 or D==1 or D&3 > 1: return False
    return conductor(D)[0] == 1
