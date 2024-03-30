from ZZlib import power

class ZZ2:
    """ Integers in Quadratic Field of discriminant D
    reference: T. Takagi
      "Lectures on Elementary Number Theory"
       section 41 (in Japanese)
    """
    D,Dm4,D4 = 0,0,0
    # D = discriminant of quadratic field
    # Dm4 = D mod 4 (= 0 or 1)
    # D4 = D/4 if D==0 (mod 4), (D-1)/4 if D==1 (mod 4)
    class Push:
        def __init__(self):
            """ save current values of D,Dm4,D4 """
            self.D = ZZ2.D
            self.Dm4 = ZZ2.Dm4
            self.D4 = ZZ2.D4

        def __del__(self):
            """ restore D,Dm4,D4 when this object is destructed """
            ZZ2.D = self.D
            ZZ2.Dm4 = self.Dm4
            ZZ2.D4 = self.D4

    @staticmethod
    def init(D, push=False):
        """ set discriminant D
        if push is True, save old values of D,Dm4,D4 to
          Push object, and return the object
        """
        if push: p = ZZ2.Push()
        Dm4 = D&3
        if Dm4 > 1:
            raise RuntimeError('D must be 0 or 1(mod 4)')
        ZZ2.D = D
        ZZ2.Dm4 = Dm4
        ZZ2.D4 = (D - Dm4)>>2
        if push: return p

    def __init__(a, x=0, y=0):
        """ set x,y components of integer a = x + y*w
        if D==0 (mod 4), w = sqrt(D)/2
        if D==1 (mod 4), w = (1 + sqrt(D))/2
        """
        a.x = int(x)
        a.y = int(y)

    def __repr__(a):
        return '(' + str(a.x) + ' ' + str(a.y) + ')'

    def __hash__(a):
        return a.x ^ a.y

    def __neg__(a): # -a
        return ZZ2(-a.x, -a.y)

    def __add__(a,b): # a+b
        if isinstance(b, ZZ2):
            return ZZ2(a.x + b.x, a.y + b.y)
        else:
            return ZZ2(a.x + b, a.y)

    def __sub__(a,b): # a-b
        if isinstance(b, ZZ2):
            return ZZ2(a.x - b.x, a.y - b.y)
        else:
            return ZZ2(a.x - b, a.y)

    def __mul__(a,b): # a*b
        if isinstance(b, ZZ2):
            s = a.x * b.x
            t = a.y * b.y
            u = (a.x + a.y)*(b.x + b.y) - s
            if a.Dm4 == 0: u -= t
            return ZZ2(s + t*a.D4, u)
        else:
            return ZZ2(b*a.x, b*a.y)

    def __truediv__(a,b): # a/b (floor)
        if isinstance(b, ZZ2):
            a *= conj(b)
            b = norm(b)
        return ZZ2(a.x//b, a.y//b)

    def __pow__(a,e): # a**e
        """ assume e>=0 """
        return power(a,e)

    def __eq__(a,b): # a==b
        if isinstance(b, ZZ2):
            return a.x == b.x and a.y == b.y
        else:
            return a.x == b and a.y == 0

    def __radd__(a,b):
        return ZZ2(b + a.x, a.y)

    def __rsub__(a,b):
        return ZZ2(b - a.x, -a.y)

    def __rmul__(a,b):
        return ZZ2(b*a.x, b*a.y)

    def __rtruediv__(a,b):
        return a.conj()*b/a.norm()

    def norm(a): # a*conj(a)
        if a.Dm4:
            return a.x*(a.x + a.y) - a.y**2*a.D4
        else: return a.x**2 - a.y**2*a.D4

    def conj(a): # conjugate of a
        """ replace sqrt(D) with -sqrt(D) """
        if a.Dm4:
            return ZZ2(a.x + a.y, -a.y)
        else: return ZZ2(a.x, -a.y)

    def isUnit(a):
        return abs(a.norm()) == 1

    def times_w(a): # a*w
        if a.Dm4:
            return ZZ2(a.y*a.D4, a.x + a.y)
        else: return ZZ2(a.y*a.D4, a.x)
