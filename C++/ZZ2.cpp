// uses NTL
//   http://www.shoup.net/ntl

#include "ZZ2.h"
#include<exception>
using namespace NTL;

ZZ ZZ2::D;// discriminant
ZZ ZZ2::D4;// D/4 or (D-1)/4 for D==0,1(mod 4)
long ZZ2::Dm4;// D mod 4

void ZZ2::init(const ZZ& D_) {// set discriminant 
    long d(D_%4);
    if(d>1) throw std::runtime_error("D must be 0 or 1 mod 4");
    sub(D4, D=D_, Dm4=d);
    D4 >>= 2;
}

ZZ2& ZZ2::operator=(const ZZ& a) { x=a; clear(y); return *this; }
ZZ2& ZZ2::operator=(long a) { x=a; clear(y); return *this; }

ZZ2& operator+=(ZZ2& b, const ZZ2& a) {// b += a
    b.x += a.x;
    b.y += a.y;
    return b;
}

ZZ2& operator-=(ZZ2& b, const ZZ2& a) {// b -= a
    b.x -= a.x;
    b.y -= a.y;
    return b;
}

ZZ2& operator*=(ZZ2& b, const ZZ2& a) { mul(b,b,a); return b; }// b*=a
ZZ2& operator*=(ZZ2& b, const ZZ& a) { mul(b,b,a); return b; }// b*=a

void clear(ZZ2& a) { clear(a.x); clear(a.y); }// a=0
void set(ZZ2& a) { set(a.x); clear(a.y); }// a=1
void set(ZZ2& a, const ZZ& x, const ZZ& y) { a.x = x; a.y = y; }// a=x+yw
void set(ZZ2& a, long x, long y) { a.x = x; a.y = y; }// a = x+yw
void set(ZZ2& a, const ZZ& x, long y) { a.x = x; a.y = y; }// a = x+yw
void set(ZZ2& a, long x, const ZZ& y) { a.x = x; a.y = y; }// a = x+yw

void norm(ZZ& x, const ZZ2& a) {// assume &x!=&a.x and &x!=&a.y
    ZZ s;
    if(ZZ2::Dm4) {
        add(x, a.x, a.y);
        x *= a.x;
    }
    else sqr(x, a.x);
    sqr(s, a.y);
    MulSubFrom(x, s, ZZ2::D4);
}

void negate(ZZ2& b, const ZZ2& a) {// b=-a
    negate(b.x, a.x);
    negate(b.y, a.y);
}

void conj(ZZ2& b, const ZZ2& a) {// replace sqrt(D) by -sqrt(D)
    if(&b!=&a) b.x = a.x;
    if(ZZ2::Dm4) b.x += a.y;
    negate(b.y, a.y);
}

void mul_w(ZZ2& b, const ZZ2& a) {// b = a*w
    if(&b==&a) {
        swap(b.x, b.y);
        if(ZZ2::Dm4) b.y += b.x;
        b.x *= ZZ2::D4;
    }
    else {
        mul(b.x, a.y, ZZ2::D4);
        b.y = a.x;
        if(ZZ2::Dm4) b.y += a.y;
    }
}

long IsUnit(const ZZ2& a) {// test if a==1
    ZZ s;
    norm(s,a);
    abs(s,s);
    return IsOne(s);
}

long IsAssoc(const ZZ2& a, const ZZ2& b) {// test if a/b is unit
    ZZ2 q;
    return div(q,a,b) && IsUnit(q);
}

void add(ZZ2& c, const ZZ2& a, const ZZ2& b) {// c=a+b
    add(c.x, a.x, b.x);
    add(c.y, a.y, b.y);
}

void sub(ZZ2& c, const ZZ2& a, const ZZ2& b) {// c=a-b
    sub(c.x, a.x, b.x);
    sub(c.y, a.y, b.y);
}

void mul(ZZ2& c, const ZZ2& a, const ZZ2& b) {// c=a*b
    ZZ s,t,u,v;
    mul(s, a.x, b.x);
    mul(t, a.y, b.y);
    add(u, a.x, a.y);
    add(v, b.x, b.y);
    mul(c.x, t, ZZ2::D4);
    mul(c.y, u, v);
    c.x += s;
    c.y -= s;
    if(ZZ2::Dm4 == 0) c.y -= t;
}

void mul(ZZ2& c, const ZZ2& a, const ZZ& b) {// c=a*b
    mul(c.x, a.x, b);
    mul(c.y, a.y, b);
}

void mul(ZZ2& c, const ZZ& b, const ZZ2& a) {// c=a*b
    mul(c.x, b, a.x);
    mul(c.y, b, a.y);
}

void sqr(ZZ2& b, const ZZ2& a) {// b = a*a
    ZZ s,t,u;
    sqr(s, a.x);
    sqr(t, a.y);
    mul(u, a.x, a.y);
    mul(b.x, t, ZZ2::D4);
    b.x += s;
    LeftShift(b.y, u, 1);
    if(ZZ2::Dm4) b.y += t;
}

long div(ZZ2& q, const ZZ2& a, const ZZ2& b) {// q = floor(a/b)
    ZZ s;
    ZZ2 c;
    norm(s,b);
    conj(c,b);
    c *= a;
    return divide(q.x, c.x, s) && divide(q.y, c.y, s);
}

long div(ZZ2& q, const ZZ2& a, const ZZ& b) {// q = floor(a/b)
    ZZ s,t;
    if(!divide(s, a.x, b) || !divide(t, a.y, b))
        return 0;
    q.x = s;
    q.y = t;
    return 1;
}

long div(const ZZ2& a, const ZZ2& b) {// test if b divides a
    ZZ s;
    ZZ2 c;
    norm(s,b);
    conj(c,b);
    c *= a;
    return divide(c.x, s) && divide(c.y, s);
}

long div(const ZZ2& a, const ZZ& b) {// test if b divides a
    return divide(a.x, b) && divide(a.y, b);
}

void power(ZZ2& b, const ZZ2& a, long n) {// b=a^n, assume n>=0
    if(n==0) { set(b); return; }
    if(&b==&a) { power(b, ZZ2(a), n); return; }
    long m(1<<(NumBits(n)-1));
    b = a;
    for(m>>=1; m; m>>=1) {
        sqr(b,b);
        if(n&m) b *= a;
    }
}


std::ostream& operator<<(std::ostream& s, const ZZ2& a) {
    s << '[' << a.x << ' ' << a.y << ']';
    return s;
}