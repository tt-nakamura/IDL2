// uses NTL
//   http://www.shoup.net/ntl

#include "BQF.h"
using namespace NTL;

void set(BQF& f, long a, long b, long c)
{ f.a = a; f.b = b; f.c = c; }

void primitive(ZZ& c, BQF& g, const BQF& f)
// c = content (c>0), g = primitive part of f
{
    GCD(c, f.a, f.b);
    GCD(c, c, f.c);
    if(sign(f.a) < 0) negate(c,c);
    div(g.a, f.a, c);
    div(g.b, f.b, c);
    div(g.c, f.c, c);
}

void discriminant(ZZ& D, const BQF& f)
// D = b^2 - 4ac
{
    sqr(D, f.b);
    MulSubFrom(D, f.a<<2, f.c);
}

void eval(ZZ& z, const BQF& f, const ZZ& x, const ZZ& y)
{
    ZZ s,t;
    sqr(s,x); s *= f.a;
    sqr(t,y); MulAddTo(s, t, f.c);
    mul(t,x,y); t *= f.b;
    add(z,s,t);
}