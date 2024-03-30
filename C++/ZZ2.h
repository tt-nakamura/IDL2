// uses NTL
//   http://www.shoup.net/ntl

#ifndef __ZZ2_h__
#define __ZZ2_h__

#include<NTL/ZZ.h>

struct ZZ2
// Quadratic Integer x + yw
// if D==0(mod 4), w = sqrt(D)/2
// if D==1(mod 4), w = (1 + sqrt(D))/2
// reference: T. Takagi
//   "Lectures on Elementary Number Theory"
//    section 41 (in Japanese)
{
    NTL::ZZ x,y;// components of x+yw
    static NTL::ZZ D;// discriminant
    static NTL::ZZ D4;// D/4 or (D-1)/4 for D==0,1(mod 4)
    static long Dm4; // D mod 4 (0 or 1)
    static void init(const NTL::ZZ& D);// set discriminant D
    // if D!=0,1 (mod 4), raise runtime_error
    static void init(long D) { init(NTL::ZZ(D)); }
    ZZ2() {;}
    // set x,y components
    ZZ2(const NTL::ZZ& a, const NTL::ZZ& b) : x(a), y(b) {;}
    ZZ2(long a, long b) : x(a), y(b) {;}
    ZZ2(const NTL::ZZ& a) : x(a) {;} // y=0
    ZZ2(long a) : x(a) {;} // y=0
    ZZ2& operator=(const NTL::ZZ& a);// x=a, y=0
    ZZ2& operator=(long a);// x=a, y=0
};

struct ZZ2Push {
    NTL::ZZ D,D4;
    long Dm4;
    ZZ2Push() : D(ZZ2::D), D4(ZZ2::D4), Dm4(ZZ2::Dm4) {;}
    // save current values of D,D4,Dm4
    ~ZZ2Push() { ZZ2::D=D; ZZ2::D4=D4; ZZ2::Dm4=Dm4; }
    // restore D,D4,Dm4 when this object is destructed
};

void clear(ZZ2& a); // a=0+0w
void set(ZZ2& a); // a=1+0w
void set(ZZ2& a, const NTL::ZZ& x, const NTL::ZZ& y);// a=x+yw
void set(ZZ2& a, long x, const NTL::ZZ& y);// a=x+yw
void set(ZZ2& a, const NTL::ZZ& x, long y);// a=x+yw
void set(ZZ2& a, long x, long y);// a=x+yw
void negate(ZZ2& b, const ZZ2& a); // b=-a
void mul_w(ZZ2& b, const ZZ2& a);// b = a*(0+1w)

void conj(ZZ2& b, const ZZ2& a);
// b = conjugate of a, i.e., replace sqrt(D) by -sqrt(D)

void norm(NTL::ZZ& n, const ZZ2& a);
// n = a*conj(a); assume &n!=&a.x and &n!=&a.y

long IsAssoc(const ZZ2& a, const ZZ2& b);
// test if a/b is unit

inline long IsZero(const ZZ2& a)// test if a==0+0w
{ return NTL::IsZero(a.x) && NTL::IsZero(a.y); }
long IsUnit(const ZZ2& a);// test if |norm(a)|==1

void add(ZZ2& c, const ZZ2& a, const ZZ2& b);// c = a+b
void sub(ZZ2& c, const ZZ2& a, const ZZ2& b);// c = a-b
void mul(ZZ2& c, const ZZ2& a, const ZZ2& b);// c = a*b
void mul(ZZ2& c, const ZZ2& a, const NTL::ZZ& b);// c = a*b
void mul(ZZ2& c, const NTL::ZZ& a, const ZZ2& b);// c = a*b
void sqr(ZZ2& b, const ZZ2& a);// b = a*a

// test if b divides a
long div(ZZ2& q, const ZZ2& a, const ZZ2& b);// q = floor(a/b)
long div(ZZ2& q, const ZZ2& a, const NTL::ZZ& b);// q = floor(a/b)
long div(const ZZ2& a, const ZZ2& b);
long div(const ZZ2& a, const NTL::ZZ& b);

void power(ZZ2& b, const ZZ2& a, long e);// b = a**e (e>=0)

inline long operator==(const ZZ2& a, const ZZ2& b)// test if a==b
{ return a.x == b.x && a.y == b.y; }
inline long operator!=(const ZZ2& a, const ZZ2& b)// test if a!=b
{ return a.x != b.x || a.y != b.y; }

ZZ2& operator+=(ZZ2& a, const ZZ2& b);// a = a+b
ZZ2& operator-=(ZZ2& a, const ZZ2& b);// a = a-b
ZZ2& operator*=(ZZ2& a, const ZZ2& b);// a = a*b
ZZ2& operator*=(ZZ2& a, const NTL::ZZ& b);// a = a*b

std::ostream& operator<<(std::ostream&, const ZZ2&);// for printing

#endif // __ZZ2_h__