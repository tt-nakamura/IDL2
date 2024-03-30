// uses NTL
//   http://www.shoup.net/ntl

#ifndef __IDL2ClassGroup_h__
#define __IDL2ClassGroup_h__

#include "IDL2.h"
#include<NTL/pair.h>

struct ICG2 : IDL2
// Ideal Class Group in Quadratic fields
{
    static NTL::ZZ amax; // Minkowski bound for a
    static void init(const NTL::ZZ& D);// set discriminant D
    // if D is not fundamental, raise rutime_error
    static void init(long D) { init(NTL::ZZ(D)); }
    static void ClassNum(NTL::ZZ& h);
    static void ClassNum(NTL::ZZ& h, const NTL::ZZ& D);
    // h = class number of new discriminant D
    // old value of D is restored on exit
};

struct ICG2Push : IDL2Push {
    NTL::ZZ amax;
    ICG2Push() : amax(ICG2::amax) {;}
    // save current values of D,D4,Dm4,S,W,W1,amax
    ~ICG2Push() { ICG2::amax=amax; }
    // restore old values when this object is destructed
};

inline long operator==(const ICG2& A, const ICG2& B)
{ return IsEquiv(A,B); }// test equivalence of A and B
inline long operator!=(const ICG2& A, const ICG2& B)
{ return !IsEquiv(A,B); }

inline long IsUnit(const ICG2& A)
{ return IsPrincipal(A); }// test principality of A

inline void inv(ICG2& B, const ICG2& A) { conj(B,A); }
// B = inverse class of A

void mul(ICG2& C, const ICG2& A, const ICG2& B);// C=A*B
void sqr(ICG2& B, const ICG2& A);// B=A*A

inline void operator*=(ICG2& B, const ICG2& A) { mul(B,B,A); }// B=B*A
                                                     
void power(ICG2& B, const ICG2& A, long n);// B=A^n (n may be n<0)

long generator(NTL::Vec<NTL::Pair<ICG2, long> >& G, long min=1);
// generator of class group
// G = vector of (generator, order) pair
// return class number = product of G[i].b
// if min==1, group representatives are
//    chosen from ideals of minimum value of a

long generator(NTL::Vec<NTL::Pair<ICG2, long> >& G,
               const NTL::ZZ& D, long min=1);
// generator of class group of new discriminant D
// old value of D is restored on exit

#endif // __IDL2ClassGroup_h__