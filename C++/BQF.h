// uses NTL
//   http://www.shoup.net/ntl

#ifndef __BQF_h__
#define __BQF_h__

#include<NTL/mat_ZZ.h>
#include "ZZ2.h"

struct BQF {// binary quadratic form ax^2 + bxy + cy^2
    NTL::ZZ a,b,c;
    BQF() {;}
    BQF(const NTL::ZZ& a_,
        const NTL::ZZ& b_,
        const NTL::ZZ& c_) : a(a_),b(b_),c(c_) {;}
    BQF(long a_,
        long b_,
        long c_) : a(a_),b(b_),c(c_) {;}
};

void set(BQF& f, long a, long b, long c);// f = (a,b,c)

void primitive(NTL::ZZ& k, BQF& g, const BQF& f);
// primitive part of f = (a,b,c)
// k = gcd(a,b,c) if a>0 else -gcd(a,b,c)
// g = (a/k, b/k, c/k) so that a/k > 0

void discriminant(NTL::ZZ& D, const BQF& f);
// D = b^2 - 4ac

void eval(NTL::ZZ& z, const BQF& f, const NTL::ZZ& x, const NTL::ZZ& y);
// z = ax^2 + bxy + cy^2
inline void eval(NTL::ZZ& z, const BQF& f, long x, long y)
{ eval(z,f, NTL::ZZ(x), NTL::ZZ(y)); }

long SolveBQE(NTL::mat_ZZ& xy, const BQF& F, const NTL::ZZ& N, long M=0);
// solve binary quadratic equation ax^2 + bxy + cy^2 = N
// F = (a,b,c), N = RHS of equation
// M = maximum number of solutions to search for.
// If M<=0, only fundamental solutions are searched.
// xy[i][0] = x component of i-th solution (0<=i<m).
// xy[i][1] = y component of i-th solution (0<=i<m).
//   where m is the number of solutions found (0<=m<=M).
// return sign (=-1,0,1) of D = b^2 - 4ac.
// note that if D>=0, number of solutions is either inf or 0.

void AssocSol(NTL::mat_ZZ& xy, const NTL::mat_ZZ& xy0,
              const BQF& f, long k=1);
// xy = associate solutions of ax^2 + bxy + cy^2 = N
// Assume b^2 - 4ac > 0.
// xy0 = output of SolveBQE, F = coeffients a,b,c
// k = exponent of unit e^k to be multiplied to
//     quadratic integer that corresponds to xy0

#endif // __BQF_h__