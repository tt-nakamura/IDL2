// uses NTL
//   http://www.shoup.net/ntl

#ifndef __IDL2Factoring_h__
#define __IDL2Factoring_h__

#include "IDL2.h"
#include "ZZFactoring.h"

void factor(NTL::Vec<NTL::Pair<IDL2, long> >& F, const IDL2& A);
// prime decomposition of ideal A
// F = vector of (prime ideal, exponent) pair
//     so that product of F[i].a**F[i].b = A

void IDL2FromNorm(NTL::Vec<IDL2>& J, const NTL::ZZ& n);
// J = vector of ideals that have given norm |n|

void mul(IDL2& A, NTL::Vec<NTL::Pair<IDL2, long> >& F);
// A = product of F[i].a**F[i].b

#endif // __IDL2Factoring_h__