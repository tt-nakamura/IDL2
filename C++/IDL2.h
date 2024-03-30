// uses NTL
//   http://www.shoup.net/ntl

#ifndef __IDL2_h__
#define __IDL2_h__

#include "ZZ2.h"

struct IDL2
// integral ideal in quadratic field aZ + bZ
// reference: T. Takagi
//   "Lectures on Elementary Number Theory"
//    section 42 (in Japanese)
{
    NTL::ZZ a;// integral basis a, a>0
    ZZ2 b; // integral basis b, 0<=b.x<a, b.y>0
    // b.y divides both a and b.x; a divides norm(b)
    static NTL::ZZ S;// floor(sqrt(D)) (used only for D>0)
    static NTL::ZZ W;// floor(w) (used only for D>0)
    static NTL::ZZ W1;// floor(-conj(w)) (used only for D>0)
    static void init(const NTL::ZZ& D);// set discriminant D
    // if D!=0,1 (mod 4) D is square, raise runtime_error
    static void init(long D) { init(NTL::ZZ(D)); }
    static long kron(const NTL::ZZ& p);// Kronecker symbol
    static long kron(long p) { return kron(NTL::ZZ(p)); }
    static long FundUnit(ZZ2& u);// get fundamenatal unit
    static long FundUnit(ZZ2& u, const NTL::ZZ& D);
    // u = fundamenatal unit for new discriminant D
    // old value of D is restored on exit
    IDL2() {;}
    IDL2(const ZZ2& a);// principal ideal (a)
    IDL2(const NTL::ZZ& a);// principal ideal (a)
    IDL2(long a);// principal ideal (a)
};

struct IDL2Push : ZZ2Push {
    NTL::ZZ S,W,W1;
    IDL2Push() : S(IDL2::S), W(IDL2::W), W1(IDL2::W1) {;}
    // save current values of D,D4,Dm4,S,W,W1
    ~IDL2Push() { IDL2::S=S; IDL2::W=W; IDL2::W1=W1; }
    // restore old values when this object is destructed
};

struct infra
// infrastructure for D>0
// reference: H. Cohen
//   "A Course in Computational Algebraic Number Theory"
//    section 5.8.1
{
    ZZ2 n;// numerator
    NTL::ZZ d;// denominator
    infra() : n(1), d(1) {;}
    void operator*=(const IDL2& A) { n *= A.b; d *= A.a; }
    void eval(ZZ2& a) { div(a,n,d); }// a = exp(distance)
};

void clear(IDL2& A);// A = zero ideal
void set(IDL2& A);// A = unit ideal
void conv(IDL2& A, const ZZ2& a);// principal ideal (a)
void conv(IDL2& A, long a);// principal ideal (a)
void conv(IDL2& A, const NTL::ZZ& a);// principal ideal (a)
void conj(IDL2& B, const IDL2& A);// conjugate, 0<=b.x<a

inline void norm(NTL::ZZ& n, const IDL2& A) { NTL::mul(n, A.a, A.b.y); }
// n = norm of ideal A; n>0

inline long IsZero(const IDL2& A) { return NTL::IsZero(A.a); }// test if A==0
inline long IsUnit(const IDL2& A) { return NTL::IsOne(A.a); }// test if A==1

void add(IDL2& C, const IDL2& A, const IDL2& B);// C = A+B
void add(IDL2& B, const IDL2& A, const ZZ2& a);// B = A+(a)
void add(IDL2& B, const IDL2& A, const NTL::ZZ& a);// B = A+(a)

// B = (a)+A
inline void add(IDL2& B, const NTL::ZZ& a, const IDL2& A) { add(B,A,a); }
inline void add(IDL2& B, const ZZ2& a, const IDL2& A) { add(B,A,a); }

void mul(IDL2& C, const IDL2& A, const IDL2& B);// C = A*B
void mul(IDL2& B, const IDL2& A, const NTL::ZZ& a);// B = A*(a)
void mul(IDL2& B, const IDL2& A, const ZZ2& a);// B = A*(a)
void sqr(IDL2& B, const IDL2& A);// B = A*A

// B = (a)*A
inline void mul(IDL2& B, const NTL::ZZ& a, const IDL2& A) { mul(B,A,a); }
inline void mul(IDL2& B, const ZZ2& a, const IDL2& A) { mul(B,A,a); }

long div(const IDL2& A, const IDL2& B);// test if B divides A
void power(IDL2& B, const IDL2& A, long n);// B = A**n (n>=0)

IDL2& operator+=(IDL2& A, const IDL2& B);// A = A+B
IDL2& operator+=(IDL2& A, const ZZ2& b);// A = A+(b)
IDL2& operator+=(IDL2& A, const NTL::ZZ& b);// A = A+(b)
IDL2& operator*=(IDL2& A, const IDL2& B);// A = A*B
IDL2& operator*=(IDL2& A, const NTL::ZZ& b);// A = A*(b)
IDL2& operator*=(IDL2& A, const ZZ2& b);// A = A*(b)

// test equality, assuming 0<=b.x<a
inline long operator==(const IDL2& A, const IDL2& B)
{ return A.a == B.a && A.b == B.b; }
inline long operator!=(const IDL2& A, const IDL2& B)
{ return A.a != B.a || A.b != B.b; }

// common factor of (a, b.x, b.y)
inline NTL::ZZ& content(IDL2& A) { return A.b.y; }
inline const NTL::ZZ& content(const IDL2& A) { return A.b.y; }

void SetPrime(IDL2& A, const NTL::ZZ& p);
// A = prime ideal above prime number p
// assume kronecker symbol (p/D) >= 0
inline void SetPrime(IDL2& A, long p) { SetPrime(A, NTL::ZZ(p)); }

void primitive(IDL2& B, const IDL2& A);
// B = primitive part of A (common factor removed)

long cfrac(IDL2& A, long red=0);
// one step of reduction algorithm
// (a.k.a. continued fraction expansion)
// assume A is primitive.
// normalize A, i.e., add mutiple of A.a to A.b such that
// if D>0 and A.a < sqrt(D), then -A.a < conj(A.b) < 0,
// else -A.a < 2*A.b.x + (D mod d) <= A.a.
// if red==1, test if A is reduced and
//   if A is reduced, return 1 with A unchanged
// else, let A.a = |norm(A.b)/A.a|, A.b = -conj(A.b),
//   let A.b.x %= A.a so that 0 <= A.b.x < A.a
// and return 0
// reference: J. Buchmann and U. Vollmer,
//   "Binary Quadratic Forms" sections 5.2, 6.1, 6.4

long cfrac(infra& i, IDL2& A, long red=0);
// after normalizing A, multiply A.b and |norm(A.b)/A.a|
// to i.n and i.d, respectively.
// if red==1, test if A is reduced and
//   if A is reduced, return 1 with A and i unchanged
// else, return (cfrac(A), False)

void reduce(IDL2&B, const IDL2& A);
// B = reduced ideal of A

long IsEquiv(const IDL2& A, const IDL2& B);
// test if A and B are equivalent ideals

long IsEquiv(ZZ2& a, const IDL2& A, const IDL2& B);
// if A and B are equivalent, output a such that A = (a)*B/norm(B)
// else a is unchanged

long IsPrincipal(const IDL2& A);
// test if A is a principal ideal

long IsPrincipal(ZZ2& a, const IDL2& A);
// if A is principal, output a such that A = (a)
// else a is unchanged

std::ostream& operator<<(std::ostream&, const IDL2&);// for printing

#endif // __IDL2_h__