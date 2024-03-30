// uses NTL
//   http://www.shoup.net/ntl

#include "IDL2.h"
#include<NTL/mat_ZZ.h>
#include<exception>
using namespace NTL;

ZZ IDL2::S;// floor(sqrt(D)) (used only for D>0)
ZZ IDL2::W;// floor(w) (used only for D>0)
ZZ IDL2::W1;// floor(-conj(w)) (used only for D>0)

void HermitNF(mat_ZZ&, const mat_ZZ&);

void IDL2::init(const ZZ& D) {// set discriminant
    if(sign(D) < 0) { ZZ2::init(D); return; }
    ZZ s,t; SqrRoot(s,D); sqr(t,s);
    if(t==D) throw std::runtime_error("D must not be square");
    ZZ2::init(D);
    add(W, S=s, ZZ2::Dm4); W >>= 1;
    sub(W1, S, ZZ2::Dm4); W1 >>= 1;// -conj(w)
}

IDL2::IDL2(const ZZ2& a) { conv(*this, a); }// principal ideal (a)
IDL2::IDL2(const ZZ& a) { conv(*this, a); }// principal ideal (a)
IDL2::IDL2(long a) { conv(*this, a); }// principal ideal (a)

void clear(IDL2& A) { clear(A.a); clear(A.b); }// A=(0)
void set(IDL2& A) { set(A.a); clear(A.b.x); set(A.b.y); }// A=(1)
void conv(IDL2& A, long a) { A.a = A.b.y = a; clear(A.b.x); }// A=(a)
void conv(IDL2& A, const ZZ& a) { A.a = A.b.y = a; clear(A.b.x); }// A=(a)

void conj(IDL2& B, const IDL2& A) {// B = conjugate of A, 0<=B.b.y<B.a
    if(&B!=&A) B=A;
    conj(B.b, B.b);
    negate(B.b, B.b);
    if(sign(B.b.x) < 0) B.b.x += B.a;
}

void conv(IDL2& A, const ZZ2& a) {// principal ideal A = (a)
    if(IsZero(a.y)) { conv(A, a.x); return; }
    ZZ2 b;
    mat_ZZ W;
    W.SetDims(2,2);
    mul_w(b,a);
    W[0][0] = a.x;
    W[0][1] = a.y;
    W[1][0] = b.x;
    W[1][1] = b.y;
    HermitNF(W,W);
    A.a = W[0][0];
    A.b.x = W[1][0];
    A.b.y = W[1][1];
}

void add(IDL2& C, const IDL2& A, const IDL2& B) {// C = A+B
    if(IsZero(A)) { if(&C!=&B) C=B; return; }
    if(IsZero(B)) { if(&C!=&A) C=A; return; }
    mat_ZZ W;
    W.SetDims(4,2);
    W[0][0] = A.a;
    W[1][0] = A.b.x;
    W[1][1] = A.b.y;
    W[2][0] = B.a;
    W[3][0] = B.b.x;
    W[3][1] = B.b.y;
    HermitNF(W,W);
    C.a = W[0][0];
    C.b.x = W[1][0];
    C.b.y = W[1][1];
}

void add(IDL2& B, const IDL2& A, const ZZ& a)
{ IDL2 C; conv(C,a); add(B,A,C); }// B = A+(a)

void add(IDL2& B, const IDL2& A, const ZZ2& a)
{ IDL2 C; conv(C,a); add(B,A,C); }// B = A+(a)

void mul(IDL2& C, const IDL2& A, const IDL2& B) {// C = A*B
    if(IsZero(A) ||
       IsZero(B)) { clear(C); return; }
    if(IsUnit(A)) { if(&C!=&B) C=B; return; }
    if(IsUnit(B)) { if(&C!=&A) C=A; return; }
    ZZ2 s;
    mat_ZZ W;
    W.SetDims(4,2);
    mul(W[0][0], A.a, B.a);
    mul(W[1][0], A.a, B.b.x);
    mul(W[1][1], A.a, B.b.y);
    mul(W[2][0], A.b.x, B.a);
    mul(W[2][1], A.b.y, B.a);
    mul(s, A.b, B.b);
    W[3][0] = s.x;
    W[3][1] = s.y;
    HermitNF(W,W);
    C.a = W[0][0];
    C.b.x = W[1][0];
    C.b.y = W[1][1];
}

void mul(IDL2& B, const IDL2& A, const ZZ2& a) {// B = A*(a)
    if(IsZero(A)) { clear(B); return; }
    if(IsUnit(A)) { conv(B,a); return; }
    if(IsZero(a.y)) { mul(B, A, a.x); return; }
    ZZ2 s;
    mat_ZZ W;
    W.SetDims(2,2);
    mul(W[0][0], A.a, a.x);
    mul(W[0][1], A.a, a.y);
    mul(s, A.b, a);
    W[1][0] = s.x;
    W[1][1] = s.y;
    HermitNF(W,W);
    B.a = W[0][0];
    B.b.x = W[1][0];
    B.b.y = W[1][1];
}

void mul(IDL2& B, const IDL2& A, const ZZ& a) {// B = A*(a)
    ZZ b; abs(b,a);
    mul(B.a, A.a, b);
    mul(B.b, A.b, b);
}

void sqr(IDL2& B, const IDL2& A) {// B = A*A
    if(IsZero(A)) { clear(B); return; }
    if(IsUnit(A)) { set(B); return; }
    ZZ2 s;
    mat_ZZ W;
    W.SetDims(3,2);
    sqr(W[0][0], A.a);
    mul(W[1][0], A.a, A.b.x);
    mul(W[1][1], A.a, A.b.y);
    sqr(s, A.b);
    W[2][0] = s.x;
    W[2][1] = s.y;
    HermitNF(W,W);
    B.a = W[0][0];
    B.b.x = W[1][0];
    B.b.y = W[1][1];
}

long div(const IDL2& A, const IDL2& B) {// test if B divides A
    IDL2 C;
    add(C,A,B);
    return B==C;
}

void power(IDL2& B, const IDL2& A, long n) {// B = A^n, assume n>=0
    if(n==0 || IsUnit(A)) { set(B); return; }
    if(&B==&A) { power(B, IDL2(A), n); return; }
    long m(1<<(NumBits(n)-1));
    B = A;
    for(m>>=1; m; m>>=1) {
        sqr(B,B);
        if(n&m) B *= A;
    }
}

IDL2& operator+=(IDL2& B, const IDL2& A) { add(B,B,A); return B; }// B += A
IDL2& operator+=(IDL2& A, const ZZ& a)   { add(A,A,a); return A; }// A += (a)
IDL2& operator+=(IDL2& A, const ZZ2& a) { add(A,A,a); return A; }// A += (a)
IDL2& operator*=(IDL2& B, const IDL2& A) { mul(B,B,A); return B; }// B *= A
IDL2& operator*=(IDL2& A, const ZZ& a)   { mul(A,A,a); return A; }// A *= (a)
IDL2& operator*=(IDL2& A, const ZZ2& a) { mul(A,A,a); return A; }// A *= (a)

long IDL2::kron(const ZZ& p) {// Kronecker symbol
    if(IsOdd(p)) return Jacobi(ZZ2::D%p, p);
    if(ZZ2::Dm4 == 0) return 0;
    return (ZZ2::D%8 > 1 ? -1:1);
}

void primitive(IDL2& B, const IDL2& A) {// B = primitive part of A
    if(&B!=&A) B=A;
    if(IsZero(A) || IsOne(A.b.y)) return;
    B.a /= A.b.y;
    B.b.x /= A.b.y;
    set(B.b.y);
}

void SetPrime(IDL2& A, const ZZ& p)
// A = prime ideal (p) for prime number p.
// assume Kronecker(D,p) >= 0.
{
    A.a = p;
    if(IsOdd(p)) {
        ZZ s;
        rem(s, ZZ2::D, p);
        SqrRootMod(s,s,p);
        if(IsOdd(s)^ZZ2::Dm4) sub(s,p,s);
        set(A.b, s>>=1, 1);
    }
    else set(A.b, ZZ2::D%8>>2, 1);
}

void normalize(IDL2& B, const IDL2& A)
// assume A is primitive and A!=0.
// private function, used only internally by cfrac
{
    ZZ n;
    if(sign(ZZ2::D) < 0 || A.a > IDL2::S) {
        LeftShift(n, A.b.x, 1);
        if(ZZ2::Dm4) n++;
        if(n > A.a) sub(B.b.x, A.b.x, A.a);
        else if(&B!=&A) B.b.x = A.b.x;
    }
    else {
        sub(B.b.x, IDL2::W1, A.b.x);
        sub(B.b.x, IDL2::W1, B.b.x%=A.a);
    }
    set(B.b.y);
    norm(n, B.b);
    div(B.a, n, A.a);
    abs(B.a, B.a);
}

long IsReduced(const ZZ& a, const ZZ& b, const ZZ& c)
// test if binary quadratic form (a,b,c) is reduced.
// assume A is normalized.
// c is used only for D<0.
// private function, used only internally by cfrac
{
    if(sign(ZZ2::D) < 0)
        return a<c || a==c && sign(b) >= 0;
    else return a-b <= IDL2::W;
}

long cfrac(IDL2& A, long red)
// continued fractio expansion.
// assume A is primitive and A!=0
{
    IDL2 B;
    normalize(B,A);
    if(red && IsReduced(A.a, B.b.x, B.a))
        return 1;
    conj(A,B);
    A.b.x %= A.a;
    return 0;
}

long cfrac(infra& i, IDL2& A, long red)
// continued fractio expansion.
// assume A is primitive and A!=0
{
    IDL2 B;
    normalize(B,A);
    if(red && IsReduced(A.a, B.b.x, B.a))
        return 1;
    i *= B;
    conj(A,B);
    A.b.x %= A.a;
    return 0;
}

void reduce(IDL2& B, const IDL2& A)
// B = reduction of primitive part of A
// assume A!=0
{
    primitive(B,A);
    while(cfrac(B,1) == 0) {;}
}

long IsEquiv(const IDL2& A, const IDL2& B)
// test if A and B are equivalent
{
    if(IsZero(A) || IsZero(B)) return 1;
    IDL2 C,D;
    reduce(C,A);
    reduce(D,B);
    if(C==D) return 1;
    if(sign(ZZ2::D) < 0) return 0;
    IDL2 E(C);
    for(cfrac(C); C!=D; cfrac(C))
        if(C==E) return 0;
    return 1;
}

long IsEquiv(ZZ2& a, const IDL2& A, const IDL2& B)
// test if A and B are equivalent
// and output a such that A = (a)*B/norm(B)
{
    IDL2 C;
    conj(C,B);
    C *= A;
    return IsPrincipal(a,C);
}

long IsPrincipal(const IDL2& A)
// test if A is principal ideal
{
    if(IsZero(A)) return 1;
    IDL2 B;
    reduce(B,A);
    if(IsUnit(B)) return 1;
    if(sign(ZZ2::D) < 0) return 0;
    IDL2 C(B);
    for(cfrac(B); !IsUnit(B); cfrac(B))
        if(B==C) return 0;
    return 1;
}

long IsPrincipal(ZZ2& a, const IDL2& A)
// test if A = (a) for some a
{
    if(IsZero(A)) { clear(a); return 1; }
    IDL2 B;
    infra i;
    primitive(B,A);
    while(cfrac(i,B,1) == 0) {;}
    if(!IsUnit(B)) {
        if(sign(ZZ2::D) < 0) return 0;
        IDL2 C(B);
        for(cfrac(i,B); !IsUnit(B); cfrac(i,B))
            if(B==C) return 0;
    }
    i.eval(a);
    a *= content(A);
    if(sign(a.y) < 0) negate(a,a);// upper half
    return 1;
}

long IDL2::FundUnit(ZZ2& u)
// u = fundamental unit, return norm(u)
{
    if(ZZ2::D < -4) { set(u); return 1; }
    if(sign(ZZ2::D) < 0) { set(u,0,1); return 1; }
    long s(1);
    IDL2 A(1);
    infra i;
    do { cfrac(i,A); s=-s; } while(!IsUnit(A));
    i.eval(u);
    return s;
}

long IDL2::FundUnit(ZZ2& u, const ZZ& D)
// set discriminant D and output u = fundamental unit.
// old value of D is restored on exit
{
    IDL2Push p;
    IDL2::init(D);
    return IDL2::FundUnit(u);
}

std::ostream& operator<<(std::ostream& s, const IDL2& A) {
    s << '[' << A.a << ' ' << A.b << ']';
    return s;
}