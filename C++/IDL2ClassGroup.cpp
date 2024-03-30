// uses NTL
//   http://www.shoup.net/ntl

#include "IDL2ClassGroup.h"
#include "GroupGenerator.h"
#include "ZZFactoring.h"
#include<exception>
#include<list>
using namespace NTL;

ZZ ICG2::amax;// Minkowski bound for a

long IsFundDisc(const NTL::ZZ& a);

void ICG2::init(const ZZ& D) {// set discriminant
    if(!IsFundDisc(D))
        throw std::runtime_error("D is not fundamental");
    IDL2::init(D);
    if(sign(D) > 0)
        RightShift(amax, IDL2::S, 1);
    else {
        mul(amax, D, -3);
        SqrRoot(amax, amax);
        amax /= 3;
    }
}

void mul(ICG2& C, const ICG2& A, const ICG2& B) {// C=A*B
    mul((IDL2&)C, (IDL2&)A, (IDL2&)B);
    reduce(C,C);
}

void sqr(ICG2& B, const ICG2& A) {// B=A*A
    sqr((IDL2&)B, (IDL2&)A);
    reduce(B,B);
}

void power(ICG2& B, const ICG2& A, long n) {// B=A^n (n may be n<0)
    if(n==0) set(B);
    else if(n<0) { ICG2 C; inv(C,A); power(B,C,-n); }
    else if(&B==&A) power(B, ICG2(A), n);
    else {
        long m((1<<(NumBits(n)-1))>>1);
        for(B=A; m; m>>=1) { sqr(B,B); if(n&m) B*=A; }
    }
}

static void ImQIClassNum(ZZ& h)
// class number of imaginary quadratic fields
// reference: H. Cohen
//  "A Course in Computational Algebraic Number Theory"
//   Algorithm 5.3.5
{
    long i,j;
    ZZ s,r,b,ac;
    Vec<ZZ> a;
    ZZ2 q;
    RightShift(s, ICG2::amax, 1);
    for(clear(h); r<=s; r++) {
        set(q,r,1);
        norm(ac,q);
        LeftShift(b,r,1);
        if(ZZ2::Dm4) b++;
        divisor(a,ac);
        for(i=0, j=a.length()-1; i<=j; i++, j--) {
            if(a[i]==b || i==j || IsZero(b)) h++;
            else if(a[i] > b) h += 2;
        }
    }
}

static void ReQIClassNum(ZZ& h)
// class number of real quadratic fields
// reference: J. Buchmann and U. Vollmer
//  "Binary Quadratic Forms" section 6.17
{
    long i,j,k;
    ZZ r,s,ac;
    Vec<ZZ> a;
    IDL2 A;
    ZZ2 &q(A.b);
    std::list<IDL2> L;
    std::list<IDL2>::iterator p;
    if(ZZ2::Dm4 == 0) set(r);
    for(; r <= IDL2::W1; r++) {
        sub(s, IDL2::W1, r);
        set(q,r,1);
        norm(ac,q);
        divisor(a,ac);
        for(i=0, j=a.length()-1; i<=j; i++, j--) {
            if(a[i] <= s) continue;
            A.a = a[i]; L.push_back(A);
            if(i==j) continue;
            A.a = a[j]; L.push_back(A);
        }
    }
    for(p = L.begin(); p != L.end(); p++)
        (*p).b.x %= (*p).a;
    for(clear(h); !L.empty(); h++) {
        for(A = L.front();; cfrac(A)) {
            for(p = L.begin(); p != L.end(); p++)
                if(*p == A) break;
            if(p == L.end()) break;
            L.erase(p);            
        }
    }
}

void ICG2::ClassNum(ZZ& h) {
    if(sign(ZZ2::D) < 0) ImQIClassNum(h);
    else ReQIClassNum(h);
}

void ICG2::ClassNum(ZZ& h, const ZZ& D) {
    ICG2Push p;// save old D
    ICG2::init(D);// set new discriminant
    ICG2::ClassNum(h);
}

long generator(Vec<Pair<ICG2, long> >& G, long min)
// generator of class group
{
    if(!ICG2::amax.SinglePrecision())
        throw std::runtime_error("|D| is too large");
    long i,j,k;
    ICG2 A,B,C;
    Vec<ICG2> P;
    PrimeSeq ps;
    // gather candidates
    while((k = ps.next()) <= ICG2::amax) {
        if(IDL2::kron(k) < 0) continue;
        SetPrime(A,k);
        reduce(A,A);
        P.append(A);
    }
    G.SetLength(0);
    if(P.length() == 0) return 1;
    k = GroupGenerator(G,P);// utilize general routine
    if(!min) return k;
    // find minimum generators
    for(i=0; i<G.length(); i++) {
        set(A);
        B = G[i].a;
        for(j=1; j<G[i].b; j++) {
            A *= B;
            if(GCD(j, G[i].b) > 1) continue;
            if(A.a < G[i].a.a) G[i].a = A;
            if(sign(ZZ2::D) < 0) continue;
            // find a minimum in cfrac cycle for D>0
            for(cfrac(C=A); (IDL2)C!=A; cfrac(C))
                if(C.a < G[i].a.a) G[i].a = C;
        }
    }
    return k;
}

long generator(Vec<Pair<ICG2, long> >& G,
               const ZZ& D, long min) {
    ICG2Push p;// save old D
    ICG2::init(D);// set new discriminant
    return generator(G, min);
}