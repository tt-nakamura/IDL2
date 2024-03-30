// uses NTL
//   http://www.shoup.net/ntl

#include "BQF.h"
#include "IDL2Factoring.h"
#include<exception>
using namespace NTL;

long SolveBQE(mat_ZZ & xy, const BQF& F, const ZZ& N, long M)
// solve binary quadratic equation ax^2 + bxy + cy^2 = N
// F = (a,b,c), n = RHS of equation,
// M = maximum number of solutions to search for.
// xy[i][0] = x component of i-th solution (0<=i<m).
// xy[i][1] = y component of i-th solution (0<=i<m).
//   where m is the number of solutions found (0<=m<=M).
// return sign (-1,0,1) of D = b^2 - 4ac.
// note that if D>=0, number of solutions is either inf or 0.
// implements the algorithm described in T. Takagi
//  "Lectures on Elementary Number Theory"
//   sections 49,52 (in Japanese)
{
    long i,j,k(0),m(1);
    ZZ D,f,d,r,s,t,n(N);
    BQF G;

    xy.SetDims(0,0);
    primitive(s,G,F);
    if(!IsOne(s) && !divide(n,n,s)) return 0;

    discriminant(D,G);
    if(sign(D) < 0 && sign(n) < 0) return -1;
    if(IsZero(D)) {// complete square
        if(sign(n) < 0) return 0;
        mul(s, G.a, n);
        SqrRoot(r,s);// s must be square
        if(sqr(r) != s) return 0;
        XGCD(d,s,t, G.a, G.b>>=1);
        if(!divide(r,r,d)) return 0;
        xy.SetDims(1,2);
        mul(xy[0][0], s, r);
        mul(xy[0][1], t, r);
        G.a /= d;
        G.b /= d;
        for(i=1; i<M; i++) {
            j = (i+1)>>1;
            if(i&1) j = -j;
            xy.SetDims(i+1, 2);
            sub(xy[i][0], xy[0][0], j*G.b);
            add(xy[i][1], xy[0][1], j*G.a);
        }
        return 0;
    }
    conductor(f,d,D);
    IDL2Push __p__;
    try { IDL2::init(d); }
    catch(std::exception) {
        // D is square (uninteresting cases)
        vec_ZZ d,v;
        mat_ZZ A;
        A.SetDims(2,2);
        v.SetLength(2);
        A[0][0] = A[0][1] = G.a;
        SqrRoot(s,D);
        add(A[1][0], G.b, s); A[1][0] >>= 1;
        sub(A[1][1], G.b, s); A[1][1] >>= 1;
        mul(s, G.a, n);
        divisor(d,s);
        for(i=0; i<d.length(); i++) {
            v[0] = d[i];
            v[1] = d[d.length()-1-i];
            solve1(t,v,A,v);
            if(!IsOne(t)) continue;
            xy.SetDims(k+1, 2);
            xy[k][0] = v[0];
            xy[k][1] = v[1];
            if(++k == M) break;
        }
        return 1;
    }
    // quadratic irrationality (interesting cases)
    long sgn;
    ZZ2 q,e,e1;
    IDL2 A(G.a);
    Vec<IDL2> J;
    if((sgn = IDL2::FundUnit(e)) < 0) sqr(e1,e);
    else e1=e;
    if(!IsOne(f)) {
        if(sign(d) > 0)
            for(q=e1; !divide(q.y, f); q*=e1) m++;
        else if(d == -4) m = 2;
        else if(d == -3) m = 3;
    }
    r = G.b;
    if(IsOdd(d)) r -= f;// r is even
    r >>= 1;
    set(q,r,f);
    A += q;
    IDL2FromNorm(J,n);
    for(i=0; i<J.length(); i++) {
        if(!IsPrincipal(q, J[i]*=A)) continue;
        norm(s,q);
        if(sign(s) != sign(n)) {
            if(sgn > 0) continue;
            else q *= e;
        }
        for(j=0; j<m; j++, q*=e1) {
            if(!divide(t, q.y, f)) continue;
            mul(s,r,t);
            sub(s, q.x, s);
            if(!divide(s, s, G.a)) continue;
            xy.SetDims(k+1, 2);
            xy[k][0] = s;
            xy[k][1] = t;
            if(++k == M) return sign(D);
        }
    }
    if((m=k)==0 || M<1 || D<-4) return sign(D);
    // search for associate solutions
    ZZ2 u,p(1),p1(1);
    IDL2::init(D);
    if(IDL2::FundUnit(e) < 0) sqr(e,e);
    conj(e1,e);
    r = G.b;
    if(IsOdd(D)) r--;
    r >>= 1;
    if(D==-4) M = min(M, 2*m); else
    if(D==-3) M = min(M, 3+m);
    for(i=0;; i++) {
        if(i&1) u = (p *= e);
        else u = (p1 *= e1);
        for(j=0; j<m; j++) {
            q.y = xy[j][1];
            mul(q.x, G.a, xy[j][0]);
            MulAddTo(q.x, q.y, r);
            q *= u;
            mul(s, q.y, r);
            sub(s, q.x, s);
            xy.SetDims(k+1, 2);
            div(xy[k][0], s, G.a);
            xy[k][1] = q.y;
            if(++k == M) return sign(D);
        }
    }
    return sign(D);
}

void AssocSol(mat_ZZ& xy, const mat_ZZ& xy0, const BQF& F, long k)
// xy = associate solutions of ax^2 + bxy + cy^2 = N
// Assume b^2 - 4ac > 0.
// xy0 = output of SolveBQE, F = coeffients a,b,c
// k = exponent of unit e^k to be multiplied to
//     quadratic integer that corresponds to xy0
{
    long i,j;
    ZZ s,r,D;
    BQF G;

    primitive(s,G,F);
    discriminant(D,G);
    if(sign(D) <= 0) return;

    ZZ2 e,q,p;
    IDL2Push __p__;
    try { IDL2::init(D); }
    catch(std::exception) { return; }
    if(IDL2::FundUnit(e) < 0) sqr(e,e);
    if(k<0) conj(e,e);
    power(p, e, abs(k));
    r = G.b;
    if(IsOdd(D)) r--;
    r >>= 1;
    if(&xy!=&xy0) xy = xy0;
    for(i=0; i<xy.NumRows(); i++) {
        q.y = xy[i][1];
        mul(q.x, G.a, xy[i][0]);
        MulAddTo(q.x, q.y, r);
        q *= p;
        mul(s, q.y, r);
        sub(s, q.x, s);
        div(xy[i][0], s, G.a);
        xy[i][1] = q.y;
    }
}