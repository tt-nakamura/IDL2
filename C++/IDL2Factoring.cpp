// uses NTL
//   http://www.shoup.net/ntl

#include "IDL2Factoring.h"
using namespace NTL;

void factor(Vec<Pair<IDL2, long> >& F, const IDL2& A)
// prime decomposition of ideal A
{
    int i,j,k(0);
    ZZ s,t;
    IDL2 B;
    Vec<Pair<ZZ, long> > g,h;
    F.SetLength(0);
    if(IsZero(A) || IsUnit(A)) return;
    primitive(B,A);
    norm(s,B);
    factor(g,s);
    factor(h,content(A));
    for(i=j=0; i<h.length(); i++) {
        if(IDL2::kron(h[i].a) == -1) {// inert
            F.SetLength(k+1);
            conv(F[k].a, h[i].a);
            F[k++].b = h[i].b;
        }
        else { if(i>j) h[j] = h[i]; j++; }
    }
    h.SetLength(j);
    i=j=0;
    while(i<g.length() || j<h.length()) {
        if(j==h.length() || i<g.length() && g[i].a < h[j].a) {
            F.SetLength(k+1);
            SetPrime(F[k].a, g[i].a);
            if(!div(B, F[k].a)) conj(F[k].a, F[k].a);
            F[k++].b = g[i++].b;
        }
        else if(IDL2::kron(h[j].a) > 0) {// split
            F.SetLength(k+2);
            SetPrime(F[k].a, h[j].a);
            conj(F[k+1].a, F[k].a);
            F[k+1].b = F[k].b = h[j].b;
            if(i<g.length() && g[i].a == h[j].a) {
                if(div(B, F[k].a)) F[k].b += g[i++].b;
                else F[k+1].b += g[i++].b;
            }
            j++; k+=2;
        }
        else {// ramify
            F.SetLength(k+1);
            SetPrime(F[k].a, h[j].a);
            F[k].b = 2*h[j].b;
            if(i<g.length() && g[i].a == h[j].a)
                F[k].b += g[i++].b;
            j++; k++;
        }
    }
}

void mul(IDL2& A, Vec<Pair<IDL2, long> >& F)
// A = product of F[i].a**F[i].b
{
    int i;
    IDL2 B;
    set(A);
    for(i=0; i<F.length(); i++) {
        power(B, F[i].a, F[i].b);
        A *= B;
    }
}

void IDL2FromNorm(Vec<IDL2>& J, const ZZ& n)
// J = vector of ideals that have given norm |n|.
// implements the algorithm described in T. Takagi
//   "Lectures on Elementary Number Theory" section 50
{
    long i,j,k,l,m;
    ZZ p;
    IDL2 A,B;
    Vec<IDL2> C;
    Vec<Pair<ZZ, long> > f;
    factor(f,n);
    J.SetLength(1);
    set(J[0]);
    for(i=0; i<f.length(); i++) {
        l = J.length();
        if((m = IDL2::kron(f[i].a)) <= 0) {
            power(p, f[i].a, f[i].b>>1);
            if(f[i].b & 1) {
                if(m<0) break;
                SetPrime(A, f[i].a);
                if(f[i].b > 1) A *= p;
                while(l) J[--l] *= A;
            }
            else while(l) J[--l] *= p;
        }
        else {
            C.SetLength(((f[i].b + 1)>>1)<<1);
            SetPrime(C[0], f[i].a);
            conj(C[1], C[0]);
            if(f[i].b > 1) {
                sqr(A, C[0]);
                sqr(B, C[1]);
                p = f[i].a;
            }
            if((f[i].b & 1) == 0) {
                C[0] = A;
                C[1] = B;
            }
            for(j=2; j<C.length(); j+=2) {
                mul(C[j],   C[j-2], A);
                mul(C[j+1], C[j-1], B);
            }
            for(j=C.length()-4; j>=0; j-=2) {
                C[j] *= p;
                C[j+1] *= p;
                p *= f[i].a;
            }
            J.SetLength(m = l*(f[i].b + 1));
            for(j=0; j<C.length(); j++)
                for(k=l-1; k>=0; k--)
                    mul(J[--m], J[k], C[j]);
            while(m) J[--m] *= p;
        }
    }
    if(i<f.length()) J.SetLength(0);
}