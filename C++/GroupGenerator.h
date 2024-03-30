// uses NTL
//   http://www.shoup.net/ntl

#ifndef __GroupGenerator_h__
#define __GroupGenerator_h__

#include<NTL/mat_ZZ.h>
#include<NTL/pair.h>
#include<list>

void SmithNF(NTL::vec_ZZ& D, NTL::mat_ZZ& U, const NTL::mat_ZZ& A);

template<class T>
long GroupGenerator(NTL::Vec<NTL::Pair<T, long> >& G,
                    const NTL::Vec<T>& P)
// generating system of commutative group
// P = candidates of generators
// G = vector of (generator, order) pair
// return order of group = product of G[i].b
// T = type of group elements
// T must suport following funcitons:
//   set(T& a) : a is set to unit element of T
//   mul(T& c, T& a, T& b) : group operation c = a*b
//   operator*=(T& b, T& a) : group operation b *= a
//   operator==(T& a, T& b) : equality of group elements
//   power(T& b, T& a, long n) : b = a^n (n may be n<0)
// reference: J. Buchmann and U. Vollmer
//   "Binary Quadratic Forms" Algorithm 9.1
{
    long i,j,k,l(0),m(0);
    T f;
    NTL::Vec<T> F,H;
    NTL::Vec<NTL::Vec<long> > U;
    NTL::Vec<NTL::Pair<T, NTL::Vec<long> > > S;
    std::list<T> L;
    typename std::list<T>::iterator p;
    // find relation basis in Hermite Normal Form (HNF)
    for(i=0; i<P.length(); i++)
        L.push_back(P[i]);
    S.SetLength(1);
    set(S[0].a);// unit element
    for(;; l++) {
        for(i=m; i<S.length(); i++)
            for(p = L.begin(); p != L.end();)
                if(*p == S[i].a) p = L.erase(p);// equality
                else p++;
        if(L.empty()) break;
        F.append(f = L.front());
        H.SetLength(0);
        k = m = S.length();
        for(;;) {
            for(i=0; i<m; i++)
                if(f == S[i].a) break;// equality
            if(i<m) break;
            H.append(f);
            f *= L.front();// group operation
        }
        U.SetLength(l+1);
        for(j=0; j<l; j++)
            U[l].append(-S[i].b[j]);
        U[l].append(H.length() + 1);
        for(i=0; i<H.length(); i++) {
            for(j=0; j<m; j++, k++) {
                S.SetLength(k+1);
                mul(S[k].a, S[j].a, H[i]);// group operation
                S[k].b = S[j].b;
                S[k].b.append(i+1);
            }
        }
        for(i=0; i<m; i++)
            S[i].b.append(0);
    }
    // convert HNF to Smith Normal Form
    NTL::mat_ZZ B,V;
    NTL::vec_ZZ d;
    B.SetDims(l,l);
    for(i=0; i<l; i++)
        for(j=0; j<=i; j++)
            B[i][j] = U[i][j];// HNF
    SmithNF(d,V,B);
    for(k=d.length(); k>0; k--)
        if(d[k-1] > 1) break;
    G.SetLength(k);
    for(i=0; i<k; i++) {
        set(G[i].a);// unit element
        for(j=0; j<l; j++) {
            conv(m, V[i][j]);
            power(f, F[j], m);// group operation
            G[i].a *= f;// group operation
        }
        conv(G[i].b, d[i]);
    }
    return S.length();
}

#endif // __GroupGenerator_h__