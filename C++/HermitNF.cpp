// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/mat_ZZ.h>
using namespace NTL;

void HermitNF(mat_ZZ& W, const mat_ZZ& A)
// A = matrix of shape (m,n) (m>=n)
// W = Hermite Normal Form A, i.e.,
//     lower triangular matrix of shape (n,n)
// reference: H. Cohen
//   "A Course in Computational Algebraic Number Theory"
//    Algorithm 2.4.5
{
    int m(A.NumRows()),n(A.NumCols());
    int i,j,k(m-1),l;
    ZZ d,s,t,u,v,b;
    if(m<n) Error("m<n in HermitNF");
    if(&W!=&A) W=A;
    for(l=n-1; l>=0; l--) {
        for(i=k-1; i>=0; i--) {
            if(IsZero(W[i][l])) continue;
            // DO NOT interchange
            // (u, W[i][l]) and (v, W[k][l]) in XGCD
            XGCD(d,u,v, W[i][l], W[k][l]);
            div(s, W[k][l], d);
            div(t, W[i][l], d);
            W[k][l] = d;
            clear(W[i][l]);
            for(j=0; j<l; j++) {
                mul(b, u, W[i][j]);
                MulAddTo(b, v, W[k][j]);
                W[i][j] *= s;
                MulSubFrom(W[i][j], t, W[k][j]);
                W[k][j] = b;
            }
        }
        if(IsZero(W[k][l])) continue;
        if(sign(W[k][l]) < 0)
           for(j=0; j<=l; j++)
               negate(W[k][j], W[k][j]);
        for(i=k+1; i<m; i++) {
            DivRem(s, W[i][l], W[i][l], W[k][l]);
            for(j=0; j<l; j++)
                MulSubFrom(W[i][j], s, W[k][j]);
        }
        k--;
    }
    if((k=m-n)==0) return;
    for(i=0; i<n; i++, k++)
        for(j=0; j<=i; j++)
            W[i][j] = W[k][j];
    W.SetDims(n,n);
}