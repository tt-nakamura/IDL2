#include "BQF.h"
using namespace NTL;

void examle(const BQF& f, const ZZ& n) {
    long i;
    ZZ s;
    mat_ZZ xy;
    i = SolveBQE(xy,f,n);
    std::cout << xy << ' ' << i << std::endl;
    for(i=0; i<xy.NumRows(); i++) {
        eval(s, f, xy[i][0], xy[i][1]);
        if(s!=n) std::cout << "wrong" << std::endl;
    }
    AssocSol(xy, xy, f, -1);
    std::cout << xy << std::endl;
    for(i=0; i<xy.NumRows(); i++) {
        eval(s, f, xy[i][0], xy[i][1]);
        if(s!=n) std::cout << "wrong" << std::endl;
    }
}

main() {
    BQF f;
    ZZ n;
    set(f, 1,1,1); n = 175; examle(f,n);// D==-3
    set(f, 2,-3,-2); n = 18; examle(f,n);// D==25
    set(f, 18, -60, 50); n = 8; examle(f,n);// D==0
    set(f, 111,-24,-391);
    eval(n,f,29,31); examle(f,n);
}