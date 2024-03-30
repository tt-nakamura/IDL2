#include "IDL2ClassGroup.h"
using namespace NTL;

main() {
    long d, d1(-1000), d2(d1-50), i;
    ZZ h;
    Vec<Pair<ICG2, long> > G;
    for(d=d1; d>=d2; d--) {
        try { ICG2::init(d); }
        catch(std::exception) { continue; }
        ICG2::ClassNum(h);
        i = generator(G);
        if(i!=h) Error("i!=h");
        std::cout << ZZ2::D << ' ';
        std::cout << h << ' ';
        std::cout << G << ' ';
        std::cout << std::endl;
    }
}