#include "IDL2Factoring.h"
using namespace NTL;

main() {
    ZZ2::init(-20);
    IDL2 A,B,C;
    Vec<Pair<IDL2, long> > F;
    SetPrime(A,7); conj(A,A);
    SetPrime(B,5);
    mul(C,A,B);
    factor(F,C);
    std::cout << F << std::endl;
    mul(A,F);
    std::cout << A << std::endl;
    std::cout << C << std::endl;
}