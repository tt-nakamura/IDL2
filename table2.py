from IDL2 import FundUnit
from IDL2ClassGroup import ICG2, ClassNum, generator

D1 = 1000
D2 = D1 + 50

for D in range(D1,D2):
    try: ICG2.init(D)
    except: continue
    e,s = FundUnit()
    h = ClassNum()
    G,h1 = generator()
    if h!=h1: print('wrong')
    print(D,h,G,e,s)
