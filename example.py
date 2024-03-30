from SolveBQE import SolveBQE, AssocSol

def f(F,x,y):
    a,b,c = F
    return a*x**2 + b*x*y + c*y**2

F = [(11,-24,-45),
     (18,-60,50),# D == 0
     (2,-3,-2)] # D == 25
n = [f(F[0],29,31), 8, 18]
N,k = 8, -3

for F,n in zip(F,n):
    xy = SolveBQE(F,n,N)
    print(xy)
    for x,y in xy:
        if f(F,x,y) != n: print('wrong')
    xy = AssocSol(xy,F,k)
    print(xy)
    for x,y in xy:
        if f(F,x,y) != n: print('wrong')
