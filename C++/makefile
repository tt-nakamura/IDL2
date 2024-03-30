NTL = -lntl -lgmp -L/usr/local/lib
OBJ = ZZ2.o IDL2.o HermitNF.o ZZFactoring.o ZZlib.o mpqs.o rho.o
BQF = BQF.o SolveBQE.o
CG = IDL2ClassGroup.o SmithNF.o

example1: example1.o $(BQF) IDL2Factoring.o $(OBJ)
	g++ example1.o $(BQF) IDL2Factoring.o $(OBJ) $(NTL)
example2: example2.o IDL2Factoring.o $(OBJ)
	g++ example2.o IDL2Factoring.o $(OBJ) $(NTL)
table1: table1.o $(CG) $(OBJ)
	g++ table1.o $(CG) $(OBJ) $(NTL)
table2: table2.o $(CG) $(OBJ)
	g++ table2.o $(CG) $(OBJ) $(NTL)