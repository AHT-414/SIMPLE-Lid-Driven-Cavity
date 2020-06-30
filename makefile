MODS = variables.o
SUBS = cavityFlow.f95 allocateArrays.f95 xLinkCoeffs.f95 yLinkCoeffs.f95 solvePressureMatrix.f95 Update.f95 solve_tridiag.f95

cavityFlow: $(MODS) $(SUBS)
		gfortran -o cavityFlow $(MODS) $(SUBS) -L/usr/local/lib -llapack -lblas

variables.o: variables.f95
		gfortran -c variables.f95
