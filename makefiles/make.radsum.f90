FC = f90
FCFLAG = -mips4 -r10000 -lfastm -O3 -64 -r8 -i8 -TENV:X=0

radsum     :     radsum.o util_sgi.o
	$(FC) $(FCFLAG) -o radsum radsum.o util_sgi.o
#
radsum.o: src/radsum.f
	$(FC) $(FCFLAG) -c src/radsum.f
#
util_sgi.o: src/util_sgi.f
	$(FC) $(FCFLAG) -c src/util_sgi.f

