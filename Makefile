OPT = -ffpe-trap=invalid,zero,overflow -O2

EXE = meteochem.exe

all: $(EXE)

$(EXE): main.f90 chemf.o opkdmain.o opkda1.o opkda2.o
	gfortran $(OPT) main.f90 chemf.o opkdmain.o opkda1.o opkda2.o -o $(EXE)

chemf.o: chemf.f90 opkdmain.o opkda1.o opkda2.o
	gfortran $(OPT) -c chemf.f90

opkda1.o: opkda1.f
	gfortran $(OPT) -std=legacy -w -c opkda1.f

opkda2.o: opkda2.f
	gfortran $(OPT) -std=legacy -c opkda2.f

opkdmain.o: opkdmain.f
	gfortran $(OPT) -std=legacy -c opkdmain.f

run:
	$(EXE)
clean:
	rm *.o, *.mod, *.exe

plot:
	python.exe vis.py


