IFC=`which ifort`
FLAGS=${MYLIB} ${LAP}

a.out : global_constant.o input_data.o potentials.o \
        sp_basis.o main.o
	${IFC} global_constant.o input_data.o potentials.o \
        sp_basis.o main.o ${FLAGS}
global_constant.o : global_constant.f90
	${IFC} ${FLAGS} -c global_constant.f90
input_data.o : input_data.f90
	${IFC} ${FLAGS} -c input_data.f90
potentials.o : potentials.f90
	${IFC} ${FLAGS} -c potentials.f90
sp_basis.o : sp_basis.f90
	${IFC} ${FLAGS} -c sp_basis.f90
main.o : main.f90
	${IFC} ${FLAGS} -c main.f90

.PHONY : clean
clean :
	$(RM) *.o *.mod

