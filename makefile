FORT=`which gfortran`
FLAGS=${MYLIB} ${LAP}

a.out : global_constant.o input_data.o potentials.o \
        sp_basis.o print_array.o main.o
	${FORT} global_constant.o input_data.o potentials.o \
        sp_basis.o print_array.o main.o ${FLAGS}
global_constant.o : global_constant.f90
	${FORT} ${FLAGS} -c global_constant.f90
input_data.o : input_data.f90
	${FORT} ${FLAGS} -c input_data.f90
potentials.o : potentials.f90
	${FORT} ${FLAGS} -c potentials.f90
print_array.o : print_array.f90
	${FORT} ${FLAGS} -c print_array.f90
sp_basis.o : sp_basis.f90
	${FORT} ${FLAGS} -c sp_basis.f90
main.o : main.f90
	${FORT} ${FLAGS} -c main.f90

.PHONY : clean
clean :
	$(RM) *.o *.mod

