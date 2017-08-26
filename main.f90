program main
  use, intrinsic :: iso_fortran_env
  use input_data
  use potentials
  use sp_basis
  use print_array
  implicit none
  integer :: L, s, i, q
  real(8), parameter :: epsr=1.0d-15
  type(inp) :: ip
  type(sp_states) :: sp
  type(potential) :: pot

  call ip%read_input("input.txt")
  call pot%potential_(ip)
  call sp%sp_states_(ip, pot)

  do q=0, 1
    if (q == 0) then
      write(output_unit,*) "neutron"
    else
      write(output_unit,*) "proton"
    end if

    do s=0, 1
      do L=s, 15
        if (.not. sp%bind_energy(-ip%V0, epsr, q, s, L)) exit
        call print_vec(sp%E)
      end do
    end do
  end do

end program
