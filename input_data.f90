module input_data
  type inp
    real(8) :: A, Z
    real(8) :: V0, av, r0, Rn, mass, Vso, rc0, Rc
    real(8) :: rmin, rmax, dr, rgrid

    contains

    procedure :: read_input
  end type

  contains

    subroutine read_input(this, fname)
      use global_constant
      implicit none
      class(inp), intent(out) :: this
      character(len=*), intent(in) :: fname

      open(7,file=fname, action="read")
      read(7,*) this%A, this%Z
      read(7,*) this%V0, this%r0, this%av
      read(7,*) this%Vso
      read(7,*) this%rmax, this%dr
      close(7)

      this%mass = mass_n
      this%rgrid = nint(this%rmax/this%dr)
      this%rmin = 1.0d-100
      this%rmax = this%rmin + dble(this%rgrid) * this%dr
      this%rc0 = this%r0
!     call setPot(V, r, a, Vls, this%A)

    end subroutine
end module
