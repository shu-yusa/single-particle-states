module sp_basis
  use input_data
  use global_constant, only : mass_n, mass_p, hbar
  use potentials

  type sp_states
    type(inp), pointer :: ip
    type(potential), pointer :: pot
    real(8), allocatable, dimension(:) :: E
    real(8), allocatable, dimension(:,:,:) :: psi
    real(8) :: c

    contains
    procedure :: sp_states_
    procedure :: bind_energy
  end type

  contains
    subroutine sp_states_(this, ip, pot)
     implicit none
      class(sp_states), intent(out) :: this
      type(inp), intent(in), target :: ip
      type(potential), intent(in), target :: pot

      this%ip => ip
      this%pot => pot
      this%c = 2.0d0 * mass_n / (hbar*hbar) * this%ip%dr * this%ip%dr

    end subroutine


    function bind_energy(this, Est, epsr, q, s, L) result(b)
      implicit none
      class(sp_states), intent(inout) :: this
      integer :: Nnode, Nstate, Nstate_max, n, i
      integer, intent(in) :: L, s, q
      real(8), intent(in) :: Est, epsr
      real(8) :: E0, Emin, Emax
      logical :: b

      Emin = Est
      Emax = 2.0d0
      Nstate = 0
      call solve_eq_Numerov(this, q, s, L, 2.0d0, Nstate_max)
      write(6,*) "L=",L,"n=",Nstate_max
      b = (Nstate_max /= 0)
      if (.not. b) return
      if (allocated(this%E)) then
        deallocate(this%E)
      end if
      allocate(this%E(Nstate_max))

      do
        this%E(Nstate+1) = 0.5d0 * Emin
        E0 = 2.0d0
        do 
          call solve_eq_Numerov(this, q, s, L, this%E(Nstate+1), Nnode)
          if (Nnode >= Nstate + 1) then
            Emax = this%E(Nstate+1)
          else if (Nnode == Nstate) then
            Emin = this%E(Nstate+1)
          end if
          this%E(Nstate+1) = 0.5d0 * (Emin + Emax)
          if (abs(E0 - this%E(Nstate+1)) < abs(this%E(Nstate+1)) * epsr) exit
          E0 = this%E(Nstate+1)
        end do
        Nstate = Nstate + 1
        if (Nstate == Nstate_max) exit
        Emin = this%E(Nstate)
        Emax = 2.0d0
      end do
    end function


    subroutine solve_eq_Numerov(this, q, s, L, E, Nnode)
      implicit none
      class(sp_states), intent(in) :: this
      logical :: flag_p, flag_n
      integer, intent(in) :: L, s, q
      integer, intent(out) :: Nnode
      integer :: i
      real(8), intent(in) :: E
      real(8) :: r,  u0, u, c, w1, w2, u1
      real(8) :: E_V0, E_V1


      call Init_Numerov(q, s, L, E, r, u0, u1, this%ip%dr)
      Nnode = 0
      flag_p = .true.
      flag_n = .false.
      E_V0 = E-this%pot%Vpot(q,s,L,this%ip%rmin)
      E_V1 = E-this%pot%Vpot(q,s,L,this%ip%rmin+this%ip%dr)

      do i=1, this%ip%rgrid+1
        r = this%ip%rmin + dble(i+1) * this%ip%dr
        w1 = 1.0d0 + this%c * (E_V0) / 12.0d0
        w2 = (2.0d0 - 5.0d0/6.0d0*this%c*(E_V1))
        E_V0 = E_V1
        E_V1 = E - this%pot%Vpot(q,s,L,r)
        u = (w2*u1 - w1*u0) / (1.0d0 + this%c*(E_V1)/12.0d0)
        u0 = u1
        u1 = u
        if (flag_n .and. u > 0.0d0) then
          flag_p = .true.
          flag_n = .false.
          Nnode = Nnode + 1
        else if (flag_p .and. u < 0.0d0) then
          flag_p = .false.
          flag_n = .true.
          Nnode = Nnode + 1
        end if
      end do

      contains
        subroutine Init_Numerov(q, s, L, E, x, y0, y, h)
          implicit none
          integer, intent(in) :: L, s, q
          real(8), intent(in) :: E, h
          real(8), intent(out) :: x, y0, y
          real(8) :: yt(2)

          x = this%ip%rmin
          y0 = this%ip%rmin ** (L+1)
          yt(1) = y0
          yt(2) = dble(L+1) * this%ip%rmin ** L
          call Runge_Kutta(this, q, s, L, E, x, yt, h)
          y = yt(1)

        end subroutine
    end subroutine

    subroutine Runge_Kutta(this, q, s, L, E, x, y, h)
      implicit none
      class(sp_states), intent(in) :: this
      integer, intent(in) :: s, L, q
      real(8), intent(in) :: x, h, E
      real(8), intent(inout) :: y(2)
      real(8), dimension(2) :: k1, k2, k3, k4, y2
      real(8) :: h2

      h2 = 0.5d0 * h
      call fct(q, s, L, E, x, y, k1)
      y2 = y + h2*k1
      call fct(q, s, L, E, x+h2, y2, k2)
      y2 = y + h2*k2
      call fct(q, s, L, E, x+h2, y2, k3)
      y2 = y + h*k3
      call fct(q, s, L, E, x+h,  y2,  k4)
      y = y + h * (k1 + 2.0d0 * (k2 + k3) + k4) / 6.0d0

      contains

        subroutine fct(q, s, L, E, r, y, f)
        implicit none
        integer, intent(in) :: L, s, q
        real(8), intent(in) :: r, E, y(2)
        real(8), intent(out) :: f(2)

        f(1) = y(2)
        f(2) = - y(1) * 2.0d0 * mass_n / (hbar*hbar) * (E-this%pot%Vpot(q,s,L,r))

        return
        end subroutine
    end subroutine

    subroutine solve_eq(this, q, s, L, E, Nnode)
      implicit none
      class(sp_states), intent(in) :: this
      logical :: flag_p, flag_n
      integer :: i
      integer, intent(in) :: L, s, q
      integer, intent(out) :: Nnode
      real(8), intent(in) :: E
      real(8) :: r, u(2)

      u(1) = this%ip%rmin ** (L+1)
      u(2) = dble(L+1) * this%ip%rmin ** L

      Nnode = 0
      flag_p = .true.
      flag_n = .false.
      do i=1, this%ip%rgrid+1
        r = this%ip%rmin + dble(i-1) * this%ip%dr
        call Runge_Kutta(this, q, s, L, E, r, u, this%ip%dr)
        if (flag_n .and. u(1) > 0.0d0) then
          flag_p = .true.
          flag_n = .false.
          Nnode = Nnode + 1
        else if (flag_p .and. u(1) < 0.0d0) then
          flag_p = .false.
          flag_n = .true.
          Nnode = Nnode + 1
        end if
      end do

    end subroutine

end module
