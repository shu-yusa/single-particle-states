module potentials
  use input_data
  use global_constant

  type potential
    type(inp), pointer :: ip
    real(8) :: V0, av, r0, Rn, Vso, rc0, Rc
    real(8) :: c1, Ze2
    
    contains

    procedure :: potential_
    procedure :: Vpot
  end type

  contains

    subroutine potential_(this, ip)
      implicit none
      class(potential), intent(out) :: this
      type(inp), intent(in), target :: ip

      this%ip => ip
      this%V0 = ip%V0
      this%av = ip%av
      this%r0 = ip%r0
      this%Vso = ip%Vso
      this%rc0 = ip%rc0

      this%Rn = this%r0 * (ip%A-1.0d0) ** (1.0d0/3.0d0)
      this%Rc = this%rc0 *  (ip%A-1.0d0) ** (1.0d0/3.0d0)

      this%c1 = hbar * hbar / (2.0d0 * mass_n)
      this%Ze2 = (this%ip%Z-1.0d0) * hbar / 137.0d0

    end subroutine

    function Vpot(this, q, s, L, r) result(V)
      implicit none
      class(potential), intent(in) :: this
      integer, intent(in) :: s, L, q
      real(8), intent(in) :: r
      real(8) :: V, er

      er = exp((r - this%Rn) / this%av)
      if (q == 0) then
        select case(s)
        case(-1)
          V = Vn(this, r, er) + Vcnt(this, L, r)
        case default
          V = Vn(this, r, er) + Vcnt(this, L, r) + Vls(this, s, L, r, er)
        end select
      else
        select case(s)
        case(-1)
          V = Vn(this, r, er) + Vcnt(this, L, r) + Vc(this, r)
        case default
          V = Vn(this, r, er) + Vcnt(this, L, r) + Vls(this, s, L, r, er) + Vc(this, r)
        end select
      end if

    end function

    function Vn(this, r, er) result(V)
      implicit none
      class(potential), intent(in) :: this
      real(8), intent(in) :: r, er
      real(8) :: V

!     V = - this%V0 / (1.0d0 + exp((r-this%Rn)/this%av))
      V = - this%V0 / (1.0d0 + er)
    end function


    function Vls(this, s, L, r, er) result(V)
      implicit none
      class(potential), intent(in) :: this
      integer, intent(in) :: s, L
      real(8), intent(in) :: r, er
      real(8) :: V, df_dr

      df_dr = - er / (this%av * (1.0d0 + er)**2)
      select case(s)
      case(0)
        V = 0.5d0 * this%Vso * df_dr / r * dble(L)
      case(1)
        if (L == 0) then
          V = 0
        else
          V = - 0.5d0 * this%Vso * df_dr / r * dble(L+1)
        end if
      case default
        stop "wrong argument"
      end select

    end function

    function Vcnt(this, L, r) result(V)
      implicit none
      class(potential), intent(in) :: this
      integer, intent(in) :: L
      real(8), intent(in) :: r
      real(8) :: V

      V = dble(L*(L+1)) * this%c1 / (r*r)

    end function

    function Vc(this, r) result(V)
      implicit none
      class(potential), intent(in) :: this
      real(8), intent(in) :: r
      real(8) :: V

      if (r > this%Rc) then
        V = this%Ze2 / r
      else
        V = this%Ze2 * (3.0d0*this%Rc*this%Rc - r*r) / (2.0d0 * this%Rc**3)
      end if

    end function

end module
