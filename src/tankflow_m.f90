module tankflow_m
    implicit none
    private
    real, parameter :: PI = 3.1416
    real, parameter :: Patmo = 101325.          ! atmospheric pressure at sea level [Pa]
    real, parameter :: GAM = 1.4                ! specific heat ratio of air
    real, parameter :: GRAVITY = 9.81           ! acceleration on earth [kg/m-s2]
    integer, parameter, public :: nmax = 100
    public :: tank_t, orifice_t
    public :: exportfcn
  
    type tank_t
        private
        real, public, dimension(nmax+1) :: z, d, vol
    contains
      procedure :: init_tk, update_tk
    end type
  
    type orifice_t
      real, public, dimension(nmax+1) :: mdot, velo, area, d
    contains
        procedure :: init_or, update_or
    end type
  
  contains
    subroutine init_tk(self, dtank, z)
        ! given as m, m
        implicit none
        class(tank_t), intent(inout) :: self
        real, intent(in) :: z, dtank
  
        self%z(1) = z
        self%d(1) = dtank
        self%vol(1) = z * 0.25 * PI * dtank ** 2
  
        !print *, "init : z, d, vol =", self%z(1), self%d(1), self%vol(1)
    end subroutine
  
    subroutine init_or(self, tank, dorifice)
        implicit none
        class(orifice_t), intent(inout) :: self
        class(tank_t), intent(in) :: tank
        real, intent(in) :: dorifice
  
        self%area(1) = 0.25 * PI * dorifice ** 2
        self%velo(1) = (2.0 * GRAVITY * tank%z(1)) ** 0.5
        self%mdot(1) = self%area(1) * self%velo(1)
        !print *, "init : velo, mdot =", self%velo(1), self%mdot(1)
        print *, ""
  end subroutine
  
    subroutine update_tk(self, orifice, it_, dt)
        implicit none
        class(tank_t), intent(inout) :: self
        type(orifice_t), intent(in) :: orifice
        real :: dt
        integer :: it_
  
        self%vol(it_) = self%vol(it_-1) - orifice%mdot(it_-1) * dt / 1.0
        self%d(it_) = self%d(it_-1)
        self%z(it_) = max(self%vol(it_) / (0.25 * PI * self%d(it_) ** 2), 0.0)
  
        !print *, "vol=", self%vol(it_), "m3, z=", self%z(it_), " m"
    end subroutine
  
    subroutine update_or(self, tank, it_)
        implicit none
        class(orifice_t), intent(inout) :: self
        class(tank_t), intent(in) :: tank
        integer :: it_
  
        self%area(it_) = self%area(it_-1)
        self%d(it_) = self%d(it_-1)
        self%velo(it_) = (2.0 * GRAVITY * tank%z(it_)) ** 0.5
        self%mdot(it_) = 1.0 * self%area(it_) * self%velo(it_)
        !print *, "velo=", self%velo(it_), " m/s, mdot=", self%mdot(it_), " kg/s"
    end subroutine
    
    function exportfcn(tank, orifice, t_) result(output_)
        implicit none
        type(tank_t), intent(in) :: tank
        type(orifice_t), intent(in) :: orifice
        real, dimension(:), allocatable, intent(in) :: t_
        real, dimension(:, :), allocatable :: output_
        integer :: n
  
        n = size(t_) - 1
        allocate(output_(n, 5))
  
        output_(1:n, 1) = t_(1:n)
        output_(1:n, 2) = tank%z(1:n)
        output_(1:n, 3) = tank%vol(1:n)
        output_(1:n, 4) = orifice%velo(1:n)
        output_(1:n, 5) = orifice%mdot(1:n)
    end function
  end module tankflow_m
  