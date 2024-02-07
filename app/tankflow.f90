program tankflow
  use tankflow_m, only : tank_t, orifice_t, exportfcn, nmax
  
  implicit none
  real, dimension(:), allocatable :: time
  real :: tmax = 500.0, dt
  integer :: i, i_current, i_target
  integer :: unit0
  real, dimension(:, :), allocatable :: output0

  type(tank_t) :: tank
  type(orifice_t) :: orifice


  dt = tmax / nmax
  print *, "maximum time = ", tmax, " sec"
  print *, "that's ", nmax, " iterations, with dt=", dt, " sec"
  print *,""

  allocate(time(nmax))
  time = [(dt * (i - 1), i=1, nmax + 1)]
  allocate(output0(nmax+1, 5))


  ! INITIAL
  i_current = 1

  call tank%init_tk(10.0, 40.0)     ! diameter, height
  call orifice%init_or(tank, 1.0)         ! diameter
  
  print *, "Initialization:"
  print *, "    tank : z, d, vol =", tank%z(1), tank%d(1), tank%vol(1)
  print *, "    iorifice : velo, mdot =", orifice%velo(1), orifice%mdot(1)


  ! ITERATE
  i_current = 2
  do while (i_current < nmax)
    call tank%update_tk(orifice, i_current, dt)
    call orifice%update_or(tank, i_current)
    i_current = i_current + 1
  end do

  !call tank%datarecord(time)

output0 = exportfcn(tank, orifice, time)

print *, "Results:"

i_target = nmax
open(newunit=unit0, file="output.txt", form='FORMATTED')
write(unit0, *) "time(s)", "z(m)", "vol(m3)", "v(m/s)", "mdot(kg/s)"
do i=1, nmax - 1
  write(unit0, *) output0(i, 1:5)
  if (output0(i, 2) > 5.0 .and. output0(i+1, 2) <= 5.0) then
    i_target = i + 1
    print *, "    the tank has drained down to z=", output0(i_target,2), " m at t=", time(i_target), " sec"
    print *, ""
  end if
end do
close(unit0)


print *, "    at the end, z=", output0(nmax,2), " m at t=", time(nmax), " sec"

call execute_command_line ("gnuplot output.plt", wait=.true.)

end program tankflow