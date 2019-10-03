!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function eos_eni(rho,pre,gamma)
  ! EOS: int. energy
  include "precision.h"
  ! I/O
  real(kind=rkind)             :: eos_eni
  real(kind=rkind), intent(in) :: rho, pre, gamma

!  print*,rho, pre, gamma

  eos_eni = pre / rho / (gamma-1.0_d)
  return
end function eos_eni 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function eos_pre(rho,ein,gamma)
  ! EOS: pressure
  include "precision.h"
  ! I/O
  real(kind=rkind)             :: eos_pre
  real(kind=rkind), intent(in) :: rho, ein, gamma

  eos_pre = ein * rho * (gamma-1.0_d)
  return
end function eos_pre

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function eos_sspe(ein,gamma)
  ! EOS: speed of sound from int. energy
  include "precision.h"
  ! I/O
  real(kind=rkind)             :: eos_sspe
  real(kind=rkind), intent(in) :: ein, gamma

  real(kind=rkind) :: ein_safe

  if (ein.lt.0.0_d) print*,'EOS_SSPE: negative int. energy:',ein
  ein_safe = max(ein,0.0_d)

  eos_sspe = sqrt(gamma*(gamma-1.0_d)*ein_safe)
  return
end function eos_sspe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function eos_sspp(rho,pre,gamma)
  ! EOS: speed of sound from density and pressure
  include "precision.h"
  ! I/O
  real(kind=rkind)             :: eos_sspp
  real(kind=rkind), intent(in) :: rho, pre, gamma

  real(kind=rkind) :: rho_safe, pre_safe

  if (rho.le.0.0_d) print*,'EOS_SSPP: too small density:',rho
  rho_safe = max(rho,1.0e-20_d)
  if (pre.lt.0.0_d) print*,'EOS_SSPP: negative pressure:',pre
  pre_safe = max(pre,0.0_d)

  eos_sspp = sqrt(gamma * pre_safe / rho_safe)
  return
end function eos_sspp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

