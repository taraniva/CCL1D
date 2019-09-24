  implicit none
  ! 
  save
  !
  ! Miminum accuracy for real numbers (No. of significant digits, exponent range)
  integer, parameter :: rkind = selected_real_kind(P=15,R=307)
  integer, parameter :: d = rkind
