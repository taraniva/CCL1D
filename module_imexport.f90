module module_imexport
	!Dependent on DATA nad GEOMETRY

	contains

	subroutine setup_test(topo,vars,par)
		use module_data
		use module_geometry

		implicit none

		include "eos_func.h"
  		include "map_func.h"
  		! I/O
  		type(topo_type), intent(inout)  :: topo
  		type(vars_type), intent(inout)  :: vars
		type(para_type), intent(inout)  :: par

		!Local
		integer :: nc, nn

		! Pointers definition
		!Nodal variables
		real(kind=rkind), dimension(:), pointer   :: massn => null()

		!Cell variables
		real(kind=rkind), dimension(:), pointer   :: gamma => null()
		real(kind=rkind), dimension(:), pointer   :: xc => null()
		real(kind=rkind), dimension(:), pointer   :: uc => null()
		real(kind=rkind), dimension(:), pointer   :: vol => null()
		real(kind=rkind), dimension(:), pointer   :: mass => null()
		real(kind=rkind), dimension(:), pointer   :: rho => null()
		real(kind=rkind), dimension(:), pointer   :: pre => null()
		real(kind=rkind), dimension(:), pointer   :: eni => null()
		real(kind=rkind), dimension(:), pointer   :: ent => null()
		real(kind=rkind), dimension(:), pointer   :: ssp => null()

		!2nd order
		!real(kind=rkind), dimension(:), pointer   :: du => null()
		!real(kind=rkind), dimension(:), pointer   :: dpre => null()

		! Pointers assignment
		massn => vars%mass_n
		gamma => vars%gamma_c
		xc	  => vars%x_c
		uc    => vars%u_c
		vol   => vars%vol_c
		mass  => vars%mass_c
		rho   => vars%rho_c
		pre   => vars%pre_c
		eni   => vars%eni_c
		ent   => vars%ent_c
		ssp   => vars%ssp_c

		nc = topo%nc
		nn = topo%nn
		  
		select case(par%prob)
			! Noh implosion test
			case("noh")
			
			! Sod shock tube
			case("sod")
		end select

	end subroutine setup_test
end module module_imexport