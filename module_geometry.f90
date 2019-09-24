module module_geometry
	contains

	subroutine calculate_geometry(topo,par,xn,vol,xc)
		use module_data
		implicit none
		include "map_func.h"

		!In/Out
		type(topo_type), intent(inout) :: topo
		type(para_type), intent(inout) :: par
		real(kind=rkind), dimension(topo%nn), intent(in)	:: xn
		real(kind=rkind), dimension(topo%nc), intent(inout)	:: vol
		real(kind=rkind), dimension(topo%nc), intent(inout) :: xc

		!Local
		integer				:: ic
		real(kind=rkind)	:: cvol, cxc, x_lef, x_rig

		print*,"ENTERING calculate_geometry"

		!Cell volume calculation (1D == length)
		do ic = 1,topo%nc
			cvol = 0.0_d
			cxc  = 0.0_d
			x_lef = 0.0_d
			x_rig = 0.0_d

			x_lef = xn(ic)
			x_rig = xn(ic+1)

			cvol = abs(x_rig - x_lef)
			cxc = x_lef + 0.5_d*cvol

			xc(ic) = cxc
			vol(ic) = cvol
		end do

		print*, "EXIT calculate_geometry"
		return
	end subroutine calculate_geometry

	!/////////////////////////////////////////////////////////////////

	

end module module_geometry