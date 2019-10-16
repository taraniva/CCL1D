module module_geometry
	!DEPENDENT ON DATA
	contains

	subroutine setup_mesh(topo,vars,par)
		use module_data
		implicit none

		! I/O
		type(topo_type), intent(inout) :: topo
		type(para_type), intent(inout) :: par
		type(vars_type), intent(inout) :: vars

		! Local
		integer :: i
		integer :: nn, nc
		real(kind=rkind) :: x_lef, x_rig, len0

		nc = topo%nc
		nn = topo%nn

		x_lef = par%xminn
		x_rig = par%xmaxn

		len0 = abs(x_rig - x_lef)/nc

		vars%x_n(1) = x_lef
		vars%x_n(nn) = x_rig

		do i = 2,nn-1
			vars%x_n(i) = vars%x_n(i-1) + len0
		end do

		!do i = 1,nn
		!	print*, vars%x_n(i)
		!end do

	end subroutine setup_mesh

	subroutine calculate_geometry(topo,vars,par,bool_err)
		use module_data
		implicit none

		!In/Out
		type(topo_type), intent(inout) :: topo
		type(para_type), intent(inout) :: par
		type(vars_type), intent(inout) :: vars
		logical,		 intent(inout) :: bool_err

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

			x_lef = vars%x_n(ic)
			x_rig = vars%x_n(ic+1)

			cvol = abs(x_rig - x_lef)
			cxc = x_lef + 0.5_d*cvol

			if(cvol.le.0.0_d) then
				print*, "VOLUME OF CELL", ic ,"IS .LE. ZERO!!!!"
				bool_err = .FALSE.
				exit
			endif

			vars%x_c(ic) = cxc
			vars%vol_c(ic) = cvol
			
		end do

		print*, "EXIT calculate_geometry"
		return
	end subroutine calculate_geometry

	!/////////////////////////////////////////////////////////////////

	

end module module_geometry