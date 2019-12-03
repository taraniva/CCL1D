module module_imexport
	!Dependent on DATA nad GEOMETRY

	implicit none

	contains

	subroutine copy_all_vars(topo, vars_from, vars_to)
		use module_data
		implicit none

		! I/O
		type(topo_type), intent(inout)  :: topo
  		type(vars_type), intent(inout)  :: vars_from
		type(vars_type), intent(inout)  :: vars_to
		
		! Local
		integer :: in, ic

		print*, " - COPYING ALL THE STUFF..."

		! Nodal
		do in=1,topo%nn
			vars_to%gamma_n(in) = vars_from%gamma_n(in)
			vars_to%x_n(in) = vars_from%x_n(in)
			vars_to%u_n(in) = vars_from%u_n(in)
			vars_to%pre_n(in) = vars_from%pre_n(in)
			vars_to%du_n(in) = vars_from%du_n(in)
			vars_to%dpre_n(in) = vars_from%dpre_n(in)
		end do

		! Cell
		do ic=1,topo%nc
			vars_to%gamma_c(ic) = vars_from%gamma_c(ic)
			vars_to%x_c(ic) = vars_from%x_c(ic)
			vars_to%u_c(ic) = vars_from%u_c(ic)
			vars_to%vol_c(ic) = vars_from%vol_c(ic)
			vars_to%mass_c(ic) = vars_from%mass_c(ic)
			vars_to%rho_c(ic) = vars_from%rho_c(ic)
			vars_to%pre_c(ic) = vars_from%pre_c(ic)
			vars_to%eni_c(ic) = vars_from%eni_c(ic)
			vars_to%ent_c(ic) = vars_from%ent_c(ic)
			vars_to%ssp_c(ic) = vars_from%ssp_c(ic)
			! Not reallly necessary
			vars_to%slope_vel(ic) = vars_from%slope_vel(ic)
			vars_to%slope_pre(ic) = vars_from%slope_pre(ic)
		end do

		print*, "STUFF COPIED SUCCESSFULLY"

	end subroutine copy_all_vars


!//////////////////////////////////////////////////////////////////////////////////////////////////

	subroutine setup_test(topo,vars,par)
		use module_data
		use module_geometry
		use module_eos

		implicit none

		!include "eos_func.h"
  		! I/O
  		type(topo_type), intent(inout)  :: topo
  		type(vars_type), intent(inout)  :: vars
		type(para_type), intent(inout)  :: par

		!Local
		integer :: ic
		integer :: nc, nn
		real(kind=rkind) :: gam0
		logical :: bool_placeholder
		! Pointers definition
		!Nodal variables
		real(kind=rkind), dimension(:), pointer   :: pren => null()

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
		pren => vars%pre_n
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
		bool_placeholder = .TRUE.
		  
		select case(par%prob)
			! Noh implosion test
			case("noh")
				par%time_fin = 0.6_d
				par%time_prn = 0.1_d

				par%xminn = 0.0_d
				par%xmaxn = 1.0_d

				call setup_mesh(topo,vars,par)
				call calculate_geometry(topo,vars,par,bool_placeholder)

				! Define cell variables
				do ic=1,nc
				end do

			
			! Sod shock tube
			case("sod")
				par%time_fin = 0.20_d
				par%time_prn = 0.05_d
				
				par%xminn = 0.0_d
				par%xmaxn = 1.0_d

				call setup_mesh(topo,vars,par)
				call calculate_geometry(topo,vars,par,bool_placeholder)

				gam0 = 1.4_d
				gamma = gam0

				do ic=1,nc
					if(xc(ic) .le. 0.5_d) then
						rho(ic) = 1.0_d
						pre(ic) = 1.0_d
						uc(ic) = 0.0_d
					elseif(xc(ic) .gt. 0.5_d) then
						rho(ic) = 0.125_d
						pre(ic) = 0.1_d
						uc(ic) = 0.0_d
					endif
					eni(ic) = eos_eni(rho(ic),pre(ic),gam0)
				end do
			case("vil")
				par%time_fin = 0.8_d
				par%time_prn = 0.1_d
				
				par%xminn = 0.0_d
				par%xmaxn = 1.0_d

				call setup_mesh(topo,vars,par)
				call calculate_geometry(topo,vars,par,bool_placeholder)

				gam0 = 3.0_d
				gamma = gam0

				do ic=1,nc
					rho(ic) = 1.0 + 0.1*sin(2*pi*xc(ic))
					pre(ic) = rho(ic)**gam0
					uc(ic)  = 0.0_d
					eni(ic) = eos_eni(rho(ic),pre(ic),gam0)
				enddo

			case default
				print*,'SETUP_VARS: unknown test ',par%prob
				stop
		end select

		do ic = 1,nc
			mass(ic) = vol(ic) * rho(ic)
			ent(ic)  = eni(ic) + 0.5_d * (uc(ic)**2)
			ssp(ic)  = eos_sspp(rho(ic),pre(ic),gam0)
		enddo

	end subroutine setup_test

!//////////////////////////////////////////////////////////////////////////////////////////////////

	subroutine read_input_file(topo,vars,par,filename)
		! Read parameters from input file
		use module_data
		implicit none

		! I/O
		type(topo_type), intent(inout)  :: topo
		type(vars_type), intent(inout)  :: vars
		type(para_type), intent(inout)  :: par
		character(len=5), intent(in)    :: filename

		! Local vars
		integer :: nn
		
		print*," - READING INPUT FILE"
	  
		open(unit=8, file=filename,access='sequential',status='old')
		read(8,*) !First line is a comment which denotes the test

		! Test problem
		read(8,*) par%prob
		print*, "PROBLEM NAME:", par%prob

		! Mesh topology
		read(8,*) topo%nc
		nn = topo%nc + 1
		topo%nn = nn
		print*, "NUMBER OF CELLS:", topo%nc
		print*, "NUMBER OF NODES:", topo%nn

		! Boundary conditions
		read(8,*) par%bc_type_lef, par%bc_val_lef
		read(8,*) par%bc_type_rig, par%bc_val_rig
		print*, "-----------------------------------------"
		print*, "LEFT BOUNDARY CONDITION:"
		print*, "Type: ", par%bc_type_lef, " || Value: ", par%bc_val_lef
		print*, "RIGHT BOUNDARY CONDITION:"
		print*, "Type: ", par%bc_type_rig, " || Value: ", par%bc_val_rig
		print*, "-----------------------------------------"

		! Control parameters
		read(8,*) par%method
		read(8,*) par%dt_init
		read(8,*) par%cfl
		read(8,*) par%max_vol_change
		read(8,*) par%max_dt_change
		read(8,*) par%limiter_type
		print*, "METHOD (", par%method, ")"
		print*, "Initial timestep length (", par%dt_init, ")"
		print*, "CFL (", par%cfl, ")"
		print*, "Max volume change (", par%max_vol_change, ")" 
		print*, "Max timestep change (", par%max_dt_change, ")"
		if(par%method.eq.'ph2o') then
			select case(par%limiter_type)
			case('baj')
				print*, "Limiter type ( Barth-Jespersen )"
			case('vkr')
				print*, "Limiter type ( Vankatarishnan )"
			case('non')
				print*, "!!! NO SLOPE LIMITER !!!"
			end select
		endif
		print*, "INITIAL PARAMETERS LOADED SUCCESSFULLY"
		close(8)
		!print*, "PRESS ENTER TO CONTINUE..."
		!read*

	end subroutine read_input_file

!//////////////////////////////////////////////////////////////////////////////////////////////////

	subroutine get_next_fnum_string(fnum,fnumstr)
		! Increase number of file and return it as a three character string with leading zeros
		implicit none
		! I/O
		integer, intent(inout)          :: fnum
		character(len=3), intent(inout) :: fnumstr
		! Local vars
		integer   :: ndig
		integer   :: i, ndig0
	  
		ndig = 3
	  
		if (fnum > 998) then
		   print*,'GET_FNUM_STRING: Too many output files'
		   print*,'                 Reusing files No. 999'
	  !     read*
		   fnum = 999
		else
		   fnum = fnum+1
		endif
		
		! Convert num to a ndig-character string with leading '0'
		ndig = 3
		write(fnumstr,'(i3)') fnum 
		fnumstr= trim(adjustl(fnumstr))
		ndig0 = ndig-len(trim(adjustl(fnumstr))) 
		do i=1,ndig0
		   fnumstr = '0' // fnumstr
		enddo
	  
		return
	  end subroutine get_next_fnum_string

	  !//////////////////////////////////////////////////////////////////////////////////////////////////

	  subroutine write_cell_vals(topo,vars,par,filenamebase)
		use module_data
		implicit none

		! I/O
		type(topo_type), intent(inout)  :: topo
		type(vars_type), intent(inout)  :: vars
		type(para_type), intent(inout)  :: par
		character(len=4), intent(in)    :: filenamebase

		! Local vars
		character(len=12) :: filename
		integer           :: ic
		real(kind=rkind), dimension(1:5) :: round

		! Pointers definition for -> vars
		real(kind=rkind), dimension(:), pointer :: xc => null()
		real(kind=rkind), dimension(:), pointer :: uc  => null()
		real(kind=rkind), dimension(:), pointer   :: rho => null()
		real(kind=rkind), dimension(:), pointer   :: pre => null()
		real(kind=rkind), dimension(:), pointer   :: eni => null()
	  
		! Pointers assignment for -> vars%
		xc  => vars%x_c
		uc   => vars%u_c
		rho  => vars%rho_c
		pre  => vars%pre_c
		eni  => vars%eni_c
	  
		
		filename = filenamebase // '_' // par%file_output_string // '.dat'
		open(10,file=filename)
	  
		do ic = 1,topo%nc   
		   ! prevent exponent overflow
		   round(1) = xc(ic)
		   round(2) = uc(ic)
		   round(3)   = rho(ic)
		   round(4)   = pre(ic)
		   round(5)   = eni(ic)
		   where (abs(round).le.1.0e-30_d) 
			  round = 0.0_d
		   end where
		   ! write it
		   write(10,917) round(1:5)
		end do
		
		close(10)
	  
	  917 format(1pe15.6,6('  ',1pe15.6))
	  
		return
	  
	  end subroutine write_cell_vals

	  !//////////////////////////////////////////////////////////////////////////////////////////////////

	  subroutine write_node_coords(topo,vars,par,filenamebase)
		! Save the mesh (node positions) to file
		use module_data
		implicit none

		! I/O
		type(topo_type), intent(inout)  :: topo
		type(vars_type), intent(inout)  :: vars
		type(para_type), intent(inout)  :: par
		character(len=4), intent(in)    :: filenamebase

		! Local vars
		character(len=12) :: filename
		integer           :: in
		real(kind=rkind), dimension(1) :: round
		! Pointers definition for -> vars
		real(kind=rkind), dimension(:), pointer :: xn => null()
	  
		! Pointers assignment for -> vars%
		xn => vars%x_n
	  
		filename = filenamebase // '_' // par%file_output_string // '.dat'
		open(10,file=filename)
	  
		do in = 1,topo%nn   
		   ! prevent exponent overflow
		   round(:) = xn(in)
		   where (abs(round).le.1.0e-30_d) 
			  round = 0.0_d
		   end where
		   ! write it
		   write(10,912) round(1)
		end do
		
		close(10)
	  
	  912 format(1pe15.6,'  ',1pe15.6)
	  
		return
	  
	  end subroutine write_node_coords

end module module_imexport