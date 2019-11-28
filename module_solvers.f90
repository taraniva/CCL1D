module module_solvers

	! DEPENDENT ON DATA, EOS

	implicit none

	contains

	subroutine nodal_solver(topo,vars,vars_np1,par)

		! Solving a Riemann problem at the cell boundary to get nodal velocity

		use module_data

		! I/O
		type(topo_type), intent(inout) :: topo
		type(para_type), intent(inout) :: par
		type(vars_type), intent(inout) :: vars
		type(vars_type), intent(inout) :: vars_np1

		! Local
		integer :: in, nn, nc
		real(kind=rkind) :: rho_lef, u_lef, p_lef, rho_rig, u_rig, p_rig, z_lef, z_rig
		real(kind=rkind) :: z_denom, u_nodal, pre_nodal


		! Pointers definition

		! Nodal variables
		real(kind=rkind), dimension(:), pointer   :: xn => null()
		real(kind=rkind), dimension(:), pointer   :: un => null()
		real(kind=rkind), dimension(:), pointer   :: pren => null()

		! Cell variables
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
		xn => vars%x_n
		un => vars%u_n
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

		in = 0

		! Loop for internal nodes
		do in = 2,nn-1
			! Left state
			rho_lef = rho(in-1)
			u_lef = uc(in-1)
			p_lef = pre(in-1)
			z_lef = rho(in-1)*ssp(in-1)

			! Right state
			rho_rig = rho(in)
			u_rig = uc(in)
			p_rig = pre(in)
			z_rig = rho(in)*ssp(in)

			! Denominator z_i + z_i+1
			z_denom = z_lef + z_rig

			! Nodal velocities
			u_nodal = z_lef*u_lef + z_rig*u_rig - p_rig + p_lef
			u_nodal = u_nodal/z_denom
			un(in) = u_nodal

			! Nodal pressure
			pre_nodal = z_lef*p_rig + z_rig*p_lef - z_lef*z_rig*(u_rig - u_lef)
			pre_nodal = pre_nodal/z_denom
			pren(in) = pre_nodal

		end do

		! Left boundary
		if (par%bc_type_lef.eq.'tra') then
			! Prescribed pressure
			pren(1) = par%bc_val_lef
			!un(1) = uc(1) - (pre(1) - pren(1))/(rho(1)*ssp(1))
			un(1) = uc(1) - (pre(1) - par%bc_val_lef)/(rho(1)*ssp(1))
		elseif(par%bc_type_lef.eq.'vel') then
			! Prescribed velocity
			un(1) = par%bc_val_lef
			pren(1) = pre(1) + (rho(1)*ssp(1))*(un(1) - uc(1))
		endif

		! Right boundary
		if (par%bc_type_rig.eq.'tra') then
			! Prescribed pressure
			pren(nn) = par%bc_val_rig
			!un(nn) = uc(nc) + (pre(nc) - pren(nn))/(rho(nc)*ssp(nc))
			un(nn) = uc(nc) + (pre(nc) - par%bc_val_rig)/(rho(nc)*ssp(nc))
		elseif(par%bc_type_rig.eq.'vel') then
			! Prescribed velocity
			un(nn) = par%bc_val_rig
			pren(nn) = pre(nc) - (rho(nc)*ssp(nc))*(un(nn) - uc(nc))
		endif

		! Update nodal positions in vars_np1
		do in=1,nn
			vars_np1%x_n(in) = xn(in) + par%dt * un(in)
			!print*, "Nodal velocity"
			!print*, un(in)
			!print*, vars_np1%x_n(in)
		end do

		print*, "EXIT NODAL SOLVER"

	end subroutine nodal_solver

	!/////////////////////////////////////////////////////////////////////////

	subroutine calculate_cell_variables(topo,vars_n,vars_np1,par,bool_success)

		use module_data
		use module_eos

		!include "eos_func.h"
		!include "precision.h"

		! I/O
		type(topo_type), intent(inout) :: topo
		type(para_type), intent(inout) :: par
		type(vars_type), intent(inout) :: vars_n
		type(vars_type), intent(inout) :: vars_np1
		logical,		 intent(inout) :: bool_success

		! Local
		integer :: ic, nc
		real(kind=rkind) :: eta_np1, eta_n

		nc = topo%nc
		ic = 0

		
		! Loop over all cells
		do ic=1,nc
			!vars_np1%pre_c(ic) = 0.0_d

			select case(par%method)
		
			case('ph1o')
				! Specific density
				eta_np1 = 0.0
				eta_n = 1.0/vars_n%rho_c(ic)
				eta_np1 = eta_n &
						+ par%dt/vars_n%mass_c(ic)*(vars_n%u_n(ic+1)-vars_n%u_n(ic))
				vars_np1%rho_c(ic) = 1.0/eta_np1
				! Cell centre speed
				vars_np1%u_c(ic) = vars_n%u_c(ic) &
										- par%dt/vars_n%mass_c(ic)*(vars_n%pre_n(ic+1)-vars_n%pre_n(ic))
				! Total energy (density of total energy)
				vars_np1%ent_c(ic) = vars_n%ent_c(ic) - par%dt/vars_n%mass_c(ic)	&
									* ( vars_n%pre_n(ic+1)*vars_n%u_n(ic+1) - vars_n%pre_n(ic)*vars_n%u_n(ic) )

			case('ph2o')
				eta_np1 = 0.0
				eta_n = 1.0/vars_n%rho_c(ic)
				eta_np1 = eta_n &
							+ par%dt/vars_n%mass_c(ic)*(vars_n%u_np12(ic+1)-vars_n%u_np12(ic))
				vars_np1%rho_c(ic) = 1.0/eta_np1
				! Cell centre speed
				vars_np1%u_c(ic) = vars_n%u_c(ic) &
									- par%dt/vars_n%mass_c(ic)*(vars_n%pre_np12(ic+1)-vars_n%pre_np12(ic))
				! Total energy (density of total energy)
				vars_np1%ent_c(ic) = vars_n%ent_c(ic) - par%dt/vars_n%mass_c(ic)	&
									* ( vars_n%pre_np12(ic+1)*vars_n%u_np12(ic+1) - vars_n%pre_np12(ic)*vars_n%u_np12(ic) )
			end select
			
			! DENSITY of Internal energy!!
			vars_np1%eni_c(ic) = vars_np1%ent_c(ic) - 0.5_d*(vars_np1%u_c(ic)**2)

			! Negative internal energy test
			if (vars_np1%eni_c(ic).le.0.0_d) then
				print*, "NEGATIVE INTERNAL ENERGY IN CELL ", ic
				bool_success = .FALSE.
				read*
			else
				
				vars_np1%rho_c(ic) = vars_np1%mass_c(ic) / vars_np1%vol_c(ic)
				vars_np1%pre_c(ic) = eos_pre(vars_np1%rho_c(ic), vars_np1%eni_c(ic), vars_np1%gamma_c(ic))

				bool_success = .TRUE.
			endif
		end do

	end subroutine calculate_cell_variables

!/////////////////////////////////////////////////////////////////////////

	subroutine calculate_energy_balance(topo,vars)
		! Used at the end of computation to evaluate the success of used method

		use module_data

		! I/O
		type(topo_type), intent(inout) :: topo
		type(vars_type), intent(inout) :: vars

		! Local
		integer :: ic,nc
		real(kind=rkind) :: e_bala

		nc = topo%nc
		e_bala = 0.0_d

		print*, " - CALCULATING TOTAL ENERGY BALANCE"

		do ic=1,nc
			e_bala = e_bala + vars%ent_c(ic)
		end do

		print*,"---------------------------------------------------"
		print*,"TOTAL ENERGY IN SYSTEM: (", e_bala, ")"
		print*,"---------------------------------------------------"

	end subroutine calculate_energy_balance

	!/////////////////////////////////////////////////////////////////////////

	subroutine calculate_cell_sound_speed(topo,vars,par)

		! USE BEFORE NODAL SOLVER!!!

		use module_data
		use module_eos

		!include "eos_func.h"

		! I/O
		type(topo_type), intent(inout) :: topo
		type(para_type), intent(inout) :: par
		type(vars_type), intent(inout) :: vars

		! Local
		integer :: ic, nc

		print*, " - CELL SOUND SPEED COMPUTATION"

		nc = topo%nc
		ic = 0

		if (par%prob.eq.'sed') then
			! Compute from internal energy
			do ic=1,nc
				vars%ssp_c(ic) = eos_sspe(vars%eni_c(ic),vars%gamma_c(ic))
			end do
		else
			do ic=1,nc
				! Compute from density and pressure
				vars%ssp_c(ic) = eos_sspp(vars%rho_c(ic),vars%pre_c(ic),vars%gamma_c(ic))

				! Compute from internal energy
				!vars%ssp_c(ic) = eos_sspe(vars%eni_c(ic),vars%gamma_c(ic))
				!print*, vars%ssp_c(ic)
			end do
		endif

		print*, "EXIT CELL SOUND SPEED COMPUTATION"
	end subroutine calculate_cell_sound_speed

	!/////////////////////////////////////////////////////////////////////////
	!						PHM 2nd Order subroutines
	!/////////////////////////////////////////////////////////////////////////

	!/////////////////////////////////////////////////////////////////////////

	subroutine compute_slopes(topo,vars_n,par)

		use module_data

		! I/O
		type(topo_type), intent(inout) :: topo
		type(para_type), intent(inout) :: par
		type(vars_type), intent(inout) :: vars_n

		! Local
		integer :: ic, nc
		real(kind=rkind) :: delta_vel, delta_pre
		real(kind=rkind) :: x_lef, x_rig, u_lef, u_rig, pre_lef, pre_rig
		real(kind=rkind) :: xc, uc, prec
		real(kind=rkind) :: lsq_denom

		nc = topo%nc

		! Loop for internal cells
		! LSQ procedure -> look at immediate left and right neighbours of a cell (C-C variables)
		do ic=2,nc-1

			delta_vel = 0.0_d
			delta_pre = 0.0_d

			x_lef = vars_n%x_c(ic-1)
			x_rig = vars_n%x_c(ic+1)
			u_lef = vars_n%u_c(ic-1)
			u_rig = vars_n%u_c(ic+1)
			pre_lef = vars_n%pre_c(ic-1)
			pre_rig = vars_n%pre_c(ic+1)

			xc = vars_n%x_c(ic)
			uc = vars_n%u_c(ic)
			prec = vars_n%pre_c(ic)

			! Denominator in LSQ (convenience)
			lsq_denom = (x_lef - xc)**2 + (x_rig - xc)**2

			! LSQ itself
			delta_vel = ((x_lef - xc)*(u_lef - uc) + (x_rig - xc)*(u_rig - uc))/lsq_denom
			delta_pre = ((x_lef - xc)*(pre_lef - prec) + (x_rig - xc)*(pre_rig - prec))/lsq_denom

			! Shove it into the vars structure
			vars_n%slope_vel(ic) = delta_vel
			vars_n%slope_pre(ic) = delta_pre
		end do 

		! Boundary cells

		!--------------------- L E F T   B D R Y ---------------------

		ic = 1

		delta_vel = 0.0_d
		delta_pre = 0.0_d

		x_lef   = vars_n%x_n(ic)
		x_rig   = vars_n%x_c(ic+1)
		pre_rig = vars_n%pre_c(ic+1)
		u_rig   = vars_n%u_c(ic+1)
		xc = vars_n%x_c(ic)
		uc = vars_n%u_c(ic)
		prec = vars_n%pre_c(ic)

		if (par%bc_type_lef.eq.'tra') then

			pre_lef = par%bc_val_lef
			u_lef = uc - (prec - pre_lef)/(vars_n%rho_c(ic)*vars_n%ssp_c(ic))	
			
		elseif (par%bc_type_lef.eq.'vel') then

			u_lef   = par%bc_val_lef
			pre_lef = prec + vars_n%rho_c(ic)*vars_n%ssp_c(ic)*(u_lef - uc)

		endif

		lsq_denom = (x_lef - xc)**2 + (x_rig - xc)**2

		! LSQ itself
		delta_vel = ((x_lef - xc)*(u_lef - uc) + (x_rig - xc)*(u_rig - uc))/lsq_denom
		delta_pre = ((x_lef - xc)*(pre_lef - prec) + (x_rig - xc)*(pre_rig - prec))/lsq_denom

		vars_n%slope_vel(ic)  = delta_vel
		vars_n%slope_pre(ic)  = delta_pre
		!---------------------  E N D   O F   L E F T   B D R Y ---------------------

		!--------------------- R I G H T   B D R Y ---------------------
		ic = nc

		delta_vel = 0.0_d
		delta_pre = 0.0_d

		x_lef   = vars_n%x_c(ic-1)
		x_rig   = vars_n%x_n(ic+1)

		xc 		= vars_n%x_c(ic)
		uc 		= vars_n%u_c(ic)
		prec 	= vars_n%pre_c(ic)

		if (par%bc_type_rig.eq.'tra') then

			pre_rig = par%bc_val_rig
			u_rig   = uc + (prec - pre_rig)/(vars_n%rho_c(ic)*vars_n%ssp_c(ic))
			
		elseif (par%bc_type_rig.eq.'vel') then

			u_rig  = par%bc_val_rig
			pre_rig = prec - vars_n%rho_c(ic)*vars_n%ssp_c(ic)*(u_rig - uc)

		endif

		lsq_denom = (x_lef - xc)**2 + (x_rig - xc)**2

		! LSQ itself
		delta_vel = ((x_lef - xc)*(u_lef - uc) + (x_rig - xc)*(u_rig - uc))/lsq_denom
		delta_pre = ((x_lef - xc)*(pre_lef - prec) + (x_rig - xc)*(pre_rig - prec))/lsq_denom

		vars_n%slope_vel(ic) = delta_vel
		vars_n%slope_pre(ic) = delta_pre

		!---------------------  E N D   O F   R I G H T   B D R Y ---------------------
		
	end subroutine compute_slopes

	!/////////////////////////////////////////////////////////////////////////

	subroutine limit_slopes(topo, vars_n, par)
	
		use module_data

		! I/O
		type(topo_type), intent(inout) :: topo
		type(para_type), intent(inout) :: par
		type(vars_type), intent(inout) :: vars_n

		! Local
		integer :: ic, nc
		real(kind=rkind) :: min_diff_vel, min_diff_pre, max_diff_vel, max_diff_pre
		real(kind=rkind) :: vel_reconst_Ln, vel_reconst_Rn, pre_reconst_Ln, pre_reconst_Rn
		real(kind=rkind) :: Phi_pre_L, Phi_pre_R, Phi_vel_L, Phi_vel_R, Phi_pre, Phi_vel
		character(len=3) :: limiter_type

		! Pointers declaration
		real(kind=rkind), dimension(:), pointer   :: xc => null()
		real(kind=rkind), dimension(:), pointer   :: xn => null()
		real(kind=rkind), dimension(:), pointer   :: uc => null()
		real(kind=rkind), dimension(:), pointer   :: prec => null()
		real(kind=rkind), dimension(:), pointer   :: delta_vel => null()
		real(kind=rkind), dimension(:), pointer   :: delta_pre => null()
		
		! Pointers assignment
		xc 				=> vars_n%x_c
		xn 				=> vars_n%x_n
		uc 				=> vars_n%u_c
		prec 			=> vars_n%pre_c
		delta_vel 		=> vars_n%slope_vel
		delta_pre 		=> vars_n%slope_pre

		! Auxiliary parameters
		nc = topo%nc
		limiter_type = par%limiter_type
		

		! Loop for internal cells
		do ic = 2,nc-1
			! Velocity difference
			min_diff_vel = min(uc(ic+1) - uc(ic), uc(ic-1) - uc(ic), 0.0_d)
			max_diff_vel = max(uc(ic+1) - uc(ic), uc(ic-1) - uc(ic), 0.0_d)

			! Pressure difference
			min_diff_pre = min(prec(ic+1) - prec(ic), prec(ic-1) - prec(ic), 0.0_d)
			max_diff_pre = max(prec(ic+1) - prec(ic), prec(ic-1) - prec(ic), 0.0_d)

			! Find reconstructed values at the cell nodes
			vel_reconst_Ln = uc(ic) + delta_vel(ic) * (xn(ic) - xc(ic))
			vel_reconst_Rn = uc(ic) + delta_vel(ic) * (xn(ic+1) - xc(ic))

			pre_reconst_Ln = prec(ic) + delta_pre(ic) * (xn(ic) - xc(ic))
			pre_reconst_Rn = prec(ic) + delta_pre(ic) * (xn(ic+1) - xc(ic))

			! Left --------------------------------------------------
			! Velocity
			if( (vel_reconst_Ln - uc(ic)).gt.0.0_d ) then
				Phi_vel_L = limiter(max_diff_vel/(vel_reconst_Ln - uc(ic)), limiter_type)
			elseif ( (vel_reconst_Ln - uc(ic)).lt.0.0_d) then
				Phi_vel_L = limiter(min_diff_vel/(vel_reconst_Ln - uc(ic)), limiter_type)
			elseif ( (vel_reconst_Ln - uc(ic)).eq.0.0_d ) then
				Phi_vel_L = 1.0_d
			endif

			! Pressure
			if( (pre_reconst_Ln - prec(ic)).gt.0.0_d ) then
				Phi_pre_L = limiter(max_diff_pre/(pre_reconst_Ln - prec(ic)), limiter_type)
			elseif ( (pre_reconst_Ln - prec(ic)).lt.0.0_d) then
				Phi_pre_L = limiter(min_diff_pre/(pre_reconst_Ln - prec(ic)), limiter_type)
			elseif ( (pre_reconst_Ln - prec(ic)).eq.0.0_d ) then
				Phi_pre_L = 1.0_d
			endif
			
			! Right --------------------------------------------------
			! Velocity
			if( (vel_reconst_Rn - uc(ic)).gt.0.0_d ) then
				Phi_vel_R = limiter(max_diff_vel/(vel_reconst_Rn - uc(ic)), limiter_type)
			elseif ( (vel_reconst_Rn - uc(ic)).lt.0.0_d) then
				Phi_vel_R = limiter(min_diff_vel/(vel_reconst_Rn - uc(ic)), limiter_type)
			elseif ( (vel_reconst_Rn - uc(ic)).eq.0.0_d ) then
				Phi_vel_R = 1.0_d
			endif

			! Pressure
			if( (pre_reconst_Rn - prec(ic)).gt.0.0_d) then
				Phi_pre_R = limiter(max_diff_pre/(pre_reconst_Rn - prec(ic)), limiter_type)
			elseif ( (pre_reconst_Rn - prec(ic)).lt.0.0_d) then
				Phi_pre_R = limiter(min_diff_pre/(pre_reconst_Rn - prec(ic)), limiter_type)
			elseif ( (pre_reconst_Rn - prec(ic)).eq.0.0_d ) then
				Phi_pre_R = 1.0_d
			endif

			! Find the actual limiter
			Phi_vel = min(Phi_vel_L, Phi_vel_R)
			Phi_pre = min(Phi_pre_L, Phi_pre_R)

			! Put into arrays
			vars_n%limcoeff_vel(ic) = Phi_vel
			vars_n%limcoeff_pre(ic) = Phi_pre
			!vars_n%limcoeff_vel(ic) = 1.0_d
			!vars_n%limcoeff_pre(ic) = 1.0_d

		end do

	end subroutine limit_slopes
	

	!/////////////////////////////////////////////////////////////////////////

	subroutine nodal_solver_reconst(topo,vars_n,vars_np1,par)
		! Calculates nodal vars from cell linear var. reconstructions

		use module_data

		! I/O
		type(topo_type), intent(inout) :: topo
		type(para_type), intent(inout) :: par
		type(vars_type), intent(inout) :: vars_n
		type(vars_type), intent(inout) :: vars_np1

		! Pointers definition
		real(kind=rkind), dimension(:), pointer   :: xn => null()
		real(kind=rkind), dimension(:), pointer   :: delta_vel => null()
		real(kind=rkind), dimension(:), pointer   :: delta_pre => null()
		real(kind=rkind), dimension(:), pointer   :: Phi_vel => null()
		real(kind=rkind), dimension(:), pointer   :: Phi_pre => null()
		real(kind=rkind), dimension(:), pointer   :: xc => null()
		real(kind=rkind), dimension(:), pointer   :: uc => null()
		real(kind=rkind), dimension(:), pointer   :: prec => null()

		! Local vars
		integer :: nn, nc, in
		real(kind=rkind) :: rho_lef, z_lef, u_lef, pre_lef
		real(kind=rkind) :: rho_rig, z_rig, u_rig, pre_rig
		real(kind=rkind) :: z_denom
		real(kind=rkind) :: u_nodal, pre_nodal

		! Pointers assignment
		xn	=> vars_n%x_n
		delta_vel => vars_n%slope_vel
		delta_pre => vars_n%slope_pre
		Phi_vel => vars_n%limcoeff_vel
		Phi_pre => vars_n%limcoeff_pre
		xc => vars_n%x_c
		uc => vars_n%u_c
		prec => vars_n%pre_c

		nc = topo%nc
		nn = topo%nn

		! Internal nodes
		do in=2,nn-1
			! Left state
			rho_lef = vars_n%rho_c(in-1)
			z_lef = vars_n%rho_c(in-1)*vars_n%ssp_c(in-1)
			u_lef = vars_n%u_c(in-1) + Phi_vel(in-1)*delta_vel(in-1)*(xn(in)-xc(in-1))
			pre_lef = vars_n%pre_c(in-1) + Phi_pre(in-1)*delta_pre(in-1)*(xn(in)-xc(in-1))

			! Right state
			rho_rig = vars_n%rho_c(in)
			z_rig = vars_n%rho_c(in)*vars_n%ssp_c(in)
			u_rig = vars_n%u_c(in) + Phi_vel(in)*delta_vel(in)*(xn(in)-xc(in))
			pre_rig = vars_n%pre_c(in) + Phi_pre(in)*delta_pre(in)*(xn(in)-xc(in))

			! "Determinant"
			z_denom = z_lef + z_rig

			! Nodal calculation
			u_nodal = z_lef*u_lef + z_rig*u_rig - pre_rig + pre_lef
			u_nodal = u_nodal / z_denom
			pre_nodal = z_lef*pre_rig + z_rig*pre_lef - z_lef*z_rig*(u_rig - u_lef)
			pre_nodal = pre_nodal / z_denom

			! Shove into the data structure
			vars_n%u_n(in) = u_nodal
			vars_n%pre_n(in) = pre_nodal
		end do

		!--------------------- L E F T   B D R Y ---------------------
		in = 1
		z_rig = vars_n%rho_c(in)*vars_n%ssp_c(in)
		! Reconstruct in the right cell
		u_rig = vars_n%u_c(in) + Phi_vel(in)*delta_vel(in)*(xn(in)-xc(in))
		pre_rig = vars_n%pre_c(in) + Phi_pre(in)*delta_pre(in)*(xn(in)-xc(in))

		if (par%bc_type_lef.eq.'tra') then
			! Prescribed pressure
			pre_nodal = par%bc_val_lef

			! ---------- Hint from 1o nodal solver ----------
			!un(1) = uc(1) - (pre(1) - pren(1))/(rho(1)*ssp(1))

			u_nodal = uc(1) - (pre_rig - par%bc_val_lef)/z_rig
		elseif(par%bc_type_lef.eq.'vel') then
			! Prescribed velocity
			u_nodal = par%bc_val_lef

			! ---------- Hint from 1o nodal solver ----------
			!pren(1) = pre(1) + (rho(1)*ssp(1))*(un(1) - uc(1))

			pre_nodal = pre_rig + z_rig*(par%bc_val_lef - u_rig)
		endif

		! Shove into the data structure
		vars_n%u_n(in) = u_nodal
		vars_n%pre_n(in) = pre_nodal
		!---------------------  E N D   O F   L E F T   B D R Y ---------------------

		!--------------------- R I G H T   B D R Y ---------------------
		in = nn

		z_lef = vars_n%rho_c(in-1)*vars_n%ssp_c(in-1)
		! Reconstruct in the right cell
		u_lef = vars_n%u_c(in-1) + Phi_vel(in-1)*delta_vel(in-1)*(xn(in)-xc(in-1))
		pre_lef = vars_n%pre_c(in-1) + Phi_pre(in-1)*delta_pre(in-1)*(xn(in)-xc(in-1))
		if (par%bc_type_rig.eq.'tra') then
			! Prescribed pressure
			pre_nodal = par%bc_val_rig

			!un(nn) = uc(nc) + (pre(nc) - pren(nn))/(rho(nc)*ssp(nc))

			u_nodal = u_lef + (pre_lef - par%bc_val_rig)/z_lef
		elseif(par%bc_type_rig.eq.'vel') then
			! Prescribed velocity
			u_nodal = par%bc_val_rig
			pre_nodal = pre_lef - z_lef*(u_nodal - u_lef)
		endif

		! Shove into the data structure
		vars_n%u_n(in) = u_nodal
		vars_n%pre_n(in) = pre_nodal
		!---------------------  E N D   O F   R I G H T   B D R Y ---------------------

		! Move the mesh
		do in=1,nn
		vars_np1%x_n(in) = vars_n%x_n(in) + par%dt*vars_n%u_n(in)
		end do

	end subroutine nodal_solver_reconst

	!/////////////////////////////////////////////////////////////////////////

	subroutine dnodal_solver(topo, vars_n, par)

		use module_data
		
		! I/O
		type(topo_type), intent(inout) :: topo
		type(para_type), intent(inout) :: par
		type(vars_type), intent(inout) :: vars_n

		! Pointers definition
		real(kind=rkind), dimension(:), pointer   :: delta_vel => null()
		real(kind=rkind), dimension(:), pointer   :: delta_pre => null()
		real(kind=rkind), dimension(:), pointer   :: Phi_vel => null()
		real(kind=rkind), dimension(:), pointer   :: Phi_pre => null()

		! Local
		integer :: in, nn
		real(kind=rkind) :: ssp_lef, delta_pre_lef, delta_vel_lef, z_lef
		real(kind=rkind) :: ssp_rig, delta_pre_rig, delta_vel_rig, z_rig
		real(kind=rkind) :: z_denom
		real(kind=rkind) :: du_nodal, dpre_nodal

		! Pointers assignment
		delta_vel => vars_n%slope_vel
		delta_pre => vars_n%slope_pre
		Phi_vel => vars_n%limcoeff_vel
		Phi_pre => vars_n%limcoeff_pre

		nn = topo%nn

		! Internal nodes
		do in=2,nn-1

			! Left state
			ssp_lef 		= vars_n%ssp_c(in-1)
			delta_pre_lef	= Phi_pre(in-1)*delta_pre(in-1)
			delta_vel_lef	= Phi_vel(in-1)*delta_vel(in-1)
			z_lef 			= vars_n%ssp_c(in-1)*vars_n%rho_c(in-1)

			! Right state
			ssp_rig			= vars_n%ssp_c(in)
			delta_pre_rig 	= Phi_pre(in)*delta_pre(in)
			delta_vel_rig	= Phi_vel(in)*delta_vel(in)
			z_rig 			= vars_n%ssp_c(in)*vars_n%rho_c(in)

			! "Determinant"
			z_denom = z_lef + z_rig

			dpre_nodal 	= ssp_rig*z_lef*(delta_pre_rig - z_rig*delta_vel_rig) & 
								- ssp_lef*z_rig*(delta_pre_lef + z_lef*delta_vel_lef)

			du_nodal 	= (-1.0_d)*( ssp_lef*(delta_pre_lef + z_lef*delta_vel_lef) &
								+ ssp_rig*(delta_pre_rig - z_rig*delta_vel_rig) )
			dpre_nodal = dpre_nodal / z_denom
			du_nodal = du_nodal / z_denom

			! Shove into the data structure
			vars_n%du_n(in) = du_nodal
			vars_n%dpre_n(in) = dpre_nodal
		end do

		!--------------------- L E F T   B D R Y ---------------------
		! 					IS THIS CORRECT???????
		in = 1
		!Right cell state
		ssp_rig			= vars_n%ssp_c(in)
		delta_pre_rig 	= Phi_pre(in)*delta_pre(in)
		delta_vel_rig	= Phi_vel(in)*delta_vel(in)
		z_rig 			= vars_n%ssp_c(in)*vars_n%rho_c(in)

		dpre_nodal = 0.0_d
		du_nodal = (ssp_rig*z_rig*delta_vel_rig - ssp_rig*delta_pre_rig)/z_rig
		
		! Shove into the data structure
		vars_n%du_n(in) = du_nodal
		vars_n%dpre_n(in) = dpre_nodal


		!---------------------  E N D   O F   L E F T   B D R Y ---------------------

		!--------------------- R I G H T   B D R Y ---------------------
		in = nn
		!Left cell state
		ssp_lef			= vars_n%ssp_c(in-1)
		delta_pre_lef 	= Phi_pre(in-1)*delta_pre(in-1)
		delta_vel_lef	= Phi_vel(in-1)*delta_vel(in-1)
		z_lef 			= vars_n%ssp_c(in-1)*vars_n%rho_c(in-1)

		dpre_nodal = 0.0_d
		du_nodal = (-1.0_d)*(ssp_lef*delta_pre_lef + ssp_lef*z_lef*delta_vel_lef)

		! Shove into the data structure
		vars_n%du_n(in) = du_nodal
		vars_n%dpre_n(in) = dpre_nodal
		
		!---------------------  E N D   O F   R I G H T   B D R Y ---------------------


		! Advance to half-time
		do in=1,nn
		vars_n%u_np12(in) = vars_n%u_n(in) + 0.5_d*par%dt*vars_n%du_n(in)
		vars_n%pre_np12(in) = vars_n%pre_n(in) + 0.5_d*par%dt*vars_n%dpre_n(in)
		end do

	
	end subroutine dnodal_solver

	!/////////////////////////////////////////////////////////////////////////


	!/////////////////////////////////////////////////////////////////////////
	!						T  I  M  E  S  T  E  P
	!/////////////////////////////////////////////////////////////////////////
	function time_step(topo,vars,par)

		use module_data

		real(kind=rkind)             :: time_step

		! I/O
		type(topo_type), intent(inout) :: topo
		type(para_type), intent(inout) :: par
		type(vars_type), intent(inout) :: vars

		! Local
		integer :: ic, nc
		real(kind=rkind) :: dt_e, C_m, C_fl, cell_rat_min, phm_1o_rat

		! Initialise loop
		nc = topo%nc
		C_fl = par%cfl
		C_m = par%max_dt_change
		phm_1o_rat = vars%vol_c(1)/vars%ssp_c(1)
		cell_rat_min = phm_1o_rat

		! Find minimal value for the PHM control ratio
		do ic=2,nc
			phm_1o_rat = vars%vol_c(ic)/vars%ssp_c(ic)
			cell_rat_min = min(cell_rat_min, phm_1o_rat)
		end do

		time_step = min(C_fl*cell_rat_min, C_m*par%dt_prev)
		!time_step = 5e-6_d
		!0.001_d

		return
	end function time_step

	!/////////////////////////////////////////////////////////////////////////
	!						A U X I L I A R Y
	!/////////////////////////////////////////////////////////////////////////

	! Vankatarishnan limiter
	function vkr_limiter(var)

		use module_data

		real(kind=rkind) :: vkr_limiter

		!I/O
		real(kind=rkind), intent(in) :: var

		!print*, "Vankatarishnan limiter"

		vkr_limiter = (var**2 + 2.0_d*var)/(var**2 + var + 2.0_d)

	end function vkr_limiter

	! Barth-Jespersen limiter
	function baj_limiter(var)

		use module_data

		real(kind=rkind) :: baj_limiter

		!I/O
		real(kind=rkind), intent(in) :: var

		!print*, "Barth-Jespersen limiter"

		baj_limiter = min(1.0_d, var)

	end function baj_limiter

	! Limiter switch function
	function limiter(var, char_limiter_type)

		use module_data

		real(kind=rkind) :: limiter

		!I/O
		real(kind=rkind), intent(in) :: var
		character(len=3), intent(in) :: char_limiter_type

		print*, "Entering limiter switch"

		select case(char_limiter_type)
		case('baj')
			limiter = baj_limiter(var)
		case('vkr')
			limiter = vkr_limiter(var)
		case('non')
			limiter = 1.0_d
		end select

	end function limiter

end module module_solvers