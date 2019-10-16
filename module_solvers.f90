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
									* (vars_n%pre_n(ic+1)*vars_n%u_n(ic+1) - vars_n%pre_n(ic)*vars_n%u_n(ic))
			
			! Internal energy 
			vars_np1%eni_c(ic) = vars_np1%ent_c(ic) - 0.5_d*vars_np1%mass_c(ic)*(vars_np1%u_c(ic)**2)

			! Negative internal energy test
			if (vars_np1%eni_c(ic).le.0.0_d) then
				print*, "NEGATIVE INTERNAL ENERGY IN CELL ", ic
				bool_success = .FALSE.
				read*
			else
				
				!vars_np1%rho_c(ic) = vars_np1%mass_c(ic) / vars_np1%vol_c(ic)
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

	!/////////////////////////////////////////////////////////////////////////////////////
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
		!0.001_d

		return
	end function time_step



end module module_solvers