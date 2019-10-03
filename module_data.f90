module module_data

	include "precision.h"

	real(kind=rkind), parameter :: pi=acos(-1.0_d)

	type para_type
		!Init
		character(len=3) 	:: prob
		integer				:: nstep
		character(len=4)	:: method
		integer				:: timestep_type

		!Boundary conditions
		character(len=3)	:: bc_type_lef, bc_type_rig
		real(kind=rkind)	:: bc_val_lef, bc_val_rig

		!Control parameters
		real(kind=rkind)	:: time_fin
		real(kind=rkind)	:: time_prn
		real(kind=rkind)	:: time
		real(kind=rkind)	:: dt
		real(kind=rkind)	:: dt_init
		real(kind=rkind)	:: dt_prev
		real(kind=rkind)	:: cfl
		real(kind=rkind)	:: max_dt_change
		real(kind=rkind)	:: max_vol_change
		real(kind=rkind)	:: xminn, xmaxn

		!Output file parameters
		integer				:: file_output_number
		character(len=3)	:: file_output_string

	end type para_type

	type topo_type

		integer :: nc, nn

	end type topo_type

	type vars_type
		!Nodal variables
		real(kind=rkind), dimension(:), pointer   :: gamma_n
		real(kind=rkind), dimension(:), pointer   :: x_n
		real(kind=rkind), dimension(:), pointer   :: u_n
		real(kind=rkind), dimension(:), pointer   :: pre_n

		!Cell variables
		real(kind=rkind), dimension(:), pointer   :: gamma_c
		real(kind=rkind), dimension(:), pointer   :: x_c
		real(kind=rkind), dimension(:), pointer   :: u_c
		real(kind=rkind), dimension(:), pointer   :: vol_c
		real(kind=rkind), dimension(:), pointer   :: mass_c
		real(kind=rkind), dimension(:), pointer   :: rho_c
		real(kind=rkind), dimension(:), pointer   :: pre_c
		real(kind=rkind), dimension(:), pointer   :: eni_c
		real(kind=rkind), dimension(:), pointer   :: ent_c
		real(kind=rkind), dimension(:), pointer   :: ssp_c
		real(kind=rkind), dimension(:), pointer   :: du_c
		real(kind=rkind), dimension(:), pointer   :: dpre_c
	end type vars_type

	contains

	subroutine alloc_vars_arrays(topo,vars)

		! I/O
		type(topo_type), intent(inout)  :: topo
		type(vars_type), intent(inout)  :: vars

		!Local
		integer nc, nn

		nc = topo%nc
		nn = topo%nn

		print*, " - ALLOCATING VARS ARRAYS"

		!Nodal
		allocate(vars%gamma_n(1:nn))
		allocate(vars%x_n(1:nn))
		allocate(vars%u_n(1:nn))
		allocate(vars%pre_n(1:nn))

		!Cell
		allocate(vars%gamma_c(1:nc))
		allocate(vars%x_c(1:nc))
		allocate(vars%u_c(1:nc))
		allocate(vars%vol_c(1:nc))
		allocate(vars%mass_c(1:nc))
		allocate(vars%rho_c(1:nc))
		allocate(vars%pre_c(1:nc))
		allocate(vars%eni_c(1:nc))
		allocate(vars%ent_c(1:nc))
		allocate(vars%ssp_c(1:nc))
		allocate(vars%du_c(1:nc))
		allocate(vars%dpre_c(1:nc))

		print*, "VARS ALLOCATED"

	end subroutine alloc_vars_arrays

	subroutine dealloc_vars_arrays(vars)

		! I/O
		type(vars_type), intent(inout)  :: vars

		!Local

		print*, " - DEALLOCATING VARS ARRAYS"

		!Nodal
		deallocate(vars%gamma_n)
		deallocate(vars%x_n)
		deallocate(vars%u_n)
		deallocate(vars%pre_n)

		!Cell
		deallocate(vars%gamma_c)
		deallocate(vars%x_c)
		deallocate(vars%u_c)
		deallocate(vars%vol_c)
		deallocate(vars%mass_c)
		deallocate(vars%rho_c)
		deallocate(vars%pre_c)
		deallocate(vars%eni_c)
		deallocate(vars%ent_c)
		deallocate(vars%ssp_c)
		deallocate(vars%du_c)
		deallocate(vars%dpre_c)

		print*, "VARS DEALLOCATED"

	end subroutine dealloc_vars_arrays

end module module_data