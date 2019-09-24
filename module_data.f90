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
		character(len=3)	:: bc_type(2)
		real(kind=rkind), dimension(1:2)	:: bc_val

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
		real(kind=rkind), dimension(:), pointer   :: mass_n

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

end module module_data