program CCL1D
	
	use module_data
	use module_geometry
	use module_imexport
	use module_solvers

	implicit none

	type(topo_type), pointer :: topo     => null()  
  	type(vars_type), pointer :: vars_n   => null()  
  	type(vars_type), pointer :: vars_np1 => null()  
  	type(para_type), pointer :: par      => null()

	integer :: pseudobool_output_new, ic
	real(kind=rkind) :: time_since_prn, dt_min
	logical :: bool_lagstep_successful

	! Memory allocation
	allocate(topo)
 	allocate(vars_n,vars_np1)
	allocate(par)

	! Read initial file
	call read_input_file(topo,vars_n,par,'INPUT')

	call alloc_vars_arrays(topo,vars_n)
	call alloc_vars_arrays(topo,vars_np1)

	call setup_test(topo,vars_n,par)

	! Starting point status
	par%file_output_number = -1  ! File counter, so that the first is indexed is 0
	call get_next_fnum_string(par%file_output_number,par%file_output_string)

	!---------------- I N I T    M A I N    L O O P ----------------
	par%nstep      = 0          ! Number of timestep
  	par%time       = 0.0_d      ! Actual time
  	time_since_prn = 0.0_d      ! Time since last output
  	par%dt         = par%dt_init ! User-provided length of first timestep
	dt_min         = par%time_fin / 1.e9_d ! mimimum allowed timestep length

	print*, "_________________________________________________________"
  	print*, ' Number of cells: ', topo%nc
  	print*, ' Number of nodes: ', topo%nn
  	print*,' Writing file: vals_',par%file_output_string,' at initial time'
  	print*,' Writing file: mesh_',par%file_output_string,' at initial time'
  	print*, "_________________________________________________________"
  	call write_node_coords(topo,vars_n,par,'mesh')
  	call write_cell_vals(topo,vars_n,par,'vals')

	! Initialise with two identical datasets
	call copy_all_vars(topo,vars_n,vars_np1)
	call calculate_energy_balance(topo,vars_n)

	do while(par%time.lt.par%time_fin)

		pseudobool_output_new = 0
		par%nstep = par%nstep + 1
		par%dt_prev = par%dt

		par%dt = time_step(topo,vars_n,par) ! <- in module_solvers.f90
		
		if (par%dt < dt_min) then
			print*,'> Time step too small: ',par%dt,' at step ',par%nstep,', time=',par%time
			exit
		endif

		if (par%time + par%dt .gt. par%time_fin) par%dt = par%time_fin-par%time

		bool_lagstep_successful = .FALSE.
		! (Repeat until lagstep successful) .or. (exit if dt is too small)
		do while (bool_lagstep_successful.eqv..FALSE.)!***************************************************

			!call lagstep()
			!*************************** P H M   1 o ***************************

			call calculate_geometry(topo,vars_n,par,bool_lagstep_successful)

			call calculate_cell_sound_speed(topo,vars_n,par)
			!call calculate_cell_sound_speed(topo,vars_np1,par)

			call nodal_solver(topo,vars_n,vars_np1,par)
			!call calculate_cell_sound_speed(topo,vars_np1,par)

			call calculate_geometry(topo,vars_np1,par,bool_lagstep_successful)

			call calculate_cell_variables(topo,vars_n,vars_np1,par,bool_lagstep_successful)

			! If failed, restart with dt/2
			if (bool_lagstep_successful.eqv..FALSE.) then
				par%dt = par%dt/2.0_d
				print*,'> Lagstep restart with dt= ',par%dt,' at step ',par%nstep,', time=',par%time
			
				! Stop if Lagr. step too short
				if ( par%dt < dt_min) then
				   print*,'> Time step too small: ',par%dt,' at step ',par%nstep,', time=',par%time
				   exit
				endif
			endif

			if (bool_lagstep_successful.eqv..FALSE.) exit 

		end do !*************************************************************************************


		! Update time
		par%time = par%time + par%dt

     	if (mod(par%nstep,1)==0) write(*,999) par%nstep, par%dt, par%time
		999  format('Finished ', i5,'. timestep, dt=', 1pe13.6,', t=', 1pe13.6)
		time_since_prn = time_since_prn + par%dt

		! Copy all vals for next step
		call copy_all_vars(topo,vars_np1,vars_n)

		! In-between records:
		if (time_since_prn .gt. par%time_prn) then

			call get_next_fnum_string(par%file_output_number,par%file_output_string)		! <- module_imexport.f90
			print*,'Results located in mesh/vals_', par%file_output_string,' at time t=', par%time
		
			call write_node_coords(topo,vars_n,par,'mesh')		! <- module_imexport.f90
			call write_cell_vals(topo,vars_n,par,'vals')	    ! <- module_imexport.f90
			
			time_since_prn = modulo( par%time , par%time_prn )
			pseudobool_output_new = 1
		endif

		print*, "_______________________________________________________________"
     	print*, "Next step"
	 	print*, "_______________________________________________________________"
	 
	 	print*, "Time step" ,par%dt
		 print*, "nstep:", par%nstep

		! do ic = 1,topo%nc
		! 	print*, vars_n%rho_c(ic)
		! 	print*, vars_n%u_c(ic)
		! 	print*, vars_n%ent_c(ic)
		! 	print*, vars_n%eni_c(ic)

		! 	print*, vars_np1%rho_c(ic)
		! 	print*, vars_np1%u_c(ic)
		! 	print*, vars_np1%ent_c(ic)
		! 	print*, vars_np1%eni_c(ic)
		! enddo
	end do

	print*, "********************************"
  	print*, "F I N A L   R E S U L T S"
	print*, "********************************"
	  

	call calculate_energy_balance(topo,vars_n)

	if (pseudobool_output_new==0) then
		call get_next_fnum_string(par%file_output_number,par%file_output_string)				! <- module_imexport.f90
		print*,' Writing files (mesh/vals)_',par%file_output_number,' at time t=',par%time
		call write_node_coords(topo,vars_n,par,'mesh')					 						! <- module_imexport.f90
		call write_cell_vals(topo,vars_n,par,'vals')											! <- module_imexport.f90
	 end if

	! Memory cleanup
	call dealloc_vars_arrays(vars_n)
	call dealloc_vars_arrays(vars_np1)


end program CCL1D