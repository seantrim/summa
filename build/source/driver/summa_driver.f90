! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

program summa_driver
  ! **** Driver program for SUMMA simulations ****

  ! * module access *
  ! data types
  USE nr_type                                                 ! variable types, etc.
  USE summa_type, only: summa1_type_dec                       ! master summa data type
  USE data_types, only: convergence_stats_type                ! convergence stats object for Newton iterations
  ! subroutines and functions: model setup
  USE summa_init, only: summa_initialize                      ! used to allocate/initialize summa data structures
  USE summa_setup, only: summa_paramSetup                     ! used to initialize parameter data structures (e.g. vegetation and soil parameters)
  USE summa_restart, only: summa_readRestart                  ! used to read restart data and reset the model state
  ! subroutines and functions: model simulation
  USE summa_forcing, only: summa_readForcing                  ! used to read forcing data
  USE summa_modelRun, only: summa_runPhysics                  ! used to run the summa physics for one time step
  USE summa_writeOutput, only: summa_writeOutputFiles         ! used to write the summa output files
  ! utility functions
  USE summa_util, only: stop_program                          ! used to stop the summa program (with errors)
  USE summa_util, only: handle_err                            ! used to process errors
  ! global data
  USE globalData, only: numtim                                ! number of model time steps
  USE globalData, only: print_step_freq

  ! OpenWQ coupling
#ifdef OPENWQ_ACTIVE
  USE summa_openwq,only:openwq_init
  USE summa_openwq,only:openwq_run_time_start
  USE summa_openwq,only:openwq_run_space_step
  USE summa_openwq,only:openwq_run_time_end
#endif

  implicit none

  ! * driver variables *
  ! define the master summa data structure
  type(summa1_type_dec), allocatable :: summa1_struc(:)
  ! define parameters for the model simulation
  integer(i4b), parameter            :: n=1                        ! number of instantiations
  ! define timing information
  integer(i4b)                       :: modelTimeStep              ! index of model time step
  ! Newton iteration stats
  type(convergence_stats_type)       :: convergence_stats          ! object for convergence stats for Newton iterations
  ! error control
  integer(i4b)                       :: err=0                      ! error code
  character(len=1024)                :: message=''                 ! error message

  ! Initialize
  call initialize_summa_driver

  ! Update
  call update_summa_driver

  ! Finalize
  call finalize_summa_driver

contains

  subroutine initialize_summa_driver
   ! *** Initial operations for SUMMA driver program ***

   ! allocate space for the master summa structure
   allocate(summa1_struc(n), stat=err)
   if (err/=0) call stop_program(1, 'problem allocating master summa structure')

   ! declare and allocate summa data structures and initialize model state to known values
   call summa_initialize(summa1_struc(n), err, message)
   call handle_err(err, message)

   ! initialize parameter data structures (e.g. vegetation and soil parameters)
   call summa_paramSetup(summa1_struc(n), err, message)
   call handle_err(err, message)

   ! read restart data and reset the model state
   call summa_readRestart(summa1_struc(n), err, message)
   call handle_err(err, message)

#ifdef OPENWQ_ACTIVE
   call openwq_init(err)
   if (err /= 0) call stop_program(1, 'Problem Initializing OpenWQ')
#endif
  end subroutine initialize_summa_driver

  subroutine update_summa_driver
   ! *** Update operations for SUMMA driver program ***

   ! initialize convergence_stats object
   call initialize_convergence_stats

   ! loop through time
   do modelTimeStep=1,numtim
 
     ! read model forcing data
     call summa_readForcing(modelTimeStep, summa1_struc(n), err, message)
     call handle_err(err, message)
 
#ifdef OPENWQ_ACTIVE
     call openwq_run_time_start(summa1_struc(n)) ! Passing state volumes to openWQ
#endif
 
     if (mod(modelTimeStep, print_step_freq) == 0) then
       print *, 'step ---> ', modelTimeStep
     end if
 
     ! run the summa physics for one time step
     call summa_runPhysics(modelTimeStep, summa1_struc(n), convergence_stats, err, message)
     call handle_err(err, message)
 
#ifdef OPENWQ_ACTIVE
     call openwq_run_space_step(summa1_struc(n)) ! Passing fluxes to openWQ
#endif
 
     ! write the model output
     call summa_writeOutputFiles(modelTimeStep, summa1_struc(n), err, message)
     call handle_err(err, message)
 
#ifdef OPENWQ_ACTIVE
     call openwq_run_time_end(summa1_struc(n))
#endif
 
   end do  ! end looping through time

   ! finalize convergence_stats object
   call finalize_convergence_stats

  end subroutine update_summa_driver

  subroutine finalize_summa_driver
   ! *** Final operations for SUMMA driver program ***
   ! successful end
   call stop_program(0, 'finished simulation successfully.')

   ! to prevent exiting before HDF5 has closed
   call sleep(2)
  end subroutine finalize_summa_driver

  subroutine initialize_convergence_stats
   ! *** Initialize convergence stats object ***
   USE globalData,only:gru_struc ! gru-hru mapping structures
   integer(i4b) :: iGRU,iHRU     ! indices for GRUs and HRUs

   ! allocate memory and initialize
   associate(nGRU => summa1_struc(n)%nGRU)
    allocate(convergence_stats%gru(1:nGRU))
    do iGRU = 1,nGRU
     allocate(convergence_stats%gru(iGRU)%hru(1:gru_struc(iGRU)%hruCount))
     do iHRU = 1,gru_struc(iGRU)%hruCount
      convergence_stats%gru(iGRU)%hru(iHRU)%high_level_step_reductions        = 0_i4b
      convergence_stats%gru(iGRU)%hru(iHRU)%low_level_step_reductions         = 0_i4b
      convergence_stats%gru(iGRU)%hru(iHRU)%low_level_step_reductions_coupled = 0_i4b
      convergence_stats%gru(iGRU)%hru(iHRU)%splitting_failures                = 0_i4b
      convergence_stats%gru(iGRU)%hru(iHRU)%splitting_failures_coupled        = 0_i4b
      convergence_stats%gru(iGRU)%hru(iHRU)%splitting_failures_coupled        = 0_i4b
      convergence_stats%gru(iGRU)%hru(iHRU)%classical_steps_coupled           = 0_i4b
      convergence_stats%gru(iGRU)%hru(iHRU)%nested_steps_coupled              = 0_i4b
     end do
    end do
   end associate

  end subroutine initialize_convergence_stats

  subroutine finalize_convergence_stats
    ! *** Finalize convergence stats object ***
    USE globalData,only:gru_struc               ! gru-hru mapping structures
    ! index variables for GRUs and HRUs
    integer(i4b) :: iGRU,iHRU                   ! indices for GRUs and HRUs
    integer(i4b) :: nGRU,nHRU                   ! counts for GRUs and HRUs
    ! parameters
    logical(lgt),parameter :: verbose = .true.      ! output for all HRUs?
    logical(lgt),parameter :: verbose_HRU = .false. ! output extra info per HRU?
    ! variables for computing sums
    integer(i4b) :: nHRU_sum
    integer(i4b) :: high_level_step_reductions_sum
    integer(i4b) :: low_level_step_reductions_sum
    integer(i4b) :: low_level_step_reductions_coupled_sum
    integer(i4b) :: splitting_failures_sum
    integer(i4b) :: splitting_failures_coupled_sum
    integer(i4b) :: classical_steps_coupled_sum
    integer(i4b) :: nested_steps_coupled_sum
      
    ! get number of GRUs (assumes n=1 for instantiations)
    nGRU = summa1_struc(n)%nGRU

    ! print convergence information
    print *, ""
    print *, "Convergence Statistics (All HRUs):"
    if (verbose_HRU) then ! verbose output per HRU
      do iGRU = 1,nGRU
        print *, "GRU=",iGRU
        do iHRU = 1,gru_struc(iGRU)%hruCount
          print *, "HRU=",iHRU
          print *, "coupled step reductions       =", convergence_stats%gru(iGRU)%hru(iHRU)%high_level_step_reductions
          print *, "substep reductions            =", convergence_stats%gru(iGRU)%hru(iHRU)%low_level_step_reductions
          print *, "monolithic substep reductions =", convergence_stats%gru(iGRU)%hru(iHRU)%low_level_step_reductions_coupled
          print *, "splitting failures            =", convergence_stats%gru(iGRU)%hru(iHRU)%splitting_failures
          print *, "monolithic failures           =", convergence_stats%gru(iGRU)%hru(iHRU)%splitting_failures_coupled
          print *, "monolithic classical steps    =", convergence_stats%gru(iGRU)%hru(iHRU)%classical_steps_coupled
          print *, "monolithic nested steps       =", convergence_stats%gru(iGRU)%hru(iHRU)%nested_steps_coupled
        end do
      end do
    else              ! streamlined ouput per HRU
      do iGRU = 1,nGRU
        print *, "GRU=",iGRU
        do iHRU = 1,gru_struc(iGRU)%hruCount
          print *, "HRU=",iHRU
          print *, "monolithic substep reductions =", convergence_stats%gru(iGRU)%hru(iHRU)%low_level_step_reductions_coupled
          print *, "monolithic failures           =", convergence_stats%gru(iGRU)%hru(iHRU)%splitting_failures_coupled
          print *, "monolithic classical steps    =", convergence_stats%gru(iGRU)%hru(iHRU)%classical_steps_coupled
          print *, "monolithic nested steps       =", convergence_stats%gru(iGRU)%hru(iHRU)%nested_steps_coupled
        end do
      end do
    end if
    print *, ""

    ! initialize summation variables
    nHRU_sum                              = 0_i4b
    high_level_step_reductions_sum        = 0_i4b 
    low_level_step_reductions_sum         = 0_i4b 
    low_level_step_reductions_coupled_sum = 0_i4b 
    splitting_failures_sum                = 0_i4b
    splitting_failures_coupled_sum        = 0_i4b
    classical_steps_coupled_sum           = 0_i4b 
    nested_steps_coupled_sum              = 0_i4b
    do iGRU = 1,nGRU
      nHRU = gru_struc(iGRU)%hruCount; nHRU_sum = nHRU_sum + nHRU
      high_level_step_reductions_sum        = high_level_step_reductions_sum        + &
                                            & sum(convergence_stats%gru(iGRU)%hru(1:nHRU)%high_level_step_reductions)
      low_level_step_reductions_sum         = low_level_step_reductions_sum         + &
                                            & sum(convergence_stats%gru(iGRU)%hru(1:nHRU)%low_level_step_reductions)
      low_level_step_reductions_coupled_sum = low_level_step_reductions_coupled_sum + &
                                            & sum(convergence_stats%gru(iGRU)%hru(1:nHRU)%low_level_step_reductions_coupled) 
      splitting_failures_sum                = splitting_failures_sum                + &
                                            & sum(convergence_stats%gru(iGRU)%hru(1:nHRU)%splitting_failures)
      splitting_failures_coupled_sum        = splitting_failures_coupled_sum        + &
                                            & sum(convergence_stats%gru(iGRU)%hru(1:nHRU)%splitting_failures_coupled)
      classical_steps_coupled_sum           = classical_steps_coupled_sum        + &
                                            & sum(convergence_stats%gru(iGRU)%hru(1:nHRU)%classical_steps_coupled)
      nested_steps_coupled_sum              = nested_steps_coupled_sum        + &
                                            & sum(convergence_stats%gru(iGRU)%hru(1:nHRU)%nested_steps_coupled)
    end do

    print *, "Convergence Statistics (Mean Per HRU):"
    print *, "coupled step reductions       =", real(high_level_step_reductions_sum,rkind)/real(nHRU_sum,rkind)
    print *, "substep reductions            =", real(low_level_step_reductions_sum,rkind)/real(nHRU_sum,rkind)
    print *, "monolithic substep reductions =", real(low_level_step_reductions_coupled_sum,rkind)/real(nHRU_sum,rkind)
    print *, "splitting failures            =", real(splitting_failures_sum,rkind)/real(nHRU_sum,rkind)
    print *, "monolithic failures           =", real(splitting_failures_coupled_sum,rkind)/real(nHRU_sum,rkind) 
    print *, "monolithic classical steps    =", real(classical_steps_coupled_sum,rkind)/real(nHRU_sum,rkind)
    print *, "monolithic nested steps       =", real(nested_steps_coupled_sum,rkind)/real(nHRU_sum,rkind)
    print *, ""
  end subroutine finalize_convergence_stats

end program summa_driver
