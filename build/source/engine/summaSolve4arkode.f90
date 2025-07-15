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
    
module summaSolve4arkode_module

 !======= Inclusions ===========
 USE, intrinsic :: iso_c_binding
 USE nrtype
 !USE type4ida
 
 ! access the global print flag
 USE globalData,only: globalPrintFlag
 
 ! access missing values
 USE globalData,only: integerMissing ! missing integer
 USE globalData,only: realMissing    ! missing real number
 
 ! access matrix information
 USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
 USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix
 USE globalData,only: ku             ! number of super-diagonal bands
 USE globalData,only: kl             ! number of sub-diagonal bands
 
 !! global metadata
 !USE globalData,only:flux_meta       ! metadata on the model fluxes
 !
 !! constants
 !USE multiconst,only: Tfreeze        ! temperature at freezing              (K)
 !
 !! provide access to indices that define elements of the data structures
 !USE var_lookup,only:iLookPROG       ! named variables for structure elements
 !USE var_lookup,only:iLookDIAG       ! named variables for structure elements
 !USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure
 !USE var_lookup,only:iLookDERIV     ! named variables for structure elements
 !USE var_lookup,only:iLookFLUX       ! named variables for structure elements
 !USE var_lookup,only:iLookPARAM      ! named variables for structure elements
 !USE var_lookup,only:iLookINDEX      ! named variables for structure elements
 !
 !! provide access to the derived types to define the data structures
 !USE data_types,only:&
 !                    var_i,        & ! data vector (i4b)
 !                    var_d,        & ! data vector (rkind)
 !                    var_ilength,  & ! data vector with variable length dimension (i4b)
 !                    var_dlength,  & ! data vector with variable length dimension (rkind)
 !                    model_options   ! defines the model decisions

 ! look-up values for the choice of groundwater parameterization
 USE mDecisions_module,only:       &
   qbaseTopmodel,                  & ! TOPMODEL-ish baseflow parameterization
   bigBucket,                      & ! a big bucket (lumped aquifer model)
   noExplicit                        ! no explicit groundwater parameterization

 ! look-up values for the choice of variable in energy equations (BE residual or IDA state variable)
 USE mDecisions_module,only:       &
   closedForm,                     & ! use temperature with closed form heat capacity
   enthalpyFormLU,                 & ! use enthalpy with soil temperature-enthalpy lookup tables
   enthalpyForm                      ! use enthalpy with soil temperature-enthalpy analytical solution
 
 ! look-up values for method used to compute derivative
 USE mDecisions_module,only:       &
   numerical,                      & ! numerical solution
   analytical                        ! analytical solution

 ! privacy
 implicit none
 public::summaSolve4arkode

contains

 ! ************************************************************************************
 ! * public subroutine summaSolve4arkode: solve My' = fE(t,y) + fI(t,y) using ARKODE (y is the state vector, y'=dy/dt)
 ! ************************************************************************************
 subroutine summaSolve4arkode
  ! SJT: the following follows an example from the SUNDIALS Git repository and will need to be adapted for SUMMA following summaSolve4ida
  ! note: https://github.com/LLNL/sundials/blob/main/examples/arkode/F2003_serial/ark_analytic_f2003.f90

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use farkode_mod                ! Fortran interface to the ARKODE
  use farkode_arkstep_mod        ! Fortran interface to the ARKStep time-stepper module
  use fnvector_serial_mod        ! Fortran interface to serial N_Vector
  use fsunmatrix_dense_mod       ! Fortran interface to dense SUNMatrix
  use fsunlinsol_dense_mod       ! Fortran interface to dense SUNLinearSolver
  use fsunadaptcontroller_soderlind_mod ! Fortran interface to Soderlind controller
  !use analytic_mod               ! ODE functions

  !======= Declarations =========
  implicit none

  ! local variables
  type(c_ptr)    :: ctx                      ! SUNDIALS context for the simulation
  real(c_double) :: tstart                   ! initial time
  real(c_double) :: tend                     ! final time
  real(c_double) :: rtol, atol               ! relative and absolute tolerance
  real(c_double) :: dtout                    ! output time interval
  real(c_double) :: tout                     ! output time
  real(c_double) :: tcur(1)                  ! current time
  integer(c_int) :: ierr                     ! error flag from C functions
  integer(c_int) :: nout                     ! number of outputs
  integer(c_int) :: outstep                  ! output loop counter

  ! SJT: fix access issues for the derived types below
  !type(N_Vector), pointer                 :: sunvec_y   ! sundials vector
  !type(SUNMatrix), pointer                :: sunmat_A   ! sundials matrix
  !type(SUNLinearSolver), pointer          :: sunls      ! sundials linear solver
  !type(SUNAdaptController), pointer       :: sunCtrl    ! time step controller
  !type(c_ptr)                             :: arkode_mem ! ARKODE memory
  !real(c_double), pointer, dimension(neq) :: yvec(:)    ! underlying vector
 
 end subroutine summaSolve4arkode

end module summaSolve4arkode_module

