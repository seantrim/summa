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
 USE type4kinsol ! reusing KINSOL data type due to overlap with ARKODE
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
 ! SJT: the following follows an example from the SUNDIALS Git repository and adapted for SUMMA following summaSolve4ida
 ! note: https://github.com/LLNL/sundials/blob/main/examples/arkode/F2003_serial/ark_analytic_f2003.f90
 subroutine summaSolve4arkode(&
                     ! dt_cur,                  & ! intent(in):    current stepsize
                     ! dt,                      & ! intent(in):    data time step
                     ! atol,                    & ! intent(in):    absolute tolerance
                     ! rtol,                    & ! intent(in):    relative tolerance
                     ! nSnow,                   & ! intent(in):    number of snow layers
                     ! nSoil,                   & ! intent(in):    number of soil layers
                     ! nLayers,                 & ! intent(in):    total number of layers
                      nState,                  & ! intent(in):    total number of state variables
                     ! ixMatrix,                & ! intent(in):    type of matrix (dense or banded)
                     ! firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                     ! computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                     ! scalarSolution,          & ! intent(in):    flag to indicate the scalar solution
                     ! computMassBalance,       & ! intent(in):    flag to compute mass balance
                     ! computNrgBalance,        & ! intent(in):    flag to compute energy balance
                     ! ! input: state vectors
                     ! stateVecInit,            & ! intent(in):    initial state vector
                     ! sMul,                    & ! intent(inout): state vector multiplier (used in the residual calculations)
                     ! dMat,                    & ! intent(inout): diagonal of the Jacobian matrix (excludes fluxes)
                     ! ! input: data structures
                     ! model_decisions,         & ! intent(in):    model decisions
                     ! lookup_data,             & ! intent(in):    lookup data
                     ! type_data,               & ! intent(in):    type of vegetation and soil
                     ! attr_data,               & ! intent(in):    spatial attributes
                     ! mpar_data,               & ! intent(in):    model parameters
                     ! forc_data,               & ! intent(in):    model forcing data
                     ! bvar_data,               & ! intent(in):    average model variables for the entire basin
                     ! prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                     ! ! input-output: data structures
                     ! indx_data,               & ! intent(inout): index data
                     ! diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                     ! flux_data,               & ! intent(inout): model fluxes for a local HRU
                     ! flux_sum,                & ! intent(inout): sum of fluxes model fluxes for a local HRU over a dt_cur
                     ! deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                     ! mLayerCmpress_sum,       & ! intent(inout): sum of compression of the soil matrix
                     ! ! output
                     ! ixSaturation,            & ! intent(inout)  index of the lowest saturated layer (NOTE: only computed on the first iteration)
                     ! idaSucceeds,             & ! intent(out):   flag to indicate if IDA successfully solved the problem in current data step
                     ! tooMuchMelt,             & ! intent(inout): lag to denote that there was too much melt
                     ! nSteps,                  & ! intent(out):   number of time steps taken in solver
                     ! stateVec,                & ! intent(out):   model state vector
                     ! stateVecPrime,           & ! intent(out):   derivative of model state vector
                     ! balance,                 & ! intent(inout): balance per state
                      err,message)               ! intent(out):   error control


  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fsundials_core_mod         ! Fortran interface to SUNContext
  use farkode_mod                ! Fortran interface to the ARKODE
  use farkode_arkstep_mod        ! Fortran interface to the ARKStep time-stepper module
  use fnvector_serial_mod        ! Fortran interface to serial N_Vector
  use fsunmatrix_dense_mod       ! Fortran interface to dense SUNMatrix
  use fsunlinsol_dense_mod       ! Fortran interface to dense SUNLinearSolver
  use fsunadaptcontroller_soderlind_mod ! Fortran interface to Soderlind controller
  use eval8summa_module,only: eval8summa4arkode ! RHS function evaluations
  !use analytic_mod               ! ODE functions

  !======= Declarations =========
  implicit none

  ! dummy variables
  ! input: model control
  !real(rkind),intent(in)          :: dt_cur                 ! current stepsize
  !real(qp),intent(in)             :: dt                     ! data time step
  !real(qp),intent(inout)          :: atol(:)                ! vector of absolute tolerances
  !real(qp),intent(inout)          :: rtol(:)                ! vector of relative tolerances
  !integer(i4b),intent(in)         :: nSnow                  ! number of snow layers
  !integer(i4b),intent(in)         :: nSoil                  ! number of soil layers
  !integer(i4b),intent(in)         :: nLayers                ! total number of layers
  integer(i4b),intent(in)         :: nState                  ! total number of state variables
  !integer(i4b),intent(in)         :: ixMatrix               ! form of matrix (dense or banded)
  !logical(lgt),intent(in)         :: firstSubStep           ! flag to indicate if we are processing the first sub-step
  !logical(lgt),intent(in)         :: computeVegFlux         ! flag to indicate if computing fluxes over vegetation
  !logical(lgt),intent(in)         :: scalarSolution         ! flag to denote if implementing the scalar solution
  !logical(lgt),intent(in)         :: computMassBalance      ! flag to compute mass balance
  !logical(lgt),intent(in)         :: computNrgBalance       ! flag to compute energy balance
  !! input: state vectors
  !real(rkind),intent(in)          :: stateVecInit(:)        ! model state vector
  !real(qp),intent(in)             :: sMul(:)                ! state vector multiplier (used in the residual calculations)
  !real(rkind), intent(inout)      :: dMat(:)                ! diagonal of the Jacobian matrix (excludes fluxes)
  !! input: data structures
  !type(model_options),intent(in)  :: model_decisions(:)     ! model decisions
  !type(zLookup),      intent(in)  :: lookup_data            ! lookup tables
  !type(var_i),        intent(in)  :: type_data              ! type of vegetation and soil
  !type(var_d),        intent(in)  :: attr_data              ! spatial attributes
  !type(var_dlength),  intent(in)  :: mpar_data              ! model parameters
  !type(var_d),        intent(in)  :: forc_data              ! model forcing data
  !type(var_dlength),  intent(in)  :: bvar_data              ! model variables for the local basin
  !type(var_dlength),  intent(in)  :: prog_data              ! prognostic variables for a local HRU
  ! ! input-output: data structures
  !type(var_ilength),intent(inout) :: indx_data              ! indices defining model states and layers
  !type(var_dlength),intent(inout) :: diag_data              ! diagnostic variables for a local HRU
  !type(var_dlength),intent(inout) :: flux_data              ! model fluxes for a local HRU
  !type(var_dlength),intent(inout) :: flux_sum               ! sum of fluxes model fluxes for a local HRU over a dt_cur
  !type(var_dlength),intent(inout) :: deriv_data             ! derivatives in model fluxes w.r.t. relevant state variables
  !real(rkind),intent(inout)       :: mLayerCmpress_sum(:)   ! sum of soil compress
  !! output: state vectors
  !integer(i4b),intent(inout)      :: ixSaturation           ! index of the lowest saturated layer
  !integer(i4b),intent(out)        :: nSteps                 ! number of time steps taken in solver
  !real(rkind),intent(inout)       :: stateVec(:)            ! model state vector (y)
  !real(rkind),intent(inout)       :: stateVecPrime(:)       ! model state vector (y')
  !logical(lgt),intent(out)        :: idaSucceeds            ! flag to indicate if IDA is successful
  !logical(lgt),intent(inout)      :: tooMuchMelt            ! flag to denote that there was too much melt
  !! output: residual terms and balances
  !real(rkind),intent(inout)       :: balance(:)             ! balance per state
  ! output: error control
  integer(i4b),intent(out)        :: err                    ! error code
  character(*),intent(out)        :: message                ! error message


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

  type(N_Vector), pointer                 :: sunvec_y   ! sundials vector
  type(SUNMatrix), pointer                :: sunmat_A   ! sundials matrix
  type(SUNLinearSolver), pointer          :: sunls      ! sundials linear solver
  type(SUNAdaptController), pointer       :: sunCtrl    ! time step controller
  type(c_ptr)                             :: arkode_mem ! ARKODE memory
  integer(c_long)                         :: neq        ! # of equations 
  real(c_double), pointer, dimension(neq) :: yvec(:)    ! underlying vector

  !======= Internals ============

  ! initialize error control
  err=0; message="summaSolve4arkode/"

  ! create the SUNDIALS context
  ierr = FSUNContext_Create(SUN_COMM_NULL, ctx)

  ! initialize ODE -- SJT: update these with SUMMA values (using dummy variables)
  tstart = 0.0d0
  tend = 10.0d0
  tcur = tstart
  tout = tstart
  dtout = 1.0d0
  nout = ceiling(tend/dtout)

  ! define # of equations
  neq = nstate

  ! create SUNDIALS N_Vector
  sunvec_y => FN_VNew_Serial(neq, ctx)
  if (.not. associated(sunvec_y)) then; err=20; message=trim(message)//'sunvec = NULL'; return; end if
  yvec => FN_VGetArrayPointer(sunvec_y)

  ! initialize solution vector
  call FN_VConst(0.0d0, sunvec_y)

  ! SJT: continue here -- RhsFn in SUMMA may be related to eval8summa (transformed from implicit form)
  ! create ARKStep memory
  arkode_mem = FARKStepCreate(c_null_funptr, c_funloc(eval8summa4arkode), tstart, sunvec_y, ctx)
  if (.not. c_associated(arkode_mem)) print *, 'ERROR: arkode_mem = NULL'
 
 end subroutine summaSolve4arkode

end module summaSolve4arkode_module

