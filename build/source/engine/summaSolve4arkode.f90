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
 USE type4ida ! reusing IDA data type due to overlap with ARKODE (prime variables not used)
 
! ! access the global print flag
! USE globalData,only: globalPrintFlag
 
! ! access missing values
! USE globalData,only: integerMissing ! missing integer
! USE globalData,only: realMissing    ! missing real number
 
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
 
 ! provide access to indices that define elements of the data structures
 !USE var_lookup,only:iLookPROG       ! named variables for structure elements
 !USE var_lookup,only:iLookDIAG       ! named variables for structure elements
 USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure
 !USE var_lookup,only:iLookDERIV     ! named variables for structure elements
 !USE var_lookup,only:iLookFLUX       ! named variables for structure elements
 !USE var_lookup,only:iLookPARAM      ! named variables for structure elements
 !USE var_lookup,only:iLookINDEX      ! named variables for structure elements
 
 ! provide access to the derived types to define the data structures
 USE data_types,only:&
                     var_i,        & ! data vector (i4b)
                     var_d,        & ! data vector (rkind)
                     var_ilength,  & ! data vector with variable length dimension (i4b)
                     var_dlength,  & ! data vector with variable length dimension (rkind)
                     model_options   ! defines the model decisions

 ! look-up values for the choice of groundwater parameterization
 USE mDecisions_module,only: qbaseTopmodel ! TOPMODEL-ish baseflow parameterization

! ! look-up values for the choice of variable in energy equations (BE residual or IDA state variable)
! USE mDecisions_module,only:       &
!   closedForm,                     & ! use temperature with closed form heat capacity
!   enthalpyFormLU,                 & ! use enthalpy with soil temperature-enthalpy lookup tables
!   enthalpyForm                      ! use enthalpy with soil temperature-enthalpy analytical solution
 
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
                      dt_cur,                  & ! intent(in):    current stepsize
                      dt,                      & ! intent(in):    data time step
                      atol,                    & ! intent(in):    absolute tolerance
                      rtol,                    & ! intent(in):    relative tolerance
                      fScale,                  & ! intent(inout): characteristic scale of the function evaluations (mixed units)
                      nSnow,                   & ! intent(in):    number of snow layers
                      nSoil,                   & ! intent(in):    number of soil layers
                      nLayers,                 & ! intent(in):    total number of layers
                      nState,                  & ! intent(in):    total number of state variables
                      ixMatrix,                & ! intent(in):    type of matrix (dense or banded)
                      firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                      computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                      scalarSolution,          & ! intent(in):    flag to indicate the scalar solution
                     ! computMassBalance,       & ! intent(in):    flag to compute mass balance
                     ! computNrgBalance,        & ! intent(in):    flag to compute energy balance
                     ! ! input: state vectors
                      stateVecInit,            & ! intent(in):    initial state vector
                      sMul,                    & ! intent(inout): state vector multiplier (used in the residual calculations)
                      dMat,                    & ! intent(inout): diagonal of the Jacobian matrix (excludes fluxes)
                     ! ! input: data structures
                      model_decisions,         & ! intent(in):    model decisions
                      lookup_data,             & ! intent(in):    lookup data
                      type_data,               & ! intent(in):    type of vegetation and soil
                      attr_data,               & ! intent(in):    spatial attributes
                      mpar_data,               & ! intent(in):    model parameters
                      forc_data,               & ! intent(in):    model forcing data
                      bvar_data,               & ! intent(in):    average model variables for the entire basin
                      prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                     ! ! input-output: data structures
                      indx_data,               & ! intent(inout): index data
                      diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                      flux_data,               & ! intent(inout): model fluxes for a local HRU
                     ! flux_sum,                & ! intent(inout): sum of fluxes model fluxes for a local HRU over a dt_cur
                      deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                     ! mLayerCmpress_sum,       & ! intent(inout): sum of compression of the soil matrix
                     ! ! output
                      ixSaturation,            & ! intent(inout)  index of the lowest saturated layer (NOTE: only computed on the first iteration)
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
   use fsunmatrix_band_mod        ! Fortran interface to banded SUNMatrix
   use fsunlinsol_band_mod        ! Fortran interface to dense SUNLinearSolver
  !use fsunadaptcontroller_soderlind_mod ! Fortran interface to Soderlind controller
   use eval8summa_module,only: eval8summa4arkode          ! RHS function evaluations
   use summaSolve4kinsol_module,only: setInitialCondition ! subroutine for setting initial condition (borrowed from KINSOL routines)
   use tol4ida_module,only:computWeight4ida               ! weight required for tolerances (borrowed from IDA routines)

   !======= Declarations =========
   implicit none

   ! dummy variables
   ! input: model control
   real(rkind),intent(in)          :: dt_cur                 ! current stepsize
   real(qp),intent(in)             :: dt                     ! data time step
   real(qp),intent(inout)          :: atol(:)                ! vector of absolute tolerances
   real(qp),intent(inout)          :: rtol(:)                ! vector of relative tolerances
   real(rkind),intent(inout)       :: fScale(:)              ! characteristic scale of the function evaluations (mixed units)
   integer(i4b),intent(in)         :: nSnow                  ! number of snow layers
   integer(i4b),intent(in)         :: nSoil                  ! number of soil layers
   integer(i4b),intent(in)         :: nLayers                ! total number of layers
   integer(i4b),intent(in)         :: nState                  ! total number of state variables
   integer(i4b),intent(in)         :: ixMatrix               ! form of matrix (dense or banded)
   logical(lgt),intent(in)         :: firstSubStep           ! flag to indicate if we are processing the first sub-step
   logical(lgt),intent(in)         :: computeVegFlux         ! flag to indicate if computing fluxes over vegetation
   logical(lgt),intent(in)         :: scalarSolution         ! flag to denote if implementing the scalar solution
   !logical(lgt),intent(in)         :: computMassBalance      ! flag to compute mass balance
   !logical(lgt),intent(in)         :: computNrgBalance       ! flag to compute energy balance
   !! input: state vectors
   real(rkind),intent(in)          :: stateVecInit(:)        ! model state vector
   real(qp),intent(in)             :: sMul(:)                ! state vector multiplier (used in the residual calculations)
   real(rkind), intent(inout)      :: dMat(:)                ! diagonal of the Jacobian matrix (excludes fluxes)
   !! input: data structures
   type(model_options),intent(in)  :: model_decisions(:)     ! model decisions
   type(zLookup),      intent(in)  :: lookup_data            ! lookup tables
   type(var_i),        intent(in)  :: type_data              ! type of vegetation and soil
   type(var_d),        intent(in)  :: attr_data              ! spatial attributes
   type(var_dlength),  intent(in)  :: mpar_data              ! model parameters
   type(var_d),        intent(in)  :: forc_data              ! model forcing data
   type(var_dlength),  intent(in)  :: bvar_data              ! model variables for the local basin
   type(var_dlength),  intent(in)  :: prog_data              ! prognostic variables for a local HRU
   ! ! input-output: data structures
   type(var_ilength),intent(inout) :: indx_data              ! indices defining model states and layers
   type(var_dlength),intent(inout) :: diag_data              ! diagnostic variables for a local HRU
   type(var_dlength),intent(inout) :: flux_data              ! model fluxes for a local HRU
   !type(var_dlength),intent(inout) :: flux_sum               ! sum of fluxes model fluxes for a local HRU over a dt_cur
   type(var_dlength),intent(inout) :: deriv_data             ! derivatives in model fluxes w.r.t. relevant state variables
   !real(rkind),intent(inout)       :: mLayerCmpress_sum(:)   ! sum of soil compress
   !! output: state vectors
   integer(i4b),intent(inout)      :: ixSaturation           ! index of the lowest saturated layer
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


   ! * local variables *
   ! ODE system variables
   integer(c_long) :: neq     ! # of equations 
   real(c_double)  :: tstart  ! initial time
   real(c_double)  :: tend    ! final time
   real(c_double)  :: dtout   ! output time interval
   real(c_double)  :: tout    ! output time
   real(c_double)  :: tcur(1) ! current time
   integer(c_int)  :: nout    ! number of outputs
   !integer(c_int) :: outstep                  ! output loop counter

   ! SUNDIALS variables
   type(c_ptr)                             :: ctx        ! SUNDIALS context for the simulation
   type(N_Vector), pointer                 :: sunvec_y   ! sundials vector
   type(SUNMatrix), pointer                :: sunmat_A   ! sundials matrix
   integer(c_long)                         :: mu, lu     ! in banded matrix mode in SUNDIALS type
   type(SUNLinearSolver), pointer          :: sunls      ! sundials linear solver
   !type(SUNAdaptController), pointer       :: sunCtrl    ! time step controller
   type(c_ptr)                             :: arkode_mem ! ARKODE memory
   real(c_double), pointer, dimension(neq) :: yvec(:)    ! underlying vector

   ! user data object
   type(data4ida), target  :: eqns_data    ! SUNDIALS user data - reusing IDA data type due to overlap

   ! option variables
   logical(lgt)   :: use_fdJac                ! flag to use finite difference Jacobian, controlled by decision fDerivMeth

   ! return variables
   logical(lgt) :: return_flag ! logical flag for control of return statements
   integer(i4b) :: retval      ! return value for SUNDIALS procedures

   !======= Internals ============

   ! initialize error control
   call initialize_error_control; if (return_flag) return

   ! initialize time variables and # of equations
   call initialize_ODE_system_values; if (return_flag) return

   ! load SUNDIALS user data (eqns_data object)
   call initialize_SUNDIALS_user_data; if (return_flag) return

   ! create the SUNDIALS context
   call initialize_SUNDIALS_context; if (return_flag) return

   ! initialize solution vector
   call initialize_SUNDIALS_solution_vector; if (return_flag) return 

   ! choose Jacobian type
   call initialize_Jacobian_type; if (return_flag) return

   ! create matrix and linear solver objects for SUNDIALS (general or banded storage)
   call initialize_SUNDIALS_matrix_objects; if (return_flag) return

   ! create ARKODE memory variable: attach user data and matrix objects 
   call initialize_ARKODE_memory; if (return_flag) return

   ! initialize tolerance vectors for ARKODE
   call initialize_ARKODE_tolerance_vectors; if (return_flag) return

   ! main solver loop
   !call update_ARKODE_solver_loop
  contains

   subroutine initialize_error_control
    ! *** initialize error control operations ***
    err=0; message="summaSolve4arkode/" ! initialize error code and message
    return_flag=.false.                 ! initialzie return flag
   end subroutine initialize_error_control

   subroutine initialize_ODE_system_values
    ! *** initialize ODE system values *** -- SJT: update these with SUMMA values (using dummy variables)
    tstart = 0.0d0
    tend = 10.0d0
    tcur = tstart
    tout = tstart
    dtout = 1.0d0
    nout = ceiling(tend/dtout)

    ! define # of equations
    neq = nState
 
   end subroutine initialize_ODE_system_values 

   subroutine initialize_SUNDIALS_user_data
    ! *** load SUNDIALS user data ***

    ! fill eqns_data which will be required later to call eval8summa4arkode
    eqns_data%dt_cur         = dt_cur
    eqns_data%dt             = dt
    eqns_data%nSnow          = nSnow
    eqns_data%nSoil          = nSoil
    eqns_data%nLayers        = nLayers
    eqns_data%nState         = nState
    eqns_data%ixMatrix       = ixMatrix
    eqns_data%firstSubStep   = firstSubStep
    eqns_data%firstFluxCall  = .false. ! already called for initial data window -- SJT: may need to be reset in ARKODE solver loop
    eqns_data%firstSplitOper = .false. ! already called for initial data window -- SJT: may need to be reset in ARKODE solver loop
    eqns_data%computeVegFlux = computeVegFlux
    eqns_data%scalarSolution = scalarSolution
    eqns_data%deriv_data     = deriv_data
    eqns_data%lookup_data    = lookup_data
    eqns_data%type_data      = type_data
    eqns_data%attr_data      = attr_data
    eqns_data%mpar_data      = mpar_data
    eqns_data%forc_data      = forc_data
    eqns_data%bvar_data      = bvar_data
    eqns_data%prog_data      = prog_data
    eqns_data%indx_data      = indx_data
    eqns_data%diag_data      = diag_data
    eqns_data%flux_data      = flux_data
    eqns_data%ixSaturation   = ixSaturation

    ! allocate space and fill
    eqns_data%fScale          = fScale          ! allocate on assignment
    eqns_data%model_decisions = model_decisions ! allocate on assignment
    eqns_data%atol = atol ! allocate on assignment
    eqns_data%rtol = rtol ! allocate on assignment
    eqns_data%sMul = sMul ! allocate on assignment
    eqns_data%dMat = dMat ! allocate on assignment

    ! allocate space for other variables -- SJT: commented out lines not required for eval8summa4arkode
    if (model_decisions(iLookDECISIONS%groundwatr)%iDecision==qbaseTopmodel) then
      allocate(eqns_data%dBaseflow_dMatric(nSoil,nSoil),stat=err)
    else
      allocate(eqns_data%dBaseflow_dMatric(0,0),stat=err)
    end if
    !allocate( eqns_data%mLayerTempPrev(nLayers) )
    !allocate( eqns_data%mLayerMatricHeadPrev(nSoil) )
    !allocate( eqns_data%mLayerTempTrial(nLayers) )
    !allocate( eqns_data%mLayerMatricHeadTrial(nSoil) )
    !allocate( eqns_data%mLayerTempPrime(nLayers) )
    !allocate( eqns_data%mLayerMatricHeadPrime(nSoil) )
    !allocate( eqns_data%mLayerVolFracWatPrime(nLayers) )
    !allocate( mLayerMatricHeadPrimePrev(nSoil) )
    !allocate( dCompress_dPsiPrev(nSoil) )
    allocate( eqns_data%fluxVec(nState) )
    allocate( eqns_data%resVec(nState) )
    allocate( eqns_data%resSink(nState) )
    !allocate( resVecPrev(nState) )

   end subroutine initialize_SUNDIALS_user_data

   subroutine initialize_SUNDIALS_context
    ! *** initialize SUNDIALS context variable ***
    retval = FSUNContext_Create(SUN_COMM_NULL, ctx)
    if (retval /= 0) then; err=20; message=trim(message)//'error in FSUNContext_Create'; return_flag=.true.; return; end if
   end subroutine initialize_SUNDIALS_context

   subroutine initialize_SUNDIALS_solution_vector 
    ! *** initialize soultion vector for SUNDIALS ***
    ! create SUNDIALS N_Vector
    sunvec_y => FN_VNew_Serial(neq, ctx)
    if (.not. associated(sunvec_y)) then; err=20; message=trim(message)//'sunvec = NULL'; return_flag=.true.; return; end if
    yvec => FN_VGetArrayPointer(sunvec_y)

    ! initialize solution vector
    call setInitialCondition(neq, stateVecInit, sunvec_y)
   end subroutine initialize_SUNDIALS_solution_vector 

   subroutine initialize_Jacobian_type
    ! *** initialize Jacobian type ***
    ! choose Jacobian type
    select case(model_decisions(iLookDECISIONS%fDerivMeth)%iDecision)
      case(numerical);  use_fdJac =.true.
      case(analytical); use_fdJac =.false.
      case default
       err=20; message=trim(message)//'expect choice numericl or analytic to calculate derivatives for Jacobian'
       return_flag=.true.; return
    end select
   end subroutine initialize_Jacobian_type

   subroutine initialize_SUNDIALS_matrix_objects
    ! *** initialize matrix and linear solver SUNDIALS objects ***
    ! define the form of the matrix
    select case(ixMatrix)
      case(ixBandMatrix)
        mu = ku; lu = kl;
        ! Create banded SUNMatrix for use in linear solves
        sunmat_A => FSUNBandMatrix(neq, mu, lu, ctx)
        if (.not. associated(sunmat_A)) then; err=20; message=trim(message)//'sunmat = NULL'; return_flag=.true.; return; end if

        ! Create banded SUNLinearSolver object
        sunls => FSUNLinSol_Band(sunvec_y, sunmat_A, ctx)
        if (.not. associated(sunls)) then; err=20; message=trim(message)//'sunls = NULL'; return_flag=.true.; return; end if

      case(ixFullMatrix)
        ! Create dense SUNMatrix for use in linear solves
        sunmat_A => FSUNDenseMatrix(neq, neq, ctx)
        if (.not. associated(sunmat_A)) then; err=20; message=trim(message)//'sunmat = NULL'; return_flag=.true.; return; end if

        ! Create dense SUNLinearSolver object
        sunls => FSUNLinSol_Dense(sunvec_y, sunmat_A, ctx)
        if (.not. associated(sunls)) then; err=20; message=trim(message)//'sunls = NULL'; return_flag=.true.; return; end if

        ! check
      case default; err=20; message=trim(message)//'error in type of matrix'; return_flag=.true.; return
    end select
   end subroutine initialize_SUNDIALS_matrix_objects

   subroutine initialize_ARKODE_memory
    ! *** initialize ARKODE memory variable ***
    ! create ARKStep memory - args: (explicit RHS, implicit RHS, start time, sunvec_y, SUNDIALS context)
    arkode_mem = FARKStepCreate(c_null_funptr, c_funloc(eval8summa4arkode), tstart, sunvec_y, ctx)
    if (.not. c_associated(arkode_mem)) then; err=20; message=trim(message)//'arkode_mem = NULL'; return_flag=.true.; return; end if

    ! Attach user data to memory
    retval = FARKodeSetUserData(arkode_mem, c_loc(eqns_data))
    if (retval /= 0) then; err=20; message=trim(message)//'error in FARKodeSetUserData'; return_flag=.true.; return; end if

    ! Attach the matrix and linear solver
    ! For the nonlinear solver, ARKODE uses a Newton SUNNonlinearSolver-- it is not necessary to create and attach it ** SJT: verify this **
    retval = FARKodeSetLinearSolver(arkode_mem, sunls, sunmat_A)
    if (retval /= 0) then; err=20; message=trim(message)//'error in FARKodeSetLinearSolver'; return_flag=.true.; return; end if
   end subroutine initialize_ARKODE_memory

   subroutine initialize_ARKODE_tolerance_vectors
    ! *** initialize ARKODE tolerance vectors ***
    ! set relative and absolute tolerance vectors (using components from eqns_data)
    ! note: reusing the tolerance formula from IDA routines
    retval = FARKodeWFtolerances(arkode_mem, c_funloc(computWeight4ida))
    if (retval /= 0) then; err=20; message=trim(message)//'error in FARKodeWFtolerances'; return_flag=.true.; return; end if
   end subroutine initialize_ARKODE_tolerance_vectors

   subroutine update_ARKODE_solver_loop
    ! *** main ARKODE solver loop ***
   end subroutine update_ARKODE_solver_loop

 end subroutine summaSolve4arkode

end module summaSolve4arkode_module

