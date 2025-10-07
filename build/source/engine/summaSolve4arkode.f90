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
 
 ! access missing values
 USE globalData,only: integerMissing ! missing integer
! USE globalData,only: realMissing    ! missing real number

 ! access numerical parameters
 USE globalData,only: verySmaller    ! a smaller number used as an additive constant to check if substantial difference among real numbers
 
 ! access matrix information
 USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
 USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix
 USE globalData,only: ku             ! number of super-diagonal bands
 USE globalData,only: kl             ! number of sub-diagonal bands
 
 !! global metadata
 USE globalData,only:flux_meta       ! metadata on the model fluxes
 
 ! constants
 USE multiconst,only: Tfreeze        ! temperature at freezing              (K)
 
 ! provide access to indices that define elements of the data structures
 USE var_lookup,only:iLookPROG       ! named variables for structure elements
 USE var_lookup,only:iLookDIAG       ! named variables for structure elements
 USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure
 !USE var_lookup,only:iLookDERIV     ! named variables for structure elements
 !USE var_lookup,only:iLookFLUX       ! named variables for structure elements
 !USE var_lookup,only:iLookPARAM      ! named variables for structure elements
 USE var_lookup,only:iLookINDEX      ! named variables for structure elements
 
 ! provide access to the derived types to define the data structures
 USE data_types,only:&
                     var_i,        & ! data vector (i4b)
                     var_d,        & ! data vector (rkind)
                     var_ilength,  & ! data vector with variable length dimension (i4b)
                     var_dlength,  & ! data vector with variable length dimension (rkind)
                     model_options   ! defines the model decisions

 ! look-up values for the choice of groundwater parameterization
 USE mDecisions_module,only: qbaseTopmodel ! TOPMODEL-ish baseflow parameterization

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
 private
 public  :: summaSolve4arkode

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
                      computMassBalance,       & ! intent(in):    flag to compute mass balance
                      computNrgBalance,        & ! intent(in):    flag to compute energy balance
                      ! input: state vectors
                      stateVecInit,            & ! intent(in):    initial state vector
                      sMul,                    & ! intent(inout): state vector multiplier (used in the residual calculations)
                      dMat,                    & ! intent(inout): diagonal of the Jacobian matrix (excludes fluxes)
                      ! input: data structures
                      model_decisions,         & ! intent(in):    model decisions
                      lookup_data,             & ! intent(in):    lookup data
                      type_data,               & ! intent(in):    type of vegetation and soil
                      attr_data,               & ! intent(in):    spatial attributes
                      mpar_data,               & ! intent(in):    model parameters
                      forc_data,               & ! intent(in):    model forcing data
                      bvar_data,               & ! intent(in):    average model variables for the entire basin
                      prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                      ! input-output: data structures
                      indx_data,               & ! intent(inout): index data
                      diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                      flux_data,               & ! intent(inout): model fluxes for a local HRU
                      flux_sum,                & ! intent(inout): sum of fluxes model fluxes for a local HRU over a dt_cur
                      deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                      mLayerCmpress_sum,       & ! intent(inout): sum of compression of the soil matrix
                      ! output
                      ixSaturation,            & ! intent(inout)  index of the lowest saturated layer (NOTE: only computed on the first iteration)
                      arkodeSucceeds,          & ! intent(out):   flag to indicate if IDA successfully solved the problem in current data step
                      tooMuchMelt,             & ! intent(inout): lag to denote that there was too much melt
                      nSteps,                  & ! intent(out):   number of time steps taken in solver
                      stateVec,                & ! intent(out):   model state vector
                      balance,                 & ! intent(inout): balance per state
                      err,message)               ! intent(out):   error control


   !======= Inclusions ===========
   use, intrinsic :: iso_c_binding
   use fsundials_core_mod                ! Fortran interface to SUNContext
   use farkode_mod                       ! Fortran interface to the ARKODE
   use farkode_arkstep_mod               ! Fortran interface to the ARKStep time-stepper module
   use fnvector_serial_mod               ! Fortran interface to serial N_Vector
   use fsunmatrix_dense_mod              ! Fortran interface to dense SUNMatrix
   use fsunlinsol_dense_mod              ! Fortran interface to dense SUNLinearSolver
   use fsunmatrix_band_mod               ! Fortran interface to banded SUNMatrix
   use fsunlinsol_band_mod               ! Fortran interface to dense SUNLinearSolver
   use fsunadaptcontroller_soderlind_mod ! Fortran interface to Soderlind controller
   use allocspace_module,only:allocLocal                  ! allocate local data structures
   use eval8summa_module,only: eval8summa4arkode          ! RHS function evaluations
   use summaSolve4kinsol_module,only: setInitialCondition ! subroutine for setting initial condition (borrowed from KINSOL routines)
   use tol4ida_module,only:computWeight4ida               ! weight required for tolerances (borrowed from IDA routines)
   use getVectorz_module,only:checkFeas                   ! check feasibility of state vector

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
   integer(i4b),intent(in)         :: nState                 ! total number of state variables
   integer(i4b),intent(in)         :: ixMatrix               ! form of matrix (dense or banded)
   logical(lgt),intent(in)         :: firstSubStep           ! flag to indicate if we are processing the first sub-step
   logical(lgt),intent(in)         :: computeVegFlux         ! flag to indicate if computing fluxes over vegetation
   logical(lgt),intent(in)         :: scalarSolution         ! flag to denote if implementing the scalar solution
   logical(lgt),intent(in)         :: computMassBalance      ! flag to compute mass balance
   logical(lgt),intent(in)         :: computNrgBalance       ! flag to compute energy balance
   ! input: state vectors
   real(rkind),intent(in)          :: stateVecInit(:)        ! model state vector
   real(qp),intent(in)             :: sMul(:)                ! state vector multiplier (used in the residual calculations)
   real(rkind), intent(inout)      :: dMat(:)                ! diagonal of the Jacobian matrix (excludes fluxes)
   ! input: data structures
   type(model_options),intent(in)  :: model_decisions(:)     ! model decisions
   type(zLookup),      intent(in)  :: lookup_data            ! lookup tables
   type(var_i),        intent(in)  :: type_data              ! type of vegetation and soil
   type(var_d),        intent(in)  :: attr_data              ! spatial attributes
   type(var_dlength),  intent(in)  :: mpar_data              ! model parameters
   type(var_d),        intent(in)  :: forc_data              ! model forcing data
   type(var_dlength),  intent(in)  :: bvar_data              ! model variables for the local basin
   type(var_dlength),  intent(in)  :: prog_data              ! prognostic variables for a local HRU
   ! input-output: data structures
   type(var_ilength),intent(inout) :: indx_data              ! indices defining model states and layers
   type(var_dlength),intent(inout) :: diag_data              ! diagnostic variables for a local HRU
   type(var_dlength),intent(inout) :: flux_data              ! model fluxes for a local HRU
   type(var_dlength),intent(inout) :: flux_sum               ! sum of fluxes model fluxes for a local HRU over a dt_cur
   type(var_dlength),intent(inout) :: deriv_data             ! derivatives in model fluxes w.r.t. relevant state variables
   real(rkind),intent(inout)       :: mLayerCmpress_sum(:)   ! sum of soil compress
   ! output: state vectors
   integer(i4b),intent(inout)      :: ixSaturation           ! index of the lowest saturated layer
   integer(i4b),intent(out)        :: nSteps                 ! number of time steps taken in solver
   real(rkind),intent(inout)       :: stateVec(:)            ! model state vector (y)
   !real(rkind),intent(inout)       :: stateVecPrime(:)       ! model state vector (y')
   logical(lgt),intent(out)        :: arkodeSucceeds         ! flag to indicate if ARKODE is successful
   logical(lgt),intent(inout)      :: tooMuchMelt            ! flag to denote that there was too much melt
   !! output: residual terms and balances
   real(rkind),intent(inout)       :: balance(:)             ! balance per state
   ! output: error control
   integer(i4b),intent(out)        :: err                    ! error code
   character(*),intent(out)        :: message                ! error message


   ! * local variables *
   ! ODE system variables
   integer(c_long) :: neq              ! # of equations 
   real(c_double)  :: tstart           ! initial time
   real(c_double)  :: tend             ! final time
   real(c_double)  :: tret(1),tretPrev ! current and previous times in data window
   real(c_double)  :: dt_last(1)       ! last time step
   real(rkind)     :: dt_diff          ! difference from previous timeste
   real(rkind)     :: dt_mult          ! multiplier for time step average values

   ! SUNDIALS variables
   type(c_ptr)                             :: ctx        ! SUNDIALS context for the simulation
   type(N_Vector), pointer                 :: sunvec_y   ! sundials vector
   type(SUNMatrix), pointer                :: sunmat_A   ! sundials matrix
   integer(c_long)                         :: mu, lu     ! in banded matrix mode in SUNDIALS type
   type(SUNLinearSolver), pointer          :: sunls      ! sundials linear solver
   type(SUNAdaptController), pointer       :: sunCtrl    ! time step controller
   type(c_ptr)                             :: arkode_mem ! ARKODE memory
   real(c_double), pointer, dimension(neq) :: yvec(:)    ! underlying vector

   ! ARKODE statistics
   integer(c_long) :: nStepsSun(1)
   integer(c_long) :: nREvals(1)
   integer(c_long) :: nLinSetups(1)
   integer(c_long) :: netFails(1)
   integer(c_int)  :: qLast(1)
   integer(c_int)  :: qCur(1)
   real(c_double)  :: hInitUsed(1)
   real(c_double)  :: hLast(1)
   real(c_double)  :: hCur(1)
   real(c_double)  :: tCur(1)

   ! user data object
   type(data4ida), target  :: eqns_data ! SUNDIALS user data - reusing IDA data type due to overlap

   ! arrays for previous internal ARKODE steps
   real(rkind),allocatable :: resVecPrev(:)         ! previous value for residuals
   type(var_dlength)       :: flux_prev             ! previous model fluxes for a local HRU
   real(rkind),allocatable :: mLayerCompressPrev(:) ! previous soil compressibility value
   !real(rkind),allocatable :: dCompress_dPsiPrev(:) ! previous derivative value soil compression

   ! option variables
   logical(lgt)   :: use_fdJac ! flag to use finite difference Jacobian, controlled by decision fDerivMeth

   ! logical flags
   logical(lgt) :: tinystep    ! if step goes below small size
   logical(lgt) :: feasible    ! feasibility flag

   ! return variables
   logical(lgt)    :: return_flag    ! logical flag for control of return statements
   integer(c_int)  :: retval,retvalr ! return values for SUNDIALS procedures

   ! indices
   integer(i4b) :: i    ! loop index
   integer(i4b) :: iVar ! loop index

   ! error messages
   character(256) :: cmessage ! error message

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

!   ! initialize root finding problem --- to be implemented (see draft routine below)
!   call initialize_ARKODE_root_finding; if (return_flag) return

   ! set controller for adaptive time step sizes
   call initialize_time_step_adaptivity_controller; if (return_flag) return

   ! set time integration scheme options
   call initialize_solver_options; if (return_flag) return

   ! main solver loop
   call update_ARKODE_solver_loop; if (return_flag) return
   print *, "summaSolve4arkode A:" ! SJT --- take out ---
   ! finalize
   call finalize_ARKODE_solver; if (return_flag) return
  contains

   subroutine initialize_error_control
    ! *** initialize error control operations ***
    err=0_i4b; message = "summaSolve4arkode/" ! initialize error code and message

    ! validate: must use enthalpy formulation for ARKODE
    associate(&
     ixNrgConserv => model_decisions(iLookDECISIONS%nrgConserv)%iDecision & ! choice of energy formulation
    &)
     if ((ixNrgConserv /= enthalpyFormLU).and.(ixNrgConserv /= enthalpyForm)) then
      cmessage="enthalpy formulation required for ARKODE"
      err=20; message=trim(message)//trim(cmessage); return_flag=.true.; return
     end if
    end associate

    return_flag    = .false.              ! initialzie return flag
    arkodeSucceeds = .true.               ! initialize ARKODE success flag
   end subroutine initialize_error_control

   subroutine initialize_ODE_system_values
    ! *** initialize ODE system values ***
    tstart  = 0._rkind ! same as IDA
    tend    = dt_cur   ! end time for solver loop
    tret(1) = tstart   ! initialize time in data window

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
    eqns_data%firstFluxCall  = .false. ! already called for initial data window
    eqns_data%firstSplitOper = .false. ! already called for initial data window
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

    ! allocate space for the to save previous fluxes
    call allocLocal(flux_meta(:),flux_prev,nSnow,nSoil,err,cmessage)
    if (err/=0) then; err=20; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if

    ! allocate space for other variables -- SJT: commented out lines not required for eval8summa4arkode
    if (model_decisions(iLookDECISIONS%groundwatr)%iDecision==qbaseTopmodel) then
      allocate(eqns_data%dBaseflow_dMatric(nSoil,nSoil),stat=err)
    else
      allocate(eqns_data%dBaseflow_dMatric(0,0),stat=err)
    end if
    !allocate( eqns_data%mLayerTempPrev(nLayers) )     ! may be required for root finding problem
    !allocate( eqns_data%mLayerMatricHeadPrev(nSoil) )
    !allocate( eqns_data%mLayerTempTrial(nLayers) )    ! may be required for root finding problem
    !allocate( eqns_data%mLayerMatricHeadTrial(nSoil) )
    !allocate( eqns_data%mLayerTempPrime(nLayers) )
    !allocate( eqns_data%mLayerMatricHeadPrime(nSoil) )
    !allocate( eqns_data%mLayerVolFracWatPrime(nLayers) )
    !allocate( mLayerMatricHeadPrimePrev(nSoil) )
    !allocate( dCompress_dPsiPrev(nSoil) )
    allocate( mLayerCompressPrev(nSoil) ) ! note: added for soil compressibility sum calculation for ARKODE (without primed variables)
    allocate( eqns_data%fluxVec(nState) )
    allocate( eqns_data%resVec(nState) )
    allocate( eqns_data%resSink(nState) )
    allocate( resVecPrev(nState) )

    ! need the following values for the first substep
    do iVar=1,size(flux_meta)  ! loop through fluxes
      flux_prev%var(iVar)%dat(:)      = 0._rkind
    end do
    eqns_data%scalarCanopyTempPrev    = prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1) ! required for ARKODE?
    !eqns_data%mLayerTempPrev(:)       = prog_data%var(iLookPROG%mLayerTemp)%dat(:)      ! may be required for root finding problem
    eqns_data%scalarCanopyTempTrial   = prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1) ! required for ARKODE?
    !eqns_data%mLayerTempTrial(:)      = prog_data%var(iLookPROG%mLayerTemp)%dat(:)      ! may be required for root finding problem
    !eqns_data%mLayerMatricHeadPrev(:) = prog_data%var(iLookPROG%mLayerMatricHead)%dat(:)
    !mLayerMatricHeadPrimePrev         = 0._rkind
    !dCompress_dPsiPrev(:)             = 0._rkind
    mLayerCompressPrev(:)             = 0._rkind  ! note: added for soil compressibility sum calculation for ARKODE (without primed variables)
    resVecPrev(:)                     = 0._rkind
    balance(:)                        = 0._rkind

   end subroutine initialize_SUNDIALS_user_data

   subroutine initialize_SUNDIALS_context
    ! *** initialize SUNDIALS context variable ***
    retval = FSUNContext_Create(SUN_COMM_NULL, ctx)
    if (retval /= 0_c_int) then; err=20_i4b; message=trim(message)//'error in FSUNContext_Create'; return_flag=.true.; return; end if
   end subroutine initialize_SUNDIALS_context

   subroutine initialize_SUNDIALS_solution_vector 
    ! *** initialize soultion vector for SUNDIALS ***
    ! create SUNDIALS N_Vector
    sunvec_y => FN_VNew_Serial(neq, ctx)
    if (.not. associated(sunvec_y)) then; err=20_i4b; message=trim(message)//'sunvec = NULL'; return_flag=.true.; return; end if
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
       err=20_i4b; message=trim(message)//'expect choice numericl or analytic to calculate derivatives for Jacobian'
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
        if (.not. associated(sunmat_A)) then; err=20_i4b; message=trim(message)//'sunmat = NULL'; return_flag=.true.; return; end if

        ! Create banded SUNLinearSolver object
        sunls => FSUNLinSol_Band(sunvec_y, sunmat_A, ctx)
        if (.not. associated(sunls)) then; err=20_i4b; message=trim(message)//'sunls = NULL'; return_flag=.true.; return; end if

      case(ixFullMatrix)
        ! Create dense SUNMatrix for use in linear solves
        sunmat_A => FSUNDenseMatrix(neq, neq, ctx)
        if (.not. associated(sunmat_A)) then; err=20_i4b; message=trim(message)//'sunmat = NULL'; return_flag=.true.; return; end if

        ! Create dense SUNLinearSolver object
        sunls => FSUNLinSol_Dense(sunvec_y, sunmat_A, ctx)
        if (.not. associated(sunls)) then; err=20_i4b; message=trim(message)//'sunls = NULL'; return_flag=.true.; return; end if

        ! check
      case default; err=20_i4b; message=trim(message)//'error in type of matrix'; return_flag=.true.; return
    end select
   end subroutine initialize_SUNDIALS_matrix_objects

   subroutine initialize_ARKODE_memory
    ! *** initialize ARKODE memory variable ***
    ! create ARKStep memory - args: (explicit RHS, implicit RHS, start time, sunvec_y, SUNDIALS context)
    arkode_mem = FARKStepCreate(c_null_funptr, c_funloc(eval8summa4arkode), tstart, sunvec_y, ctx)
    if (.not. c_associated(arkode_mem)) then; err=20_i4b; message=trim(message)//'arkode_mem = NULL'; return_flag=.true.; return; end if

    ! Attach user data to memory
    retval = FARKodeSetUserData(arkode_mem, c_loc(eqns_data))
    if (retval /= 0_c_int) then; err=20_i4b; message=trim(message)//'error in FARKodeSetUserData'; return_flag=.true.; return; end if

    ! Attach the matrix and linear solver
    ! For the nonlinear solver, ARKODE uses a Newton SUNNonlinearSolver-- it is not necessary to create and attach it
    retval = FARKodeSetLinearSolver(arkode_mem, sunls, sunmat_A)
    if (retval /= 0_c_int) then; err=20_i4b; message=trim(message)//'error in FARKodeSetLinearSolver'; return_flag=.true.; return; end if
   end subroutine initialize_ARKODE_memory

   subroutine initialize_ARKODE_tolerance_vectors
    ! *** initialize ARKODE tolerance vectors ***
    ! set relative and absolute tolerance vectors (using components from eqns_data)
    ! note: reusing the tolerance formula from IDA routines
    retval = FARKodeWFtolerances(arkode_mem, c_funloc(computWeight4ida))
    if (retval /= 0_c_int) then; err=20_i4b; message=trim(message)//'error in FARKodeWFtolerances'; return_flag=.true.; return; end if
   end subroutine initialize_ARKODE_tolerance_vectors

   !subroutine initialize_ARKODE_root_finding ! to be implemented
   ! ! *** initialize ARKODE root finding problem ***

   ! ! initialize rootfinding problem and allocate space, counting roots
   ! if(detect_events)then
   !   nRoot = 0
   !   if(ixVegNrg/=integerMissing) nRoot = nRoot+1
   !   if(nSnow>0)then
   !     do i = 1,nSnow
   !       if(ixSnowOnlyNrg(i)/=integerMissing) nRoot = nRoot+1
   !     enddo
   !   endif
   !   if(nSoil>0)then
   !     do i = 1,nSoil
   !       if(ixSoilOnlyHyd(i)/=integerMissing) nRoot = nRoot+1
   !       if(ixSoilOnlyNrg(i)/=integerMissing) nRoot = nRoot+1
   !     enddo
   !   endif
   !   allocate( rootsfound(nRoot) )
   !   allocate( rootdir(nRoot) )
   !   rootdir = 0
   !   retval = FIDARootInit(ida_mem, nRoot, c_funloc(layerDisCont4ida))
   !   if (retval /= 0) then; err=20; message=trim(message)//'error in FIDARootInit'; return; endif
   ! else ! will not use, allocate at something
   !   nRoot = 1
   !   allocate( rootsfound(nRoot) )
   !   allocate( rootdir(nRoot) )
   ! endif

   !end subroutine initialize_ARKODE_root_finding 

   subroutine initialize_time_step_adaptivity_controller
    ! *** initialize time step adaptivity controller for ARKODE ***
    sunCtrl => FSUNAdaptController_ImpGus(ctx)
    if (.not. associated(sunCtrl)) then
      err=20_i4b; message=trim(message)//'error: sunCtrl = NULL'; return_flag=.true.; return
    end if
    retval = FARKodeSetAdaptController(arkode_mem, sunCtrl)
    if (retval /= 0_c_int) then; err=20_i4b; message=trim(message)//'error in FARKodeSetAdaptController'; return_flag=.true.; return; end if
   end subroutine initialize_time_step_adaptivity_controller

   subroutine initialize_solver_options
    ! *** set ARKODE solver options ***
    ! note: implicit methods are used due to NULL input for the explicit RHS function in FARKStepCreate call in initialize_ARKODE_memory
    logical(lgt), parameter  :: use_Butcher_tableau = .false. ! flag controlling use of specified Butcher tableau (else use order parameter)
    character(:),allocatable :: method          ! string for ARKODE Butcher tableau
    integer(c_int),parameter :: order = 3_c_int ! order of time integration scheme if not using Butcher tableau (2 <= order <= 5)

    if (use_Butcher_tableau) then ! specify a built-in ARKODE Butcher tableau
      method = "ARKODE_SDIRK_2_1_2"
      retval = FARKStepSetTableName(arkode_mem, method, "ARKODE_ERK_NONE")
      if (retval /= 0_c_int) then; err=20_i4b; message=trim(message)//'error in FARKStepSetTableName'; return_flag=.true.; return; end if
    else                          ! use default method with specified order
      retval = FARKStepSetOrder(arkode_mem, order)
      if (retval /= 0_c_int) then; err=20_i4b; message=trim(message)//'error in FARKStepSetOrder'; return_flag=.true.; return; end if
    end if
   end subroutine initialize_solver_options

   subroutine update_ARKODE_solver_loop
    ! *** main ARKODE solver loop ***

    ! Enforce the solver to stop at end of the time step
    retval = FARKodeSetStopTime(arkode_mem, dt_cur)
    if (retval /= 0_c_int) then; err=20_i4b; message=trim(message)//'error in FARKodeSetStopTime'; return_flag=.true.; return; end if

    ! the following is based on the looping strategy from summaSolve4ida
    tinystep = .false.
    tret(1)  = tstart  ! initial time
    tretPrev = tret(1)
    nSteps   = 0_i4b   ! initialize number of time steps taken in solver

    do while(tret(1) < dt_cur)

      ! SJT: need to set up ARKODE root finding before implementing this block
      ! ! call this at beginning of step to reduce root bouncing (only looking in one direction)
      ! if(detect_events .and. .not.tinystep)then
      !   call find_rootdir(eqns_data, rootdir)
      !   retval = FIDASetRootDirection(ida_mem, rootdir)
      !   if (retval /= 0_c_int) then; err=20_i4b; message=trim(message)//'error in FIDASetRootDirection'; return_flag=.true.; return; end if
      ! endif

      eqns_data%firstFluxCall  = .false. ! already called for initial data window
      eqns_data%firstSplitOper = .false. ! already called for initial data window

      ! call ARKodeEvolve, advance solver just one internal step
      retvalr = FARKodeEvolve(arkode_mem, dt_cur, sunvec_y, tret, ARK_ONE_STEP)
      ! early return if ARKodeEvolve failed
      if( retvalr < 0_c_int )then ! all failures are captured using negative return values
        arkodeSucceeds = .false.
        ! fail from summa problem
        if (eqns_data%err /= 0_i4b) then; message=trim(message)//trim(eqns_data%message); return_flag=.true.; return; end if
        ! fail from solver problem                  
        call getARKodeEvolveMessage
        message=trim(message)//trim(cmessage)
        ! note: the following step to handle exceeding the max # of steps may be implemented in the future
        !if (retvalr == ARK_TOO_MUCH_WORK) err = -20_i4b ! exit and reduce the data window time in varSubStep (not implemented) 
        exit
      end if

      ! loop through non-missing energy state variables in the snow domain to see if need to merge
      tooMuchMelt = .false.
      associate(&
        ixSnowOnlyNrg => eqns_data%indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat & ! intent(in): indices for energy states in the snow subdomain
      &)
        do concurrent (i=1:nSnow,ixSnowOnlyNrg(i) /= integerMissing)
          if (model_decisions(iLookDECISIONS%nrgConserv)%iDecision /= closedForm) then ! using enthalpy as state variable
            if (stateVec(ixSnowOnlyNrg(i)) > 0._rkind) tooMuchMelt = .true. ! need to merge
          else
            if (stateVec(ixSnowOnlyNrg(i)) > Tfreeze)  tooMuchMelt = .true. ! need to merge
          end if
        end do
      end associate
      if (tooMuchMelt) exit

      ! get the last stepsize and difference from previous end time, not necessarily the same
      retval = FARKodeGetLastStep(arkode_mem, dt_last)
      if (retval /= 0_c_int) then; err=20_i4b; message=trim(message)//'error in FARKodeGetLastStep'; return_flag=.true.; return; end if
      dt_diff = tret(1) - tretPrev
      nSteps = nSteps + 1_i4b ! number of time steps taken in solver

      ! possible that vegetation water may go a bit negative because of discontinous canopy wetting derivatives, so check and correct
      associate(&
        ixVegHyd => eqns_data%indx_data%var(iLookINDEX%ixVegHyd)%dat(1) & ! intent(in): index of canopy hydrology state variable (mass)
      &)
        if (ixVegHyd /= integerMissing) then 
          if (stateVec(ixVegHyd) < 0._rkind .and. stateVec(ixVegHyd)>= -verySmaller*1.e3_rkind) stateVec(ixVegHyd) = 0._rkind ! set to zero
        end if
      end associate

      ! check the feasibility of the solution
      feasible=.true.
      call checkFeas(&
                     ! input
                     stateVec,                                             & ! intent(in):    model state vector (mixed units)
                     eqns_data%mpar_data,                                  & ! intent(in):    model parameters
                     eqns_data%prog_data,                                  & ! intent(in):    model prognostic variables for a local HRU
                     eqns_data%indx_data,                                  & ! intent(in):    indices defining model states and layers
                     model_decisions(iLookDECISIONS%nrgConserv)%iDecision.ne.closedForm, & ! intent(in): flag to indicate if we are using enthalpy as state variable
                     ! output: feasibility
                     feasible,                                             & ! intent(inout):   flag to denote the feasibility of the solution
                     ! output: error control
                     err,cmessage)                                           ! intent(out):   error control

      ! early return for non-feasible solutions, right now will just fail if goes infeasible
      if (.not.feasible) then
        arkodeSucceeds = .false.
        message=trim(message)//trim(cmessage)//'non-feasible' ! err=0 is already set, could make this a warning and reduce the data window time in varSubStep
        exit
      end if

      ! sum of fluxes smoothed over the time step, average from instantaneous values
      if (nSteps > 1_i4b) then
        dt_mult = dt_diff/2._rkind ! results in trapezoidal integration rule
      else ! first step no averaging
        dt_mult = dt_diff
      end if

      do iVar=1,size(flux_meta)
        flux_sum%var(iVar)%dat(:) = flux_sum%var(iVar)%dat(:) + ( eqns_data%flux_data%var(iVar)%dat(:) &
                                  & + flux_prev%var(iVar)%dat(:) ) * dt_mult
      end do
      ! note: soil compressibility not computed using primed variables in ARKODE
      ! here is the IDA version with primed variables
      !mLayerCmpress_sum(:) = mLayerCmpress_sum(:) + ( eqns_data%deriv_data%var(iLookDERIV%dCompress_dPsi)%dat(:) * eqns_data%mLayerMatricHeadPrime(:) &
      !                     & + dCompress_dPsiPrev(:) * mLayerMatricHeadPrimePrev(:) ) * dt_mult
      ! here is the ARKODE version without primed variables (for the time integral of mLayerCompress)
      mLayerCmpress_sum(:) = mLayerCmpress_sum(:) + ( eqns_data%diag_data%var(iLookDIAG%mLayerCompress)%dat(:) &
                           & + mLayerCompressPrev(:) ) * dt_mult

      ! ----
      ! * compute energy balance, from residuals *
      !------------------------
      associate(&
        ixCasNrg      => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)      ,& ! intent(in): index of canopy air space energy state variable
        ixVegNrg      => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)      ,& ! intent(in): index of canopy energy state variable
        ixSnowSoilNrg => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat    ,& ! intent(in): indices for energy states in the snow+soil subdomain
        nSnowSoilNrg  => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)  & ! intent(in): number of energy state variables in the snow+soil domain
      &)
        if (computNrgBalance) then
          ! compute energy balance mean, resVec is the instantaneous residual vector from the solver
          if (ixCasNrg/=integerMissing) balance(ixCasNrg) = balance(ixCasNrg) + ( eqns_data%resVec(ixCasNrg) + resVecPrev(ixCasNrg) )*dt_mult/dt
          if (ixVegNrg/=integerMissing) balance(ixVegNrg) = balance(ixVegNrg) + ( eqns_data%resVec(ixVegNrg) + resVecPrev(ixVegNrg) )*dt_mult/dt
          if (nSnowSoilNrg > 0) then
            do concurrent (i=1:nLayers,ixSnowSoilNrg(i)/=integerMissing)
              balance(ixSnowSoilNrg(i)) = balance(ixSnowSoilNrg(i)) + ( eqns_data%resVec(ixSnowSoilNrg(i)) + resVecPrev(ixSnowSoilNrg(i)) )*dt_mult/dt
            end do
          end if
        end if
      end associate

      ! ----
      ! * compute mass balance, from residuals *
      !------------------------   
      associate(&
        ixAqWat       => indx_data%var(iLookINDEX%ixAqWat)%dat(1)       ,&      ! intent(in): index of water storage in the aquifer
        ixVegHyd      => eqns_data%indx_data%var(iLookINDEX%ixVegHyd)%dat(1), & ! intent(in): index of canopy hydrology state variable (mass)
        ixSnowSoilHyd => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat    ,&      ! intent(in): indices for hydrology states in the snow+soil subdomain
        nSnowSoilHyd  => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)  &      ! intent(in): number of hydrology variables in the snow+soil domain
      &)
        if (computMassBalance) then   
    
          ! compute mass balance mean, resVec is the instantaneous residual vector from the solver
          if (ixVegHyd/=integerMissing) balance(ixVegHyd) = balance(ixVegHyd) + ( eqns_data%resVec(ixVegHyd) + resVecPrev(ixVegHyd) )*dt_mult/dt
          if (nSnowSoilHyd>0) then    
            do concurrent (i=1:nLayers,ixSnowSoilHyd(i)/=integerMissing) 
              balance(ixSnowSoilHyd(i)) = balance(ixSnowSoilHyd(i)) + ( eqns_data%resVec(ixSnowSoilHyd(i)) + resVecPrev(ixSnowSoilHyd(i)) )*dt_mult/dt
            end do
          end if
          if (ixAqWat/=integerMissing) balance(ixAqWat) = balance(ixAqWat) + ( eqns_data%resVec(ixAqWat) + resVecPrev(ixAqWat) )*dt_mult/dt
        end if
      end associate


      ! save required quantities for next step
      eqns_data%scalarCanopyTempPrev     = eqns_data%scalarCanopyTempTrial ! required for ARKODE?
      !eqns_data%mLayerTempPrev(:)       = eqns_data%mLayerTempTrial(:)   ! may be required for root finding problem 
      !eqns_data%mLayerMatricHeadPrev(:) = eqns_data%mLayerMatricHeadTrial(:)
      !mLayerMatricHeadPrimePrev(:)      = eqns_data%mLayerMatricHeadPrime(:)
      !dCompress_dPsiPrev(:)             = eqns_data%deriv_data%var(iLookDERIV%dCompress_dPsi)%dat(:)
      mLayerCompressPrev(:)              = eqns_data%diag_data%var(iLookDIAG%mLayerCompress)%dat(:)
      tretPrev                           = tret(1)
      resVecPrev(:)                      = eqns_data%resVec(:)
      flux_prev                          = eqns_data%flux_data

      ! SJT: need to set up ARKODE root finding before implementing this block
      !! Restart for where vegetation and layers cross freezing point
      !if(detect_events)then
      !  if (retvalr .eq. IDA_ROOT_RETURN) then ! IDASolve succeeded and found one or more roots at tret(1)
      !    ! rootsfound[i]= +1 indicates that gi is increasing, -1 g[i] decreasing, 0 no root
      !    !retval = FIDAGetRootInfo(ida_mem, rootsfound)
      !    !if (retval < 0) then; err=20; message=trim(message)//'error in FIDAGetRootInfo'; return; endif
      !    !print '(a,f15.7,2x,17(i2,2x))', "time, rootsfound[] = ", tret(1), rootsfound
      !    ! Reininitialize solver for running after discontinuity and restart
      !    retval = FIDAReInit(ida_mem, tret(1), sunvec_y, sunvec_yp)
      !    if (retval /= 0) then; err=20; message=trim(message)//'error in FIDAReInit'; return; endif
      !    if(dt_last(1) < 0.1_rkind)then ! don't keep calling if step is small (more accurate with this tiny but getting hung up)
      !      retval = FIDARootInit(ida_mem, 0, c_funloc(layerDisCont4ida))
      !      tinystep = .true.
      !    else
      !      retval = FIDARootInit(ida_mem, nRoot, c_funloc(layerDisCont4ida))
      !      tinystep = .false.
      !    endif
      !    if (retval /= 0) then; err=20; message=trim(message)//'error in FIDARootInit'; return; endif
      !  endif
      !endif

    end do

   end subroutine update_ARKODE_solver_loop

   subroutine finalize_ARKODE_solver
    ! *** Finalize operations for ARKODE solver ***

    ! interface ARKODE user data object to summaSolve4arkode variables 
    if (arkodeSucceeds) then
      ! copy to output data
      diag_data     = eqns_data%diag_data
      flux_data     = eqns_data%flux_data
      deriv_data    = eqns_data%deriv_data
      ixSaturation  = eqns_data%ixSaturation
      indx_data%var(iLookINDEX%numberFluxCalc)%dat(1) = eqns_data%indx_data%var(iLookINDEX%numberFluxCalc)%dat(1) ! only number of flux calculations changes in indx_data
      err           = eqns_data%err
      message       = eqns_data%message
    end if

    ! free memory
    deallocate( eqns_data%model_decisions)
    deallocate( eqns_data%sMul )
    deallocate( eqns_data%dMat )
    deallocate( eqns_data%dBaseflow_dMatric )
    !deallocate( eqns_data%mLayerTempPrev )          ! may be required for root finding problem
    !deallocate( eqns_data%mLayerMatricHeadPrev )
    !deallocate( eqns_data%mLayerTempTrial )         ! may be required for root finding problem
    !deallocate( eqns_data%mLayerMatricHeadTrial )
    !deallocate( eqns_data%mLayerTempPrime )
    !deallocate( eqns_data%mLayerMatricHeadPrime )
    !deallocate( eqns_data%mLayerVolFracWatPrime )
    !deallocate( mLayerMatricHeadPrimePrev )
    !deallocate( dCompress_dPsiPrev )
    deallocate( eqns_data%resVec )
    deallocate( eqns_data%resSink )
    !deallocate( rootsfound ) ! may be required for root finding problem
    !deallocate( rootdir )    ! may be required for root finding problem

    ! Get Stats from ARKODE
    retval = FARKodeGetStepStats(arkode_mem, nStepsSun, hInitUsed, hLast, hCur, tCur)
    if (retval /= 0_c_int) then; err=20_i4b; message=trim(message)//'error in FARKodeGetStepStats'; return_flag=.true.; return; end if
    retval = FARKodeGetNumRhsEvals(arkode_mem, 1_c_int, nREvals) ! args = (arkode_mem,partition,nREvals) -- partition=1 for implicit RHS
    if (retval /= 0_c_int) then; err=20_i4b; message=trim(message)//'error in FARKodeGetNumRhsEvals'; return_flag=.true.; return; end if
    retval = FARKodeGetNumErrTestFails(arkode_mem, netFails)
    if (retval /= 0_c_int) then; err=20_i4b; message=trim(message)//'error in FARKodeGetNumErrTestFails'; return_flag=.true.; return; end if
    retval = FARKodeGetNumLinSolvSetups(arkode_mem, nLinSetups)
    if (retval /= 0_c_int) then; err=20_i4b; message=trim(message)//'error in FARKodeGetNumLinSolvSetups'; return_flag=.true.; return; end if

    diag_data%var(iLookDIAG%numSteps)%dat(1) = nStepsSun(1)
    diag_data%var(iLookDIAG%numResEvals)%dat(1) = nREvals(1)
    diag_data%var(iLookDIAG%numLinSolvSetups)%dat(1) = nLinSetups(1)
    diag_data%var(iLookDIAG%numErrTestFails)%dat(1) = netFails(1)
    !diag_data%var(iLookDIAG%kLast)%dat(1) = qLast(1) ! IDA only -- for variable order
    !diag_data%var(iLookDIAG%kCur)%dat(1) = qCur(1)   ! IDA only -- for variable order
    diag_data%var(iLookDIAG%hInitUsed)%dat(1) = hInitUsed(1)
    diag_data%var(iLookDIAG%hLast)%dat(1) = hLast(1)
    diag_data%var(iLookDIAG%hCur)%dat(1) = hCur(1)
    diag_data%var(iLookDIAG%tCur)%dat(1) = tCur(1)

    call FARKodeFree(arkode_mem)
    retval = FSUNLinSolFree(sunls)
    if (retval /= 0_c_int) then; err=20_i4b; message=trim(message)//'unable to free the linear solver'; return_flag=.true.; return; end if
    call FSUNMatDestroy(sunmat_A)
    call FN_VDestroy(sunvec_y)
    !call FN_VDestroy(sunvec_yp) ! IDA only
    retval = FSUNContext_Free(ctx)
    if (retval /= 0_c_int) then; err=20_i4b; message=trim(message)//'unable to free the SUNDIALS context'; return_flag=.true.; return; end if

   end subroutine finalize_ARKODE_solver

   subroutine getARKodeEvolveMessage
    ! *** Get FARKodeEvolve error message from return value *** 
    ! note: https://sundials.readthedocs.io/en/latest/arkode/Usage/User_callable.html#arkode-solver-function 
 
    if (retvalr == ARK_SUCCESS) then
     cmessage = "ARKodeEvolve successful"
    else if (retvalr == ARK_ROOT_RETURN) then
     cmessage = "succeeded and found one or more roots"
    else if (retvalr == ARK_TSTOP_RETURN) then
     cmessage = "succeeded and returned at tstop"
    else if (retvalr == ARK_MEM_NULL) then
     cmessage = "arkode_mem was null"
    else if (retvalr == ARK_NO_MALLOC) then
     cmessage = "arkode_mem was not allocated"
    else if (retvalr == ARK_ILL_INPUT) then
     cmessage = "invalid input"
    else if (retvalr == ARK_TOO_MUCH_WORK) then
     cmessage = "the solver took mxstep internal steps but could not reach tout"
    else if (retvalr == ARK_TOO_MUCH_ACC) then
     cmessage = "the solver could not satisfy the accuracy demanded by the user for some internal step"
    else if (retvalr == ARK_ERR_FAILURE) then
     cmessage = "error test failures occurred either too many times (ark_maxnef) during one internal time step or occurred with |h|=hmin"
    else if (retvalr == ARK_CONV_FAILURE) then
     cmessage = "either convergence test failures occurred too many times (ark_maxncf) during one internal time step or occurred with |h|=hmin"
    else if (retvalr == ARK_LINIT_FAIL) then
     cmessage = "the linear solver’s initialization function failed"
    else if (retvalr == ARK_LSETUP_FAIL) then
     cmessage = "the linear solver’s setup routine failed in an unrecoverable manner"
    else if (retvalr == ARK_LSOLVE_FAIL) then
     cmessage = "the linear solver’s solve routine failed in an unrecoverable manner"
    else if (retvalr == ARK_MASSINIT_FAIL) then
     cmessage = "the mass matrix solver’s initialization function failed"
    else if (retvalr == ARK_MASSSETUP_FAIL) then
     cmessage = "the mass matrix solver’s setup routine failed"
    else if (retvalr == ARK_MASSSOLVE_FAIL) then
     cmessage = "the mass matrix solver’s solve routine failed"
    else if (retvalr == ARK_VECTOROP_ERR) then
     cmessage = "a vector operation error occurred"
    else if (retvalr == ARK_DOMEIG_FAIL) then
     cmessage = "the dominant eigenvalue function failed -- it is either not provided or returns an illegal value"
    else if (retvalr == ARK_MAX_STAGE_LIMIT_FAIL) then
     cmessage = "stepper failed to achieve stable results -- either reduce the step size or increase the stage_max_limit"
    else
     cmessage = "unknown return value from ARKodeEvolve"
    end if
   end subroutine getARKodeEvolveMessage

 end subroutine summaSolve4arkode


end module summaSolve4arkode_module

