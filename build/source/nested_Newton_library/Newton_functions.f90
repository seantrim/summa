module Newton_functions
 use, intrinsic :: iso_fortran_env, only: stdout=>output_unit ! for default output
 ! nested Newton solver modules
 use kind_params,only: i4b,r8b ! kind parameters
 use Richards,only : Richards_obj ! Richards test problem
 ! SUMMA modules (for access to constant data and procedures)
 use nr_type,only: rkind,qp,lgt ! SUMMA's kind parameters (i4b is already used in kind_params module)
 use eval8summa_module, only: eval8summa,imposeConstraints           ! SUMMA's eval8summa and imposeConstraints routines
 use computJacob_module,only: computJacob                            ! SUMMA's computJacob routine 
 use summaSolv4homegrown_module,only: checkConv ! SUMMA's checkConv procedures
 use data_types,only: in_type_computJacob,out_type_computJacob ! objects for SUMMA's computJacob routine
 use data_types,only: in_type_summaSolv4homegrown,&           ! objects for SUMMA's summaSolv4homegrown routine
                     &io_type_summaSolv4homegrown,&
                     &out_type_summaSolv4homegrown 
 use data_types,only: model_options           ! type for SUMMA's model decision structure
 use data_types,only: var_ilength,var_dlength ! derived types for SUMMA data structures
 use data_types,only: var_i,var_d             ! derived types for SUMMA data vectors
 use data_types,only: zLookup                 ! derived type for SUMMA lookup tables
 use var_lookup,only: iLookDECISIONS          ! named variables for elements of the SUMMA decision structure
 use var_lookup,only: iLookINDEX              ! named variables for SUMMA structure elements
 use mDecisions_module,only:qbaseTopmodel     ! SUMMA groundwater parameterization model decision
 use stateFilter_module,only: split_select_type
 implicit none
 private

 public :: f_obj_base,f_obj_inputs,f_obj_type

 ! option parameters
 integer(i4b),parameter,public :: LAPACK_expert = 0_i4b,LAPACK_standard = 1_i4b ! linear_system_solver options
 integer(i4b),parameter,public :: LS_C = 0_i4b, LS_I = 1_i4b, LS_O = 2_i4b ! line_search_option options
 integer(i4b),parameter,public :: silent = 0_i4b,minimal = 1_i4b,production = 2_i4b,verbose = 3_i4b,debug = 4_i4b ! output options
 integer(i4b),parameter,public :: custom = 0_i4b, custom_strict = 1_i4b, strict= 2_i4b,&
                                & custom_predictive = 3_i4b, predictive = 4_i4b ! convergence criterion options
 
 ! ***** Parent Type ***** !
 type :: f_obj_base
   ! ** Default data components used by the Newton solvers ** !
   logical      :: banded            ! flag for banded Jacobians
   logical      :: nested            ! flag for nested algorithm
   logical      :: inner             ! flag to indicate the execution of inner iterations
   logical      :: converged         ! flag to indicate that the obtained solution meets the convergence criterion
   logical      :: constraints       ! flag to indicate that constraints are to be applied between outer/classical iterations
   logical      :: constraints_inner ! flag to indicate that constraints are to be applied between inner iterations
   logical      :: refinement        ! flag to indicate that refinement is to be applied following outer/classical iterations
   logical      :: scaling           ! flag to indicate that user-specified scaling is to be applied for linear systems
   logical      :: f_error           ! flag to indicate an error in the evaluation of f, f1, or f2
   logical      :: f_eval_flag       ! flag to indicate that the total non-linear function vector is to be computed
   logical      :: f1_eval_flag      ! flag to indicate that the non-linear function 1 vector is to be computed
   logical      :: f2_eval_flag      ! flag to indicate that the non-linear function 2 vector is to be computed
   logical      :: J_eval_flag       ! flag to indicate that the total Jacobian is to be computed
   logical      :: J1_eval_flag      ! flag to indicate that Jacobian 1 is to be computed
   logical      :: J2_eval_flag      ! flag to indicate that Jacobian 2 is to be computed
   logical      :: evaluate_B        ! flag to indicate if we are evaluating the RHS vector of the Newton iteration equations
   logical      :: dual              ! flag to indicate if dual method is used (.true. to switch selection of f1 and f2)
   logical      :: LAPACK_error      ! flag to indicate an LAPACK error (details given in solver warnings -- otherwise error will be silent)
   logical      :: dynamic           ! flag to indicate for dynamic selection of classical or nested iterations
   logical      :: dynamic_classical ! flag to indicate if in classical phase of dynamic selection mode
   logical      :: dynamic_revert    ! flag to indicate if initial condition should revert to original vector for nested portion of dynamic mode
   integer(i4b) :: subdiag,superdiag ! # of subdiagonals and superdiagonals for banded Jacobians
   integer(i4b) :: n                 ! vector size
   integer(i4b) :: nrow              ! # of matrix rows (adapts to storage type)
   integer(i4b) :: nrow_banded       ! # of matrix rows for banded storage
   integer(i4b) :: k,l               ! indices for classical/outer and inner iterations
   integer(i4b) :: kmax,lmax         ! max # of classical/outer and inner iterations
   integer(i4b) :: kmax_classical    ! max # of classical iterations for dynamic mode
   integer(i4b) :: lmax_loop         ! variable max # of inner iterations (initially lmax but may change depending on observed solution convergence)
   integer(i4b) :: kcount,lcount     ! total # of classical/outer and inner iterations
   integer(i4b) :: LDA,LDAF,LDX,LDB  ! leading dimensions of A, AF, X, and B LAPACK arrays
   integer(i4b) :: KL,KU             ! # of subdiagonals and superdiagonals for LAPACK
   integer(i4b),allocatable :: IPIV(:),IWORK(:)                              ! LAPACK arrays
   real(r8b),allocatable    :: WORK(:),RA(:),CA(:),B(:,:),X(:,:),AF(:,:)            ! LAPACK arrays
   real(r8b)                :: FERR(1:1),BERR(1:1)        ! forward and backward error estimates (single right-hand side assumed)
   real(r8b),pointer        :: x0(:)                      ! guess vector for vector algorithms -- must be associated with vector allocated in external program 
   real(r8b),allocatable    :: xk(:),xkp1(:)              ! intermediate root estimates for classical iterations
   real(r8b),allocatable    :: xk0(:),xkp1l(:),xkp1lp1(:) ! intermediate root estimates for nested iterations
   real(r8b),allocatable    :: xk_0(:),xk_1(:) ! solutions used for computing convergence order for dynamic mode
   real(r8b),allocatable    :: J(:,:)        ! total Jacobian
   real(r8b),allocatable    :: J1(:,:)       ! Jacobian 1
   real(r8b),allocatable    :: J2(:,:)       ! Jacobian 2
   !real(r8b),allocatable    :: Jdiff(:,:)    ! difference Jacobian
   real(r8b),allocatable    :: f_vec(:)      ! total non-linear function evaluation
   real(r8b),allocatable    :: f1_vec(:)     ! non-linear function evaluation 1
   real(r8b),allocatable    :: f2_vec(:)     ! non-linear function evaluation 2
   real(r8b),allocatable    :: initial_solution(:) ! intial solution vector for line search
   real(r8b),allocatable    :: updated_solution(:) ! updated solution vector for line search
   real(r8b),allocatable    :: p_scaled(:)         ! search direction (scaled) for line search
   real(r8b),allocatable    :: grad_L(:)           ! gradient of objective function L for line search
   real(r8b),allocatable    :: f_temp(:)           ! temporary storage of f1 for 'L' scheme for line search
   real(r8b),allocatable    :: vector(:)           ! array for temporary vector output 
   logical,allocatable      :: accept(:)     ! logical mask for accepting guess vector entries for switch to nested iterations in dynamic mode
   real(r8b)                :: tol,tol_inner ! tolerance for classical/outer and inner iterations
   real(r8b)                :: order_min     ! min convergence order for classical iterations in dynamic mode
   real(r8b)                :: R_work(-1:1)  ! work array for max residual computations 
   real(r8b)                :: R(-1:1)       ! max residual computed for outer/classical iterations j-1, j, and j+1 (estimated)  
   real(r8b)                :: R_inner(-1:1) ! max residual computed for inner iterations j-1, j, and j+1 (estimated) 
   real(r8b),allocatable    :: R_vec(:)      ! residual vector
   integer(i4b)             :: convergence          ! string for convergence control option for outer/classical iterations
   integer(i4b)             :: convergence_inner    ! string for convergence control option for inner iterations
   integer(i4b)             :: linear_system_solver ! option for selecting solver for linear systems
   integer(i4b)             :: line_search_option   ! line search option/scheme
   ! solver output
   integer(i4b) :: output      ! solver output control option
   integer(i4b) :: unit = stdout ! file unit number for solver output (default is standard output)
   logical      :: out_debug   ! output flag for debugging
   logical      :: out_detail  ! output flag for details
   logical      :: out_basic   ! output flag for basic information
   logical      :: out_warning ! output flag for warnings
   logical      :: out_error   ! output flag for errors
  contains
   ! procedures used prior to calling the solver
   procedure :: set_defaults    => f_set_defaults    ! set default options 
   procedure :: allocate_memory => f_allocate_memory ! allocate array data components 
   procedure :: set_tolerance   => f_set_tolerance   ! set tolerances and iteration count maximums
   procedure :: solver_output   => f_solver_output   ! set solver output options
   procedure :: matrix_vector_product                ! compute matrix-vector product using matmul or BLAS
 end type f_obj_base

 type,extends(f_obj_base) :: f_obj_inputs
   ! * SUMMA data *
   !type(model_options),allocatable :: model_decisions(:) ! model decisions
   type(model_options),pointer :: model_decisions(:) => null() ! model decisions

   type(zLookup)    ,pointer :: lookup_data => null() ! lookup tables
   !type(var_dlength) :: flux_init                    ! model fluxes at the start of the time step
   type(var_i)      ,pointer :: type_data => null()   ! type of vegetation and soil
   type(var_d)      ,pointer :: attr_data => null()   ! spatial attributes
   type(var_d)      ,pointer :: forc_data => null()   ! model forcing data
   type(var_dlength),pointer :: mpar_data => null()   ! model parameters
   type(var_dlength),pointer :: bvar_data => null()   ! model variables for the local basin


   type(var_ilength),pointer :: indx_data => null()    ! indices defining model states and layers
   type(var_dlength),pointer :: prog_data => null()    ! prognostic variables for a local HRU
   type(var_dlength),pointer :: diag_data => null()    ! diagnostic variables for a local HRU
   type(var_dlength),pointer :: flux_data => null()    ! temporary flux variables for a local HRU
   type(var_dlength),pointer :: deriv_data => null()   ! derivatives in model fluxes w.r.t. relevant state variables
   !real(rkind),allocatable   :: dBaseflow_dMatric(:,:) ! derivative in baseflow w.r.t. matric head (s-1)
   real(rkind),allocatable   :: dBaseflow_dWat(:,:)    ! derivative in baseflow w.r.t. water content (s-1)
   real(rkind),allocatable   :: dBaseflow_dTk(:,:)     ! derivative in baseflow w.r.t. temperature (s-1)
   real(rkind),pointer       :: dMat(:) => null()      ! diagonal matrix (excludes flux derivatives) 

   ! * summaSolve4homegrown (SS4HG) objects *
   ! classical and outer iterations
   type(in_type_summaSolv4homegrown) ,pointer :: in_SS4HG  => null() ! SS4HG input object: model control variables and previous function evaluation
   type(io_type_summaSolv4homegrown) ,pointer :: io_SS4HG  => null() ! SS4HG io object: model control variables and previous function evaluation
   type(out_type_summaSolv4homegrown),pointer :: out_SS4HG => null() ! SS4HG output object: model control variables and previous function evaluation

   ! additional variables for eval8summa call
   logical(lgt),pointer    :: firstSplitOper => null() ! flag to indicate if we are processing the first flux call in a splitting operation
   real(rkind),pointer     :: fScale(:)      => null() ! characteristic scale of the function evaluations (mixed units)
   real(rkind),pointer     :: xScale(:)      => null() ! characteristic scale of the state vector (mixed units)
   real(qp),pointer        :: sMul(:)        => null() ! NOTE: qp  ! multiplier for state vector for the residual calculations
   logical(lgt),pointer    :: feasible       => null() ! feasibility flag
   real(rkind),pointer     :: fluxVec0(:)    => null() ! flux vector (mixed units)
   real(rkind),pointer     :: fRHS(:)        => null() ! RHS function for ARKODE
   real(rkind),pointer     :: rAdd(:)        => null() ! additional terms in the residual vector
   real(qp),pointer        :: resVec(:)      => null() ! NOTE: qp  ! residual vector 

   ! scaled arrays
   real(rkind),allocatable :: rVecScaled(:) ! scaled residual
   real(rkind),allocatable :: aJacScaled(:,:) ! scaled Jacobian

   ! scalars
   integer(i4b)             :: nLeadDim ! lead dimension of SUMMA LAPACK arrays

   ! variables to handle state type non-linear function decompositions
   integer(i4b)             :: nSubset1,nSubset2
   logical(lgt),allocatable :: stateMask1(:),stateMask2(:)  
   logical :: f1_mass_flag,f1_energy_flag
   logical :: f2_mass_flag,f2_energy_flag

  contains
   ! ** routines that point to external sources ** !
   ! note: - these procedures are not directly called in the solver
   !       - however, these procedures may be called within procedures that are called in the solver

   ! * Interfaces for SUMMA procedures  *
   procedure :: SUMMA_eval8summa
   procedure :: SUMMA_computJacob

 end type f_obj_inputs

 type,extends(f_obj_inputs) :: f_obj_type
   real(r8b) :: L0 ! initial line search objective function value
  contains
   ! *** these procedures take the procedures from f_obj_inputs type as input *** !
   ! vector routines
   procedure, non_overridable :: f_vec_eval  => f_SUMMA_vec  ! solver -- f
   procedure, non_overridable :: f1_vec_eval => f_f1_SUMMA_vec_full   ! solver -- f and f1
   procedure, non_overridable :: f2_vec_eval => f_f2_SUMMA_vec_full   ! solver -- f and f2
   procedure, non_overridable :: f1_vec_only_eval => f1_SUMMA_vec_full   ! solver --- f1 only
   procedure, non_overridable :: f2_vec_only_eval => f2_SUMMA_vec_full   ! solver --- f2 only
   procedure, non_overridable :: f1_f2_vec_eval => f_f1_f2_SUMMA_vec_full   ! solver -- f, f1, and f2
   procedure, non_overridable :: J_eval  => J_SUMMA_vec       ! solver -- J
   procedure, non_overridable :: J1_eval => J1_SUMMA_vec_full ! solver -- J1
   procedure, non_overridable :: J2_eval => J2_SUMMA_vec_full ! solver -- J2
   !procedure :: J1_J2_eval => J_J1_J2_SUMMA_vec_full ! solver J, J1, and J2
   procedure, non_overridable :: apply_constraints  => SUMMA_imposeConstraints
   procedure, non_overridable :: apply_nested_line_search => SUMMA_nested_line_search
   procedure, non_overridable :: line_search_objective => SUMMA_line_search_objective
   procedure, non_overridable :: custom_convergence => SUMMA_check_convergence_flag !SUMMA_checkConv  
   procedure, non_overridable :: custom_scaling     => SUMMA_scaling  
   procedure, non_overridable :: custom_descaling   => SUMMA_descaling  
   !procedure :: get_mass_energy_masks => get_SUMMA_mass_energy_masks
   procedure, non_overridable :: get_f1_f2_flags => get_SUMMA_f1_f2_flags
   procedure, non_overridable :: f_state_SUMMA_vec_full

 end type f_obj_type

contains

!!!!!!!!!! ******************* User defined functions below ******************* !!!!!!!!!!

 ! **** Utilities **** !

 subroutine f_set_defaults(f_obj)
  ! ** set default values for options in f_obj_base class **
  class(f_obj_base),intent(inout) :: f_obj

   f_obj % banded            = .false. ! flag for banded Jacobians
   f_obj % nested            = .false. ! flag for nested algorithm
   f_obj % dynamic           = .false. ! flag for dynamic switching between classical and nested regimes
   f_obj % dynamic_classical = .false. ! flag for dynamic switching between classical and nested regimes (classical iteration portion)
   f_obj % dual              = .false. ! flag to use dual method for selecting f1 and f2
   f_obj % constraints       = .false. ! flag to indicate that constraints are to be applied between outer/classical iterations
   f_obj % constraints_inner = .false. ! flag to indicate that constraints are to be applied between inner iterations
   f_obj % refinement        = .false. ! flag to indicate that refinement is to be applied following outer/classical iterations
   f_obj % scaling           = .false. ! flag to indicate that user-specified scaling is to be applied for linear systems
   f_obj % f_eval_flag       = .true.  ! flag to indicate that the total non-linear function vector is to be computed
   f_obj % f1_eval_flag      = .true.  ! flag to indicate that the non-linear function 1 vector is to be computed
   f_obj % f2_eval_flag      = .true.  ! flag to indicate that the non-linear function 2 vector is to be computed
   f_obj % J_eval_flag       = .true.  ! flag to indicate that the total Jacobian is to be computed
   f_obj % J1_eval_flag      = .true.  ! flag to indicate that Jacobian 1 is to be computed
   f_obj % J2_eval_flag      = .true.  ! flag to indicate that Jacobian 2 is to be computed

   f_obj % kmax        = 100_i4b ! max # of classical/outer iterations
   f_obj % lmax        = 100_i4b ! max # inner iterations
   
   f_obj % tol         = 1.e-8   ! tolerance for classical/outer iterations
   f_obj % tol_inner   = 1.e-8   ! tolerance for inner iterations

   f_obj % linear_system_solver = LAPACK_standard          ! string for control of linear system solver
   f_obj % convergence          = strict                   ! string for convergence criterion method for solver
   f_obj % convergence_inner    = strict                   ! string for convergence criterion method for solver
   f_obj % output               = production               ! string for solver output control option
   
   f_obj % unit                 = stdout                   ! file unit number for solver output (this is also the initial value on declaration of f_obj)

 end subroutine f_set_defaults
 
 subroutine f_allocate_memory(f_obj)
  ! ** allocate array data components for f_obj_base class **
  class(f_obj_base),intent(inout) :: f_obj

  ! allocate solution and function arrays
  associate(n => f_obj % n)
   allocate(f_obj % f_vec(1:n))                    ! total non-linear function vector
   if (f_obj % nested) then
    allocate(f_obj % xk0(1:n),f_obj % xkp1l(1:n),f_obj % xkp1lp1(1:n)) ! intermediate root estimates for nested iterations
    allocate(f_obj % f1_vec(1:n),f_obj % f2_vec(1:n))                  ! non-linear functions vectors 1 and 2 
    if (f_obj % dynamic) then
     !allocate(f_obj % xk1(1:n),f_obj % xk2(1:n)) ! solutions used to compute convergence order
     allocate(f_obj % xk(1:n),f_obj % xkp1(1:n))    ! intermediate root estimates for classical iterations
     allocate(f_obj % accept(1:n)) ! logical mask for acceptance of guess vector entries for switch to nested iterations
    end if
   else
    allocate(f_obj % xk(1:n),f_obj % xkp1(1:n))    ! intermediate root estimates for classical iterations
   end if
   allocate(f_obj % xk_0(1:n),f_obj % xk_1(1:n)) ! solutions used to compute convergence order
   allocate(f_obj % R_vec(1:n)) ! residual vector
   if (f_obj % refinement) then ! arrays for line search
    allocate(f_obj % initial_solution(1:n),f_obj % updated_solution(1:n),f_obj % p_scaled(1:n),f_obj % grad_L(1:n),&
            &f_obj % f_temp(1:n))          
   end if
   allocate(f_obj % vector(1:n)) ! output vector (e.g., for matrix-vector products)
  end associate

  ! * allocate LAPACK arrays *

  ! LAPACK parameters independent of matrix storage type
  f_obj % LDX=f_obj % n; f_obj % LDB=f_obj % n ! leading dimensions for RHS arrays

  ! allocate memory and set LAPACK parameters for choice of matrix storage
  allocate(f_obj % B(1:f_obj % n,1:1))                            ! RHS vector (and solution following solver call)
  if (f_obj % banded) then ! banded storage
   f_obj % KL = f_obj % subdiag; f_obj % KU = f_obj % superdiag
   f_obj % LDA = f_obj % KL + f_obj % KU + 1_i4b; f_obj % LDAF = f_obj % LDA + f_obj % KL
   allocate(f_obj % AF(1:f_obj % LDAF,1:f_obj % n)) ! storing LU factors requires an additional f_obj % subdiag rows
   allocate(f_obj % IPIV(1:f_obj % n))                            ! pivot index vector
   if (f_obj % linear_system_solver .eq. LAPACK_expert) then
    allocate(f_obj % IWORK(1:f_obj % n))                          ! work integer array
    allocate(f_obj % WORK(1:3_i4b*f_obj % n))
    allocate(f_obj % RA(1:f_obj % n),f_obj % CA(1:f_obj % n))     ! row and column scale factors for A
    allocate(f_obj % X(1:f_obj % n,1:1))                          ! solution to original (unscaled) system
   end if
  else ! full matrix storage
   f_obj % LDA = f_obj % n; f_obj % LDAF = f_obj % n
   allocate(f_obj % AF(1:f_obj % n,1:f_obj % n))
   allocate(f_obj % IPIV(1:f_obj % n))                            ! pivot index vector
   if (f_obj % linear_system_solver .eq. LAPACK_expert) then
    allocate(f_obj % IWORK(1:f_obj % n))                          ! work integer array
    allocate(f_obj % WORK(1:4_i4b*f_obj % n))
    allocate(f_obj % RA(1:f_obj % n),f_obj % CA(1:f_obj % n))     ! row and column scale factors for A
    allocate(f_obj % X(1:f_obj % n,1:1))                          ! solution to original (unscaled) system
   end if
  end if

  ! allocate Jacobian arrays (and initialize to zero)
  if (f_obj % banded) then ! banded storage
    !f_obj % nrow_banded = f_obj % subdiag + f_obj % superdiag + 1_i4b
    f_obj % nrow_banded = f_obj % LDAF ! 2*KL + KU + 1 (should work for both standard and expert LAPACK solvers, although the input matrix must be loaded differently for each)
    f_obj % nrow = f_obj % nrow_banded
  else
    f_obj % nrow_banded = f_obj % n
    f_obj % nrow        = f_obj % n
  end if
  
  allocate(f_obj % J(1:f_obj % nrow,1:f_obj % n),source=0._r8b) ! SJT: available for classical and nested iterations (for testing -- not needed for nested)
  if (f_obj % nested) then
   !allocate(f_obj % J1(1:f_obj % nrow,1:f_obj % n),f_obj % J2(1:f_obj % nrow,1:f_obj % n),&
   !        &f_obj % Jdiff(1:f_obj % nrow,1:f_obj % n),source=0._r8b)
   allocate(f_obj % J1(1:f_obj % nrow,1:f_obj % n),f_obj % J2(1:f_obj % nrow,1:f_obj % n),source=0._r8b)
  end if

 end subroutine f_allocate_memory

 subroutine f_solver_output(f_obj,method)
  ! ** set output control for solver **
  class(f_obj_base),intent(inout)  :: f_obj
  integer(i4b),intent(in)          :: method
 
  if (method.eq.silent) then
   f_obj % out_debug   = .false. ! output flag for debugging
   f_obj % out_detail  = .false. ! output flag for details
   f_obj % out_basic   = .false. ! output flag for basic information
   f_obj % out_warning = .false. ! output flag for warnings
   f_obj % out_error   = .true.  ! output flag for errors
  else if (method.eq.minimal) then
   f_obj % out_debug   = .false. ! output flag for debugging
   f_obj % out_detail  = .false. ! output flag for details
   f_obj % out_basic   = .true.  ! output flag for basic information
   f_obj % out_warning = .false. ! output flag for warnings
   f_obj % out_error   = .true.  ! output flag for errors
  else if (method.eq.production) then
   f_obj % out_debug   = .false. ! output flag for debugging
   f_obj % out_detail  = .false. ! output flag for details
   f_obj % out_basic   = .true.  ! output flag for basic information
   f_obj % out_warning = .true.  ! output flag for warnings
   f_obj % out_error   = .true.  ! output flag for errors
  else if (method.eq.verbose) then
   f_obj % out_debug   = .false. ! output flag for debugging
   f_obj % out_detail  = .true.  ! output flag for details
   f_obj % out_basic   = .true.  ! output flag for basic information
   f_obj % out_warning = .true.  ! output flag for warnings
   f_obj % out_error   = .true.  ! output flag for errors
  else if (method.eq.debug) then
   f_obj % out_debug   = .true. ! output flag for debugging
   f_obj % out_detail  = .true. ! output flag for details
   f_obj % out_basic   = .true. ! output flag for basic information
   f_obj % out_warning = .true. ! output flag for warnings
   f_obj % out_error   = .true. ! output flag for errors
  else
   if (f_obj % out_error) then
    write(f_obj % unit,'(a65)') "Error in f_solver_output: method argument not currently supported"
   end if
   stop
  end if
 end subroutine f_solver_output

 ! **** Numerics **** !

 subroutine f_set_tolerance(f_obj,method,tol,kmax) ! ***** note: not currently in use (needs updates, including change to integer method options) *****
  ! ** set tolerance for f_obj_base class **
  class(f_obj_base),intent(inout) :: f_obj
  character(*),intent(in)         :: method
  real(r8b),intent(in)            :: tol  ! relative tolerance for the classical/outer solution 
  integer(i4b),intent(in)         :: kmax ! max # of classical/outer iterations

  ! assign values for classical/outer iterations from inputs
  f_obj % tol  = tol
  f_obj % kmax = kmax

  ! assign values for inner iterations if needed
  if (f_obj % nested) then ! nested Newton iterations
   if (method.eq.'strict') then ! inner and outer iterations have the same tolerance
    f_obj % tol_inner = tol
    f_obj % lmax      = kmax
   else if (method.eq.'balanced') then ! inner iterations have a higher tolerance than outer iterations
    f_obj % tol_inner = sqrt(tol)
    f_obj % lmax      = kmax/2_i4b
   else if (method.eq.'fast') then
    f_obj % tol_inner = tol**(1._r8b/6._r8b) 
    f_obj % lmax      = kmax/6_i4b
   else if (method.eq.'minimal') then ! only one inner iteration
    f_obj % tol_inner = 10._r8b
    f_obj % lmax      = 1_i4b
   else
    if (f_obj % out_error) then
     write(f_obj % unit,'(a65)') "Error in f_set_tolerance: method argument not currently supported"
    end if
    stop
   end if
  else ! default values for classical iterations
   f_obj % tol_inner = 0._r8b
   f_obj % lmax      = 0_i4b
  end if
 end subroutine f_set_tolerance

! subroutine f_initial_guess(f_obj,method) ! no longer used: x0 is now interfaced to nested Newton object using a pointer
!  ! ** initial guess strategy for time-dependent algorithms for f_obj_base class **
!  ! note: it may be possible to add filtering techniques for the initial guess to improve efficiency
!  class(f_obj_base),intent(inout) :: f_obj
!  character(*),intent(in)         :: method
!
!  ! Note: avoid unintentional reallocation of object components (use array slices for assignment statements)
!  if (method.eq.'previous') then 
!   f_obj % x0(:) = f_obj % x1(:) ! initial guess -- solution from previous time step
!  else
!   if (f_obj % out_error) then
!    write(f_obj % unit,'(a66)') "Error in f_initial_guess: method argument not currently supported."
!   end if
!   stop
!  end if
! end subroutine f_initial_guess

 subroutine matrix_vector_product(f_obj,A,x,y)
  ! *** Compute matrix vector product y=A*x ***
  ! input
  class(f_obj_base),intent(inout) :: f_obj              ! class object containing solver options (intent(inout) so that y may be a component of f_obj)
  real(r8b),intent(in),contiguous :: A(:,:)             ! input matrix 
  real(r8b),intent(in),contiguous :: x(:)               ! input vector
    
  ! output
  real(r8b),intent(out),contiguous :: y(:)              ! product vector
    
  ! local variables
  character(1),parameter :: TRANS='N'                ! option for matrix transposition
  integer(i4b),parameter :: INCX=1_i4b, INCY=1_i4b   ! increment for elements of x and y vectors
  real(r8b),parameter    :: ALPHA=1._r8b,BETA=0._r8b ! scalars used in LAPACK solvers

  if (f_obj % banded) then ! banded storage
   ! note: we are passing the banded matrix A by reference to it's first element to avoid temporary array copies from the compiler
   associate(KL => f_obj % subdiag,KU => f_obj % superdiag)
    if (f_obj % linear_system_solver == LAPACK_standard) then
     !call DGBMV(TRANS,f_obj % n,f_obj % n,KL,KU,ALPHA,A(KL+1:f_obj % LDAF,:),f_obj % LDA,x,INCX,BETA,y,INCY) ! BLAS -- DGBMV uses different banded storage scheme compared to standard solver (works but creates temporary arrays)
     call DGBMV(TRANS,f_obj % n,f_obj % n,KL,KU,ALPHA,A(KL+1,1),f_obj % LDAF,x,INCX,BETA,y,INCY) ! BLAS -- pass A by reference because A uses difference banded storage scheme compared to standard solver (avoids temporary arrays)
    else if (f_obj % linear_system_solver == LAPACK_expert) then
     !call DGBMV(TRANS,f_obj % n,f_obj % n,KL,KU,ALPHA,A,LDA,x,INCX,BETA,y,INCY) ! BLAS -- DGBMV uses same banded storage scheme as expert solver
     call DGBMV(TRANS,f_obj % n,f_obj % n,KL,KU,ALPHA,A(1,1),f_obj % LDA,x,INCX,BETA,y,INCY) ! BLAS -- DGBMV uses same banded storage scheme as expert solver (pass A by refeence to be consistent with LAPACK_standard option above)
    end if
   end associate
  else ! full matrix storage
   associate(LDA => f_obj % LDA)
    call DGEMV(TRANS,f_obj % n,f_obj % n,ALPHA,A,LDA,x,INCX,BETA,y,INCY) ! BLAS
   end associate
  end if
 end subroutine matrix_vector_product

 !! ******************************* SUMMA procedures below ******************************* !!

 !subroutine SUMMA_get_scaled_Jacobian(f_obj,J,aJacScaled) ! no longer needed -- would need to take new storage scheme for J into account (standard LAPACK solver)
 ! ! ** Get scaled SUMMA Jacobian from nested Newton solver Jacobian **
 ! use matrixOper_module,  only: scaleMatrices
 ! ! arguments
 ! type(f_obj_type),intent(inout) :: f_obj ! nested Newton object
 ! !real(r8b),intent(in)           :: J(1:f_obj % nrow,1:f_obj % n)      ! nested Newton solver Jacobian
 ! real(r8b),intent(in)           :: J(:,:)                              ! nested Newton solver Jacobian
 ! !real(rkind),intent(out)        :: aJacScaled(f_obj % in_SS4HG % nLeadDim,f_obj % in_SS4HG % nState) ! scaled SUMMA Jacobian matrix
 ! real(rkind),intent(out)        :: aJacScaled(:,:) ! scaled SUMMA Jacobian matrix

 ! ! local
 ! real(rkind)    :: aJac(f_obj % in_SS4HG % nLeadDim,f_obj % in_SS4HG % nState) ! SUMMA Jacobian matrix (descaled)
 ! integer(i4b)   :: nBands   ! SUMMA's leading dimension for banded Jacobians
 ! integer(i4b)   :: err      ! SUMMA error code
 ! character(256) :: cmessage ! error message from SUMMA

 !   ! get SUMMA Jacobian from solver Jacobian
 !   if (f_obj % banded) then ! banded storage
 !    associate(nrow_banded => f_obj % nrow_banded, n => f_obj % n, subdiag => f_obj % subdiag)
 !     nBands=nrow_banded+subdiag
 !     aJac(1:subdiag,1:n) = 0._rkind
 !     aJac(subdiag+1:nBands,1:n) = J(1:nrow_banded,1:n) ! SUMMA's aJac has extra storage rows
 !    end associate
 !   else ! full matrix storage
 !    aJac(:,:) = J(:,:)
 !   end if

 !   ! Scale Jacobian
 !   associate(ixMatrix => f_obj % in_SS4HG % ixMatrix, nState => f_obj % in_SS4HG % nState)
 !    call scaleMatrices(ixMatrix,nState,aJac,f_obj % fScale,f_obj % xScale,aJacScaled,err,cmessage) ! matches solve_linear_system
 !   end associate
 !   if (err/=0) then
 !    if (f_obj % out_error) then
 !     write(f_obj % unit,*) "Error in SUMMA_get_scaled_Jacobian: scaleMatrices message="//trim(cmessage); stop
 !    end if
 !   end if

 !end subroutine SUMMA_get_scaled_Jacobian

 subroutine SUMMA_nested_line_search(f_obj,option,nested_algorithm,p)
  ! ** nested Newton line search **
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj            ! nested Newton object
  integer(i4b)     ,intent(in)    :: option           ! line search scheme option
  logical          ,intent(in)    :: nested_algorithm ! flag for nested algorithm (takes dynamic Newton iteration type selection mode into account)
  real(r8b)        ,intent(in),contiguous :: p(:) ! search direction

  ! local
  real(r8b) :: L0,L1 ! objective function values
  real(r8b) :: L1_prev ! objective function value from previous line search iteration
  real(r8b)            :: m ! local slope
  real(r8b), parameter :: c=1.e-4_r8b   ! objective function check control parameter
  real(r8b)            :: c_m ! c times m
  !real(r8b), parameter :: tao=0.5e0_r8b ! step reduction control parameter
  !real(r8b), parameter :: m_tol=0.1_r8b !100._r8b*epsilon(1._r8b)
  real(r8b)            :: alpha      ! step size
  real(r8b)            :: alpha_temp ! step size (temporary)
  real(r8b)            :: alpha_prev ! step size (from previous line search iteration)
  real(r8b)            :: rhs1,rhs2,aCoef,bCoef,disc ! constants for cubic interpolant
  logical   :: do_line_search
  logical   :: converged ! checkConv convergence flag
  integer(i4b), parameter :: i_max = 5_i4b ! max number of line search iterations
  integer(i4b)   :: i,j      ! loop index
  integer(i4b)   :: err      ! SUMMA error code
  character(256) :: cmessage ! error message from SUMMA
  logical, parameter :: debug_output=.false.

  ! working options: LS_I for inner iterations, LS_O for last inner iteration (outer scheme), LS_C for classical (lmax=0)

  ! initial solutions and option validation
  if (option == LS_I) then ! inner case
   f_obj % initial_solution(:) = f_obj % xkp1l(:) ! previous inner iteration
  else if ((option == LS_O).or.(option == LS_C)) then ! last inner iteration (LS_O) or classical (LS_C)
   !if (f_obj % nested) then
   if (nested_algorithm) then
    f_obj % initial_solution(:) = f_obj % xk0(:)   ! previous outer iteration
   else
    f_obj % initial_solution(:) = f_obj % xk(:)   ! previous outer iteration
   end if
  else
   print *, "Error in SUMMA_nested_line_search: option is not supported"
   stop
  end if

  ! get initial objective function (scaled)
  if (option == LS_O) then
   f_obj % f_temp(:) = f_obj % f1_vec(:) ! save value of f1 for scaled residual calculations in objective function procedure
   call f_obj % line_search_objective(nested_algorithm,.true.,.true.,option,f_obj % initial_solution,f_obj % f_temp,L0) ! need to compute when switching to outer scheme
  else if (f_obj % evaluate_B) then
   if (option == LS_I) then
    call f_obj % line_search_objective(nested_algorithm,.false.,.false.,option,f_obj % initial_solution,f_obj % f_temp,L0) ! can reuse f, J, and rVecScaled values
   end if
  else
   L0 = f_obj % L0 ! reuse from systemSolv or previous Newton iteration (classical and inner schemes)
  end if

  ! compute gradient of objective function (scaled)
  ! note: uses scaled Jacobian from LAPACK system (J for classical, Jdiff=J1-J2 for nested)
  call SUMMA_computeGradient(f_obj,f_obj % aJacScaled,f_obj % rVecScaled,f_obj % grad_L)

  ! compute local slope (use scaled values)
  ! note: descaled search direction p computed from LAPACK solve
  do concurrent (i = 1:f_obj % n)
   f_obj % p_scaled(i) = p(i)/f_obj % xScale(i)
  end do
  m = dot_product(f_obj % grad_L,f_obj % p_scaled) ! confirmed against homegrown line search
  c_m = c*m ! control parameter times local slope

  ! check that local slope is negative (needed to reduce the line search objective function)
  if (m < 0._rkind) then
   do_line_search = .true.
  !else if ((0._r8b <= m).and.(m <= m_tol)) then ! non-negative slope with allowance for round-off error
  ! do_line_search = .false. ! skip line search (there would be no improvement anyway)
  !else ! non-negative but exceeding tolerance
  ! print *, "m=",m
  ! print *, "option=",option
  ! print *, "Error in SUMMA_nested_line_search: initial slope is non-negative"
  ! stop
  else
   do_line_search = .false. ! skip line search (there would be no improvement anyway)
  end if

  if (debug_output) then
   print *, "NN Line Search:"
   print *, "option=",option
   print *, "m=",m
   print *, "L0=",L0
  end if
  
  alpha = 1._r8b ! initialize line search step size
  line_search: do i=1_i4b,i_max+1_i4b
   do concurrent (j = 1:f_obj % n)
    f_obj % updated_solution(j) = f_obj % initial_solution(j) + alpha*p(j)
   end do

   ! impose constraints (calls SUMMA's imposeConstraints routine)
   call f_obj % apply_constraints(f_obj % initial_solution,f_obj % updated_solution)

   ! compute objective function
   call f_obj % line_search_objective(nested_algorithm,.true.,.true.,option,f_obj % updated_solution,f_obj % f_temp,L1)

   !! check SUMMA's feasibility flag ------------------ turn this into a recoverable error
   !if (.not.(f_obj % feasible)) then
   ! print *, "Error in SUMMA_nested_line_search: not feasible"
   ! stop
   !end if

   ! get convergence flag
   !if (option == LS_I) then
   ! converged = .false.
   !else
    converged = SUMMA_checkConv(f_obj,p,f_obj % updated_solution)
   !end if   

   if (debug_output) then
    print *, "i=",i
    print *, "alpha=",alpha
    print *, "L1=",L1
    print *, "L0 + alpha*c*m=",L0 + alpha*c*m
   end if

   ! check convergence
   if (converged) then
    call f_and_J_values ! obtain f1, f2, J1, and J2 values needed for next nested Newton iteration
    return 
   end if
   
   ! exit early if not doing the line search (only alpha=1.0 solution is used)
   if (.not.do_line_search) then
    call f_and_J_values ! obtain f1, f2, J1, and J2 values needed for next nested Newton iteration
    return
   end if

   ! check if the objective function is accepted using the Armijo-Goldstein Criterion
   !if (L1 <= L0 + alpha*c*m) then
   if (L1 <= L0 + alpha*c_m) then
    call f_and_J_values ! obtain f1, f2, J1, and J2 values needed for next nested Newton iteration
    return
   end if

   ! * adjust the step size in preparation for next line search iteration *
   !alpha = alpha * tao ! basic reduction by a constant factor

   if (i == 1_i4b) then ! first backtrack: use quadratic
    alpha_temp = -m / ( 2._r8b*(L1 - L0 - m) )
    if (alpha_temp > 0.5_r8b*alpha) alpha_temp = 0.5_r8b*alpha

   else ! subsequent backtracks: use cubic
    ! define rhs
    rhs1 = L1      - L0 - alpha     *m
    rhs2 = L1_prev - L0 - alpha_prev*m

    ! define coefficients
    aCoef = (rhs1/(alpha**2_i4b) - rhs2/(alpha_prev**2_i4b))/(alpha - alpha_prev)
    bCoef = (-alpha_prev*rhs1/(alpha**2_i4b) + alpha*rhs2/(alpha_prev**2_i4b)) / (alpha - alpha_prev)

    if (aCoef == 0._r8b) then ! check if a quadratic
     alpha_temp = -m/(2._r8b*bCoef)

    else ! calculate cubic

     ! only allow real roots of the cubic 
     disc = bCoef**2_i4b - 3._r8b*aCoef*m ! discriminant?
     if (disc < 0._r8b) then
      alpha_temp = 0.5_r8b*alpha
     else
      alpha_temp = (-bCoef + sqrt(disc))/(3._r8b*aCoef)
     end if

    end if

     ! constrain to <= 0.5*alpha
     if (alpha_temp > 0.5_r8b*alpha) alpha_temp=0.5_r8b*alpha

   end if

   ! save results
   alpha_prev = alpha
   L1_prev = L1

   ! constrain lambda and finalize
   alpha = max(alpha_temp, 0.1_r8b*alpha)


   ! if stopping criterion not reached within the max # of line search iterations, use full Newton step with constraints imposed
   ! note: an extra loop iteration is performed to get the constrained full Newton step solution  --- may be able to reduce expense by saving this solution earlier
   if (i == i_max) then
    do_line_search = .false.
    alpha = 1._r8b
   end if

  end do line_search

  ! should not be able to reach this point
  print *, "Error in SUMMA_nested_line_search: no return criteria were triggered"

 contains
 
  subroutine f_and_J_values
   ! ** post-processing to obtain function and Jacobian values needed for next nested Newton iteration **
   logical, parameter :: periodic_J1 = .false. ! periodic evaluation of J1 during inner loop?
   logical :: evaluate_J1 ! flag for only evaluating J1 on every other inner iteration

   ! update solution in nested Newton algorithm
   !if (f_obj % nested) then
   if (nested_algorithm) then
    f_obj % xkp1lp1(:) = f_obj % updated_solution(:) ! apply updated solution (all cases -- to be used in case of early loop exit)
   else
    f_obj % xkp1(:)    = f_obj % updated_solution(:) ! apply updated solution (all cases -- to be used in case of early loop exit)
   end if

   if (.not.converged) then ! if outer/classical iterations not converged, prep for next Newton iteration (if applicable)

    ! obtain remaining function and Jacobian variables needed for next Newton iteration
    if (option == LS_C) then

     !if (f_obj % nested) then
     if (nested_algorithm) then
      if (f_obj % k < f_obj % kmax) then ! not required for last outer iteration
       ! have f -- need f1, J1, f2, and J2
       !call filter_SUMMA_f(.false.,f_obj % stateMask1,f_obj % f_vec,f_obj % f1_vec) ! get f1 from total f - now obtained directly from eval8summa
       !call filter_SUMMA_f(.true.,f_obj % stateMask2,f_obj % f_vec,f_obj % f2_vec) ! get f2 from total f - now obtained directly from eval8summa
       !call f_obj % J1_J2_eval(f_obj % updated_solution) ! get J1 and J2 based on previous eval8summa call (used to compute f)
       call f_obj % J1_eval(f_obj % updated_solution) ! get J1 based on previous eval8summa call (used to compute f)
       call f_obj % J2_eval(f_obj % updated_solution) ! get J2 based on previous eval8summa call (used to compute f)
       f_obj % L0 = L1 ! store previous objective function value
      end if
     else
      if (f_obj % k < f_obj % kmax_classical) then ! not required for last outer iteration
       ! have f -- need J
       call f_obj % J_eval(f_obj % updated_solution)
       f_obj % L0 = L1 ! store previous objective function value
      end if 
     end if

    else if (option == LS_I) then
     ! have f and f1 -- need J1 (f2 and J2 don't change)

     ! SJT: testing computing J1 every other inner iteration
     if (periodic_J1) then
      ! check if f_obj % l is even
      if (mod(f_obj % l,2_i4b) == 0_i4b) then ! if l is even
       evaluate_J1 = .false. ! reuse previous value for next odd l value
      else ! if l is odd
       evaluate_J1 = .true. ! evaluate J1 to be used for next even l value
      end if
      if (evaluate_J1) call f_obj % J1_eval(f_obj % updated_solution)   
     else ! evaluate J1 on every inner iteration
      call f_obj % J1_eval(f_obj % updated_solution)   ! get J1 based on eval8summa call for total f 
     end if 

     f_obj % L0 = L1 ! store previous inner scheme objective function value (does not apply if switching to outer line search scheme)

    else if (option == LS_O) then

     if (f_obj % k < f_obj % kmax) then ! not required for last outer iteration
      ! have f and f2 -- need J2, f1, J1
      !call filter_SUMMA_f(.false.,f_obj % stateMask1,f_obj % f_vec,f_obj % f1_vec) ! get f1 from total f - now obtained directly from eval8summa
      !call f_obj % J1_J2_eval(f_obj % updated_solution) ! get J1 and J2 based on previous eval8summa call (used to compute f2)
      call f_obj % J1_eval(f_obj % updated_solution) ! get J1 based on previous eval8summa call (used to compute f2)
      call f_obj % J2_eval(f_obj % updated_solution) ! get J2 based on previous eval8summa call (used to compute f2)
     end if

    else
      print *, "Error in SUMMA_nested_line_search: option is not supported"; stop
    end if

   end if

  end subroutine f_and_J_values

 end subroutine SUMMA_nested_line_search

 subroutine SUMMA_line_search_objective(f_obj,nested_algorithm,evaluate_f,evaluate_rVecScaled,option,solution,f_temp,L)
  ! ** compute line search objective function for SUMMA **
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj ! nested Newton object
  logical          ,intent(in)    :: nested_algorithm ! flag for nested algorithm (takes dynamic Newton iteration type selection mode into account)
  logical          ,intent(in)    :: evaluate_f ! perform evaluations for f, f1, or f2? (if not, use stored values)
  logical          ,intent(in)    :: evaluate_rVecScaled ! perform evaluations for rVecScaled (if not, use stored values)
  integer(i4b)     ,intent(in)    :: option ! line search option
  real(r8b)        ,intent(in),contiguous :: solution(:) ! updated solution vector
  real(r8b)        ,intent(in),contiguous :: f_temp(:) ! storage vector for f1 for 'L' scheme
  real(r8b)        ,intent(out)   :: L ! objective function value
  ! local
  integer(i4b) :: i ! loop index

  if (option == LS_C) then ! classical case
   !if (evaluate_f) call f_obj % f_vec_eval(solution) ! update total f
   if (evaluate_f) then
    if (nested_algorithm) then
     call f_obj % f1_f2_vec_eval(solution) ! update total f (also f1 and f2 for next nested iteration)
    else
     call f_obj % f_vec_eval(solution) ! update total f
    end if
   end if
   !if (evaluate_rVecScaled) f_obj % rVecScaled(:) = f_obj % fScale(:) * f_obj % f_vec ! now obtained from f_vec_eval (e.g., eval8summa)
   L=f_obj % out_SS4HG % fNew ! scaled
  else if (option == LS_I) then ! inner case
   if (evaluate_f) call f_obj % f1_vec_eval(solution) ! update f1 (and f)
   if (evaluate_rVecScaled) then
    do concurrent (i = 1:f_obj % n)
     f_obj % vector(i) = solution(i) - f_obj % xk0(i)
    end do
    call f_obj % matrix_vector_product(f_obj % J2,f_obj % vector,f_obj % rVecScaled) 
    do concurrent (i = 1:f_obj % n)
     f_obj % rVecScaled(i) = f_obj % fScale(i) * ( f_obj % f1_vec(i) &
                         & - ( f_obj % f2_vec(i) + f_obj % rVecScaled(i) )&
                         & )
    end do
    !f_obj % rVecScaled(:) = f_obj % fScale(:) * ( f_obj % f1_vec(:) &
    !                    & - ( f_obj % f2_vec(:) + f_obj % matrix_vector_product(f_obj % J2,solution - f_obj % xk0) )&
    !                    & )
   end if
   L = 0.5_r8b*dot_product(f_obj % rVecScaled,f_obj % rVecScaled)
  else if (option == LS_O) then ! last inner iteration case
   ! OG
   !if (evaluate_f) call f_obj % f2_vec_eval(solution) ! update f2 (and f which is used for checkConv)
   !if (evaluate_rVecScaled) then 
   ! f_obj % rVecScaled(:) = f_obj % fScale(:) * ( f_obj % f1_vec(:) &
   !                     & + f_obj % matrix_vector_product(f_obj % J1,solution - f_obj % xkp1l) - f_obj % f2_vec(:)&
   !                     & )
   !end if
   ! ***** f1 is overwritten after f1_f2_vec_eval call, throwing off remaining LS iterations *****
   if (f_obj % k < f_obj % kmax) then ! not required for last outer iteration
    !if (evaluate_rVecScaled) then ! before overwriting f1, save contribution to scaled residual vector 
    ! rVec_temp(:) = f_obj % f1_vec(:)
    !end if
    if (evaluate_f) call f_obj % f1_f2_vec_eval(solution) ! update f2 (also f which is used for checkConv and f1 for next iteration)
    if (evaluate_rVecScaled) then ! f1 stored within f_temp because f1_vec is overwritten
     do concurrent (i = 1:f_obj % n)
      f_obj % vector(i) = solution(i) - f_obj % xkp1l(i)
     end do
     call f_obj % matrix_vector_product(f_obj % J1,f_obj % vector,f_obj % rVecScaled)
     do concurrent (i = 1:f_obj % n)
      f_obj % rVecScaled(i) = f_obj % fScale(i) * ( f_temp(i) &
                          & + f_obj % rVecScaled(i) - f_obj % f2_vec(i)&
                          & )
     end do
     !f_obj % rVecScaled(:) = f_obj % fScale(:) * ( f_temp(:) &
     !                    & + f_obj % matrix_vector_product(f_obj % J1,solution - f_obj % xkp1l) - f_obj % f2_vec(:)&
     !                    & )
    end if
   else ! don't need f1 for last outer iteration
    if (evaluate_f) call f_obj % f2_vec_eval(solution) ! update f2 (and f which is used for checkConv)
    if (evaluate_rVecScaled) then 
     do concurrent (i = 1:f_obj % n)
      f_obj % vector(i) = solution(i) - f_obj % xkp1l(i)
     end do
     call f_obj % matrix_vector_product(f_obj % J1,f_obj % vector,f_obj % rVecScaled)
     do concurrent (i = 1:f_obj % n)
      f_obj % rVecScaled(i) = f_obj % fScale(i) * ( f_obj % f1_vec(i) &
                          & + f_obj % rVecScaled(i) - f_obj % f2_vec(i)&
                          & )
     end do
     !f_obj % rVecScaled(:) = f_obj % fScale(:) * ( f_obj % f1_vec(:) &
     !                    & + f_obj % matrix_vector_product(f_obj % J1,solution - f_obj % xkp1l) - f_obj % f2_vec(:)&
     !                    & )
    end if
   end if
   L = 0.5_r8b*dot_product(f_obj % rVecScaled,f_obj % rVecScaled)
  else 
   print *, "Error in SUMMA_line_search_objective: option is not supported"; stop
  end if
 end subroutine SUMMA_line_search_objective

 subroutine SUMMA_computeGradient(f_obj,aJacScaled,rVecScaled,gradScaled)
  ! ** compute product of scaled residual and scaled SUMMA Jacobian **
  use matrixOper_module, only: computGradient
  ! arguments
  type(f_obj_type),intent(in) :: f_obj ! nested Newton object
  !real(rkind),intent(in) :: aJacScaled(f_obj % in_SS4HG % nLeadDim,f_obj % in_SS4HG % nState) ! scaled SUMMA Jacobian matrix
  real(rkind),intent(in),contiguous :: aJacScaled(:,:) ! scaled SUMMA Jacobian matrix
  real(r8b),intent(in)  ,contiguous :: rVecScaled(:) ! gradient of objective function L
  real(r8b),intent(out) ,contiguous :: gradScaled(:) ! gradient of objective function L

  ! local
  integer(i4b)   :: err      ! SUMMA error code
  character(256) :: cmessage ! error message from SUMMA

  associate(ixMatrix => f_obj % in_SS4HG % ixMatrix, nState => f_obj % in_SS4HG % nState)
   call computGradient(ixMatrix,nState,aJacScaled,rVecScaled,gradScaled,err,cmessage)
  end associate
  if (err/=0) then
   print *, "Error in SUMMA_computeGradient: "//trim(cmessage)
   stop
  end if
 end subroutine SUMMA_computeGradient

 subroutine SUMMA_scaling(f_obj,B)
  ! ** apply scaling from SUMMA's fScale and xScale vectors to matrix and RHS for LAPACK **
  use matrixOper_module,  only: scaleMatrices
  ! input
  class(f_obj_type),intent(inout) :: f_obj ! nested Newton object

  ! input-output
  real(r8b),intent(inout),contiguous :: B(:,:) ! right-hand side vector

  ! local
  integer(i4b)   :: i
  integer(i4b)   :: err
  character(256) :: cmessage

  ! get scaled variables (accoring to SUMMA's fScale and xScale vectors)
  ! note: need to match scaling applied in solve_linear_system subroutine in summaSolve4homegrown

  ! if computing RHS vector, scale in preparation for LAPACK
  if (f_obj % evaluate_B) then
   do concurrent (i = 1:f_obj % n)
    B(i,1) = f_obj % fScale(i) * B(i,1)
    f_obj % rVecScaled(i) = -B(i,1) ! store for reuse
   end do
  end if

  associate(&
   ixMatrix => f_obj % in_SS4HG % ixMatrix , & ! type of matrix (full or band diagonal)
   nState   => f_obj % in_SS4HG % nState   , & ! number of state variables in the current subset
   fScale   => f_obj % fScale              , & 
   xScale   => f_obj % xScale                & 
  &)
   call scaleMatrices(ixMatrix,nState,f_obj % AF,fScale,xScale,f_obj % aJacScaled,err,cmessage)
  end associate
  if (err/=0) then
   if (f_obj % out_error) then
    write(f_obj % unit,*) "Error in SUMMA_scaling: scaleMatrices message="//trim(cmessage); stop
   end if
  end if

  f_obj % AF(:,:) = f_obj % aJacScaled(:,:) ! load AF matrix for LAPACK
 end subroutine SUMMA_scaling

 subroutine SUMMA_descaling(f_obj,B)
  ! ** apply descaling from SUMMA's xScale vector to solution for LAPACK **
  ! input
  class(f_obj_type),intent(in) :: f_obj ! nested Newton object

  ! input-output
  real(r8b),intent(inout),contiguous :: B(:,:) ! solution side vector

  ! local
  integer(i4b)   :: i

  do concurrent (i = 1:f_obj % n)
   B(i,1) = B(i,1) * f_obj % xScale(i)
  end do
  
 end subroutine SUMMA_descaling

 function SUMMA_check_convergence_flag(f_obj) result(converged)
  ! ** check convergence flag from out_SS4HG object found during Newton step refinement **
  ! input
  class(f_obj_type),intent(in)   :: f_obj

  ! output
  logical :: converged

  converged = f_obj % out_SS4HG % converged

 end function SUMMA_check_convergence_flag

 function SUMMA_checkConv(f_obj,step,xvec1) result(converged)
  ! ** interface for SUMMA's checkConv subroutine **
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  !real(r8b),intent(in)            :: step(1:f_obj % n)  ! Newton step (iteration increment)
  !real(r8b),intent(in)            :: xvec1(1:f_obj % n) ! updated solution vector
  real(r8b),intent(in),contiguous  :: step(:)  ! Newton step (iteration increment)
  real(r8b),intent(in),contiguous  :: xvec1(:) ! updated solution vector

  ! output
  logical :: converged

  ! local variables
  integer(i4b) :: mSoil             ! number of soil layers in the solution vector

  ! get the number of soil layers in the solution vector
  mSoil = size(f_obj % indx_data % var(iLookINDEX % ixMatOnly) % dat)

  converged = checkConv(mSoil,f_obj % in_SS4HG,f_obj % mpar_data,f_obj % indx_data,f_obj % prog_data,&
                       &f_obj % f_vec,step,xvec1,f_obj % out_SS4HG)

 end function SUMMA_checkConv

 subroutine SUMMA_imposeConstraints(f_obj,xvec0,xvec1)
  ! ** interface for SUMMA's imposeConstraints subroutine **
  class(f_obj_type),intent(inout)    :: f_obj
  !real(r8b),intent(in)              :: xvec0(1:f_obj % n) ! previous guess vector
  !real(r8b),intent(inout)           :: xvec1(1:f_obj % n) ! current guess vector
  real(r8b),intent(in)   ,contiguous :: xvec0(:) ! previous guess vector
  real(r8b),intent(inout),contiguous :: xvec1(:) ! current guess vector

  ! increment the proposed iteration for simple error control if needed
  associate(&
   err => f_obj % out_SS4HG % err, message => f_obj % out_SS4HG % message & ! SUMMA error code and message 
  &)
     call imposeConstraints(f_obj % model_decisions,f_obj % indx_data,f_obj % prog_data,f_obj % mpar_data,& ! data structures
                           &xvec1,xvec0,&                                                                   ! state variables
                           & f_obj % in_SS4HG % nState, f_obj % in_SS4HG % nSoil,f_obj % in_SS4HG % nSnow,& ! layer variables
                           & message, err)                                                                  ! error control
     if (err /= 0) then
      if (f_obj % out_error) then
       write(f_obj % unit,*) "Error in SUMMA_imposeConstraints: imposeConstraints message="//trim(message); stop
      end if
     end if
  end associate
 end subroutine SUMMA_imposeConstraints

 subroutine SUMMA_eval8summa(f_obj,xvec)
  ! ** interface for SUMMA's eval8summa subroutine **
  ! compute SUMMA derivative values and residual vector
  ! note: - eval8summa was not refactored to use object arguments
  !       - objects for summaSolve4homegrown were reused where possible
  class(f_obj_inputs),intent(inout) :: f_obj
  !real(r8b),intent(in),contiguous   :: xvec(:) ! current guess
  real(r8b),intent(in)   :: xvec(:) ! current guess
  logical,parameter :: mass_flag=.true.,energy_flag=.true.

  ! update
  associate(&
   stateVecTrial => xvec & ! current guess for state vector
  &)
   call eval8summa(&
                    ! input: model control
                    f_obj % in_SS4HG % dt_cur,         & ! intent(in):    current stepsize
                    f_obj % in_SS4HG % dt,             & ! intent(in):    length of the entire time step (seconds) for drainage pond rate
                    f_obj % in_SS4HG % nSnow,          & ! intent(in):    number of snow layers
                    f_obj % in_SS4HG % nSoil,          & ! intent(in):    number of soil layers
                    f_obj % in_SS4HG % nLayers,        & ! intent(in):    number of layers
                    f_obj % in_SS4HG % nState,         & ! intent(in):    number of state variables in the current subset
                    .false.,                           & ! intent(in):    not inside Sundials solver
                    f_obj % in_SS4HG % firstSubStep,   & ! intent(in):    flag to indicate if we are processing the first sub-step
                    f_obj % io_SS4HG % firstFluxCall,  & ! intent(inout): flag to indicate if we are processing the first flux call
                    .false.,                           & ! intent(in):    flag to indicate if we are processing the first flux call in a splitting operation (.false. based on usage of eval8summa in summaSolve4homegrown)
                    f_obj % in_SS4HG % computeVegFlux, & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                    f_obj % in_SS4HG % scalarSolution, & ! intent(in):    flag to indicate the scalar solution
                    mass_flag,                         & ! intent(in):    flag to compute mass terms
                    energy_flag,                       & ! intent(in):    flag to compute energy terms
                    .true.,.false.,.false.,.false.,.false., & ! intent(in):    flag to compute f, f1, and f2 for nested Newton (classical iterations assumed)
                    .true.,.true.,                     & ! intent(in):    flags to compute mass and energy Jacobian terms
                    ! input: state vectors
                    stateVecTrial,                   & ! intent(in):    model state vector
                    f_obj % fScale,                  & ! intent(in):    characteristic scale of the function evaluations
                    f_obj % sMul,                    & ! intent(inout): state vector multiplier (used in the residual calculations)
                    ! input: data structures
                    f_obj % model_decisions,         & ! intent(in):    model decisions
                    f_obj % lookup_data,             & ! intent(in):    lookup tables
                    f_obj % type_data,               & ! intent(in):    type of vegetation and soil
                    f_obj % attr_data,               & ! intent(in):    spatial attributes
                    f_obj % mpar_data,               & ! intent(in):    model parameters
                    f_obj % forc_data,               & ! intent(in):    model forcing data
                    f_obj % bvar_data,               & ! intent(in):    average model variables for the entire basin
                    f_obj % prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                    ! input-output: data structures
                    f_obj % indx_data,               & ! intent(inout): index data
                    f_obj % diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                    f_obj % flux_data,               & ! intent(inout): model fluxes for a local HRU (initial flux structure)
                    f_obj % deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                    ! input-output: baseflow
                    f_obj % io_SS4HG % ixSaturation, & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                    !f_obj % dBaseflow_dMatric,       & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                    f_obj % dBaseflow_dWat,          & ! intent(out):   derivative in baseflow w.r.t. water content (s-1)
                    f_obj % dBaseflow_dTk,           & ! intent(out):   derivative in baseflow w.r.t. temperature (s-1)
                    ! output
                    f_obj % feasible,                & ! intent(out):   flag to denote the feasibility of the solution
                    f_obj % fluxVec0,                & ! intent(out):   flux vector
                    f_obj % f_vec,                   & ! intent(inout): f vector
                    f_obj % f1_vec,                  & ! intent(inout): f1 vector
                    f_obj % f2_vec,                  & ! intent(inout): f2 vector
                    f_obj % fRHS,                    & ! intent(out):   RHS function for ARKODE
                    f_obj % rAdd,                    & ! intent(out):   additional (sink) terms on the RHS of the state equation
                    f_obj % resVec,                  & ! intent(out):   residual vector
                    f_obj % rVecScaled,              & ! intent(out):   scaled residual vector
                    f_obj % out_SS4HG % fNew,        & ! intent(out):   function evaluation
                    f_obj % out_SS4HG % err,         & ! intent(out): error code
                    f_obj % out_SS4HG % message)       ! intent(out): error message (note: eval8summa uses "cmessage" instead)
  end associate
  
  ! finalize
  associate(err => f_obj % out_SS4HG % err, message => f_obj % out_SS4HG % message) 
   if (err /= 0) then
    f_obj % f_error = .true.
    if (f_obj % out_warning) then
     write(f_obj % unit,*) "Error SUMMA_eval8summa: eval8summa message="//trim(message)
    end if
   end if
  end associate
 end subroutine SUMMA_eval8summa

 subroutine SUMMA_computJacob(f_obj,mass_flag,energy_flag,&
                             &indx_data,diag_data,flux_data,deriv_data,&
                             &dMat,dBaseflow_dWat,dBaseflow_dTk,&
                             &aJac)
  ! ** Interface for SUMMA's computJacob subroutine **
  ! arguments
  class(f_obj_inputs),intent(inout) :: f_obj
  logical          ,intent(in)      :: mass_flag,energy_flag  ! flags for evaluating mass and energy Jacobians
  type(var_ilength),intent(in)      :: indx_data              ! indices defining model states and layers for selected split 
  type(var_dlength),intent(in)      :: diag_data              ! diagnostic variables for a local HRU
  type(var_dlength),intent(in)      :: flux_data              ! flux data
  type(var_dlength),intent(in)      :: deriv_data             ! derivative data
  !real(rkind)      ,intent(in) ,contiguous :: dMat(:)          ! diagonal matrix (no flux derivatives) for split
  !real(rkind)      ,intent(in) ,contiguous :: dBaseflow_dWat(:,:)    ! derivative in baseflow w.r.t. water content (s-1)
  !real(rkind)      ,intent(in) ,contiguous :: dBaseflow_dTk(:,:)     ! derivative in baseflow w.r.t. temperature (s-1)
  !real(rkind)      ,intent(out),contiguous :: aJac(:,:) ! SUMMA's unscaled Jacobian matrix
  real(rkind)      ,intent(in)  :: dMat(:)          ! diagonal matrix (no flux derivatives) for split
  real(rkind)      ,intent(in)  :: dBaseflow_dWat(:,:)    ! derivative in baseflow w.r.t. water content (s-1)
  real(rkind)      ,intent(in)  :: dBaseflow_dTk(:,:)     ! derivative in baseflow w.r.t. temperature (s-1)
  real(rkind)      ,intent(out) :: aJac(:,:) ! SUMMA's unscaled Jacobian matrix

  ! local variables
  type(in_type_computJacob)  :: in_computJacob  ! computJacob input object
  type(out_type_computJacob) :: out_computJacob ! computJacob output object  

  ! initialize
  ! *** Transfer data to in_computJacob class object from local variables in summaSolve4homegrown ***
  associate(&
   ixGroundwater  => f_obj % model_decisions(iLookDECISIONS%groundwatr)%iDecision,&  ! intent(in): [i4b] groundwater parameterization
   dt_cur         => f_obj % in_SS4HG % dt_cur         ,& ! intent(in): current stepsize
   nSnow          => f_obj % in_SS4HG % nSnow          ,& ! intent(in): number of snow layers
   nSoil          => f_obj % in_SS4HG % nSoil          ,& ! intent(in): number of soil layers
   nLayers        => f_obj % in_SS4HG % nLayers        ,& ! intent(in): total number of layers
   ixRichards     => f_obj % model_decisions(iLookDECISIONS%f_Richards)%iDecision,&  ! intent(in): form of Richards' equation
   ixMatrix       => f_obj % in_SS4HG % ixMatrix       ,& ! intent(in): type of matrix (full or band diagonal)
   computeVegFlux => f_obj % in_SS4HG % computeVegFlux  & ! intent(in): flag to indicate if computing fluxes over vegetation
  &)   
   call in_computJacob % initialize(dt_cur,nSnow,nSoil,nLayers,computeVegFlux,(ixGroundwater==qbaseTopmodel),ixRichards,&
                                   &ixMatrix,mass_flag,energy_flag)
  end associate 

   ! update
   associate(&
    prog_data         => f_obj % prog_data&         ! prognostic variables for a local HRU
   &)
    call computJacob(in_computJacob,indx_data,prog_data,diag_data,deriv_data,dBaseflow_dWat,dBaseflow_dTk,dMat,&
                    &aJac,out_computJacob)
   end associate

  ! finalize
  ! *** Transfer data from out_computJacob class object to local variables in summaSolve4homegrown ***
  ! note: "message" used for out_SS4HG data component but "cmessage" used within summaSolve4homegrown subroutine
  associate(err => f_obj % out_SS4HG % err, cmessage => f_obj % out_SS4HG % message) 
   call out_computJacob % finalize(err,cmessage)
   if (err /= 0) then
    if (f_obj % out_error) then
     write(f_obj % unit,*) "Error in SUMMA_computJacob: computJacob message="//trim(cmessage); stop
    end if
   end if
  end associate

 end subroutine SUMMA_computJacob

 subroutine get_SUMMA_f1_f2_flags(f_obj)
  ! *** Compute flags for f1 and f2 evaluations for mass and energy state variables in SUMMA ***
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj

  ! determine logical flags for mass and energy terms
  if (f_obj % nested) then ! nested iterations
   if (f_obj % dual) then ! assign masks for f1 (energy) and f2 (mass)
    f_obj % f1_mass_flag = .false.; f_obj % f1_energy_flag = .true.
    f_obj % f2_mass_flag = .true.; f_obj % f2_energy_flag = .false.
   else ! assign masks for f1 (mass) and f2 (energy)
    f_obj % f1_mass_flag = .true.; f_obj % f1_energy_flag = .false.
    f_obj % f2_mass_flag = .false.; f_obj % f2_energy_flag = .true.
   end if
  else ! classical iterations
    f_obj % f1_mass_flag = .false.; f_obj % f1_energy_flag = .false.
    f_obj % f2_mass_flag = .false.; f_obj % f2_energy_flag = .false.
  end if

 end subroutine get_SUMMA_f1_f2_flags

 subroutine get_SUMMA_mass_energy_masks(f_obj)
  ! *** Compute masks for mass and energy state variables for SUMMA ***
  ! ** NOTE: the fully-coupled solution method in SUMMA's opSplittin is assumed **
  use stateFilter_module,only: fullyCoupled,stateTypeSplit
  use stateFilter_module,only: massSplit,nrgSplit
  use stateFilter_module,only: fullDomain,subDomain
  use stateFilter_module,only: vector,scalar

  ! arguments
  class(f_obj_type),intent(inout) :: f_obj

  ! local variables
  type(split_select_type)         :: split_select      ! split select object
  character(LEN=256)              :: message           ! total error message
  character(LEN=256)              :: cmessage          ! error message of downwind routine
  integer(i4b)                    :: err               ! error code of downwind routine
  logical(lgt)                    :: return_flag
  !logical(lgt),parameter          :: dual = .true. !.false. = f1->mass, f2->energy, .true. = f1->energy, f2->mass

  ! * initialize operations for split_select object *

  !associate(nstate => f_obj % in_SS4HG % nState)
  associate(nstate => f_obj % n)
   ! initialize total # of state variables
   split_select % nState = nState 

   ! allocate data components
   allocate(split_select % stateMask(1:nState)) ! allocate split_select components
  end associate

  ! use split_select_type object to specify the desired split
  ! NOTE: we are computing the energy state mask and negating to find the mass state mask (to include pressure head state variables)
  !split_select % iSplit =                      ! iteration counter for split_select_loop (not used)
  split_select % ixCoupling = stateTypeSplit    ! splitting is used
  split_select % iStateTypeSplit = nrgSplit     ! state variable type
  split_select % ixStateThenDomain = fullDomain ! do not split the domain into sub-domains 
  !split_select % iDomainSplit =                ! only used for sub-domain splitting
  split_select % ixSolution = vector            ! vector split (not scalar)
  !split_select % iStateSplit =                 ! only used for scalar splits

  ! apply steps similar to initialize_split from opSplitting to generate logical masks
  ! note: from update_stateMask in opSplittin

  ! compute stateMask and nSubset (in split_select object) for the selected split
  call split_select % get_stateMask(f_obj % indx_data,err,cmessage,message,return_flag)
  if (return_flag) then
    if (f_obj % out_error) then
     write(f_obj % unit,*) "Error in f_state_SUMMA_vec: stateFilter message="//trim(cmessage); stop
    end if
  end if

  ! assign masks
  if (f_obj % dual) then ! assign masks for f1 (energy) and f2 (mass)
   f_obj % stateMask2 = .not.(split_select % stateMask(:)) ! negate energy mask to find mass mask --- allocate on assignment
   f_obj % stateMask1 = split_select % stateMask           ! no transformation --- allocate on assignment
  else ! assign masks for f1 (mass) and f2 (energy)
   f_obj % stateMask1 = .not.(split_select % stateMask(:)) ! negate energy mask to find mass mask --- allocate on assignment
   f_obj % stateMask2 = split_select % stateMask           ! no transformation --- allocate on assignment
  end if

  ! get counts for mass and energy splits
  if (f_obj % dual) then
   f_obj % nSubset2 = split_select % nState - split_select % nSubset ! transform to get count for mass split
   f_obj % nSubset1 = split_select % nSubset                         ! no transformation for energy split
  else
   f_obj % nSubset1 = split_select % nState - split_select % nSubset ! transform to get count for mass split
   f_obj % nSubset2 = split_select % nSubset                         ! no transformation for energy split
  end if

  if ((f_obj % nSubset1 == 0_i4b).or.(f_obj % nSubset2 == 0_i4b)) then
    if (f_obj % out_error) then
     write(f_obj % unit,*) "Error in get_SUMMA_mass_energy_masks: empty stateMask detected"
     stop
    end if
  end if
 end subroutine get_SUMMA_mass_energy_masks

 subroutine filter_SUMMA_f(negative,stateMask,f_total,f_filter)
  ! filter total f from SUMMA into f1 or f2 for use in nested Newton solver (use appropriate stateMask)
  ! arguments
  logical         ,intent(in)    :: negative      ! apply a negative sign to f_filter?
  logical(lgt)    ,intent(in)    :: stateMask(:)  ! logical mask for filtering
  real(r8b)       ,intent(in) ,contiguous :: f_total(:)  ! total f in nested Newton solver
  real(r8b)       ,intent(out),contiguous :: f_filter(:) ! filtered f in nested Newton solver 

  ! assign non-zero function values based on logical mask
  f_filter(:)=0._r8b
  if (negative) then
   f_filter(:)=merge(-f_total,f_filter,stateMask)
  else
   f_filter(:)=merge(f_total,f_filter,stateMask)
  end if
 end subroutine filter_SUMMA_f

 subroutine filter_SUMMA_Jacobian(f_obj,negative,stateMask,aJac,J_total,J_filter)
  ! filter total Jacobian from SUMMA into J1 or J2 for use in nested Newton solver (use appropriate stateMask)
  ! arguments
  type(f_obj_type),intent(inout) :: f_obj         ! nested Newton object
  logical         ,intent(in)    :: negative      ! apply a negative sign to J_filter?
  logical(lgt)    ,intent(in)   ,contiguous :: stateMask(:)  ! logical mask for filtering
  real(rkind)     ,intent(inout),contiguous :: aJac(:,:)     ! total Jacobian from SUMMA
  real(r8b)       ,intent(out)  ,contiguous :: J_total(:,:)  ! total Jacobian in nested Newton solver storage scheme
  real(r8b)       ,intent(out)  ,contiguous :: J_filter(:,:) ! filtered Jacobian in nested Newton solver storage scheme 

  ! local variables
  integer(i4b) :: i,j,k  ! loop indices
  integer(i4b) :: nBands ! # of bands for banded storage

  ! store Jacobian used in solver
  if (f_obj % banded) then ! banded storage
   associate(nrow_banded => f_obj % nrow_banded, n => f_obj % n, subdiag => f_obj % subdiag, superdiag => f_obj % superdiag)
    nBands=nrow_banded+subdiag ! number of non-zero bands
    J_total(1:nrow_banded,1:n) = aJac(subdiag+1:nBands,1:n) ! aJac has extra storage rows (store total Jacobian)
    J_filter(1:nrow_banded,1:n) = 0._r8b ! initialize -- may not be needed
    if (negative) then
     do j=1,n ! column index for dense and banded storage
      do i=max(1,j-superdiag),min(n,j+subdiag) ! row index for dense storage
       k = nrow_banded+i-j ! row index for LAPACK banded storage (nrow_banded = subdiag+superdiag+1)
       aJac(k,j) = merge(-aJac(k,j),0._rkind,stateMask(i)) ! zero the elements that are not included in J2
      end do
     end do
    else
     do j=1,n ! column index for dense and banded storage
      do i=max(1,j-superdiag),min(n,j+subdiag) ! row index for dense storage
       k = nrow_banded+i-j ! row index for LAPACK banded storage (nrow_banded = subdiag+superdiag+1)
       aJac(k,j) = merge(aJac(k,j),0._rkind,stateMask(i)) ! zero the elements that are not included in J1
      end do
     end do
    end if
    J_filter(1:nrow_banded,1:n) = aJac(subdiag+1:nBands,1:n) ! aJac has extra storage rows
   end associate
  else ! full matrix storage
   J_total(:,:) = aJac(:,:) ! store total Jacobian (for Newton step refinement)
   J_filter(:,:) = 0._r8b   ! initialize
   if (negative) then
    do i=1,f_obj % n
     J_filter(:,i) = merge(-aJac(:,i),J_filter(:,i),stateMask(:)) ! negative sign
    end do
   else
    do i=1,f_obj % n
     J_filter(:,i) = merge(aJac(:,i),J_filter(:,i),stateMask(:))
    end do
   end if
  end if
 end subroutine filter_SUMMA_Jacobian

 subroutine f_f1_SUMMA_vec_full(f_obj,xvec) ! actually evaluates f1 depending on stateMask1 (not necessarily mass)
  ! *** Compute mass non-linear function --- use fully-coupled eval8summa call and filter results ***
  ! evaluate f anf f1 based on stateMask1
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in),contiguous :: xvec(:) ! current guess -- contiguous attribute avoids a temporary array copy

  ! local
  logical,parameter               :: mass_flag = .true.,energy_flag = .true. ! flags to compute mass and energy terms

  ! note: data structures and variables for f1 are initialized in systemSolv

  call f_obj % f_state_SUMMA_vec_full(&
               &mass_flag,energy_flag,.true.,.false.,xvec,&
               &f_obj % indx_data,f_obj % diag_data,f_obj % flux_data,f_obj % deriv_data,f_obj % sMul,&
               &f_obj % dBaseflow_dWat,f_obj % dBaseflow_dTk,f_obj % resVec)

 end subroutine f_f1_SUMMA_vec_full

 subroutine J1_SUMMA_vec_full(f_obj,xvec)
  ! *** Compute Jacobian J1 for f1 non-linear function ***
  ! NOTE: assumes appropriate eval8summa call has already been made to get the fluxes
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  !real(r8b),intent(in),contiguous :: xvec(:) ! current guess (needed for interface)
  real(r8b),intent(in) :: xvec(:) ! current guess (needed for interface)

  call f_obj % SUMMA_computJacob(f_obj % f1_mass_flag,f_obj % f1_energy_flag,&
               &f_obj % indx_data,f_obj % diag_data,f_obj % flux_data,f_obj % deriv_data,&
               &f_obj % dMat,f_obj % dBaseflow_dWat,f_obj % dBaseflow_dTk,&
               &f_obj % J1)

 end subroutine J1_SUMMA_vec_full

 subroutine f_f2_SUMMA_vec_full(f_obj,xvec)
  ! *** Compute energy non-linear function --- use fully-coupled eval8summa call and filter results ***
  ! evaluates f2 and f using stateMask2
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in),contiguous :: xvec(:) ! current guess -- contiguous attribute avoids a temporary array copy

  ! local
  logical,parameter               :: mass_flag = .true.,energy_flag = .true. ! flags to compute mass and energy terms

  ! note: data structures and variables for f1 are initialized in systemSolv

  call f_obj % f_state_SUMMA_vec_full(&
               &mass_flag,energy_flag,.false.,.true.,xvec,&
               &f_obj % indx_data,f_obj % diag_data,f_obj % flux_data,f_obj % deriv_data,f_obj % sMul,&
               &f_obj % dBaseflow_dWat,f_obj % dBaseflow_dTk,f_obj % resVec)

 end subroutine f_f2_SUMMA_vec_full

 subroutine f1_SUMMA_vec_full(f_obj,xvec)
  ! *** Compute f1 non-linear function --- use fully-coupled eval8summa call with mass and energy logical flags ***
  ! evaluates f1 using f1_mass_flag and f1_energy_flag
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in),contiguous :: xvec(:) ! current guess -- contiguous attribute avoids a temporary array copy

  call f_obj % f_state_SUMMA_vec_full(&
               &f_obj % f1_mass_flag,f_obj % f1_energy_flag,.true.,.false.,xvec,&
               &f_obj % indx_data,f_obj % diag_data,f_obj % flux_data,f_obj % deriv_data,f_obj % sMul,&
               &f_obj % dBaseflow_dWat,f_obj % dBaseflow_dTk,f_obj % resVec)

 end subroutine f1_SUMMA_vec_full

 subroutine f2_SUMMA_vec_full(f_obj,xvec)
  ! *** Compute f2 non-linear function --- use fully-coupled eval8summa call with mass and energy logical flags ***
  ! evaluates f2 using f2_mass_flag and f2_energy_flag
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in),contiguous :: xvec(:) ! current guess -- contiguous attribute avoids a temporary array copy

  call f_obj % f_state_SUMMA_vec_full(&
               &f_obj % f2_mass_flag,f_obj % f2_energy_flag,.false.,.true.,xvec,&
               &f_obj % indx_data,f_obj % diag_data,f_obj % flux_data,f_obj % deriv_data,f_obj % sMul,&
               &f_obj % dBaseflow_dWat,f_obj % dBaseflow_dTk,f_obj % resVec)

 end subroutine f2_SUMMA_vec_full

 subroutine f_f1_f2_SUMMA_vec_full(f_obj,xvec)
  ! *** Compute mass and energy non-linear functions --- use fully-coupled eval8summa call and filter results ***
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in),contiguous :: xvec(:) ! current guess -- contiguous attribute avoids a temporary array copy

  ! local
  logical,parameter               :: mass_flag = .true.,energy_flag = .true. ! flags to compute mass and energy terms

  ! note: data structures and variables for f1 are initialized in systemSolv

  call f_obj % f_state_SUMMA_vec_full(&
               &mass_flag,energy_flag,.true.,.true.,xvec,&
               &f_obj % indx_data,f_obj % diag_data,f_obj % flux_data,f_obj % deriv_data,f_obj % sMul,&
               &f_obj % dBaseflow_dWat,f_obj % dBaseflow_dTk,f_obj % resVec)

 end subroutine f_f1_f2_SUMMA_vec_full   ! solver

 subroutine J2_SUMMA_vec_full(f_obj,xvec)
  ! *** Compute Jacobian for energy non-linear function --- use fully-coupled computJacob call and filter results ***
  ! evaluate J2 and J using stateMask2
  ! NOTE: assumes appropriate eval8summa call has already been made to get the fluxes
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  !real(r8b),intent(in),contiguous :: xvec(:) ! current guess (needed for interface)
  real(r8b),intent(in) :: xvec(:) ! current guess (needed for interface)

  call f_obj % SUMMA_computJacob(f_obj % f2_mass_flag,f_obj % f2_energy_flag,&
               &f_obj % indx_data,f_obj % diag_data,f_obj % flux_data,f_obj % deriv_data,&
               &f_obj % dMat,f_obj % dBaseflow_dWat,f_obj % dBaseflow_dTk,&
               &f_obj % J2)

  call flip_sign(f_obj % nrow, f_obj % n,f_obj % J2) ! apply negative sign to J2 entries (so that J = J1 - J2)

 contains

  subroutine flip_sign(nrow,n,array)
   ! ** optimized routine to flip array signs **
   ! arguments
   integer(i4b),intent(in)            :: nrow,n
   real(r8b),intent(inout),contiguous :: array(:,:) 
   !local
   integer(i4b) :: i,j

   do j = 1,n
    do concurrent (i = 1:nrow) ! permits the compiler to optimize aggressively (e.g., simd vectorization)
     array(i,j) = -array(i,j)
    end do
   end do

  end subroutine flip_sign

 end subroutine J2_SUMMA_vec_full

! subroutine J_J1_J2_SUMMA_vec_full(f_obj,xvec) ! works if logical state masks are computed first, but no longer used
!  ! *** Compute Jacobian for mass and energy non-linear functions --- use fully-coupled computJacob call and filter results ***
!  ! NOTE: assumes appropriate eval8summa call has already been made to get the fluxes
!  ! arguments
!  class(f_obj_type),intent(inout) :: f_obj
!  real(r8b),intent(in)            :: xvec(:) ! current guess (needed for interface)
!
!  ! local
!  real(rkind)  :: aJac(f_obj % in_SS4HG % nLeadDim,f_obj % in_SS4HG % nState) ! SUMMA's unscaled Jacobian matrix
!  real(rkind)  :: aJac2(f_obj % in_SS4HG % nLeadDim,f_obj % in_SS4HG % nState) ! SUMMA's unscaled Jacobian matrix
!  integer(i4b) :: i,j,k ! loop indices
!  integer(i4b) :: nBands ! # of bands for banded storage
!
!  call f_obj % SUMMA_computJacob(.true.,.true.,&
!               &f_obj % indx_data,f_obj % diag_data,f_obj % flux_data,f_obj % deriv_data,&
!               &f_obj % dMat,f_obj % dBaseflow_dMatric,&
!               &aJac)
!
!  aJac2=aJac ! save total SUMMA Jacobian because it is filtered on output after calls to filter_SUMMA_Jacobian
!
!  ! get nested Newton solver Jacobian J1
!  call filter_SUMMA_Jacobian(f_obj,.false.,f_obj % stateMask1,aJac,f_obj % J,f_obj % J1)
!
!  ! get nested Newton solver Jacobian J2
!  call filter_SUMMA_Jacobian(f_obj,.true.,f_obj % stateMask2,aJac2,f_obj % J,f_obj % J2) ! negative sign applied
!
! end subroutine J_J1_J2_SUMMA_vec_full

 subroutine f_state_SUMMA_vec_full(f_obj,mass_flag,energy_flag,f1_flag,f2_flag,xvec,&
                                  &indx_data,diag_data,flux_data,deriv_data,sMul,&
                                  &dBaseflow_dWat,dBaseflow_dTk,resVec)
  ! *** Compute SUMMA's vector non-linear function for mass or energy state variables -- uses fully-coupled eval8summa call ***
  ! ** NOTE: the fully-coupled solution method in SUMMA's opSplittin is assumed **

  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  logical,intent(in)              :: mass_flag,energy_flag ! flags to compute mass and energy terms 
  logical,intent(in)              :: f1_flag,f2_flag ! flags to f1 and f2 
  real(r8b),intent(in),contiguous :: xvec(:) ! current guess

  type(var_ilength),intent(inout) :: indx_data            ! indices defining model states and layers for selected split 
  type(var_dlength),intent(inout) :: diag_data            ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout) :: flux_data            ! flux data
  type(var_dlength),intent(inout) :: deriv_data           ! derivative data
  !real(qp)         ,intent(inout),contiguous :: sMul(:)              ! state vector multipliers
  !real(rkind)      ,intent(out)  ,contiguous :: dBaseflow_dWat(:,:)  ! derivative in baseflow w.r.t. soil water characteristic
  !real(rkind)      ,intent(out)  ,contiguous :: dBaseflow_dTk(:,:)   ! derivative in baseflow w.r.t. temperature (m s-1 K-1)
  !real(qp)         ,intent(out)  ,contiguous :: resVec(:)            ! residual vector
  real(qp)         ,intent(inout) :: sMul(:)              ! state vector multipliers
  real(rkind)      ,intent(out)   :: dBaseflow_dWat(:,:)  ! derivative in baseflow w.r.t. soil water characteristic
  real(rkind)      ,intent(out)   :: dBaseflow_dTk(:,:)   ! derivative in baseflow w.r.t. temperature (m s-1 K-1)
  real(qp)         ,intent(out)   :: resVec(:)            ! residual vector

  ! local
  logical :: f1_mass,f1_energy,f2_mass,f2_energy

  ! determine flags for function evaluations (takes choice of decomposition and input flags for f1 and f2 evaluations into account)
  f1_mass   = f_obj % f1_mass_flag.and.f1_flag
  f1_energy = f_obj % f1_energy_flag.and.f1_flag
  f2_mass   = f_obj % f2_mass_flag.and.f2_flag
  f2_energy = f_obj % f2_energy_flag.and.f2_flag

  ! evaluate residual vector for mass split
  call eval8summa(&
                   ! input: model control
                   f_obj % in_SS4HG % dt_cur,         & ! intent(in):    current stepsize
                   f_obj % in_SS4HG % dt,             & ! intent(in):    length of the entire time step (seconds) for drainage pond rate
                   f_obj % in_SS4HG % nSnow,          & ! intent(in):    number of snow layers
                   f_obj % in_SS4HG % nSoil,          & ! intent(in):    number of soil layers
                   f_obj % in_SS4HG % nLayers,        & ! intent(in):    number of layers
                   f_obj % n,                         & ! intent(in):    number of state variables in the current subset
                   .false.,                           & ! intent(in):    not inside Sundials solver
                   f_obj % in_SS4HG % firstSubStep,   & ! intent(in):    flag to indicate if we are processing the first sub-step
                   f_obj % io_SS4HG % firstFluxCall,  & ! intent(inout): flag to indicate if we are processing the first flux call
                   .false.,                           & ! intent(in):    flag to indicate if we are processing the first flux call in a splitting operation (.false. based on usage of eval8summa in summaSolve4homegrown)
                   f_obj % in_SS4HG % computeVegFlux, & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                   f_obj % in_SS4HG % scalarSolution, & ! intent(in):    flag to indicate the scalar solution
                   mass_flag,                         & ! intent(in):    flag to compute mass terms
                   energy_flag,                       & ! intent(in):    flag to compute energy terms
                   .true.,f1_mass,f1_energy,f2_mass,f2_energy, & ! intent(in): flag to compute f, f1, and f2 for nested Newton
                   f1_mass.or.f2_mass,f1_energy.or.f2_energy,  & ! intent(in): flags to compute mass and energy Jacobian terms
                   !.true.,.true.,  & ! intent(in): flags to compute mass and energy Jacobian terms -- SJT: used to create reference output
                   ! input: state vectors
                   xvec,                            & ! intent(in):    model state vector
                   f_obj % fScale,                  & ! intent(in):    characteristic scale of the function evaluations
                   sMul,                            & ! intent(inout): state vector multiplier (used in the residual calculations)
                   ! input: data structures
                   f_obj % model_decisions,         & ! intent(in):    model decisions
                   f_obj % lookup_data,             & ! intent(in):    lookup tables
                   f_obj % type_data,               & ! intent(in):    type of vegetation and soil
                   f_obj % attr_data,               & ! intent(in):    spatial attributes
                   f_obj % mpar_data,               & ! intent(in):    model parameters
                   f_obj % forc_data,               & ! intent(in):    model forcing data
                   f_obj % bvar_data,               & ! intent(in):    average model variables for the entire basin
                   f_obj % prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                   ! input-output: data structures
                   indx_data,                       & ! intent(inout): index data
                   diag_data,                       & ! intent(inout): model diagnostic variables for a local HRU
                   flux_data,                       & ! intent(inout): model fluxes for a local HRU (initial flux structure)
                   deriv_data,                      & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                   ! input-output: baseflow
                   f_obj % io_SS4HG % ixSaturation, & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                   !dBaseflow_dMatric,               & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                   dBaseflow_dWat,                  & ! intent(out):   derivative in baseflow w.r.t. water content (s-1)
                   dBaseflow_dTk,                   & ! intent(out):   derivative in baseflow w.r.t. temperature (s-1)
                   ! output
                   f_obj % feasible,                & ! intent(out):   flag to denote the feasibility of the solution
                   f_obj % fluxVec0,                & ! intent(out):   flux vector
                   f_obj % f_vec,                   & ! intent(inout): f vector
                   f_obj % f1_vec,                  & ! intent(inout): f1 vector
                   f_obj % f2_vec,                  & ! intent(inout): f2 vector
                   f_obj % fRHS,                    & ! intent(out):   RHS function for ARKODE
                   f_obj % rAdd,                    & ! intent(out):   additional (sink) terms on the RHS of the state equation
                   resVec,                          & ! intent(out):   residual vector
                   f_obj % rVecScaled,              & ! intent(out):   scaled residual vector
                   f_obj % out_SS4HG % fNew,        & ! intent(out):   function evaluation
                   f_obj % out_SS4HG % err,         & ! intent(out): error code
                   f_obj % out_SS4HG % message)       ! intent(out): error message (note: eval8summa uses "cmessage" instead)

  ! finalize
  associate(err => f_obj % out_SS4HG % err, message => f_obj % out_SS4HG % message) 
   if (err /= 0) then
    f_obj % f_error = .true.
    if (f_obj % out_warning) then
     write(f_obj % unit,*) "Error f_state_SUMMA_vec_full: eval8summa message="//trim(message)
    end if
   end if
  end associate

 end subroutine f_state_SUMMA_vec_full

 subroutine f_SUMMA_vec(f_obj,xvec)
  ! *** Compute SUMMA's vector non-linear function ***
  class(f_obj_type),intent(inout) :: f_obj
  !real(r8b),intent(in),contiguous :: xvec(:)  ! current guess
  real(r8b),intent(in) :: xvec(:)  ! current guess

  ! compute SUMMA residual (taken to be the non-linear function) based on current guess
  ! note: - eval8summa may contain extraneous computations not needed for the residual
  !       - perhaps introducing logical flags in eval8summa to isolate the required operations would boost efficiency 
  call f_obj % SUMMA_eval8summa(xvec)

  !f_obj % f_vec(:) = real(f_obj % resVec(:),r8b) ! now directly obtained from computResid
  
 end subroutine f_SUMMA_vec

 subroutine J_SUMMA_vec(f_obj,xvec)
  ! ** Compute SUMMA's Jacobian **
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  !real(r8b),intent(in),contiguous :: xvec(:) ! current guess
  real(r8b),intent(in) :: xvec(:) ! current guess
  ! local variables
  logical,parameter :: mass_flag = .true.,energy_flag = .true.

  ! compute derivatives based on current guess
  !call f_obj % SUMMA_eval8summa(xvec) ! not required if f_SUMMA_vec(f_obj,xvec) has already been called

  ! assemble Jacobian using the computed derivatives
  call f_obj % SUMMA_computJacob(mass_flag,energy_flag,&
                                &f_obj % indx_data,f_obj % diag_data,f_obj % flux_data,f_obj % deriv_data,&
                                &f_obj % dMat,f_obj % dBaseflow_dWat,f_obj % dBaseflow_dTk,&
                                &f_obj % J)

 end subroutine J_SUMMA_vec

end module Newton_functions
