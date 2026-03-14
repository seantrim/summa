module Newton_functions
 ! nested Newton solver modules
 use kind_params,only: i4b,r8b ! kind parameters
 use Richards,only : Richards_obj ! Richards test problem
 ! SUMMA modules (for access to constant data and procedures)
 use nrtype,only: rkind,qp,lgt ! SUMMA's kind parameters (i4b is already used in kind_params module)
 use eval8summa_module, only: eval8summa,imposeConstraints           ! SUMMA's eval8summa and imposeConstraints routines
 use computJacob_module,only: computJacob                            ! SUMMA's computJacob routine 
 use summaSolve4homegrown_module,only: refine_Newton_step,checkConv ! SUMMA's refine_Newton_step and checkConv procedures
 use data_types,only: in_type_computJacob,out_type_computJacob ! objects for SUMMA's computJacob routine
 use data_types,only: in_type_summaSolve4homegrown,&           ! objects for SUMMA's summaSolve4homegrown routine
                     &io_type_summaSolve4homegrown,&
                     &out_type_summaSolve4homegrown 
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

 ! ***** Parent Type ***** !
 type, public :: f_obj_base
   ! ** Default data components used by the Newton solvers ** !
   logical      :: banded            ! flag for banded Jacobians
   logical      :: nested            ! flag for nested algorithm
   logical      :: inner             ! flag to indicate the execution of inner iterations
   logical      :: converged         ! flag to indicate that the obtained solution meets the convergence criterion
   logical      :: constraints       ! flag to indicate that constraints are to be applied between outer/classical iterations
   logical      :: constraints_inner ! flag to indicate that constraints are to be applied between inner iterations
   logical      :: refinement        ! flag to indicate that refinement is to be applied following outer/classical iterations
   logical      :: refinement_inner  ! flag to indicate that refinement is to be applied following inner iterations
   logical      :: scaling           ! flag to indicate that user-specified scaling is to be applied for linear systems
   logical      :: f_eval_flag       ! flag to indicate that the total non-linear function vector is to be computed
   logical      :: f1_eval_flag      ! flag to indicate that the non-linear function 1 vector is to be computed
   logical      :: f2_eval_flag      ! flag to indicate that the non-linear function 2 vector is to be computed
   logical      :: J_eval_flag       ! flag to indicate that the total Jacobian is to be computed
   logical      :: J1_eval_flag      ! flag to indicate that Jacobian 1 is to be computed
   logical      :: J2_eval_flag      ! flag to indicate that Jacobian 2 is to be computed
   integer(i4b) :: subdiag,superdiag ! # of subdiagonals and superdiagonals for banded Jacobians
   integer(i4b) :: n                 ! vector size
   integer(i4b) :: nrow              ! # of matrix rows (adapts to storage type)
   integer(i4b) :: nrow_banded       ! # of matrix rows for banded storage
   integer(i4b) :: kmax,lmax         ! max # of classical/outer and inner iterations
   integer(i4b) :: kcount,lcount     ! total # of classical/outer and inner iterations
   integer(i4b) :: LDA,LDAF,LDX,LDB  ! leading dimensions of A, AF, X, and B LAPACK arrays
   integer(i4b) :: KL,KU             ! # of subdiagonals and superdiagonals for LAPACK
   real(r8b),allocatable    :: WORK(:),AF(:,:)            ! LAPACK arrays
   real(r8b),allocatable    :: x0(:),x1(:)                ! initial and final root estimates for vector algorithms
   real(r8b),allocatable    :: xk(:),xkp1(:)              ! intermediate root estimates for classical iterations
   real(r8b),allocatable    :: xk0(:),xkp1l(:),xkp1lp1(:) ! intermediate root estimates for nested iterations
   real(r8b),allocatable    :: J(:,:)        ! total Jacobian
   real(r8b),allocatable    :: J1(:,:)       ! Jacobian 1
   real(r8b),allocatable    :: J2(:,:)       ! Jacobian 2
   real(r8b),allocatable    :: Jdiff(:,:)    ! difference Jacobian
   real(r8b),allocatable    :: f_vec(:)      ! total non-linear function evaluation
   real(r8b),allocatable    :: f1_vec(:)     ! non-linear function evaluation 1
   real(r8b),allocatable    :: f2_vec(:)     ! non-linear function evaluation 2
   real(r8b)                :: tol,tol_inner ! tolerance for classical/outer and inner iterations
   real(r8b)                :: R(-1:1)       ! max residual computed for iterations j-1, j, and j+1 (estimated)  
   real(r8b)                :: R_inner(-1:1) ! exact max residual computed for iterations j-1, j, and j+1 (estimated) 
   character(:),allocatable :: convergence          ! string for convergence control option for outer/classical iterations
   character(:),allocatable :: convergence_inner    ! string for convergence control option for inner iterations
   character(:),allocatable :: linear_system_solver ! string for selecting solver for linear systems
   character(:),allocatable :: matrix_vector        ! string for selecting method for matrix-vector products
   ! solver output
   character(:),allocatable :: output ! string for solver output control option
   integer(i4b) :: unit        ! file unit number for solver output
   logical      :: out_debug   ! output flag for debugging
   logical      :: out_detail  ! output flag for details
   logical      :: out_basic   ! output flag for basic information
   logical      :: out_warning ! output flag for warnings
   logical      :: out_error   ! output flag for errors
  contains
   ! procedures used prior to calling the solver
   procedure :: set_defaults    => f_set_defaults    ! set default options 
   procedure :: allocate_memory => f_allocate_memory ! allocate array data components 
   procedure :: initial_guess   => f_initial_guess   ! apply initial guess strategy
   procedure :: set_tolerance   => f_set_tolerance   ! set tolerances and iteration count maximums
   procedure :: solver_output   => f_solver_output   ! set solver output options
 end type f_obj_base

 type,extends(f_obj_base),public :: f_obj_inputs
   ! * SUMMA data *
   type(model_options),allocatable :: model_decisions(:) ! model decisions

   type(zLookup)     :: lookup_data                  ! lookup tables
   type(var_dlength) :: flux_init                    ! model fluxes at the start of the time step
   type(var_i)       :: type_data                    ! type of vegetation and soil
   type(var_d)       :: attr_data                    ! spatial attributes
   type(var_d)       :: forc_data                    ! model forcing data
   type(var_dlength) :: mpar_data                    ! model parameters
   type(var_dlength) :: bvar_data                    ! model variables for the local basin


   type(var_ilength) :: indx_data,indx_data1,indx_data2    ! indices defining model states and layers
   type(var_dlength) :: prog_data,prog_temp                ! prognostic variables for a local HRU
   type(var_dlength) :: diag_data,diag_data1,diag_data2    ! diagnostic variables for a local HRU
   type(var_dlength) :: flux_data,flux_data1,flux_data2    ! temporary flux variables for a local HRU
   type(var_dlength) :: deriv_data,deriv_data1,deriv_data2 ! derivatives in model fluxes w.r.t. relevant state variables
   real(rkind),allocatable :: dBaseflow_dMatric(:,:),dBaseflow_dMatric1(:,:),dBaseflow_dMatric2(:,:) ! derivative in baseflow w.r.t. matric head (s-1)
   real(rkind),allocatable :: dMat(:),dMat1(:),dMat2(:)    ! diagonal matrix (excludes flux derivatives) 

   ! * summaSolve4homegrown (SS4HG) objects *
   ! classical and outer iterations
   type(in_type_summaSolve4homegrown)  :: in_SS4HG   ! SS4HG input object: model control variables and previous function evaluation
   type(io_type_summaSolve4homegrown)  :: io_SS4HG   ! SS4HG io object: model control variables and previous function evaluation
   type(out_type_summaSolve4homegrown) :: out_SS4HG  ! SS4HG output object: model control variables and previous function evaluation

   ! inner iterations
   type(in_type_summaSolve4homegrown)  :: in_SS4HG_inner  ! SS4HG input object: model control variables and previous function evaluation
   type(io_type_summaSolve4homegrown)  :: io_SS4HG_inner  ! SS4HG io object: model control variables and previous function evaluation
   type(out_type_summaSolve4homegrown) :: out_SS4HG_inner ! SS4HG output object: model control variables and previous function evaluation

   ! additional variables for eval8summa call
   logical(lgt)            :: firstSplitOper         ! flag to indicate if we are processing the first flux call in a splitting operation
   real(rkind),allocatable :: fScale(:)              ! characteristic scale of the function evaluations (mixed units)
   real(rkind),allocatable :: xScale(:)              ! characteristic scale of the state vector (mixed units)
   real(qp),allocatable    :: sMul(:),sMul1(:),sMul2(:) ! NOTE: qp  ! multiplier for state vector for the residual calculations
   logical(lgt) :: feasible                          ! feasibility flag
   real(rkind),allocatable :: fluxVec0(:)            ! flux vector (mixed units)
   real(rkind),allocatable :: fRHS(:)                ! RHS function for ARKODE
   real(rkind),allocatable :: rAdd(:)                ! additional terms in the residual vector
   real(qp),allocatable    :: resVec(:)  ! NOTE: qp  ! residual vector 

   ! scaled arrays
   real(rkind),allocatable :: rVecScaled(:)   ! scaled residual
   real(rkind),allocatable :: aJacScaled(:,:) ! scaled Jacobian

   ! variables to handle state type non-linear function decompositions
   integer(i4b)             :: nLeadDim1,nLeadDim2
   integer(i4b)             :: nSubset1,nSubset2
   logical(lgt),allocatable :: stateMask1(:),stateMask2(:)  

  contains
   ! ** routines that point to external sources ** !
   ! note: - these procedures are not directly called in the solver
   !       - however, these procedures may be called within procedures that are called in the solver

   ! * Interfaces for SUMMA procedures  *
   procedure :: SUMMA_eval8summa
   procedure :: SUMMA_computJacob

!   ! f=space minus time
!   ! scalar input routines
!   procedure :: f1 => f_Rich_space
!   procedure :: f2 => f_Rich_time
!   procedure :: df1dx_element => df_Rich_dh_element_space
!   procedure :: df2dx_element => df_Rich_dh_element_time

   ! f=time minus space
   ! scalar input routines
   procedure :: f1 => f_Rich_time
   procedure :: f2 => f_Rich_space
   procedure :: df1dx_element => df_Rich_dh_element_time
   procedure :: df2dx_element => df_Rich_dh_element_space
 end type f_obj_inputs

 type,extends(f_obj_inputs),public :: f_obj_type
  contains
   ! *** these procedures take the procedures from f_obj_inputs type as input *** !
   ! vector routines
   procedure :: f_vec_eval  => f_SUMMA_vec  ! solver
   !procedure :: f_vec_eval => f_diff_vec   ! solver (note: f_diff_vec requires nested iterations to be activated)
   !procedure :: f1_vec_eval => f1_SUMMA_vec ! solver -- trivial split
   !procedure :: f2_vec_eval => f2_zero_vec  ! solver -- trivial split
   procedure :: f1_vec_eval => f_mass_SUMMA_vec_full   ! solver
   procedure :: f2_vec_eval => f_energy_SUMMA_vec_full ! solver
   procedure :: dfdx_vec  => dfdx_diff_vec 
   procedure :: df1dx_vec => df1_Rich_dh_vec
   procedure :: df2dx_vec => df2_Rich_dh_vec
   procedure :: J_eval  => Jacobian_f_SUMMA_vec  ! solver
   !procedure :: J1_eval => Jacobian_f1_SUMMA_vec ! solver -- trivial split
   !procedure :: J2_eval => Jacobian_f2_zero_vec  ! solver -- trivial split
   procedure :: J1_eval => Jacobian_f_mass_SUMMA_vec_full   ! solver
   procedure :: J2_eval => Jacobian_f_energy_SUMMA_vec_full ! solver
   procedure :: J1_J2_eval => Jacobian_f_mass_energy_SUMMA_vec_full ! solver
   procedure :: apply_constraints  => SUMMA_imposeConstraints
   procedure :: apply_refinement_classical   => SUMMA_refine_Newton_step_classical
   procedure :: apply_refinement_inner       => SUMMA_refine_Newton_step_inner
   procedure :: apply_refinement_outer       => SUMMA_refine_Newton_step_outer
   procedure :: apply_nested_line_search     => SUMMA_nested_line_search
   procedure :: line_search_objective => SUMMA_line_search_objective
   procedure :: custom_convergence => SUMMA_check_convergence_flag !SUMMA_checkConv  
   procedure :: custom_scaling     => SUMMA_scaling  
   procedure :: custom_descaling   => SUMMA_descaling  
   procedure :: Jacobian_f1_SUMMA_vec_numerical ! SJT: testing ----- take out -----
   procedure :: Jacobian_f2_SUMMA_vec_numerical ! SJT: testing ----- take out -----
   procedure :: Jacobian_f_SUMMA_vec_numerical  ! SJT: testing ----- take out -----
   procedure :: get_mass_energy_masks => get_SUMMA_mass_energy_masks
   procedure :: f_state_SUMMA_vec_full
   procedure :: f_mass_SUMMA_vec_full
   procedure :: f_energy_SUMMA_vec_full
   procedure :: Jacobian_f_mass_SUMMA_vec_full
   procedure :: Jacobian_f_energy_SUMMA_vec_full

   ! scalar routines
   procedure :: f     => f_diff 
   procedure :: dfdx  => dfdx_diff 
   procedure :: df1dx => df1_Rich_dh 
   procedure :: df2dx => df2_Rich_dh 
 end type f_obj_type

contains

!!!!!!!!!! ******************* User defined functions below ******************* !!!!!!!!!!

 ! **** Utilities **** !

 subroutine f_set_defaults(f_obj)
  ! ** set default values for options in f_obj_base class **
  use, intrinsic :: iso_fortran_env, only: stdout=>output_unit ! for default output
  class(f_obj_base),intent(inout) :: f_obj

   f_obj % banded            = .false. ! flag for banded Jacobians
   f_obj % nested            = .false. ! flag for nested algorithm
   f_obj % constraints       = .false. ! flag to indicate that constraints are to be applied between outer/classical iterations
   f_obj % constraints_inner = .false. ! flag to indicate that constraints are to be applied between inner iterations
   f_obj % refinement        = .false. ! flag to indicate that refinement is to be applied following outer/classical iterations
   f_obj % refinement_inner  = .false. ! flag to indicate that refinement is to be applied following inner iterations
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

   f_obj % linear_system_solver = "LAPACK_expert"          ! string for control of linear system solver
   f_obj % matrix_vector        = "BLAS"                   ! string for selecting matrix-vector product method
   f_obj % convergence          = "strict"                 ! string for convergence criterion method for solver
   f_obj % convergence_inner    = "strict"                 ! string for convergence criterion method for solver
   f_obj % output               = "production"             ! string for solver output control option
   
   f_obj % unit                 = stdout                   ! file unit number for solver output

 end subroutine f_set_defaults
 
 subroutine f_allocate_memory(f_obj)
  ! ** allocate array data components for f_obj_base class **
  class(f_obj_base),intent(inout) :: f_obj

  ! allocate solution and function arrays
  associate(n => f_obj % n)
   allocate(f_obj % x0(1:n),f_obj % x1(1:n))          ! initial and final root estimates
   allocate(f_obj % f_vec(1:n))                       ! total non-linear function vector
   if (f_obj % nested) then
    allocate(f_obj % xk0(1:n),f_obj % xkp1l(1:n),f_obj % xkp1lp1(1:n)) ! intermediate root estimates for nested iterations
    allocate(f_obj % f1_vec(1:n),f_obj % f2_vec(1:n))                  ! non-linear functions vectors 1 and 2 
   else
    allocate(f_obj % xk(1:n),f_obj % xkp1(1:n))                        ! intermediate root estimates for classical iterations
   end if
  end associate

  ! * allocate LAPACK arrays *

  ! LAPACK parameters independent of matrix storage type
  f_obj % LDX=f_obj % n; f_obj % LDB=f_obj % n ! leading dimensions for RHS arrays

  ! allocate memory and set LAPACK parameters for choice of matrix storage
  if (f_obj % banded) then ! banded storage
   f_obj % KL = f_obj % subdiag; f_obj % KU = f_obj % superdiag
   f_obj % LDA = f_obj % KL + f_obj % KU + 1_i4b; f_obj % LDAF = f_obj % LDA + f_obj % KL
   allocate(f_obj % AF(1:f_obj % LDAF,1:f_obj % n)) ! storing LU factors requires an additional f_obj % subdiag rows
   if (f_obj % linear_system_solver .eq. "LAPACK_expert") allocate(f_obj % WORK(1:3_i4b*f_obj % n))
  else ! full matrix storage
   f_obj % LDA = f_obj % n; f_obj % LDAF = f_obj % n
   allocate(f_obj % AF(1:f_obj % n,1:f_obj % n))
   if (f_obj % linear_system_solver .eq. "LAPACK_expert") allocate(f_obj % WORK(1:4_i4b*f_obj % n))
  end if

  ! allocate Jacobian arrays (and initialize to zero)
  if (f_obj % banded) then ! banded storage
    f_obj % nrow_banded = f_obj % subdiag + f_obj % superdiag + 1_i4b
    f_obj % nrow = f_obj % nrow_banded
  else
    f_obj % nrow_banded = f_obj % n
    f_obj % nrow        = f_obj % n
  end if
  
  allocate(f_obj % J(1:f_obj % nrow,1:f_obj % n),source=0._r8b) ! SJT: available for classical and nested iterations (for testing -- not needed for nested)
  if (f_obj % nested) then
   allocate(f_obj % J1(1:f_obj % nrow,1:f_obj % n),f_obj % J2(1:f_obj % nrow,1:f_obj % n),&
           &f_obj % Jdiff(1:f_obj % nrow,1:f_obj % n),source=0._r8b)
  end if

 end subroutine f_allocate_memory

 subroutine f_solver_output(f_obj,method,unit)
  ! ** set output control for solver **
  use, intrinsic :: iso_fortran_env, only: stdout=>output_unit ! for default output
  class(f_obj_base),intent(inout)  :: f_obj
  character(*),intent(in)          :: method
  integer(i4b),optional,intent(in) :: unit

  ! set file unit for solver output - default is standard output
  if (present(unit)) then
   f_obj % unit = unit 
  else
   f_obj % unit = stdout 
  end if  
 
  if (method.eq.'debug') then
   f_obj % out_debug   = .true. ! output flag for debugging
   f_obj % out_detail  = .true. ! output flag for details
   f_obj % out_basic   = .true. ! output flag for basic information
   f_obj % out_warning = .true. ! output flag for warnings
   f_obj % out_error   = .true. ! output flag for errors
  else if (method.eq.'verbose') then
   f_obj % out_debug   = .false. ! output flag for debugging
   f_obj % out_detail  = .true.  ! output flag for details
   f_obj % out_basic   = .true.  ! output flag for basic information
   f_obj % out_warning = .true.  ! output flag for warnings
   f_obj % out_error   = .true.  ! output flag for errors
  else if (method.eq.'production') then
   f_obj % out_debug   = .false. ! output flag for debugging
   f_obj % out_detail  = .false. ! output flag for details
   f_obj % out_basic   = .true.  ! output flag for basic information
   f_obj % out_warning = .true.  ! output flag for warnings
   f_obj % out_error   = .true.  ! output flag for errors
  else if (method.eq.'minimal') then
   f_obj % out_debug   = .false. ! output flag for debugging
   f_obj % out_detail  = .false. ! output flag for details
   f_obj % out_basic   = .true.  ! output flag for basic information
   f_obj % out_warning = .false. ! output flag for warnings
   f_obj % out_error   = .true.  ! output flag for errors
  else if (method.eq.'silent') then
   f_obj % out_debug   = .false. ! output flag for debugging
   f_obj % out_detail  = .false. ! output flag for details
   f_obj % out_basic   = .false. ! output flag for basic information
   f_obj % out_warning = .false. ! output flag for warnings
   f_obj % out_error   = .true.  ! output flag for errors
  else
   if (f_obj % out_error) then
    write(f_obj % unit,'(a65)') "Error in f_solver_output: method argument not currently supported"
   end if
   stop
  end if
 end subroutine f_solver_output

 ! **** Numerics **** !

 subroutine f_set_tolerance(f_obj,method,tol,kmax)
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

 subroutine f_initial_guess(f_obj,method)
  ! ** initial guess strategy for time-dependent algorithms for f_obj_base class **
  ! note: it may be possible to add filtering techniques for the initial guess to improve efficiency
  class(f_obj_base),intent(inout) :: f_obj
  character(*),intent(in)         :: method

  ! Note: avoid unintentional reallocation of object components (use array slices for assignment statements)
  if (method.eq.'previous') then 
   f_obj % x0 = f_obj % x1 ! initial guess -- solution from previous time step
  else
   if (f_obj % out_error) then
    write(f_obj % unit,'(a66)') "Error in f_initial_guess: method argument not currently supported."
   end if
   stop
  end if
 end subroutine f_initial_guess

 
 ! **** Richards Problem **** !

 real(r8b) function f_Rich_space(f_obj,x) result(f_space)
  ! ** space terms for discrete Richards' equation **
  class(f_obj_inputs),intent(in) :: f_obj
  real(r8b),intent(in)         :: x ! current guess

  associate(i => Richards_obj % i)
    Richards_obj % h(i) = x ! populate h array with current guess 
    f_space=Richards_obj % KT()-Richards_obj % S()
  end associate
 end function f_Rich_space

 real(r8b) function f_Rich_time(f_obj,x) result(f_time)
  ! ** time terms for discrete Richards' equation **
  class(f_obj_inputs),intent(in) :: f_obj
  real(r8b),intent(in)         :: x ! current guess

  associate(i => Richards_obj % i)
    Richards_obj % h(i) = x ! populate h array with current guess 
    f_time=Richards_obj % CT()
  end associate
 end function f_Rich_time

 real(r8b) function df_Rich_dh_element_space(f_obj,x,j) result(dfdh_element_space)
  class(f_obj_inputs),intent(in) :: f_obj
  real(r8b),intent(in)         :: x ! current guess
  integer(i4b),intent(in)      :: j

  associate(i => Richards_obj % i)
    Richards_obj % h(i) = x ! populate h array with current guess 
    dfdh_element_space=Richards_obj % dKTdh(j)
  end associate
 end function df_Rich_dh_element_space

 real(r8b) function df_Rich_dh_element_time(f_obj,x,j) result(dfdh_element_time)
  class(f_obj_inputs),intent(in) :: f_obj
  real(r8b),intent(in)         :: x ! current guess
  integer(i4b),intent(in)      :: j

  associate(i => Richards_obj % i)
    Richards_obj % h(i) = x ! populate h array with current guess 
    dfdh_element_time=Richards_obj % dCTdh(j)
  end associate
 end function df_Rich_dh_element_time

!!!!!!!!!!!!!!! ****************** Functions that adapt to specified scalar functions below ****************** !!!!!!!!!!!!!!! 

 function Jacobian_f_Rich_vec(f_obj,xvec) result(J)
  class(f_obj_type),intent(in) :: f_obj
  real(r8b),intent(in)         :: xvec(1:f_obj % n) ! current guess
  real(r8b),allocatable        :: J(:,:)
  !real(r8b)                    :: J(1:f_obj % n,1:f_obj % n) ! original Jacobian (full matrix storage)
  integer(i4b)                 :: icol,irow
  integer(i4b)                 :: nrow_banded ! # of rows for LAPACK banded matrix storage

  if (f_obj % banded) then ! banded storage
   associate(n => f_obj % n, subdiag => f_obj % subdiag, superdiag => f_obj % superdiag)
    nrow_banded=subdiag+superdiag+1
    allocate(J(1:nrow_banded,1:n))
    do icol=1,n
     do irow=max(1,icol-superdiag),min(n,icol+subdiag)
      Richards_obj % i = irow
      J(superdiag+1+irow-icol,icol)=f_obj % dfdx_vec(xvec,icol)
     end do
    end do 
   end associate
  else ! full matrix storage
   associate(n => f_obj % n)
    allocate(J(1:n,1:n))
    do icol=1,n
     do irow=1,n
      Richards_obj % i = irow
      J(irow,icol)=f_obj % dfdx_vec(xvec,icol)
     end do
    end do 
   end associate
  end if
 end function Jacobian_f_Rich_vec

 subroutine Jacobian_f1_Rich_vec(f_obj,xvec)
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)            :: xvec(1:f_obj % n) ! current guess
  integer(i4b)                    :: icol,irow

  if (f_obj % banded) then ! banded storage
   associate(n => f_obj % n, subdiag => f_obj % subdiag, superdiag => f_obj % superdiag)
    do icol=1,n
     do irow=max(1,icol-superdiag),min(n,icol+subdiag)
      Richards_obj % i = irow
      f_obj % J1(superdiag+1+irow-icol,icol)=f_obj % df1dx_vec(xvec,icol)
     end do
    end do 
   end associate
  else ! full matrix storage
   associate(n => f_obj % n)
    do icol=1,n
     do irow=1,n
      Richards_obj % i = irow
      f_obj % J1(irow,icol)=f_obj % df1dx_vec(xvec,icol)
     end do
    end do 
   end associate
  end if
 end subroutine Jacobian_f1_Rich_vec

 subroutine Jacobian_f2_Rich_vec(f_obj,xvec)
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)            :: xvec(1:f_obj % n) ! current guess
  integer(i4b)                    :: icol,irow

  if (f_obj % banded) then ! banded storage
   associate(n => f_obj % n, subdiag => f_obj % subdiag, superdiag => f_obj % superdiag)
    do icol=1,n
     do irow=max(1,icol-superdiag),min(n,icol+subdiag)
      Richards_obj % i = irow
      f_obj % J2(superdiag+1+irow-icol,icol)=f_obj % df2dx_vec(xvec,icol)
     end do
    end do 
   end associate
  else ! full matrix storage
   associate(n => f_obj % n)
    do icol=1,n
     do irow=1,n
      Richards_obj % i = irow
      f_obj % J2(irow,icol)=f_obj % df2dx_vec(xvec,icol)
     end do
    end do 
   end associate
  end if
 end subroutine Jacobian_f2_Rich_vec

 subroutine f_diff_vec(f_obj,xvec)
  ! *** form non-linear vector function using the decomposition ***
  ! note: f1_vec and f2_vec components currently only allocated for nested iterations
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)            :: xvec(1:f_obj % n)  ! current guess

  if (f_obj % f1_eval_flag) call f_obj % f1_vec_eval(xvec)
  if (f_obj % f2_eval_flag) call f_obj % f2_vec_eval(xvec)

  f_obj % f_vec = f_obj % f1_vec(:) - f_obj % f2_vec(:)
 end subroutine f_diff_vec

 subroutine f1_Rich_vec(f_obj,xvec)
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)            :: xvec(1:f_obj % n) ! current guess
  integer(i4b)                    :: i

  associate(n => f_obj % n)
   do i=1,n ! interior grid points
    Richards_obj % i = i
    ! populate h array with current guess on ith stencil
    if (i.ne.1) Richards_obj % h(i-1) = xvec(i-1) ! BC 
                Richards_obj % h(i) = xvec(i) 
    if (i.ne.n) Richards_obj % h(i+1) = xvec(i+1) ! BC 
    f_obj % f1_vec(i) = f_obj % f1(xvec(i))
   end do
  end associate
 end subroutine f1_Rich_vec

 subroutine f2_Rich_vec(f_obj,xvec)
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)            :: xvec(1:f_obj % n) ! current guess
  integer(i4b)                    :: i

  associate(n => f_obj % n)
   do i=1,n ! interior grid points
    Richards_obj % i = i
    ! populate h array with current guess on ith stencil
    if (i.ne.1) Richards_obj % h(i-1) = xvec(i-1) ! BC 
                Richards_obj % h(i) = xvec(i) 
    if (i.ne.n) Richards_obj % h(i+1) = xvec(i+1) ! BC 
    f_obj % f2_vec(i) = f_obj % f2(xvec(i))
   end do
  end associate
 end subroutine f2_Rich_vec

 real(r8b) function f_diff(f_obj,x) result(f)
  ! ** complete scalar non-linear function from Jordan decomposition **
  class(f_obj_type),intent(in) :: f_obj
  real(r8b),intent(in)         :: x ! current guess

  f=f_obj % f1(x)-f_obj % f2(x)
 end function f_diff

 real(r8b) function dfdx_diff(f_obj,x) result(dfdx)
  class(f_obj_type),intent(in) :: f_obj
  real(r8b),intent(in)         :: x ! current guess

  dfdx=f_obj % df1dx(x)-f_obj % df2dx(x)
 end function dfdx_diff

 real(r8b) function dfdx_diff_vec(f_obj,x,j) result(dfdx_vec)
  class(f_obj_type),intent(in) :: f_obj
  real(r8b),intent(in)         :: x(1:f_obj % n) ! current guess
  integer(i4b),intent(in)      :: j   
  dfdx_vec=f_obj % df1dx_vec(x,j)-f_obj % df2dx_vec(x,j)
 end function dfdx_diff_vec

 real(r8b) function df1_Rich_dh_vec(f_obj,x,j) result(df1dh)
  class(f_obj_type),intent(in) :: f_obj
  real(r8b),intent(in)         :: x(1:f_obj % n) ! current guess
  integer(i4b),intent(in)      :: j   

  associate(i => Richards_obj % i)
   ! populate h array with current guess
   if (i.ne.1)                 Richards_obj % h(i-1) = x(i-1)
                               Richards_obj % h(i) = x(i)     
   if (i.ne.Richards_obj % nz) Richards_obj % h(i+1) = x(i+1)
   ! evaluate derivative of f WRT x(j)
   df1dh=f_obj % df1dx_element(x(i),j)
  end associate
 end function df1_Rich_dh_vec

 real(r8b) function df2_Rich_dh_vec(f_obj,x,j) result(df2dh)
  class(f_obj_type),intent(in) :: f_obj
  real(r8b),intent(in)         :: x(1:f_obj % n) ! current guess
  integer(i4b),intent(in)      :: j   

  associate(i => Richards_obj % i)
   ! populate h array with current guess
   if (i.ne.1)                 Richards_obj % h(i-1) = x(i-1)
                               Richards_obj % h(i) = x(i)     
   if (i.ne.Richards_obj % nz) Richards_obj % h(i+1) = x(i+1)
   ! evaluate derivative of f WRT x(j)
   df2dh=f_obj % df2dx_element(x(i),j)
  end associate
 end function df2_Rich_dh_vec

 real(r8b) function df1_Rich_dh(f_obj,x) result(df1dh)
  class(f_obj_type),intent(in) :: f_obj
  real(r8b),intent(in)         :: x ! current guess

  associate(i => Richards_obj % i)
    df1dh=f_obj % df1dx_element(x,i) ! index j = i for scalar problem
  end associate
 end function df1_Rich_dh

 real(r8b) function df2_Rich_dh(f_obj,x) result(df2dh)
  class(f_obj_type),intent(in) :: f_obj
  real(r8b),intent(in)         :: x ! current guess

  associate(i => Richards_obj % i)
    df2dh=f_obj % df2dx_element(x,i) ! index j = i for scalar problem
  end associate
 end function df2_Rich_dh


 !! ******************************* SUMMA procedures below ******************************* !!
 subroutine SUMMA_refine_Newton_step_classical(f_obj,J,xvec0,xStep,xvec1)
  ! ** interface to SUMMA's refine_Newton_step subroutine **
  use matrixOper_module,  only: scaleMatrices
  ! object
  class(f_obj_type),intent(inout)   :: f_obj
  ! input
  real(r8b),intent(in)              :: J(1:f_obj % nrow,1:f_obj % n) ! nested Newton solver Jacobian matrix
  real(r8b),intent(in)              :: xvec0(1:f_obj % n) ! previous guess vector
  real(r8b),intent(in)              :: xStep(1:f_obj % n) ! unrefined Newton step
  ! input-output
  real(r8b),intent(inout)           :: xvec1(1:f_obj % n) ! current guess vector

  call SUMMA_refine_Newton_step(f_obj,J,xvec0,xStep,xvec1)

  ! store non-linear function vector for next Newton iteration
  f_obj % f_vec(:) = real(f_obj % resVec(:),r8b)

  ! update function value for line search
  f_obj % in_SS4HG % fOld = f_obj % out_SS4HG % fNew

 end subroutine SUMMA_refine_Newton_step_classical

 subroutine SUMMA_refine_Newton_step_inner(f_obj,J,xvec0,xStep,xvec1)
  ! ** interface to SUMMA's refine_Newton_step subroutine **
  use matrixOper_module,  only: scaleMatrices
  ! object
  class(f_obj_type),intent(inout)   :: f_obj
  ! input
  real(r8b),intent(in)              :: J(1:f_obj % nrow,1:f_obj % n) ! nested Newton solver Jacobian matrix
  real(r8b),intent(in)              :: xvec0(1:f_obj % n) ! previous guess vector
  real(r8b),intent(in)              :: xStep(1:f_obj % n) ! unrefined Newton step
  ! input-output
  real(r8b),intent(inout)           :: xvec1(1:f_obj % n) ! current guess vector
  ! local
  logical,parameter :: trivial_decomposition = .false.

  ! update function value for line search (based on previous evaluation of f1, f2, or f)
  f_obj % in_SS4HG % fOld = f_obj % out_SS4HG % fNew

  call SUMMA_refine_Newton_step(f_obj,J,xvec0,xStep,xvec1)

  ! store non-linear function vector for next Newton iteration
  if (trivial_decomposition) then
   f_obj % f1_vec(:) = real(f_obj % resVec(:),r8b) ! trivial decomposition (f2=0)
  end if

 end subroutine SUMMA_refine_Newton_step_inner

 subroutine SUMMA_refine_Newton_step_outer(f_obj,J,xvec0,xStep,xvec1)
  ! ** interface to SUMMA's refine_Newton_step subroutine **
  use matrixOper_module,  only: scaleMatrices
  ! object
  class(f_obj_type),intent(inout)   :: f_obj
  ! input
  real(r8b),intent(in)              :: J(1:f_obj % nrow,1:f_obj % n) ! nested Newton solver Jacobian matrix
  real(r8b),intent(in)              :: xvec0(1:f_obj % n) ! previous guess vector
  real(r8b),intent(in)              :: xStep(1:f_obj % n) ! unrefined Newton step
  ! input-output
  real(r8b),intent(inout)           :: xvec1(1:f_obj % n) ! current guess vector
  ! local
  logical,parameter :: trivial_decomposition = .false.

  ! update function value for line search (based on previous evaluation of f1, f2, or f)
  f_obj % in_SS4HG % fOld = f_obj % out_SS4HG % fNew

  call SUMMA_refine_Newton_step(f_obj,J,xvec0,xStep,xvec1)

  ! store non-linear function vector for next Newton iteration
  if (trivial_decomposition) then
   f_obj % f2_vec(:) = real(f_obj % resVec(:),r8b) ! trivial decomposition (f1=0)
  else
   if (.not.f_obj % f1_eval_flag) call f_obj % f1_vec_eval(xvec1) ! non-trivial decomposition
   if (.not.f_obj % f2_eval_flag) call f_obj % f2_vec_eval(xvec1) ! non-trivial decomposition
  end if

 end subroutine SUMMA_refine_Newton_step_outer

 subroutine SUMMA_refine_Newton_step(f_obj,J,xvec0,xStep,xvec1)
  ! ** interface to SUMMA's refine_Newton_step subroutine **
  use matrixOper_module,  only: scaleMatrices
  ! object
  class(f_obj_type),intent(inout)   :: f_obj
  ! input
  real(r8b),intent(in)              :: J(1:f_obj % nrow,1:f_obj % n) ! nested Newton solver Jacobian matrix
  real(r8b),intent(in)              :: xvec0(1:f_obj % n) ! previous guess vector
  real(r8b),intent(in)              :: xStep(1:f_obj % n) ! unrefined Newton step
  ! input-output
  real(r8b),intent(inout)           :: xvec1(1:f_obj % n) ! current guess vector
  ! local
  integer(i4b) :: nBands ! SUMMA's leading dimension for banded Jacobians
  integer(i4b) :: mSoil  ! number of soil layers in the solution vector
  real(rkind)  :: aJac(f_obj % in_SS4HG % nLeadDim,f_obj % in_SS4HG % nState) ! SUMMA Jacobian matrix
  real(rkind)  :: newtStepScaled(1:f_obj % in_SS4HG % nState)                 ! full newton step (scaled)
  real(rkind)  :: stateVecTrial(1:f_obj % in_SS4HG % nState)                  ! unrefined guess
  real(rkind)  :: stateVecNew(1:f_obj % in_SS4HG % nState)                    ! refined guess
  logical(lgt)   :: return_flag
  integer(i4b)   :: err
  character(256) :: cmessage
  
  if ((f_obj % scaling).and.(.not.f_obj % nested)) then ! reuse scaled arrays already computed

   newtStepScaled = xStep(:) ! assume xStep is already scaled on input

  else ! scale arrays

   ! * Get Scaled Jacobians *
   if (f_obj % nested) then ! nested iterations

    ! get SUMMA Jacobian from solver Jacobian for non-linear function f1
    if (f_obj % banded) then ! banded storage
     associate(nrow_banded => f_obj % nrow_banded, n => f_obj % n, subdiag => f_obj % subdiag)
      nBands=nrow_banded+subdiag
      aJac(1:subdiag,1:n) = 0._rkind
      aJac(subdiag+1:nBands,1:n) = f_obj % J1(1:nrow_banded,1:n) ! SUMMA's aJac has extra storage rows
     end associate
    else ! full matrix storage
     aJac(:,:) = f_obj % J1(:,:)
    end if

    ! Scale Jacobian for f1
    associate(ixMatrix => f_obj % in_SS4HG % ixMatrix, nState => f_obj % in_SS4HG % nState)
     call scaleMatrices(ixMatrix,nState,aJac,f_obj % fScale,f_obj % xScale,f_obj % in_SS4HG % aJac1Scaled,err,cmessage) ! matches solve_linear_system
    end associate
    if (err/=0) then
     if (f_obj % out_error) then
      write(f_obj % unit,*) "Error in SUMMA_refine_Newton_step: scaleMatrices message="//trim(cmessage); stop
     end if
    end if

    ! get SUMMA Jacobian from solver Jacobian for non-linear function f2
    if (f_obj % banded) then ! banded storage
     associate(nrow_banded => f_obj % nrow_banded, n => f_obj % n, subdiag => f_obj % subdiag)
      nBands=nrow_banded+subdiag
      aJac(1:subdiag,1:n) = 0._rkind
      aJac(subdiag+1:nBands,1:n) = f_obj % J2(1:nrow_banded,1:n) ! SUMMA's aJac has extra storage rows
     end associate
    else ! full matrix storage
     aJac(:,:) = f_obj % J2(:,:)
    end if

    ! Scale Jacobian for f2
    associate(ixMatrix => f_obj % in_SS4HG % ixMatrix, nState => f_obj % in_SS4HG % nState)
     call scaleMatrices(ixMatrix,nState,aJac,f_obj % fScale,f_obj % xScale,f_obj % in_SS4HG % aJac2Scaled,err,cmessage) ! matches solve_linear_system
    end associate
    if (err/=0) then
     if (f_obj % out_error) then
      write(f_obj % unit,*) "Error in SUMMA_refine_Newton_step: scaleMatrices message="//trim(cmessage); stop
     end if
    end if

   else ! classical iterations

    ! get SUMMA Jacobian from solver Jacobian for total non-linear function f
    if (f_obj % banded) then ! banded storage
     associate(nrow_banded => f_obj % nrow_banded, n => f_obj % n, subdiag => f_obj % subdiag)
      nBands=nrow_banded+subdiag
      aJac(1:subdiag,1:n) = 0._rkind
      aJac(subdiag+1:nBands,1:n) = J(1:nrow_banded,1:n) ! SUMMA's aJac has extra storage rows
     end associate
    else ! full matrix storage
     aJac(:,:) = J(:,:)
    end if

    ! scale Jacobian for f 
    associate(ixMatrix => f_obj % in_SS4HG % ixMatrix, nState => f_obj % in_SS4HG % nState)
     call scaleMatrices(ixMatrix,nState,aJac,f_obj % fScale,f_obj % xScale,f_obj % aJacScaled,err,cmessage) ! matches solve_linear_system
    end associate
    if (err/=0) then
     if (f_obj % out_error) then
      write(f_obj % unit,*) "Error in SUMMA_refine_Newton_step: scaleMatrices message="//trim(cmessage); stop
     end if
    end if

   end if

   ! * Get Other Scaled Variables *
   ! note: need to match scaling applied in solve_linear_system subroutine in summaSolve4homegrown
   if (f_obj % nested) then ! if computing the Newton step
    if (f_obj % scaling) then
     ! assume xvec1 is already scaled on input (e.g., from LAPACK solution for nested iterations)
     newtStepScaled = xvec1(:) - xvec0(:)/f_obj % xScale(:) 
    else 
     newtStepScaled(:) = (xvec1(:) - xvec0(:)) / f_obj % xScale(:) ! get scaled Newton step (consistent with scaling for aJacScaled and rVecScaled)
    end if
    f_obj % in_SS4HG % xk0 = f_obj % xk0
    f_obj % in_SS4HG % xkp1l = f_obj % xkp1l
    f_obj % in_SS4HG % xkp1lp1 = f_obj % xkp1lp1
    f_obj % io_SS4HG % f_vec(:)  = f_obj % f_vec(:)  ! initialize function values for line search
    f_obj % io_SS4HG % f1_vec(:) = f_obj % f1_vec(:)
    f_obj % io_SS4HG % f2_vec(:) = f_obj % f2_vec(:)
    !f_obj % rVecScaled(:) = f_obj % fScale(:) * f_obj % f_vec(:) ! step refinement uses updated total f
    f_obj % rVecScaled(:) = f_obj % fScale(:) * (f_obj % f1_vec(:)-f_obj % f2_vec(:)) ! uses constant f2 during inner iterations
    f_obj % in_SS4HG % fOld = 0.5_rkind*dot_product(f_obj % rVecScaled,f_obj % rVecScaled)
    f_obj % in_SS4HG % inner = f_obj % inner ! which nested iteration loop?
   else ! if Newton step is provided on input
    newtStepScaled(:) = xStep(:) / f_obj % xScale(:)             ! get scaled Newton step (consistent with scaling for aJacScaled and rVecScaled)
    f_obj % rVecScaled(:) = f_obj % fScale(:) * f_obj % f_vec(:) ! matches solve_linear_system
   end if

  end if

  ! get the number of soil layers in the solution vector
  mSoil = size(f_obj % indx_data % var(iLookINDEX % ixMatOnly) % dat)

  ! set unrefined guess
  stateVecTrial(:) = xvec0(:) 

  associate(&
   ! input
   in_SS4HG       => f_obj % in_SS4HG , & 
   fScale         => f_obj % fScale   , & 
   xScale         => f_obj % xScale   , & 
   ! input: SUMMA data structures
   model_decisions => f_obj % model_decisions , &
   lookup_data     => f_obj % lookup_data     , &
   type_data       => f_obj % type_data       , &
   attr_data       => f_obj % attr_data       , &
   mpar_data       => f_obj % mpar_data       , &
   forc_data       => f_obj % forc_data       , &
   bvar_data       => f_obj % bvar_data       , &
   prog_data       => f_obj % prog_data       , &
   ! input-output
   sMul              => f_obj % sMul              , &
   io_SS4HG          => f_obj % io_SS4HG          , &
   indx_data         => f_obj % indx_data         , & 
   diag_data         => f_obj % diag_data         , &
   flux_data         => f_obj % flux_data         , & 
   deriv_data        => f_obj % deriv_data        , &
   dBaseflow_dMatric => f_obj % dBaseflow_dMatric , &
   ! output
   fluxVecNew      => f_obj % fluxVec0      , &
   resSinkNew      => f_obj % rAdd          , &
   resVecNew       => f_obj % resVec        , &
   out_SS4HG       => f_obj % out_SS4HG       &  
  &)
    call refine_Newton_step(in_SS4HG,mSoil,stateVecTrial,newtStepScaled,f_obj%aJacScaled,f_obj%rVecScaled,fScale,xScale,& ! input
                           &model_decisions,lookup_data,type_data,attr_data,mpar_data,forc_data,bvar_data,prog_data,&     ! input
                           &sMul,io_SS4HG,indx_data,diag_data,flux_data,deriv_data,dBaseflow_dMatric,&                    ! input-output
                           &stateVecNew,fluxVecNew,resSinkNew,resVecNew,out_SS4HG,return_flag)                            ! output
  end associate

  ! check for errors in refine_Newton_step call
  if (return_flag) then
   if (f_obj % out_error) then
    write(f_obj % unit,*) "Error in SUMMA_refine_Newton_step: refine_Newton_step message="//trim(f_obj % out_SS4HG % message); stop
   end if
  end if

  ! update inner and outer iteration solutions based on line search
  if (f_obj % nested) then
   associate(nState => f_obj % in_SS4HG % nState)
    ! update inner iteration values
    f_obj % xkp1lp1(:) = f_obj % io_SS4HG % stateVecNewNested(1:nState) ! solution vector
    f_obj % f1_vec(:) = f_obj % io_SS4HG % f1_vec(:)                    ! non-linear function f1
    call filter_SUMMA_Jacobian(f_obj,.false.,f_obj % stateMask1,f_obj % io_SS4HG % aJac1,f_obj % J,f_obj % J1) ! transform aJac1 into J1
    ! update outer iteration values
    f_obj % xk0(:) = f_obj % io_SS4HG % stateVecNewNested(nState+1:2*nState) ! solution vector
    f_obj % f2_vec(:) = f_obj % io_SS4HG % f2_vec(:)                    ! non-linear function f2
    call filter_SUMMA_Jacobian(f_obj,.true.,f_obj % stateMask2,f_obj % io_SS4HG % aJac2,f_obj % J,f_obj % J2) ! transform aJac2 into J2 (negative sign applied)
   end associate
  end if

  ! store refined guess
  xvec1(:) = stateVecNew(:)

 end subroutine SUMMA_refine_Newton_step

 subroutine SUMMA_get_scaled_Jacobian(f_obj,J,aJacScaled)
  ! ** Get scaled SUMMA Jacobian from nested Newton solver Jacobian **
  use matrixOper_module,  only: scaleMatrices
  ! arguments
  type(f_obj_type),intent(inout) :: f_obj ! nested Newton object
  real(r8b),intent(in)           :: J(1:f_obj % nrow,1:f_obj % n)                              ! nested Newton solver Jacobian
  real(rkind),intent(out)        :: aJacScaled(f_obj % in_SS4HG % nLeadDim,f_obj % in_SS4HG % nState) ! scaled SUMMA Jacobian matrix

  ! local
  real(rkind)    :: aJac(f_obj % in_SS4HG % nLeadDim,f_obj % in_SS4HG % nState) ! SUMMA Jacobian matrix (descaled)
  integer(i4b)   :: nBands   ! SUMMA's leading dimension for banded Jacobians
  integer(i4b)   :: err      ! SUMMA error code
  character(256) :: cmessage ! error message from SUMMA

    ! get SUMMA Jacobian from solver Jacobian
    if (f_obj % banded) then ! banded storage
     associate(nrow_banded => f_obj % nrow_banded, n => f_obj % n, subdiag => f_obj % subdiag)
      nBands=nrow_banded+subdiag
      aJac(1:subdiag,1:n) = 0._rkind
      aJac(subdiag+1:nBands,1:n) = J(1:nrow_banded,1:n) ! SUMMA's aJac has extra storage rows
     end associate
    else ! full matrix storage
     aJac(:,:) = J(:,:)
    end if

    ! Scale Jacobian
    associate(ixMatrix => f_obj % in_SS4HG % ixMatrix, nState => f_obj % in_SS4HG % nState)
     call scaleMatrices(ixMatrix,nState,aJac,f_obj % fScale,f_obj % xScale,aJacScaled,err,cmessage) ! matches solve_linear_system
    end associate
    if (err/=0) then
     if (f_obj % out_error) then
      write(f_obj % unit,*) "Error in SUMMA_get_scaled_Jacobian: scaleMatrices message="//trim(cmessage); stop
     end if
    end if

 end subroutine SUMMA_get_scaled_Jacobian

 subroutine SUMMA_nested_line_search(f_obj,option)
  ! ** nested Newton line search **
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj ! nested Newton object
  character(1)     ,intent(in)    :: option

  ! local
  real(r8b) :: L0,L1 ! objective function values
  real(r8b) :: L1_prev ! objective function value from previous line search iteration
  real(r8b) :: initial_solution(1:f_obj % n) ! intial solution vector
  real(r8b) :: updated_solution(1:f_obj % n) ! updated solution vector
  real(r8b) :: p(1:f_obj % n) ! search direction
  real(r8b) :: grad_L(1:f_obj % n) ! gradient of objective function L
  real(rkind) :: aJacScaled(f_obj % in_SS4HG % nLeadDim,f_obj % in_SS4HG % nState) ! scaled SUMMA Jacobian matrix
  real(r8b)            :: m ! local slope
  real(r8b), parameter :: c=1.e-4_r8b!1.e-4_r8b   ! objective function check control parameter
  real(r8b), parameter :: tao=0.5e0_r8b ! step reduction control parameter
  real(r8b), parameter :: m_tol=100._r8b*epsilon(1._r8b)
  real(r8b)            :: alpha ! step size
  real(r8b)            :: alpha_temp ! step size (temporary)
  real(r8b)            :: alpha_prev ! step size (from previous line search iteration)
  real(r8b)            :: rhs1,rhs2,aCoef,bCoef,disc ! constants for cubic interpolant
  logical   :: do_line_search
  logical   :: converged ! checkConv convergence flag
  integer(i4b), parameter :: i_max = 10_i4b!5_i4b ! max number of line search iterations
  integer(i4b)   :: i ! loop index
  integer(i4b)   :: err      ! SUMMA error code
  character(256) :: cmessage ! error message from SUMMA
  logical, parameter :: debug_output=.false.

  ! initial solutions
  if ((option == 'N').or.(option == 'I')) then ! nested or inner cases
   initial_solution(:) = f_obj % xkp1l(:) ! previous inner iteration
  else if (option == 'O') then ! outer case
   initial_solution(:) = f_obj % xk0(:)   ! previous outer iteration
  end if

  ! **** initial setup operations to be moved outside of nested Newton solver loop ****
  ! compute initial objective function (scaled)
  call f_obj % line_search_objective(.true.,option,initial_solution,L0)

  ! compute initial Jacobian (scaled)
  if (option == 'N') then ! nested case
   f_obj % J(:,:) = f_obj % J1(:,:) - f_obj % J2(:,:) ! solver Jacobian
   call SUMMA_get_scaled_Jacobian(f_obj,f_obj % J,aJacScaled) ! get scaled SUMMA Jacobian
  else if (option == 'I') then ! inner case
   f_obj % J(:,:) = f_obj % J1(:,:) - f_obj % J2(:,:) ! solver Jacobian
   call SUMMA_get_scaled_Jacobian(f_obj,f_obj % J,aJacScaled) ! get scaled SUMMA Jacobian
  else if (option == 'O') then ! outer case
   call f_obj % J_eval(initial_solution)
   !f_obj % J(:,:) = f_obj % J1(:,:) - f_obj % J2(:,:) ! solver Jacobian
   call SUMMA_get_scaled_Jacobian(f_obj,f_obj % J,aJacScaled) ! get scaled SUMMA Jacobian  
  else
   print *, "Error in SUMMA_nested_line_search: option is not supported"
   stop
  end if
  ! **** end initial setup operations to be moved outside of nested Newton solver loop ****


  ! compute gradient of objective function (scaled)
  !if ((option == 'N').or.(option == 'I').or.(option == 'O')) then ! all cases?
   call SUMMA_computeGradient(f_obj,aJacScaled,f_obj % rVecScaled,grad_L)
  !end if

  ! compute search direction
  if ((option == 'N').or.(option == 'I')) then ! nested or inner cases
   p(:)=f_obj % xkp1lp1 - f_obj % xkp1l ! inner Newton step
  else if (option == 'O') then ! outer case
   p(:)=f_obj % xkp1lp1 - f_obj % xk0   ! outer Newton step
  end if


  ! compute local slope (use scaled values)
  m = dot_product(grad_L,p(:)/f_obj % xScale(:)) ! confirmed against homegrown line search

  ! check that local slope is negative (needed to reduce the line search objective function)
  if (m < 0._rkind) then
   do_line_search = .true.
  else if ((0._r8b <= m).and.(m <= m_tol)) then ! non-negative slope with allowance for round-off error
   do_line_search = .false. ! skip line search (there would be no improvement anyway)
  else ! non-negative but exceeding tolerance
   print *, "m=",m
   print *, "option=",option
   print *, "Error in SUMMA_nested_line_search: initial slope is non-negative"
   stop
  end if

  if (debug_output) then
   print *, "NN Line Search:"
   print *, "m=",m
   print *, "L0=",L0
  end if
  
  alpha = 1._r8b ! initialize line search step size
  line_search: do i=1_i4b,i_max+1_i4b
   updated_solution(:) = initial_solution(:) + alpha*p(:)

   ! impose constraints (calls SUMMA's imposeConstraints routine)
   call f_obj % apply_constraints(initial_solution,updated_solution)

   ! compute objective function
   call f_obj % line_search_objective(.true.,option,updated_solution,L1)

   ! update variables in nested Newton algorithm
   if (option == 'N') then ! nested case
    ! f1 quantities
    call filter_SUMMA_f(.false.,f_obj % stateMask1,f_obj % f_vec,f_obj % f1_vec) ! get f1 from total f
    call f_obj % J1_eval(updated_solution)   ! get J1 based on eval8summa call for total f 
    f_obj % xkp1lp1(:) = updated_solution(:) ! apply updated solution
    ! f2 quantities
    call filter_SUMMA_f(.true.,f_obj % stateMask2,f_obj % f_vec,f_obj % f2_vec) ! get f2 from total f
    call f_obj % J2_eval(updated_solution)   ! get J2 based on eval8summa call for total f 
    f_obj % xk0(:) = updated_solution(:)     ! apply updated solution
   else if (option == 'I') then ! inner case
    !call f_obj % J1_eval(updated_solution)   ! get J1 based on eval8summa call for f1 
    f_obj % xkp1lp1(:) = updated_solution(:) ! apply updated solution
   else if (option == 'O') then ! outer case
    !call filter_SUMMA_f(.true.,f_obj % stateMask2,f_obj % f_vec,f_obj % f2_vec) ! get f2 from total f -------- testing
    !call f_obj % J2_eval(updated_solution)   ! get J2 based on eval8summa call for f2 
    f_obj % xkp1lp1(:) = updated_solution(:) ! apply updated solution
   end if

   ! check SUMMA's feasibility flag ------------------ turn this into a recoverable error
   if (.not.(f_obj % feasible)) then
    print *, "Error in SUMMA_nested_line_search: not feasible"
    stop
   end if

   ! get convergence flag
   converged = SUMMA_checkConv(f_obj,p,f_obj % xkp1lp1)
   
   if (debug_output) then
    print *, "i=",i
    print *, "alpha=",alpha
    print *, "L1=",L1
    print *, "L0 + alpha*c*m=",L0 + alpha*c*m
   end if

   ! check convergence
   if (converged) then
    f_obj % xkp1lp1(:) = updated_solution(:) ! accept updated solution and exit
    return 
   end if
   
   ! exit early if not doing the line search (only alpha=1.0 solution is used)
   if (.not.do_line_search) return

   ! check if the objective function is accepted using the Armijo-Goldstein Criterion
   if (L1 <= L0 + alpha*c*m) return

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

 end subroutine SUMMA_nested_line_search

 subroutine SUMMA_line_search_objective(f_obj,evaluate_f,option,solution,L)
  ! ** compute line search objective function for SUMMA **
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj ! nested Newton object
  logical         ,intent(in)    :: evaluate_f ! perform evaluations for f, f1, or f2? (if not, use stored values)
  character(1)    ,intent(in)    :: option ! line search option
  real(r8b)       ,intent(in)    :: solution(1:f_obj % n) ! updated solution vector
  real(r8b)       ,intent(out)   :: L ! objective function value

  ! local
  real(rkind) :: aJacScaled(f_obj % in_SS4HG % nLeadDim,f_obj % in_SS4HG % nState) ! scaled SUMMA Jacobian matrix

  if (option == 'N') then ! nested case
   if (evaluate_f) call f_obj % f_vec_eval(solution) ! update total f
   L=f_obj % out_SS4HG % fNew ! scaled
  else if (option == 'I') then ! inner case
   if (evaluate_f) call f_obj % f1_vec_eval(solution) ! update f1
   if (f_obj % banded) then
    print *, "Error in SUMMA_line_search_objective: banded Jacobians not implemented"; stop
   else
    f_obj % rVecScaled(:) = f_obj % fScale(:) * ( &
                        & f_obj % f1_vec(:) - ( f_obj % f2_vec(:) + matmul(f_obj % J2,solution - f_obj % xk0) )&
                        & )
   end if
   L = 0.5_rkind*dot_product(f_obj % rVecScaled,f_obj % rVecScaled)
  else if (option == 'O') then ! outer case
   !if (evaluate_f) call f_obj % f2_vec_eval(solution) ! update f2
   !if (f_obj % banded) then
   ! print *, "Error in SUMMA_line_search_objective: banded Jacobians not implemented"; stop
   !else
   ! f_obj % rVecScaled(:) = f_obj % fScale(:) * ( &
   !                     & f_obj % f1_vec(:) + matmul(f_obj % J1,solution - f_obj % xkp1l) - f_obj % f2_vec(:) &
   !                     & )
   !end if
   if (evaluate_f) call f_obj % f_vec_eval(solution) ! update total f
   L=f_obj % out_SS4HG % fNew ! scaled
   L = 0.5_rkind*dot_product(f_obj % rVecScaled,f_obj % rVecScaled)
  else 
   print *, "Error in SUMMA_line_search_objective: option is not supported"; stop
  end if
 end subroutine SUMMA_line_search_objective

 subroutine SUMMA_computeGradient(f_obj,aJacScaled,rVecScaled,gradScaled)
  ! ** compute product of scaled residual and scaled SUMMA Jacobian **
  use matrixOper_module, only: computeGradient
  ! arguments
  type(f_obj_type),intent(in) :: f_obj ! nested Newton object
  real(rkind),intent(in) :: aJacScaled(f_obj % in_SS4HG % nLeadDim,f_obj % in_SS4HG % nState) ! scaled SUMMA Jacobian matrix
  real(r8b),intent(in)   :: rVecScaled(1:f_obj % n) ! gradient of objective function L
  real(r8b),intent(out)  :: gradScaled(1:f_obj % n) ! gradient of objective function L

  ! local
  integer(i4b)   :: err      ! SUMMA error code
  character(256) :: cmessage ! error message from SUMMA

  associate(ixMatrix => f_obj % in_SS4HG % ixMatrix, nState => f_obj % in_SS4HG % nState)
   call computeGradient(ixMatrix,nState,aJacScaled,rVecScaled,gradScaled,err,cmessage)
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
  real(r8b),intent(inout) :: B(1:f_obj % n,1:1) ! right-hand side vector

  ! local
  integer(i4b)   :: err
  character(256) :: cmessage

  ! get scaled variables (accoring to SUMMA's fScale and xScale vectors)
  ! note: need to match scaling applied in solve_linear_system subroutine in summaSolve4homegrown
  B(:,1) = f_obj % fScale(:) * B(:,1) ! matches solve_linear_system
  if (f_obj % nested) then ! nested iterations
   f_obj % rVecScaled(:) = (f_obj % f1_vec(:) - f_obj % f2_vec(:)) * f_obj % fScale(:) ! compute scaled residual for Newton step refinement
  else ! classical iterations
   f_obj % rVecScaled(:) = -B(:,1) ! save scaled residual for reuse in Newton step refinement
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
  real(r8b),intent(inout) :: B(1:f_obj % n,1:1) ! solution side vector

  B(:,1) = B(:,1) * f_obj % xScale(:)
  
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
  real(r8b),intent(in)            :: step(1:f_obj % n)  ! Newton step (iteration increment)
  real(r8b),intent(in)            :: xvec1(1:f_obj % n) ! updated solution vector

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
  class(f_obj_type),intent(inout)   :: f_obj
  real(r8b),intent(in)              :: xvec0(1:f_obj % n) ! previous guess vector
  real(r8b),intent(inout)           :: xvec1(1:f_obj % n) ! current guess vector

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
  real(r8b),intent(in)              :: xvec(1:f_obj % n) ! current guess

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
                    f_obj % dBaseflow_dMatric,       & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                    ! output
                    f_obj % feasible,                & ! intent(out):   flag to denote the feasibility of the solution
                    f_obj % fluxVec0,                & ! intent(out):   flux vector
                    f_obj % fRHS,                    & ! intent(out):   RHS function for ARKODE
                    f_obj % rAdd,                    & ! intent(out):   additional (sink) terms on the RHS of the state equation
                    f_obj % resVec,                  & ! intent(out):   residual vector
                    f_obj % out_SS4HG % fNew,        & ! intent(out):   function evaluation
                    f_obj % out_SS4HG % err,         & ! intent(out): error code
                    f_obj % out_SS4HG % message)       ! intent(out): error message (note: eval8summa uses "cmessage" instead)
  end associate

  ! finalize
  associate(err => f_obj % out_SS4HG % err, message => f_obj % out_SS4HG % message) 
   if (err /= 0) then
    if (f_obj % out_error) then
     write(f_obj % unit,*) "Error in SUMMA_eval8summa: eval8summa message="//trim(message); stop
    end if
   end if
  end associate
 end subroutine SUMMA_eval8summa

 subroutine SUMMA_computJacob(f_obj,&
                             &indx_data,diag_data,flux_data,deriv_data,&
                             &dMat,dBaseflow_dMatric,&
                             &aJac)
  ! ** Interface for SUMMA's computJacob subroutine **
  ! arguments
  class(f_obj_inputs),intent(inout) :: f_obj
  type(var_ilength),intent(in)      :: indx_data              ! indices defining model states and layers for selected split 
  type(var_dlength),intent(in)      :: diag_data              ! diagnostic variables for a local HRU
  type(var_dlength),intent(in)      :: flux_data              ! flux data
  type(var_dlength),intent(in)      :: deriv_data             ! derivative data
  real(rkind)      ,intent(in)      :: dMat(:)          ! diagonal matrix (no flux derivatives) for split
  real(rkind)      ,intent(in)      :: dBaseflow_dMatric(:,:) ! derivative in baseflow w.r.t. matric head (s-1)
  real(rkind)      ,intent(out)     :: aJac(:,:) ! SUMMA's unscaled Jacobian matrix

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
   ixMatrix       => f_obj % in_SS4HG % ixMatrix       ,& ! intent(in): type of matrix (full or band diagonal)
   computeVegFlux => f_obj % in_SS4HG % computeVegFlux  & ! intent(in): flag to indicate if computing fluxes over vegetation
  &)   
   call in_computJacob % initialize(dt_cur,nSnow,nSoil,nLayers,computeVegFlux,(ixGroundwater==qbaseTopmodel),ixMatrix)
  end associate 

   ! update
   associate(&
    prog_data         => f_obj % prog_data&         ! prognostic variables for a local HRU
   &)
    call computJacob(in_computJacob,indx_data,prog_data,diag_data,deriv_data,dBaseflow_dMatric,dMat,aJac,out_computJacob)
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
  logical(lgt),parameter          :: dual = .false.

  ! * initialize operations for split_select object *

  associate(nstate => f_obj % in_SS4HG % nState)
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
  if (dual) then ! assign masks for f1 (energy) and f2 (mass)
   f_obj % stateMask2 = .not.(split_select % stateMask(:)) ! negate energy mask to find mass mask --- allocate on assignment
   f_obj % stateMask1 = split_select % stateMask           ! no transformation --- allocate on assignment
  else ! assign masks for f1 (mass) and f2 (energy)
   f_obj % stateMask1 = .not.(split_select % stateMask(:)) ! negate energy mask to find mass mask --- allocate on assignment
   f_obj % stateMask2 = split_select % stateMask           ! no transformation --- allocate on assignment
  end if

  ! get counts for mass and energy splits
  if (dual) then
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
  real(r8b)       ,intent(in)    :: f_total(:)  ! total f in nested Newton solver
  real(r8b)       ,intent(out)   :: f_filter(:) ! filtered f in nested Newton solver 

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
  logical(lgt)    ,intent(in)    :: stateMask(:)  ! logical mask for filtering
  real(rkind)     ,intent(inout) :: aJac(:,:)     ! total Jacobian from SUMMA
  real(r8b)       ,intent(out)   :: J_total(:,:)  ! total Jacobian in nested Newton solver storage scheme
  real(r8b)       ,intent(out)   :: J_filter(:,:) ! filtered Jacobian in nested Newton solver storage scheme 

  ! local variables
  integer(i4b) :: i,j,k  ! loop indices
  integer(i4b) :: nBands ! # of bands for banded storage

  ! store Jacobian used in solver
  if (f_obj % banded) then ! banded storage
   associate(nrow_banded => f_obj % nrow_banded, n => f_obj % n, subdiag => f_obj % subdiag, superdiag => f_obj % superdiag)
    nBands=nrow_banded+subdiag ! number of non-zero bands
    J_total(1:nrow_banded,1:n) = aJac(subdiag+1:nBands,1:n) ! aJac has extra storage rows (store total Jacobian)
    if (negative) then
     do j=1,n ! column index for dense and banded storage
      do i=max(1,j-superdiag),min(n,j+subdiag) ! row index for dense storage
       k = nrow_banded+i-j ! row index for LAPACK banded storage (nrow_banded = subdiag+superdiag+1)
       aJac(k,j) = merge(-aJac(k,j),0._rkind,stateMask(i)) ! zero the elements that are not included in J1
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
    J_filter(1:nrow_banded,1:n) =  aJac(subdiag+1:nBands,1:n) ! aJac has extra storage rows
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

 subroutine f_mass_SUMMA_vec_full(f_obj,xvec)
  ! *** Compute mass non-linear function --- use fully-coupled eval8summa call and filter results ***
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)            :: xvec(1:f_obj % n) ! current guess

  ! local
  real(qp)                        :: resVec(1:f_obj % n) ! residual vector for split

  ! note: data structures and variables for f1 are initialized in systemSolv

  call f_obj % f_state_SUMMA_vec_full(&
               &xvec,&
               &f_obj % indx_data,f_obj % diag_data,f_obj % flux_data,f_obj % deriv_data,f_obj % sMul,&
               &f_obj % dBaseflow_dMatric,resVec)

  ! store total non-linear function
  f_obj % f_vec(:) = real(resVec(:),r8b)

  ! assign non-zero function values based on logical mask
  call filter_SUMMA_f(.false.,f_obj % stateMask1,f_obj % f_vec,f_obj % f1_vec)

 end subroutine f_mass_SUMMA_vec_full

 subroutine Jacobian_f_mass_SUMMA_vec_full(f_obj,xvec)
  ! *** Compute Jacobian for mass non-linear function --- use fully-coupled computJacob call and filter results ***
  ! NOTE: assumes appropriate eval8summa call has already been made to get the fluxes
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)            :: xvec(1:f_obj % n) ! current guess (needed for interface)

  ! local
  real(rkind)  :: aJac(f_obj % in_SS4HG % nLeadDim,f_obj % in_SS4HG % nState) ! SUMMA's unscaled Jacobian matrix
  integer(i4b) :: i,j,k  ! loop indices
  integer(i4b) :: nBands ! # of bands for banded storage

  call f_obj % SUMMA_computJacob(&
               &f_obj % indx_data,f_obj % diag_data,f_obj % flux_data,f_obj % deriv_data,&
               &f_obj % dMat,f_obj % dBaseflow_dMatric,&
               &aJac)

  ! get nested Newton solver Jacobian J1
  call filter_SUMMA_Jacobian(f_obj,.false.,f_obj % stateMask1,aJac,f_obj % J,f_obj % J1)

 end subroutine Jacobian_f_mass_SUMMA_vec_full

 subroutine f_energy_SUMMA_vec_full(f_obj,xvec)
  ! *** Compute energy non-linear function --- use fully-coupled eval8summa call and filter results ***
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)            :: xvec(1:f_obj % n) ! current guess

  ! local
  real(qp)                        :: resVec(1:f_obj % n) ! residual vector for split

  ! note: data structures and variables for f1 are initialized in systemSolv

  call f_obj % f_state_SUMMA_vec_full(&
               &xvec,&
               &f_obj % indx_data,f_obj % diag_data,f_obj % flux_data,f_obj % deriv_data,f_obj % sMul,&
               &f_obj % dBaseflow_dMatric,resVec)

  ! store total non-linear function
  f_obj % f_vec(:) = real(resVec(:),r8b)

  ! assign non-zero function values based on logical mask
  call filter_SUMMA_f(.true.,f_obj % stateMask2,f_obj % f_vec,f_obj % f2_vec)

 end subroutine f_energy_SUMMA_vec_full

 subroutine Jacobian_f_energy_SUMMA_vec_full(f_obj,xvec)
  ! *** Compute Jacobian for energy non-linear function --- use fully-coupled computJacob call and filter results ***
  ! NOTE: assumes appropriate eval8summa call has already been made to get the fluxes
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)            :: xvec(1:f_obj % n) ! current guess (needed for interface)

  ! local
  real(rkind)  :: aJac(f_obj % in_SS4HG % nLeadDim,f_obj % in_SS4HG % nState) ! SUMMA's unscaled Jacobian matrix
  integer(i4b) :: i,j,k ! loop indices
  integer(i4b) :: nBands ! # of bands for banded storage

  call f_obj % SUMMA_computJacob(&
               &f_obj % indx_data,f_obj % diag_data,f_obj % flux_data,f_obj % deriv_data,&
               &f_obj % dMat,f_obj % dBaseflow_dMatric,&
               &aJac)

  ! get nested Newton solver Jacobian J2
  call filter_SUMMA_Jacobian(f_obj,.true.,f_obj % stateMask2,aJac,f_obj % J,f_obj % J2) ! negative sign applied

 end subroutine Jacobian_f_energy_SUMMA_vec_full

 subroutine Jacobian_f_mass_energy_SUMMA_vec_full(f_obj,xvec)
  ! *** Compute Jacobian for mass and energy non-linear functions --- use fully-coupled computJacob call and filter results ***
  ! NOTE: assumes appropriate eval8summa call has already been made to get the fluxes
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)            :: xvec(1:f_obj % n) ! current guess (needed for interface)

  ! local
  real(rkind)  :: aJac(f_obj % in_SS4HG % nLeadDim,f_obj % in_SS4HG % nState) ! SUMMA's unscaled Jacobian matrix
  integer(i4b) :: i,j,k ! loop indices
  integer(i4b) :: nBands ! # of bands for banded storage

  call f_obj % SUMMA_computJacob(&
               &f_obj % indx_data,f_obj % diag_data,f_obj % flux_data,f_obj % deriv_data,&
               &f_obj % dMat,f_obj % dBaseflow_dMatric,&
               &aJac)

  ! get nested Newton solver Jacobian J1
  call filter_SUMMA_Jacobian(f_obj,.false.,f_obj % stateMask1,aJac,f_obj % J,f_obj % J1)

  ! get nested Newton solver Jacobian J2
  call filter_SUMMA_Jacobian(f_obj,.true.,f_obj % stateMask2,aJac,f_obj % J,f_obj % J2) ! negative sign applied

 end subroutine Jacobian_f_mass_energy_SUMMA_vec_full

 subroutine f_state_SUMMA_vec_full(f_obj,xvec,&
                                  &indx_data,diag_data,flux_data,deriv_data,sMul,&
                                  &dBaseflow_dMatric,resVec)
  ! *** Compute SUMMA's vector non-linear function for mass or energy state variables -- uses fully-coupled eval8summa call ***
  ! ** NOTE: the fully-coupled solution method in SUMMA's opSplittin is assumed **

  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)            :: xvec(1:f_obj % n) ! current guess

  type(var_ilength),intent(inout) :: indx_data            ! indices defining model states and layers for selected split 
  type(var_dlength),intent(inout) :: diag_data            ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout) :: flux_data            ! flux data
  type(var_dlength),intent(inout) :: deriv_data           ! derivative data
  real(qp)         ,intent(inout) :: sMul(:)              ! state vector multipliers
  real(rkind)      ,intent(out)   :: dBaseflow_dMatric(:,:) ! baseflow derivative matrix w.r.t pressure head
  real(qp)         ,intent(out)   :: resVec(1:f_obj % n)    ! residual vector

  ! local
  real(rkind) :: fluxVec0(1:f_obj % n) ! flux vector
  real(rkind) :: fRHS(1:f_obj % n)     ! RHS function for ARKODE
  real(rkind) :: rAdd(1:f_obj % n)     ! additional (sink) terms on the RHS of the state equation

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
                   dBaseflow_dMatric,               & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                   ! output
                   f_obj % feasible,                & ! intent(out):   flag to denote the feasibility of the solution
                   fluxVec0,                        & ! intent(out):   flux vector
                   fRHS,                            & ! intent(out):   RHS function for ARKODE
                   rAdd,                            & ! intent(out):   additional (sink) terms on the RHS of the state equation
                   resVec,                          & ! intent(out):   residual vector
                   f_obj % out_SS4HG % fNew,        & ! intent(out):   function evaluation
                   f_obj % out_SS4HG % err,         & ! intent(out): error code
                   f_obj % out_SS4HG % message)       ! intent(out): error message (note: eval8summa uses "cmessage" instead)

  ! finalize
  associate(err => f_obj % out_SS4HG % err, message => f_obj % out_SS4HG % message) 
   if (err /= 0) then
    if (f_obj % out_error) then
     write(f_obj % unit,*) "Error f_state_SUMMA_vec: eval8summa message="//trim(message); stop
    end if
   end if
  end associate

 end subroutine f_state_SUMMA_vec_full


 subroutine f_mass_SUMMA_vec(f_obj,xvec)
  ! SJT: -------- incomplete and may not include contributions due to entire input state vector --------
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)            :: xvec(1:f_obj % n) ! current guess

  ! local
  logical,parameter               :: mass_state_type = .true. ! perform transformations from energy split to mass split 
  real(rkind),allocatable         :: resVec_split(:) ! residual vector for split
  integer(i4b)                    :: nBands          ! # of bands for SUMMA's Jacobian

  ! initialize data structures (ensure that fully-coupled structures are not overwritten)
  f_obj % indx_data1  = f_obj % indx_data 
  f_obj % diag_data1  = f_obj % diag_data 
  f_obj % flux_data1  = f_obj % flux_data 
  f_obj % deriv_data1 = f_obj % deriv_data 
  f_obj % dBaseflow_dMatric1 = f_obj % dBaseflow_dMatric ! allocate

  ! compute residual for split
  call f_state_SUMMA_vec(f_obj,xvec,mass_state_type,&
                        &f_obj % indx_data1,f_obj % diag_data1,f_obj % flux_data1,f_obj % deriv_data1,&
                        &f_obj % nSubset1,f_obj % stateMask1,f_obj % dMat1,f_obj % dBaseflow_dMatric1,resVec_split)

  ! remaining non-linear function values are zero
  f_obj % f1_vec(:) = 0._r8b
  f_obj % f1_vec(:) = unpack(resVec_split,f_obj % stateMask1,f_obj % f1_vec)

  ! get leading dimension for Jacobians
  if (f_obj % banded) then
   nBands=f_obj % nrow_banded + f_obj% subdiag
   f_obj % nLeadDim1 = nBands  
  else
   f_obj % nLeadDim1 = f_obj % nSubset1 
  end if

  !!!! SJT: start test block ---- take out ----
  !print *, "mass_state_type =",mass_state_type
  !print *, split_select % nState
  !print *, f_obj % nLeadDim1,f_obj % nSubset1,size(f_obj % dMat1)
  !print *, f_obj % stateMask1(:)
  !print *, resVec_split
  !print *, f_obj % f1_vec
  !print *, sum(resVec_split)
  !!!! SJT: end test block ---- take out ----
 end subroutine f_mass_SUMMA_vec

 subroutine f_energy_SUMMA_vec(f_obj,xvec)
  ! SJT: -------- incomplete and may not include contributions due to entire input state vector --------
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)            :: xvec(1:f_obj % n) ! current guess

  ! local
  logical,parameter               :: mass_state_type = .false. ! perform transformations from energy split to mass split 
  real(rkind),allocatable         :: resVec_split(:) ! residual vector for split
  integer(i4b)                    :: nBands          ! # of bands for SUMMA's Jacobian

  ! initialize data structures (ensure that fully-coupled structures are not overwritten)
  f_obj % indx_data2  = f_obj % indx_data
  f_obj % diag_data2  = f_obj % diag_data 
  f_obj % flux_data2  = f_obj % flux_data 
  f_obj % deriv_data2 = f_obj % deriv_data 
  f_obj % dBaseflow_dMatric2 = f_obj % dBaseflow_dMatric ! allocate

  ! compute residual for split
  call f_state_SUMMA_vec(f_obj,xvec,mass_state_type,&
                        &f_obj % indx_data2,f_obj % diag_data2,f_obj % flux_data2,f_obj % deriv_data2,&
                        &f_obj % nSubset2,f_obj % stateMask2,f_obj % dMat2,f_obj % dBaseflow_dMatric2,resVec_split)

  ! remaining non-linear function values are zero
  f_obj % f2_vec(:) = 0._r8b
  f_obj % f2_vec(:) = unpack(-resVec_split,f_obj % stateMask2,f_obj % f2_vec) ! note: sign change for f2 so that f=f1-f2

  ! get leading dimension for Jacobians
  if (f_obj % banded) then
   nBands=f_obj % nrow_banded + f_obj% subdiag
   f_obj % nLeadDim2 = nBands  
  else
   f_obj % nLeadDim2 = f_obj % nSubset2 
  end if

  !!!! SJT: start test block ---- take out ----
  !print *, "mass_state_type =",mass_state_type
  !print *, split_select % nState
  !print *, split_select % nSubset
  !print *, f_obj % stateMask2(:)
  !print *, resVec_split
  !print *, f_obj % f2_vec
  !print *, sum(resVec_split)
  !!!! SJT: end test block ---- take out ----
 end subroutine f_energy_SUMMA_vec

 subroutine f_state_SUMMA_vec(f_obj,xvec,mass_state_type,&
                             &indx_data,diag_data,flux_data,deriv_data,&
                             &nSubset,stateMask,dMat_split,dBaseflow_dMatric,resVec_split)
  ! *** Compute SUMMA's vector non-linear function for mass or energy state variables ***
  ! ** NOTE: the fully-coupled solution method in SUMMA's opSplittin is assumed **
  ! SJT: -------- incomplete and may not include contributions due to entire input state vector --------
  use stateFilter_module,only: fullyCoupled,stateTypeSplit
  use stateFilter_module,only: massSplit,nrgSplit
  use stateFilter_module,only: fullDomain,subDomain
  use stateFilter_module,only: vector,scalar
  use indexState_module ,only: indexSplit                             ! get state indices from stateMask
  use data_types        ,only: in_type_indexSplit,out_type_indexSplit ! argument objects for indexSplit
  use getVectorz_module ,only: popStateVec                            ! populate the state vector
  use getVectorz_module ,only: getScaling                             ! scale factors for residual and solution vectors
  use updateVars_module ,only: updateProg                             ! update prognostic (state) variable data structure
  use mDecisions_module ,only: closedForm                             ! use temperature with closed form heat capacity
  use mDecisions_module ,only: enthalpyFormLU                         ! use look up tables for soil enthalpy
  use mDecisions_module ,only: ida                                    ! use IDA solver

  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b)        ,intent(in)    :: xvec(1:f_obj % n)    ! current guess
  logical          ,intent(in)    :: mass_state_type      ! perform transformations from energy split to mass split 
  type(var_ilength),intent(inout) :: indx_data            ! indices defining model states and layers for selected split 
  type(var_dlength),intent(inout) :: diag_data            ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout) :: flux_data            ! flux data
  type(var_dlength),intent(inout) :: deriv_data           ! derivative data
  integer(i4b)            ,intent(out) :: nSubset           ! # of state variables in split
  logical(lgt),allocatable,intent(out) :: stateMask(:)           ! logical mask array for split
  real(rkind),allocatable ,intent(out) :: dMat_split(:)          ! diagonal matrix (no flux derivatives) for split
  real(rkind)             ,intent(out) :: dBaseflow_dMatric(:,:) ! derivative in baseflow w.r.t. matric head (s-1)
  real(rkind),allocatable ,intent(out) :: resVec_split(:)        ! residual vector for split

  ! local variables
  type(split_select_type)         :: split_select      ! split select object
  type(in_type_indexSplit)        :: in_indexSplit     ! indexSplit arguments
  type(out_type_indexSplit)       :: out_indexSplit
  real(rkind)                     :: stateVecPrime(1:f_obj % n) ! trial state vector for full solution (prime variables -- not used here)
  real(rkind)                     :: untappedMelt(1:f_obj % n)  ! untapped melt energy
  real(rkind),allocatable         :: stateVecTrial(:)  ! trial state vector for split
  real(rkind),allocatable         :: fScale_split(:)   ! residual vector scale factors for split
  real(rkind),allocatable         :: xScale_split(:)   ! solution vector scale factors for split
  real(qp)   ,allocatable         :: sMul_split(:)     ! state vector multipliers for split
  real(rkind),allocatable         :: fluxVec0_split(:) ! flux vector for split
  real(rkind),allocatable         :: fRHS_split(:)     ! RHS function for ARKODE for split
  real(rkind),allocatable         :: rAdd_split(:)     ! additional (sink) terms on the RHS of the state equation for split
  logical(lgt)                    :: enthalpyStateVec  ! flag to use enthalpy as a state variable (ida)
  logical(lgt)                    :: computeEnthTemp   ! flag to use enthalpy temperature
  logical(lgt)                    :: use_lookup        ! flag to use enthalpy lookup tables
  logical(lgt)                    :: waterBalanceError,nrgFluxModified ! output flags for updateProg
  real(rkind)                     :: balance(1:f_obj % n) ! balance error from updateProg
  character(LEN=256)              :: message           ! total error message
  character(LEN=256)              :: cmessage          ! error message of downwind routine
  integer(i4b)                    :: err               ! error code of downwind routine
  logical(lgt)                    :: return_flag

  logical(lgt) :: firstFluxCall

  ! * updateProg *
  ! put xvec state variable values into prog_temp object
 
  f_obj % prog_temp = f_obj % prog_data ! initialize
  untappedMelt  = 0._rkind  ! set untapped melt energy to zero (matches systemSolv)
  stateVecPrime = 0._rkind  ! state vector (primed variables -- not used)
  associate(&
   nSnow          => f_obj % in_SS4HG % nSnow          ,& ! intent(in): number of snow layers
   nSoil          => f_obj % in_SS4HG % nSoil          ,& ! intent(in): number of soil layers
   nLayers        => f_obj % in_SS4HG % nLayers        ,& ! intent(in): total number of layers
   ixNumericalMethod => f_obj % model_decisions(iLookDECISIONS%num_method)%iDecision,& ! intent(in): [i4b] choice of numerical solver
   ixNrgConserv      => f_obj % model_decisions(iLookDECISIONS%nrgConserv)%iDecision,& ! intent(in): [i4b] choice of variable in either energy backward Euler residual or IDA state variable
   doAdjustTemp      => mass_state_type, & ! flag to adjust temperature (same behaviour as homegrown mass split from opSplittin assumed) 
   computeVegFlux    => f_obj % in_SS4HG % computeVegFlux,  & ! intent(in): flag to indicate if computing fluxes over vegetation
   computMassBalance => .false., &
   computNrgBalance  => .false.  &
  &)

   ! compute flags based on solver choices (following usage in varSubstep)
   enthalpyStateVec = (ixNrgConserv .ne. closedForm .and. ixNumericalMethod==ida) ! enthalpy as state variable (ida -- matches usage in varSubstep)
   computeEnthTemp = ((ixNrgConserv .ne. closedForm .or. computNrgBalance) .and. ixNumericalMethod .ne. ida) ! use enthTemp to conserve energy or compute energy balance
   use_lookup = (ixNrgConserv==enthalpyFormLU) ! use lookup tables for soil enthalpy instead of analytical solution

   call updateProg(f_obj % in_SS4HG % dt_cur,nSnow,nSoil,nLayers,untappedMelt,xVec,stateVecPrime,& ! input: states
                  &doAdjustTemp,computeVegFlux,computMassBalance,computNrgBalance,computeEnthTemp,enthalpyStateVec,use_lookup,& ! input: model control
                  &f_obj % model_decisions,f_obj % lookup_data,&
                  &f_obj % mpar_data,indx_data,flux_data,f_obj % prog_temp,diag_data,deriv_data,    & ! input-output: data structures
                  &f_obj % fluxVec0,f_obj % resVec,balance,waterBalanceError,nrgFluxModified,err,message) ! input-output: balances, flags, and error control
  end associate

  ! * initialize operations for split_select object *

  associate(nstate => f_obj % in_SS4HG % nState)
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

  ! apply steps similar to initialize_split from opSplitting to generate logical masks (probably skip save/restore operations)
  ! note: from update_stateMask in opSplittin

  ! compute stateMask and nSubset (in split_select object) for the selected split
  call split_select % get_stateMask(indx_data,err,cmessage,message,return_flag)
  if (return_flag) then
    if (f_obj % out_error) then
     write(f_obj % unit,*) "Error in f_state_SUMMA_vec: stateFilter message="//trim(cmessage); stop
    end if
  end if

  ! transform variables for energy split into mass split
  if (mass_state_type) then
   stateMask = .not.(split_select % stateMask(:))       ! negate energy mask to find mass mask --- allocate on assignment
   split_select % stateMask(:) = stateMask(:)           ! update split_select object in case of future use
   split_select % nSubset = split_select % nState - split_select % nSubset ! count for new stateMask
  else
   stateMask = split_select % stateMask ! no transformation --- allocate on assignment
  end if
  nSubset = split_select % nSubset ! for argument list

  ! * indexSplit *
  associate(&
   nSnow          => f_obj % in_SS4HG % nSnow          ,& ! intent(in): number of snow layers
   nSoil          => f_obj % in_SS4HG % nSoil          ,& ! intent(in): number of soil layers
   nLayers        => f_obj % in_SS4HG % nLayers         & ! intent(in): total number of layers
  &)   
   call in_indexSplit % initialize(nSnow,nSoil,nLayers,split_select % nSubset)
  end associate
  call indexSplit(in_indexSplit,stateMask,indx_data,out_indexSplit) ! update indx_data based on stateMask
  call out_indexSplit % finalize(err,cmessage)
  if (err/=0_i4b) then
    if (f_obj % out_error) then
     write(f_obj % unit,*) "Error in f_state_SUMMA_vec: indexSplit message="//trim(cmessage); stop
    end if
  end if


  ! call eval8summa to get non-linear function values for mass state type
  ! update
  associate(&
   nState            => split_select % nSubset                                      ,& ! # of state variables in split
   ixNumericalMethod => f_obj % model_decisions(iLookDECISIONS%num_method)%iDecision,& ! intent(in): [i4b] choice of numerical solver
   ixNrgConserv      => f_obj % model_decisions(iLookDECISIONS%nrgConserv)%iDecision & ! intent(in): [i4b] choice of variable in either energy backward Euler residual or IDA state variable
  &)

   ! allocate arrays for split
   allocate(stateVecTrial(1:nState))  ! state vector
   allocate(fScale_split(1:nState) )  ! residual scale factors
   allocate(xScale_split(1:nState) )  ! solution scale factors
   allocate(dMat_split(1:nState)   )  ! diagonal matrix (no flux derivatives)
   allocate(sMul_split(1:nState)   )  ! state vector multipliers
 
   allocate(fluxVec0_split(1:nState)) ! flux vector
   allocate(fRHS_split(1:nState)    ) ! RHS function for ARKODE
   allocate(rAdd_split(1:nState)    ) ! additional (sink) terms on the RHS of the state equation
   allocate(resVec_split(1:nState)  ) ! residual vector

   enthalpyStateVec = (ixNrgConserv .ne. closedForm .and. ixNumericalMethod==ida) ! enthalpy as state variable (ida -- matches usage in varSubstep)

   ! initialize state vectors
   call popStateVec(&
                   ! input
                   nState,             & ! intent(in):  number of desired state variables
                   enthalpyStateVec,   & ! intent(in):  flag to use enthalpy as a state variable
                   f_obj % prog_temp,  & ! intent(in):  model prognostic variables for a local HRU
                   diag_data,          & ! intent(in):  model diagnostic variables for a local HRU
                   indx_data,          & ! intent(in):  indices defining model states and layers
                   ! output
                   stateVecTrial,      & ! intent(out): initial model state vector (mixed units)
                   err,cmessage)         ! intent(out): error control
   if (err/=0_i4b) then
     if (f_obj % out_error) then
      write(f_obj % unit,*) "Error in f_state_SUMMA_vec: popStateVec message="//trim(cmessage); stop
     end if
   end if

!   !!!! SJT: testing -- take dependency on xvec input argument into account
!   stateVecTrial = pack(xvec,stateMask)
!   !!!! SJT: end testing

   ! compute scale factors
   call getScaling(diag_data,indx_data,fScale_split,xScale_split,sMul_split,dMat_split,err,cmessage)     
   if (err/=0_i4b) then
     if (f_obj % out_error) then
      write(f_obj % unit,*) "Error in f_state_SUMMA_vec: getScaling message="//trim(cmessage); stop
     end if
   end if

   ! evaluate residual vector for mass split
   firstFluxCall = .true. ! may not be needed
   call eval8summa(&
                    ! input: model control
                    f_obj % in_SS4HG % dt_cur,         & ! intent(in):    current stepsize
                    f_obj % in_SS4HG % dt,             & ! intent(in):    length of the entire time step (seconds) for drainage pond rate
                    f_obj % in_SS4HG % nSnow,          & ! intent(in):    number of snow layers
                    f_obj % in_SS4HG % nSoil,          & ! intent(in):    number of soil layers
                    f_obj % in_SS4HG % nLayers,        & ! intent(in):    number of layers
                    nState,                            & ! intent(in):    number of state variables in the current subset
                    .false.,                           & ! intent(in):    not inside Sundials solver
                    f_obj % in_SS4HG % firstSubStep,   & ! intent(in):    flag to indicate if we are processing the first sub-step
                    firstFluxCall,&!f_obj % io_SS4HG % firstFluxCall,  & ! intent(inout): flag to indicate if we are processing the first flux call
                    .true.,&!.false.,                           & ! intent(in):    flag to indicate if we are processing the first flux call in a splitting operation (.false. based on usage of eval8summa in summaSolve4homegrown)
                    f_obj % in_SS4HG % computeVegFlux, & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                    f_obj % in_SS4HG % scalarSolution, & ! intent(in):    flag to indicate the scalar solution
                    ! input: state vectors
                    stateVecTrial,                   & ! intent(in):    model state vector
                    fScale_split,                    & ! intent(in):    characteristic scale of the function evaluations
                    sMul_split,                      & ! intent(inout): state vector multiplier (used in the residual calculations)
                    ! input: data structures
                    f_obj % model_decisions,         & ! intent(in):    model decisions
                    f_obj % lookup_data,             & ! intent(in):    lookup tables
                    f_obj % type_data,               & ! intent(in):    type of vegetation and soil
                    f_obj % attr_data,               & ! intent(in):    spatial attributes
                    f_obj % mpar_data,               & ! intent(in):    model parameters
                    f_obj % forc_data,               & ! intent(in):    model forcing data
                    f_obj % bvar_data,               & ! intent(in):    average model variables for the entire basin
                    f_obj % prog_temp,               & ! intent(in):    model prognostic variables for a local HRU
                    ! input-output: data structures
                    indx_data,                       & ! intent(inout): index data
                    diag_data,                       & ! intent(inout): model diagnostic variables for a local HRU
                    flux_data,                       & ! intent(inout): model fluxes for a local HRU (initial flux structure)
                    deriv_data,                      & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                    ! input-output: baseflow
                    f_obj % io_SS4HG % ixSaturation, & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                    dBaseflow_dMatric,               & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                    ! output
                    f_obj % feasible,                & ! intent(out):   flag to denote the feasibility of the solution
                    fluxVec0_split,                  & ! intent(out):   flux vector
                    fRHS_split,                      & ! intent(out):   RHS function for ARKODE
                    rAdd_split,                      & ! intent(out):   additional (sink) terms on the RHS of the state equation
                    resVec_split,                    & ! intent(out):   residual vector
                    f_obj % out_SS4HG % fNew,        & ! intent(out):   function evaluation
                    f_obj % out_SS4HG % err,         & ! intent(out): error code
                    f_obj % out_SS4HG % message)       ! intent(out): error message (note: eval8summa uses "cmessage" instead)
  end associate

  ! finalize
  associate(err => f_obj % out_SS4HG % err, message => f_obj % out_SS4HG % message) 
   if (err /= 0) then
    if (f_obj % out_error) then
     write(f_obj % unit,*) "Error f_state_SUMMA_vec: eval8summa message="//trim(message); stop
    end if
   end if
  end associate

 end subroutine f_state_SUMMA_vec

 subroutine f_SUMMA_vec(f_obj,xvec)
  ! *** Compute SUMMA's vector non-linear function ***
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)            :: xvec(1:f_obj % n)  ! current guess

  ! compute SUMMA residual (taken to be the non-linear function) based on current guess
  ! note: - eval8summa may contain extraneous computations not needed for the residual
  !       - perhaps introducing logical flags in eval8summa to isolate the required operations would boost efficiency 
  call f_obj % SUMMA_eval8summa(xvec)

  f_obj % f_vec(:) = real(f_obj % resVec(:),r8b)
  
 end subroutine f_SUMMA_vec

 subroutine f1_SUMMA_vec(f_obj,xvec)
  ! *** Compute SUMMA's vector non-linear function ***
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)            :: xvec(1:f_obj % n)  ! current guess

  ! compute SUMMA residual (taken to be the non-linear function) based on current guess
  ! note: - eval8summa may contain extraneous computations not needed for the residual
  !       - perhaps introducing logical flags in eval8summa to isolate the required operations would boost efficiency 
  call f_obj % SUMMA_eval8summa(xvec)

  f_obj % f1_vec(:) = real(f_obj % resVec(:),r8b)
  
 end subroutine f1_SUMMA_vec

 subroutine f2_zero_vec(f_obj,xvec)
  ! *** Compute zero vector non-linear function ***
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)            :: xvec(1:f_obj % n)  ! current guess -- needed for argument interface

  f_obj % f2_vec(:) = 0._rkind
  
 end subroutine f2_zero_vec

 subroutine Jacobian_f_SUMMA_vec(f_obj,xvec)
  ! ** Compute SUMMA's Jacobian **
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)            :: xvec(1:f_obj % n) ! current guess
  ! local variables
  integer(i4b) :: nBands                                                      ! SUMMA's leading dimension for banded Jacobians
  real(rkind)  :: aJac(f_obj % in_SS4HG % nLeadDim,f_obj % in_SS4HG % nState) ! SUMMA's unscaled Jacobian matrix

  ! compute derivatives based on current guess
  !call f_obj % SUMMA_eval8summa(xvec) ! not required if f_SUMMA_vec(f_obj,xvec) has already been called

  ! assemble Jacobian using the computed derivatives
  call f_obj % SUMMA_computJacob(f_obj % indx_data,f_obj % diag_data,f_obj % flux_data,f_obj % deriv_data,&
                                &f_obj % dMat,f_obj % dBaseflow_dMatric,&
                                &aJac)

  ! store Jacobian used in solver
  if (f_obj % banded) then ! banded storage
   associate(nrow_banded => f_obj % nrow_banded, n => f_obj % n, subdiag => f_obj % subdiag)
    nBands=nrow_banded+subdiag
    f_obj % J(1:nrow_banded,1:n) = aJac(subdiag+1:nBands,1:n) ! aJac has extra storage rows
   end associate
  else ! full matrix storage
   f_obj % J(:,:) = aJac(:,:)
  end if

 end subroutine Jacobian_f_SUMMA_vec

 subroutine Jacobian_f1_SUMMA_vec(f_obj,xvec)
  ! ** Compute SUMMA's Jacobian **
  ! Note: for trivial decomposition with f2=0 such that f1=f
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)            :: xvec(1:f_obj % n) ! current guess
  ! local variables
  integer(i4b) :: nBands                                                      ! SUMMA's leading dimension for banded Jacobians
  real(rkind)  :: aJac(f_obj % in_SS4HG % nLeadDim,f_obj % in_SS4HG % nState) ! SUMMA's unscaled Jacobian matrix

  ! compute derivatives based on current guess
  !call f_obj % SUMMA_eval8summa(xvec) ! not required if f_SUMMA_vec(f_obj,xvec) has already been called

  ! assemble Jacobian using the computed derivatives
  call f_obj % SUMMA_computJacob(f_obj % indx_data,f_obj % diag_data,f_obj % flux_data,f_obj % deriv_data,&
                                &f_obj % dMat,f_obj % dBaseflow_dMatric,&
                                &aJac)

  ! store Jacobian used in solver
  if (f_obj % banded) then ! banded storage
   associate(nrow_banded => f_obj % nrow_banded, n => f_obj % n, subdiag => f_obj % subdiag)
    nBands=nrow_banded+subdiag
    f_obj % J1(1:nrow_banded,1:n) = aJac(subdiag+1:nBands,1:n) ! aJac has extra storage rows
   end associate
  else ! full matrix storage
   f_obj % J1(:,:) = aJac(:,:)
  end if

 end subroutine Jacobian_f1_SUMMA_vec

 subroutine Jacobian_f2_zero_vec(f_obj,xvec)
  ! ** Compute zero Jacobian **
  ! solver variables
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)            :: xvec(1:f_obj % n) ! current guess -- does not impact zero Jacobian but needed for argument interface

  ! store Jacobian used in solver
  f_obj % J2(:,:) = 0._rkind

 end subroutine Jacobian_f2_zero_vec

 subroutine Jacobian_f_mass_SUMMA_vec(f_obj,xvec)
  ! ** Compute SUMMA's Jacobian for mass **
  ! SJT: -------- incomplete and may not include contributions due to entire input state vector --------
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)            :: xvec(1:f_obj % n) ! current guess
  ! local variables
  integer(i4b) :: nBands                                                      ! SUMMA's leading dimension for banded Jacobians
  real(rkind)  :: aJac(f_obj % nLeadDim1,f_obj % nSubset1) ! SUMMA's unscaled Jacobian matrix

  ! compute derivatives based on current guess
  !call f_obj % SUMMA_eval8summa(xvec) ! not required if f_SUMMA_vec(f_obj,xvec) has already been called

  ! assemble Jacobian using the computed derivatives
  call f_obj % SUMMA_computJacob(f_obj % indx_data1,f_obj % diag_data1,f_obj % flux_data1,f_obj % deriv_data1,&
                                &f_obj % dMat1,f_obj % dBaseflow_dMatric1,&
                                &aJac)

  !! SJT: testing --- take out ---
  print *, "sum(aJac_mass)=",sum(aJac) 

! SJT: update the following taking proper sizes into account
! SJT: n probably refers to nSubset from the split (found from size(dMat) in computJacob)
!  ! store Jacobian used in solver
!  if (f_obj % banded) then ! banded storage
!   associate(nrow_banded => f_obj % nrow_banded, n => f_obj % n, subdiag => f_obj % subdiag)
!    nBands=nrow_banded+subdiag
!    f_obj % J1(1:nrow_banded,1:n) = aJac(subdiag+1:nBands,1:n) ! aJac has extra storage rows
!   end associate
!  else ! full matrix storage
!   f_obj % J1(:,:) = aJac(:,:)
!  end if

 end subroutine Jacobian_f_mass_SUMMA_vec

 subroutine Jacobian_f_energy_SUMMA_vec(f_obj,xvec)
  ! ** Compute SUMMA's Jacobian for mass **
  ! SJT: -------- incomplete and may not include contributions due to entire input state vector --------
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)            :: xvec(1:f_obj % n) ! current guess
  ! local variables
  integer(i4b) :: nBands                                                      ! SUMMA's leading dimension for banded Jacobians
  real(rkind)  :: aJac(f_obj % nLeadDim2,f_obj % nSubset2) ! SUMMA's unscaled Jacobian matrix

  ! compute derivatives based on current guess
  !call f_obj % SUMMA_eval8summa(xvec) ! not required if f_SUMMA_vec(f_obj,xvec) has already been called

  ! assemble Jacobian using the computed derivatives
  call f_obj % SUMMA_computJacob(f_obj % indx_data2,f_obj % diag_data2,f_obj % flux_data2,f_obj % deriv_data2,&
                                &f_obj % dMat2,f_obj % dBaseflow_dMatric2,&
                                &aJac)

  !! SJT: testing --- take out ---
  print *, "sum(aJac_energy)=",sum(aJac) 

! SJT: update the following taking proper sizes into account
! SJT: n probably refers to nSubset from the split (found from size(dMat) in computJacob)
!  ! store Jacobian used in solver
!  if (f_obj % banded) then ! banded storage
!   associate(nrow_banded => f_obj % nrow_banded, n => f_obj % n, subdiag => f_obj % subdiag)
!    nBands=nrow_banded+subdiag
!    f_obj % J2(1:nrow_banded,1:n) = aJac(subdiag+1:nBands,1:n) ! aJac has extra storage rows
!   end associate
!  else ! full matrix storage
!   f_obj % J2(:,:) = aJac(:,:)
!  end if

 end subroutine Jacobian_f_energy_SUMMA_vec

 subroutine Jacobian_f1_SUMMA_vec_numerical(f_obj,xvec)
  ! ** Compute SUMMA's Jacobian for mass **
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)            :: xvec(1:f_obj % n) ! current guess
  ! local variables
  integer(i4b) :: j ! loop index
  real(r8b) :: f1_j(1:f_obj % n)
  real(r8b) :: xvec_jp1(1:f_obj % n)
  real(r8b) :: delta ! step size for FD approximations

  ! store initial function vector (evaluated at xvec_j = xvec)
  f1_j(:) = f_obj % f1_vec(:)

  do j=1,f_obj % n

   ! set test x vector (perturb jth entry)
   xvec_jp1(:) = xvec(:)         ! note: xvec_j = xvec
   xvec_jp1(j) = 1.000001_r8b*xvec(j)+0.000001_r8b ! note: xvec_j = xvec -- avoid division by zero
   delta = xvec_jp1(j) - xvec(j)

   call f_obj % f1_vec_eval(xvec_jp1)

   !print *, "A0",xvec_jp1
   !print *, "B0",xvec
   !print *, "A",f_obj % f1_vec
   !print *, "B",f1_j
   f_obj % J1(:,j) = (f_obj % f1_vec(:) - f1_j(:)) / delta ! note: f1_jp1 = f_obj % f1_vec

  end do

  ! restore initial function vector (evaluated at xvec_j = xvec)
  f_obj % f1_vec(:) = f1_j(:)

  print *, "sum(J1_num)",sum(f_obj % J1)
 end subroutine Jacobian_f1_SUMMA_vec_numerical

 subroutine Jacobian_f2_SUMMA_vec_numerical(f_obj,xvec)
  ! ** Compute SUMMA's Jacobian for mass **
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)            :: xvec(1:f_obj % n) ! current guess
  ! local variables
  integer(i4b) :: j ! loop index
  real(r8b) :: f2_j(1:f_obj % n)
  real(r8b) :: xvec_jp1(1:f_obj % n)
  real(r8b) :: delta ! step size for FD approximations

  ! store initial function vector (evaluated at xvec_j = xvec)
  f2_j(:) = f_obj % f2_vec(:)

  do j=1,f_obj % n

   ! set test x vector (perturb jth entry)
   xvec_jp1(:) = xvec(:)         ! note: xvec_j = xvec
   xvec_jp1(j) = 1.000001_r8b*xvec(j)+0.000001_r8b ! note: xvec_j = xvec -- avoid division by zero
   delta = xvec_jp1(j) - xvec(j)

   call f_obj % f2_vec_eval(xvec_jp1)

   !print *, "A0",xvec_jp1
   !print *, "B0",xvec
   !print *, "A",f_obj % f2_vec
   !print *, "B",f2_j
   f_obj % J2(:,j) = (f_obj % f2_vec(:) - f2_j(:)) / delta ! note: f2_jp1 = f_obj % f2_vec

  end do

  ! restore initial function vector (evaluated at xvec_j = xvec)
  f_obj % f2_vec(:) = f2_j(:)

  print *, "sum(J2_num)",sum(f_obj % J2)
 end subroutine Jacobian_f2_SUMMA_vec_numerical

 subroutine Jacobian_f_SUMMA_vec_numerical(f_obj,xvec)
  ! ** Compute SUMMA's Jacobian for mass **
  ! arguments
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)            :: xvec(1:f_obj % n) ! current guess
  ! local variables
  integer(i4b) :: j ! loop index
  real(r8b) :: f_j(1:f_obj % n)
  real(r8b) :: xvec_jp1(1:f_obj % n)
  real(r8b) :: delta ! step size for FD approximations

  if (.not.allocated(f_obj % J)) allocate(f_obj % J(1:f_obj % n,1:f_obj % n))

  ! store initial function vector (evaluated at xvec_j = xvec)
  f_j(:) = f_obj % f_vec(:)

  do j=1,f_obj % n

   ! set test x vector (perturb jth entry)
   xvec_jp1(:) = xvec(:)         ! note: xvec_j = xvec
   xvec_jp1(j) = 1.000001_r8b*xvec(j)+0.000001_r8b ! note: xvec_j = xvec -- avoid division by zero
   delta = xvec_jp1(j) - xvec(j)

   call f_obj % f_vec_eval(xvec_jp1) ! call f_SUMMA_vec

   !print *, "A0",xvec_jp1
   !print *, "B0",xvec
   !print *, "A",f_obj % f2_vec
   !print *, "B",f2_j
   f_obj % J(:,j) = (f_obj % f_vec(:) - f_j(:)) / delta ! note: f_jp1 = f_obj % f_vec

  end do

  ! restore initial function vector (evaluated at xvec_j = xvec)
  f_obj % f_vec(:) = f_j(:)

  print *, "sum(J_num)",sum(f_obj % J)
 end subroutine Jacobian_f_SUMMA_vec_numerical

end module Newton_functions
