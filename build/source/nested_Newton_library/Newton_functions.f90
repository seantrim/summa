module Newton_functions
 ! nested Newton solver modules
 use kind_params,only: i4b,r8b ! kind parameters
 use Richards,only : Richards_obj ! Richards test problem
 ! SUMMA modules (for access to constant data and procedures)
 use nrtype,only: rkind,qp,lgt ! SUMMA's kind parameters (i4b is already used in kind_params module)
 use eval8summa_module, only: eval8summa                       ! SUMMA's eval8summa routine
 use computJacob_module,only: computJacob                      ! SUMMA's computJacob routine 
 use data_types,only: in_type_computJacob,out_type_computJacob ! objects for SUMMA's computJacob routine
 use data_types,only: in_type_summaSolve4homegrown,& ! objects for SUMMA's summaSolve4homegrown routine
                     &io_type_summaSolve4homegrown,&
                     &out_type_summaSolve4homegrown 
 use data_types,only: model_options           ! type for SUMMA's model decision structure
 use data_types,only: var_ilength,var_dlength ! derived types for SUMMA data structures
 use data_types,only: var_i,var_d             ! derived types for SUMMA data vectors
 use data_types,only: zLookup                 ! derived type for SUMMA lookup tables
 use var_lookup,only: iLookDECISIONS          ! named variables for elements of the SUMMA decision structure
 use mDecisions_module,only:qbaseTopmodel     ! SUMMA groundwater parameterization model decision
 implicit none
 private

 ! ***** Parent Type ***** !
 type, public :: f_obj_base
   ! ** Default data components used by the Newton solvers ** !
   logical      :: banded  ! flag for banded Jacobians
   logical      :: nested  ! flag for nested algorithm
   logical      :: verbose ! flag for full output
   logical      :: inner   ! flag to indicate the execution of inner iterations
   integer(i4b) :: subdiag,superdiag ! # of subdiagonals and superdiagonals for banded Jacobians
   integer(i4b) :: n                 ! vector size
   integer(i4b) :: kmax,lmax         ! max # of classical/outer and inner iterations
   integer(i4b) :: kcount,lcount     ! total # of classical/outer and inner iterations
   integer(i4b) :: unit              ! file unit number for solver output
   real(r8b),allocatable    :: x0(:),x1(:)   ! initial and final root estimates for vector algorithms
   real(r8b)                :: tol,tol_inner ! tolerance for classical/outer and inner iterations
   real(r8b)                :: R(-1:1)       ! max residual computed for iterations j-1, j, and j+1 (estimated)  
   real(r8b)                :: R_inner(-1:1) ! exact max residual computed for iterations j-1, j, and j+1 (estimated) 
   character(:),allocatable :: convergence   ! string for convergence control option
  contains
   ! procedures used prior to calling the solver
   procedure :: allocate_memory => f_allocate_memory ! allocate array data components 
   procedure :: initial_guess   => f_initial_guess   ! apply initial guess strategy
   procedure :: set_tolerance   => f_set_tolerance   ! set tolerances and iteration count maximums
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


   type(var_ilength) :: indx_data                    ! indices defining model states and layers
   type(var_dlength) :: prog_data                    ! prognostic variables for a local HRU
   type(var_dlength) :: diag_data                    ! diagnostic variables for a local HRU
   type(var_dlength) :: deriv_data                   ! derivatives in model fluxes w.r.t. relevant state variables
   real(rkind),allocatable :: dBaseflow_dMatric(:,:) ! derivative in baseflow w.r.t. matric head (s-1)
   real(rkind),allocatable :: dMat(:)                ! diagonal matrix (excludes flux derivatives) 

   type(in_type_summaSolve4homegrown)  :: in_SS4HG   ! summaSolve4homegrown input object: model control variables and previous function evaluation
   type(io_type_summaSolve4homegrown)  :: io_SS4HG   ! summaSolve4homegrown io object: model control variables and previous function evaluation
   type(out_type_summaSolve4homegrown) :: out_SS4HG  ! summaSolve4homegrown output object: model control variables and previous function evaluation

   ! additional variables for eval8summa call
   logical(lgt)            :: firstSplitOper         ! flag to indicate if we are processing the first flux call in a splitting operation
   real(rkind),allocatable :: fScale(:)              ! characteristic scale of the function evaluations (mixed units)
   real(qp),allocatable    :: sMul(:)    ! NOTE: qp  ! multiplier for state vector for the residual calculations
   logical(lgt) :: feasible                          ! feasibility flag
   real(rkind),allocatable :: fluxVec0(:)            ! flux vector (mixed units)
   real(rkind),allocatable :: fRHS(:)                ! RHS function for ARKODE
   real(rkind),allocatable :: rAdd(:)                ! additional terms in the residual vector
   real(qp),allocatable    :: resVec(:)  ! NOTE: qp  ! residual vector 
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
   procedure :: f_vec  => f_SUMMA_vec ! solver
   !procedure :: f_vec  => f_diff_vec  ! solver
   procedure :: f1_vec => f1_Rich_vec ! solver
   procedure :: f2_vec => f2_Rich_vec ! solver
   procedure :: dfdx_vec  => dfdx_diff_vec 
   procedure :: df1dx_vec => df1_Rich_dh_vec
   procedure :: df2dx_vec => df2_Rich_dh_vec
   procedure :: J  => Jacobian_f_SUMMA_vec  ! solver
   !procedure :: J  => Jacobian_f_Rich_vec  ! solver
   procedure :: J1 => Jacobian_f1_Rich_vec ! solver
   procedure :: J2 => Jacobian_f2_Rich_vec ! solver
   
   ! scalar routines
   procedure :: f     => f_diff 
   procedure :: dfdx  => dfdx_diff 
   procedure :: df1dx => df1_Rich_dh 
   procedure :: df2dx => df2_Rich_dh 
 end type f_obj_type

contains

!!!!!!!!!! ******************* User defined functions below ******************* !!!!!!!!!!

 ! **** Utilities **** !
 
 subroutine f_allocate_memory(f_obj)
  ! ** allocate array data components for f_obj_base class **
  class(f_obj_base),intent(inout) :: f_obj

  associate(n => f_obj % n)
   allocate(f_obj % x0(1:n),f_obj % x1(1:n))
  end associate
 end subroutine f_allocate_memory

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
    write(f_obj % unit,'(a65)') "Error in f_set_tolerance: method argument not currently supported"
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
   write(f_obj % unit,'(a66)') "Error in f_initial_guess: method argument not currently supported."
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

 function Jacobian_f1_Rich_vec(f_obj,xvec) result(J1)
  class(f_obj_type),intent(in) :: f_obj
  real(r8b),intent(in)         :: xvec(1:f_obj % n) ! current guess
  real(r8b),allocatable        :: J1(:,:)
  !real(r8b)                    :: J1(1:f_obj % n,1:f_obj % n)
  integer(i4b)                 :: icol,irow
  integer(i4b)                 :: nrow_banded ! # of rows for LAPACK banded matrix storage

  if (f_obj % banded) then ! banded storage
   associate(n => f_obj % n, subdiag => f_obj % subdiag, superdiag => f_obj % superdiag)
    nrow_banded=subdiag+superdiag+1
    allocate(J1(1:nrow_banded,1:n))
    do icol=1,n
     do irow=max(1,icol-superdiag),min(n,icol+subdiag)
      Richards_obj % i = irow
      J1(superdiag+1+irow-icol,icol)=f_obj % df1dx_vec(xvec,icol)
     end do
    end do 
   end associate
  else ! full matrix storage
   associate(n => f_obj % n)
    allocate(J1(1:n,1:n))
    do icol=1,n
     do irow=1,n
      Richards_obj % i = irow
      J1(irow,icol)=f_obj % df1dx_vec(xvec,icol)
     end do
    end do 
   end associate
  end if
 end function Jacobian_f1_Rich_vec

 function Jacobian_f2_Rich_vec(f_obj,xvec) result(J2)
  class(f_obj_type),intent(in) :: f_obj
  real(r8b),intent(in)         :: xvec(1:f_obj % n) ! current guess
  real(r8b),allocatable        :: J2(:,:)
  !real(r8b)                    :: J2(1:f_obj % n,1:f_obj % n)
  integer(i4b)                 :: icol,irow
  integer(i4b)                 :: nrow_banded ! # of rows for LAPACK banded matrix storage

  if (f_obj % banded) then ! banded storage
   associate(n => f_obj % n, subdiag => f_obj % subdiag, superdiag => f_obj % superdiag)
    nrow_banded=subdiag+superdiag+1
    allocate(J2(1:nrow_banded,1:n))
    do icol=1,n
     do irow=max(1,icol-superdiag),min(n,icol+subdiag)
      Richards_obj % i = irow
      J2(superdiag+1+irow-icol,icol)=f_obj % df2dx_vec(xvec,icol)
     end do
    end do 
   end associate
  else ! full matrix storage
   associate(n => f_obj % n)
    allocate(J2(1:n,1:n))
    do icol=1,n
     do irow=1,n
      Richards_obj % i = irow
      J2(irow,icol)=f_obj % df2dx_vec(xvec,icol)
     end do
    end do 
   end associate
  end if
 end function Jacobian_f2_Rich_vec

 function f_diff_vec(f_obj,xvec) result(f_vec)
  ! *** form non-linear vector function using the Jordan decomposition ***
  class(f_obj_type),intent(in) :: f_obj
  real(r8b),intent(in)         :: xvec(1:f_obj % n) ! current guess
  real(r8b)                    :: f_vec(1:f_obj % n) ! non-linear function vector

  f_vec=f_obj % f1_vec(xvec)- f_obj % f2_vec(xvec)
 end function f_diff_vec

 function f1_Rich_vec(f_obj,xvec) result(f1_vec)
  class(f_obj_type),intent(in) :: f_obj
  real(r8b),intent(in)         :: xvec(1:f_obj % n) ! current guess
  real(r8b)                    :: f1_vec(1:f_obj % n) ! non-linear function vector
  integer(i4b)                 :: i

  associate(n => f_obj % n)
   do i=1,n ! interior grid points
    Richards_obj % i = i
    ! populate h array with current guess on ith stencil
    if (i.ne.1) Richards_obj % h(i-1) = xvec(i-1) ! BC 
                Richards_obj % h(i) = xvec(i) 
    if (i.ne.n) Richards_obj % h(i+1) = xvec(i+1) ! BC 
    f1_vec(i)=f_obj % f1(xvec(i))
   end do
  end associate
 end function f1_Rich_vec

 function f2_Rich_vec(f_obj,xvec) result(f2_vec)
  class(f_obj_type),intent(in) :: f_obj
  real(r8b),intent(in)         :: xvec(1:f_obj % n) ! current guess
  real(r8b)                    :: f2_vec(1:f_obj % n) ! non-linear function vector
  integer(i4b)                 :: i

  associate(n => f_obj % n)
   do i=1,n ! interior grid points
    Richards_obj % i = i
    ! populate h array with current guess on ith stencil
    if (i.ne.1) Richards_obj % h(i-1) = xvec(i-1) ! BC 
                Richards_obj % h(i) = xvec(i) 
    if (i.ne.n) Richards_obj % h(i+1) = xvec(i+1) ! BC 
    f2_vec(i)=f_obj % f2(xvec(i))
   end do
  end associate
 end function f2_Rich_vec

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

 subroutine SUMMA_eval8summa(f_obj,xvec)
  ! ** interface for SUMMA's eval8summa subroutine **
  ! compute SUMMA derivative values and residual vector
  ! note: - eval8summa was not refactored to use object arguments
  !       - objects for summaSolve4homegrown were reused where possible
  class(f_obj_inputs),intent(inout) :: f_obj
  real(r8b),intent(in)         :: xvec(1:f_obj % n) ! current guess

  ! update
  associate(&
   stateVecTrial => xvec & ! current guess for state vector
  &)
   call eval8summa(&
                    ! input: model control
                    f_obj % in_SS4HG % dt_cur,                  & ! intent(in):    current stepsize
                    f_obj % in_SS4HG % dt,                      & ! intent(in):    length of the entire time step (seconds) for drainage pond rate
                    f_obj % in_SS4HG % nSnow,                   & ! intent(in):    number of snow layers
                    f_obj % in_SS4HG % nSoil,                   & ! intent(in):    number of soil layers
                    f_obj % in_SS4HG % nLayers,                 & ! intent(in):    number of layers
                    f_obj % in_SS4HG % nState,                  & ! intent(in):    number of state variables in the current subset
                    .false.,                 & ! intent(in):    not inside Sundials solver
                    f_obj % in_SS4HG % firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                    f_obj % io_SS4HG % firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                    f_obj % firstSplitOper,  & ! intent(in):    flag to indicate if we are processing the first flux call in a splitting operation
                    f_obj % in_SS4HG % computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                    f_obj % in_SS4HG % scalarSolution,          & ! intent(in):    flag to indicate the scalar solution
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
                    f_obj % flux_init,               & ! intent(inout): model fluxes for a local HRU (initial flux structure)
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
                    f_obj % in_SS4HG % fOld,         & ! intent(out):   function evaluation
                    f_obj % out_SS4HG % err,         & ! intent(out): error code
                    f_obj % out_SS4HG % message)       ! intent(out): error message (note: eval8summa uses "cmessage" instead)
  end associate

  ! finalize
  ! note: "message" used for out_SS4HG data component but "cmessage" used within summaSolve4homegrown subroutine
  associate(err => f_obj % out_SS4HG % err, cmessage => f_obj % out_SS4HG % message) 
   if (err /= 0) then
    write(f_obj % unit,*) "Error in SUMMA_eval8summa: eval8summa message="//trim(cmessage); stop
   end if
  end associate
 end subroutine SUMMA_eval8summa

 subroutine SUMMA_computJacob(f_obj,J)
  ! ** Interface for SUMMA's computJacob subroutine **
  ! solver variables
  class(f_obj_inputs),intent(inout) :: f_obj
  real(r8b),allocatable,intent(out) :: J(:,:)
  integer(i4b)                 :: nrow_banded ! # of rows for LAPACK banded matrix storage
  ! SUMMA variables
  type(in_type_computJacob)    :: in_computJacob  ! computJacob input object
  type(out_type_computJacob)   :: out_computJacob ! computJacob output object  


  ! memory allocation for Jacobian 
  if (f_obj % banded) then ! banded storage
   associate(n => f_obj % n, subdiag => f_obj % subdiag, superdiag => f_obj % superdiag)
    nrow_banded=subdiag+superdiag+1
    allocate(J(1:nrow_banded,1:n))
   end associate
  else ! full matrix storage
   associate(n => f_obj % n)
    allocate(J(1:n,1:n))
   end associate
  end if

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
    indx_data         => f_obj % indx_data,&         ! indices defining model states and layers
    prog_data         => f_obj % prog_data,&         ! prognostic variables for a local HRU
    diag_data         => f_obj % diag_data,&         ! diagnostic variables for a local HRU
    deriv_data        => f_obj % deriv_data,&        ! derivatives in model fluxes w.r.t. relevant state variables
    dBaseflow_dMatric => f_obj % dBaseflow_dMatric,& ! derivative in baseflow w.r.t. matric head (s-1)
    dMat              => f_obj % dMat,&              ! diagonal matrix (excludes flux derivatives) 
    aJac              => J &                         ! Jacobian
   &)
    call computJacob(in_computJacob,indx_data,prog_data,diag_data,deriv_data,dBaseflow_dMatric,dMat,aJac,out_computJacob)
   end associate

  ! finalize
  ! *** Transfer data from out_computJacob class object to local variables in summaSolve4homegrown ***
  ! note: "message" used for out_SS4HG data component but "cmessage" used within summaSolve4homegrown subroutine
  associate(err => f_obj % out_SS4HG % err, cmessage => f_obj % out_SS4HG % message) 
   call out_computJacob % finalize(err,cmessage)
   if (err /= 0) then
    write(f_obj % unit,*) "Error in Jacobian_f_SUMMA_vec: computJacob message="//trim(cmessage); stop
   end if
  end associate

 end subroutine SUMMA_computJacob

 function f_SUMMA_vec(f_obj,xvec) result(f_vec)
  ! *** Compute SUMMA's vector non-linear function ***
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)         :: xvec(1:f_obj % n) ! current guess
  real(r8b)                    :: f_vec(1:f_obj % n) ! non-linear function vector

  ! compute SUMMA residual (taken to be the non-linear function) based on current guess
  ! note: - eval8summa may contain extraneous computations not needed for the residual
  !       - perhaps introducing logical flags in eval8summa to isolate the required operations would boost efficiency 
  call f_obj % SUMMA_eval8summa(xvec)

  f_vec=real(f_obj % resVec,r8b)
 end function f_SUMMA_vec

 function Jacobian_f_SUMMA_vec(f_obj,xvec) result(J)
  ! ** Compute SUMMA's Jacobian **
  ! solver variables
  class(f_obj_type),intent(inout) :: f_obj
  real(r8b),intent(in)         :: xvec(1:f_obj % n) ! current guess
  real(r8b),allocatable        :: J(:,:)

  ! compute derivatives based on current guess
  call f_obj % SUMMA_eval8summa(xvec)

  ! assemble Jacobian using the computed derivatives
  call f_obj % SUMMA_computJacob(J)

 end function Jacobian_f_SUMMA_vec

end module Newton_functions
