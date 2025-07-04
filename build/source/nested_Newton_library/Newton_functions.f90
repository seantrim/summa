module Newton_functions
 use kind_params,only: i4b,r8b
 use Richards,only : Richards_obj ! Richards test problem
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

 type,extends(f_obj_base),public :: f_obj_input_functions
  contains
   ! ** routines that point to external sources ** !
   ! note: - these procedures are not directly called in the solver
   !       - however, these procedures may be called within procedures that are called in the solver

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
 end type f_obj_input_functions

 type,extends(f_obj_input_functions),public :: f_obj_type
  contains
   ! *** these procedures take the procedures from f_obj_input_functions type as input *** !
   ! vector routines
   procedure :: f_vec  => f_diff_vec  ! solver
   procedure :: f1_vec => f1_Rich_vec ! solver
   procedure :: f2_vec => f2_Rich_vec ! solver
   procedure :: dfdx_vec  => dfdx_diff_vec 
   procedure :: df1dx_vec => df1_Rich_dh_vec
   procedure :: df2dx_vec => df2_Rich_dh_vec
   procedure :: J  => Jacobian_f_Rich_vec  ! solver
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
   print *, "Error in f_initial_guess: method argument not currently supported."
  end if
 end subroutine f_initial_guess

 
 ! **** Richards Problem **** !

 real(r8b) function f_Rich_space(f_obj,x) result(f_space)
  ! ** space terms for discrete Richards' equation **
  class(f_obj_input_functions),intent(in) :: f_obj
  real(r8b),intent(in)         :: x ! current guess

  associate(i => Richards_obj % i)
    Richards_obj % h(i) = x ! populate h array with current guess 
    f_space=Richards_obj % KT()-Richards_obj % S()
  end associate
 end function f_Rich_space

 real(r8b) function f_Rich_time(f_obj,x) result(f_time)
  ! ** time terms for discrete Richards' equation **
  class(f_obj_input_functions),intent(in) :: f_obj
  real(r8b),intent(in)         :: x ! current guess

  associate(i => Richards_obj % i)
    Richards_obj % h(i) = x ! populate h array with current guess 
    f_time=Richards_obj % CT()
  end associate
 end function f_Rich_time

 real(r8b) function df_Rich_dh_element_space(f_obj,x,j) result(dfdh_element_space)
  class(f_obj_input_functions),intent(in) :: f_obj
  real(r8b),intent(in)         :: x ! current guess
  integer(i4b),intent(in)      :: j

  associate(i => Richards_obj % i)
    Richards_obj % h(i) = x ! populate h array with current guess 
    dfdh_element_space=Richards_obj % dKTdh(j)
  end associate
 end function df_Rich_dh_element_space

 real(r8b) function df_Rich_dh_element_time(f_obj,x,j) result(dfdh_element_time)
  class(f_obj_input_functions),intent(in) :: f_obj
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
  real(r8b)                    :: x ! current guess
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
  real(r8b)                    :: x ! current guess
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
  real(r8b)                    :: x ! current guess
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
  ! *** form vector objective function using the Jordan decomposition ***
  class(f_obj_type),intent(in) :: f_obj
  real(r8b),intent(in)         :: xvec(1:f_obj % n) ! current guess
  real(r8b)                    :: f_vec(1:f_obj % n) ! objective function vector
  real(r8b)                    :: x
  integer(i4b)                 :: i

  f_vec=f_obj % f1_vec(xvec)- f_obj % f2_vec(xvec)
 end function f_diff_vec

 function f1_Rich_vec(f_obj,xvec) result(f1_vec)
  class(f_obj_type),intent(in) :: f_obj
  real(r8b),intent(in)         :: xvec(1:f_obj % n) ! current guess
  real(r8b)                    :: f1_vec(1:f_obj % n) ! objective function vector
  real(r8b)                    :: x
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
  real(r8b)                    :: f2_vec(1:f_obj % n) ! objective function vector
  real(r8b)                    :: x
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
  ! ** complete scalar objective function from Jordan decomposition **
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

end module Newton_functions
