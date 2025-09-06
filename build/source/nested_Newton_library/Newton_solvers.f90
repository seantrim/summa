module Newton_solvers
 ! ****** Classical and Nested Newton Solvers ******
 ! Note: objective functions and derivatives are specified via classes in the Newton_functions module
 use kind_params, only: r8b,i4b
 use Newton_functions, only: f_obj_type
 implicit none
 private
 public :: Newton_solve

contains

 subroutine Newton_solve(f_obj)
  type(f_obj_type),intent(inout) :: f_obj 

  if (f_obj % n .gt. 0_i4b) then
   if (f_obj % nested) then
    call nested_Newton_vector(f_obj)   
   else
    call Newton_vector(f_obj)
   end if
  else
   if (f_obj % out_error) then
    write(f_obj % unit,*) "Error in Newton_solve: problem size is not valid."
   end if
   stop
  end if
 end subroutine Newton_solve

 subroutine Newton_vector(f_obj)
  ! Newton solver for vector problems
  type(f_obj_type),intent(inout) :: f_obj 
  real(r8b) :: xk(1:f_obj % n),xkp1(1:f_obj % n) ! x^k and x^{k+1}
  real(r8b),allocatable :: J(:,:)                ! Jacobian
  real(r8b) :: fval(1:f_obj % n)                 ! for storing the current function value
  real(r8b) :: final_mean                        ! mean of final solution vector
  real(r8b) :: R_est                             ! estimated max relative difference in solution between iterations
  integer(i4b) :: k                              ! iteration counter
  logical :: exit_flag                           ! exit flag
  ! LAPACK Variables
  real(r8b),allocatable :: A(:,:)                ! input and result matrix
  real(r8b) :: B(1:f_obj % n)                    ! right-hand side / solution vector
  integer(i4b) :: M                              ! # of rows/columns for linear system
  integer(i4b) :: nrow_banded                    ! # of rows for banded storage

  ! initialize convergence flag
  f_obj % converged = .false.

  ! allocate memory for choice of Jacobian storage
  if (f_obj % banded) then ! banded storage
   nrow_banded=f_obj % subdiag + f_obj % superdiag + 1
   allocate(J(1:nrow_banded,1:f_obj % n),A(1:nrow_banded,1:f_obj % n))
  else ! full matrix storage
   allocate(J(1:f_obj % n,1:f_obj % n),A(1:f_obj % n,1:f_obj % n))
  end if

  ! initialize LAPACK parameters
  M=f_obj % n 

  f_obj % inner = .false. ! classical iterations only
  exit_flag=.false.
  xk=f_obj % x0 ! initialize
  do k=0,f_obj % kmax
   fval=f_obj % f_vec(xk)
   J=f_obj % J(xk) ! compute Jacobian

   ! begin LAPACK operations
   A=J     ! initialize matrix used by LAPACK
   B=-fval ! initialize right-side vector used by LAPACK
   call linear_solve(f_obj,M,A,B,f_obj % tol) ! Solve Ax=B -- x stored in B on output -- M is the # of rows/columns of A

   xkp1=xk+B ! update guess

   call check_residual_vector(f_obj,k,xkp1,xk,f_obj % tol,R_est,exit_flag)
   if (f_obj % out_detail) write(f_obj % unit,'(i4,3(g23.15))') k,sum(xk)/f_obj % n,f_obj % R(0),R_est
   if (exit_flag) then ! exit loop if convergence criterion is met
    f_obj % converged = .true.
    exit
   end if

   if (f_obj % constraints) call f_obj % apply_constraints(xk,xkp1) ! apply constraints without interfering with the convergence criterion
   xk=xkp1 ! prep for next iteration - can probably evaluate in place
  end do
  ! final output
  final_mean=sum(xkp1)/f_obj % n 
  if (exit_flag.eqv..true.) then
   if (f_obj % out_detail) write(f_obj % unit,'(i4,3(g23.15))') k+1,final_mean,f_obj % R(1)
   f_obj % kcount=k+1
  else
   f_obj % kcount=k
  end if

  if (k.gt.f_obj % kmax) then
   if (f_obj % out_warning) then
    write(f_obj % unit,*) "Warning - classical Newton solver has reached the maximum number of iterations&
                          & - accuracy may not be sufficient."
   end if
  end if

  f_obj % x1=xkp1
  !write(f_obj % unit,*) "Mean Solution=",final_mean
  if (f_obj % out_basic) write(f_obj % unit,*) "Convergence Error=",f_obj % R(1)

 end subroutine Newton_vector

 subroutine nested_Newton_vector(f_obj)
  ! Newton solver
  type(f_obj_type),intent(inout) :: f_obj 
  real(r8b) :: xk0(1:f_obj % n),xkp1l(1:f_obj % n),xkp1lp1(1:f_obj % n) ! x^k, x^{k+1,l} and x^{k+1,l+1} 
  real(r8b),allocatable :: Jdiff(:,:),J1(:,:),J2(:,:)
  real(r8b) :: f1val(1:f_obj % n),f2val(1:f_obj % n)  ! for storing the current function value
  real(r8b) :: final_mean                        ! mean value of final solution vector
  real(r8b) :: R_est                             ! estimated max relative difference in solution between iterations
  integer(i4b) :: k,l                            ! iteration counters
  integer(i4b) :: l_total                        ! total number of inner iterations
  logical :: exit_outer,exit_inner 
  ! LAPACK Variables
  real(r8b),allocatable :: A(:,:)                ! input and result matrix
  real(r8b) :: B(1:f_obj % n)                    ! right-hand side / solution vector
  integer(i4b) :: M                              ! # of rows/columns for linear system
  integer(i4b) :: nrow_banded                    ! # of rows for banded storage

  ! initialize convergence flag
  f_obj % converged = .false.

  ! allocate memory for choice of Jacobian storage
  if (f_obj % banded) then ! banded storage
   nrow_banded=f_obj % subdiag + f_obj % superdiag + 1
   allocate(Jdiff(1:nrow_banded,1:f_obj % n),J1(1:nrow_banded,1:f_obj % n),J2(1:nrow_banded,1:f_obj % n))
   allocate(A(1:nrow_banded,1:f_obj % n))
  else ! full matrix storage
   allocate(Jdiff(1:f_obj % n,1:f_obj % n),J1(1:f_obj % n,1:f_obj % n),J2(1:f_obj % n,1:f_obj % n))
   allocate(A(1:f_obj % n,1:f_obj % n))
  end if

  ! initialize LAPACK parameters
  M=f_obj % n
 
  l_total=0
  exit_outer=.false.
  xk0=f_obj % x0 ! initial guess
  outer: do k=0,f_obj % kmax
   J2=f_obj % J2(xk0) ! compute Jacobian
   exit_inner=.false.
   xkp1l=xk0 !initial guess for inner iterations
   f_obj % inner=.true. ! inner iterations for next loop
   inner: do l=0,f_obj % lmax ! inner iterations
    J1=f_obj % J1(xkp1l) ! compute Jacobian
    Jdiff=J1-J2
    f1val = f_obj % f1_vec(xkp1l)
    f2val = f_obj % f2_vec(xk0)

    ! begin LAPACK operations
    A=Jdiff ! initialize matrix used by LAPACK
    ! initialize right-side vector used by LAPACK
    B=f2val-matrix_vector_product(f_obj,M,J2,xk0)-f1val+matrix_vector_product(f_obj,M,J1,xkp1l) 
!    B=f2val-matmul(J2,xk0)-f1val+matmul(J1,xkp1l) ! initialize right-side vector used by LAPACK
    call linear_solve(f_obj,M,A,B,f_obj % tol) ! Solve Ax=B -- x stored in B on output -- M is the # of rows/columns of A
    xkp1lp1=B ! update guess

    call check_residual_vector(f_obj,l,xkp1lp1,xkp1l,f_obj % tol_inner,R_est,exit_inner)
    ! print exact convergence error
    if (f_obj % out_detail) write(f_obj % unit,'(a2,i4,3(g23.15))') "  ",l,sum(xkp1l)/f_obj % n,f_obj % R_inner(0),R_est
    if (exit_inner) exit inner
    xkp1l=xkp1lp1 ! set up next inner iteration
   end do inner
   if (l.gt.f_obj % lmax) then
    if (f_obj % out_warning) then
     write(f_obj % unit,*) "Warning - nested Newton solver has reached the maximum number of inner iterations&
                           & - accuracy may not be sufficient."
    end if
   end if
   ! inner iteration counts 
   if (exit_inner) then
    ! final output for inner iterations
    if (f_obj % out_detail) write(f_obj % unit,'(a2,i4,2(g23.15))') "  ",l+1,sum(xkp1lp1)/f_obj % n,f_obj % R_inner(1) 
    l_total=l_total+(l+1)
   else
    l_total=l_total+l
   end if

   f_obj % inner=.false.
   call check_residual_vector(f_obj,k,xkp1lp1,xk0,f_obj % tol,R_est,exit_outer)
   if (f_obj % out_detail) then ! convergence error info for iteration k
    write(f_obj % unit,'(i4,3(g23.15))') k,sum(xk0)/f_obj % n,f_obj % R(0),R_est 
   end if
   if (exit_outer) then ! exit loop if convergence criterion is met
    f_obj % converged = .true.
    exit outer
   end if

   if (f_obj % constraints) call f_obj % apply_constraints(xk0,xkp1lp1) ! apply constraints without interfering with the convergence criterion
   xk0=xkp1lp1
  end do outer

  ! outer iterations counts 
  final_mean=sum(xkp1lp1)/f_obj % n
  if (exit_outer) then
   ! final convergence error for outer iterations (if early loop exit occurred)
   if (f_obj % out_detail) then
    write(f_obj % unit,'(i4,2(g23.15))') k+1,final_mean,f_obj % R(1) ! mean of final solution 
   end if
   f_obj % kcount = k+1
  else
   f_obj % kcount = k
  end if

  if (k.gt.f_obj % kmax) then
   if (f_obj % out_warning) then
    write(f_obj % unit,*) "Warning - nested Newton solver has reached the maximum number of outer iterations&
                          & - accuracy may not be sufficient."
   end if
  end if

  f_obj % x1=xkp1lp1
  !write(f_obj % unit,*) "Mean Solution=",final_mean
  if (f_obj % out_basic) write(f_obj % unit,*) "Convergence Error=",f_obj % R(1)
  f_obj % lcount = l_total
  if (f_obj % out_detail) then
   write(f_obj % unit,*) "# of outer iterations=",f_obj % kcount
   write(f_obj % unit,*) "# of inner iterations=",l_total
  end if
 end subroutine nested_Newton_vector

 subroutine check_residual_vector(f_obj,iteration,xkp1,xk,tol,R_est,exit_flag)
  ! *** Check residual vector for potential loop exit ***
  type(f_obj_type),intent(inout) :: f_obj 
  integer(i4b),intent(in) :: iteration  ! interation count
  real(r8b),intent(in)    :: xkp1(1:f_obj % n)  ! current root estimate
  real(r8b),intent(in)    :: xk(1:f_obj % n)    ! previous root estimate
  real(r8b),intent(in)    :: tol        ! tolerance
  logical,intent(inout)   :: exit_flag  ! exit flag
  real(r8b),intent(out)   :: R_est      ! estimated R for current iteration (computed in the previous call)
  ! local variables
  real(r8b)               :: R(-1:1)    ! maximum residual array (two previous exact values and prediction for next iteration)
  integer(i4b)            :: i                  ! index for residual vector
  real(r8b)               :: R_vec(1:f_obj % n) ! residual vector
  real(r8b)               :: b                  ! exponent used for convergence error estimation 

  do i=1,f_obj % n
   if (xk(i).ne.0._r8b) then
    R_vec(i)=abs((xkp1(i)-xk(i))/xk(i))
   else if (xkp1(i).ne.0._r8b) then
    R_vec(i)=abs(xkp1(i)-xk(i)) ! avoid residuals of unity (since xk(i) equals zero)
   else
    R_vec(i)=0._r8b ! both xk and xkp1 are zero -- set the residual to zero
   end if
  end do

  ! store previous residuals
  if (iteration.gt.0) then
   if (f_obj % inner) then
    R_est=f_obj % R_inner(1) ! store previous estimate for reference
    f_obj % R_inner(-1)=f_obj % R_inner(0)
    R(-1)=f_obj % R_inner(-1) ! exact residual for iteration-1
   else
    R_est=f_obj % R(1) ! store previous estimate for reference
    f_obj % R(-1)=f_obj % R(0)
    R(-1)=f_obj % R(-1)       ! exact residual for iteration-1
   end if
  end if

  R(0)=maxval(R_vec) ! actual worst case residual for input iteration
  if (f_obj % convergence.eq.'strict') then ! strict estimate
   R(1)=R(0) ! estimated residual for iteration+1
  else if (f_obj % convergence.eq.'predictive') then
   if ((iteration.eq.0)) then ! initial prediction is conservative due to lack of information
    R(1)=R(0) ! estimated residual for iteration+1   
   else ! compute prediction based on power function
    b=log10(R(0)/R(-1)) ! exponent
    R(1)=R(0)*10**b     ! power function -- estimated residual for iteration+1
   end if
  else ! method not valid
   if (f_obj % out_error) then
    write(f_obj % unit,*) "Error in check_residual_vector: method argument not currently supported."
   end if
   stop
  end if

  if (f_obj % inner) then ! inner iterations
   f_obj % R_inner(0) = R(0) ! store exact residual for current iteration
   f_obj % R_inner(1) = R(1) ! store estimated residual for iteration+1
  else                    ! outer/classical iterations
   f_obj % R(0) = R(0)       ! store exact residual for current iteration
   f_obj % R(1) = R(1)       ! store estimated residual for iteration+1
  end if
  if (iteration.eq.0) R_est=R(1) ! initialize R_est for iteration zero

  ! check exact error from current iteration and estimated error for next iteration
  if ((R(0).lt.tol).or.(R(1).lt.tol)) then
   exit_flag=.true.; return  ! set exit flag if criterion is satisfied
  end if
 end subroutine check_residual_vector

 function matrix_vector_product(f_obj,M,A,x) result(y)
  ! *** Compute matrix vector product y=A*x ***
  ! input
  type(f_obj_type),intent(in) :: f_obj        ! class object containing solver options
  integer(i4b),intent(in) :: M                ! # of rows/columns for linear system
  real(r8b),allocatable,intent(in) :: A(:,:)  ! input matrix (result is LU factorization of scaled matrix)
  real(r8b),intent(in) :: x(1:M)              ! input vector

  ! output
  real(r8b) :: y(1:M)                         ! product vector

  ! local variables
  character(1),parameter :: TRANS='N'                ! option for matrix transposition
  integer(i4b) :: N                                  ! # of columns of A
  integer(i4b) :: KL,KU                              ! # of subdiagonals and superdiagonals of A (banded storage)
  integer(i4b) :: LDA                                ! first dimension of A
  integer(i4b),parameter :: INCX=1_i4b, INCY=1_i4b   ! increment for elements of x and y vectors
  real(r8b),parameter    :: ALPHA=1._r8b,BETA=0._r8b ! scalars used in LAPACK solvers
  
  ! set general LAPACK parameters
  N=M
 
  ! set LAPACK parameters for choice of matrix storage
  if (f_obj % banded) then ! banded storage
   KL=f_obj % subdiag; KU=f_obj % superdiag 
   LDA=KL + KU + 1
  else ! full matrix storage
   LDA=M
  end if
  
  if (f_obj % banded) then ! banded storage
   call DGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,x,INCX,BETA,y,INCY) ! BLAS
  else ! full matrix storage
   !y=matmul(A,x)
   call DGEMV(TRANS,M,N,ALPHA,A,LDA,x,INCX,BETA,y,INCY) ! BLAS
  end if
 end function matrix_vector_product

 subroutine linear_solve(f_obj,M,A,B,tol)
  ! *** Solve Ax=B -- x stored in B on output -- M is the # of rows/columns of A *** 
  type(f_obj_type),intent(in) :: f_obj           ! class object containing solver options
  ! LAPACK Variables
  integer(i4b),intent(in) :: M                   ! # of rows/columns for linear system
  real(r8b),allocatable,intent(inout) :: A(:,:)  ! input and result matrix (result is LU factorization of scaled matrix)
  !real(r8b),intent(inout) :: A(1:M,1:M)          ! input and result matrix (result is LU factorization of scaled matrix)
  real(r8b),intent(inout) :: B(1:M)              ! right-hand side / solution vector
  real(r8b),intent(in)    :: tol                 ! tolerance value used by the calling routine
  ! local variables
  character(1) :: FACT                           ! option for matrix factoring
  character(1) :: TRANS                          ! option for matrix transposition
  character(1) :: EQUED                          ! specifies equilibration type
  integer(i4b) :: N                              ! # of matrix columns
  integer(i4b) :: KL,KU                          ! # of subdiagonals and superdiagonals
  integer(i4b) :: LDA,LDAF,LDX,LDB               ! leading dimensions of A, AF, X, and B arrays
  integer(i4b) :: NRHS                           ! # of right-hand sides in B vector
  integer(i4b) :: INFO                           ! error code
  integer(i4b) :: IPIV(1:M)                      ! pivot index vector
  integer(i4b) :: IWORK(1:M)                     ! work integer array
  real(r8b) :: RA(1:M),CA(1:M)                   ! row and column scale factors for A
  real(r8b) :: RCOND                             ! estimate of condition number reciprocal
  real(r8b),allocatable :: AF(:,:)               ! output matrix (result is LU factorization of scaled matrix)
  !real(r8b) :: AF(1:M,1:M)                       ! output matrix (result is LU factorization of scaled matrix)
  real(r8b) :: X(1:M)                            ! solution to original (unscaled) system
  real(r8b) :: FERR(1:1),BERR(1:1)               ! forward and backward error estimates (single right-hand side assumed)
  real(r8b),allocatable :: WORK(:)               ! work array (reciprocal pivot growth factor in work(1) on exit)
  !real(r8b) :: WORK(1:4_i4b*M)                   ! work array (reciprocal pivot growth factor in work(1) on exit)
  ! testing variables
  logical, parameter :: test=.false.

  ! allocate memory for choice of matrix storage
  if (f_obj % banded) then ! banded storage
   KL=f_obj % subdiag; KU=f_obj % superdiag 
   LDA=KL + KU + 1; LDAF=LDA+KL
   allocate(AF(1:LDAF,1:M),WORK(1:3_i4b*M)) ! storing LU factors requires an additional f_obj % subdiag rows
  else ! full matrix storage
   LDA=M; LDAF=M
   allocate(AF(1:M,1:M),WORK(1:4_i4b*M))
  end if

  ! LAPACK Parameters
  FACT='E'  ! equilibrate matrix prior to factoring
  TRANS='N' ! assume no transposition
  EQUED='N' ! assume no initial equilibration
  N=M; LDX=M; LDB=M ! assume a square system and use arrays of minimum size
  NRHS=1 ! assume a single right-hand side vector

  ! begin LAPACK operations
  ! Use expert solver with scaling and iterative refinement
  if (f_obj % banded) then ! banded matrix storage
   call DGBSVX(FACT,TRANS,N,KL,KU,NRHS,A,LDA,AF,LDAF,IPIV,EQUED,RA,CA,B,LDB,X,LDX,RCOND,FERR,BERR,WORK,IWORK,INFO)
  else ! full matrix storage
   call DGESVX(FACT,TRANS,N,NRHS,A,LDA,AF,LDAF,IPIV,EQUED,RA,CA,B,LDB,X,LDX,RCOND,FERR,BERR,WORK,IWORK,INFO)
  end if

  if (INFO.ne.0) then
   if (INFO.eq.(N+1_i4b)) then
    if (f_obj % out_warning) then
     write(f_obj % unit,*) "LAPACK Warning: RCOND=",RCOND,"may be too low for an accurate solution."
    end if
   else
    if (f_obj % out_error) then
     write(f_obj % unit,*) "LAPACK Error: DGESVX exited with an error code of",info,"."; stop
    end if
   end if
  end if
  if (test) then
   write(f_obj % unit,*) "LAPACK Test:",RCOND,FERR,BERR ! print error information for testing
  else if ((FERR(1).gt.tol).or.(BERR(1).gt.tol)) then
   if (f_obj % out_warning) then
    write(f_obj % unit,*) "LAPACK Warning -- tolerance not met:",RCOND,FERR,BERR ! print error information if tolerance is not met
   end if
  end if
  B=X ! put solution in output vector

 end subroutine linear_solve

end module Newton_solvers
