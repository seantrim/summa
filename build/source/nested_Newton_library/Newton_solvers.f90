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
  type(f_obj_type),intent(inout) :: f_obj        ! nested Newton object 
  real(r8b)    :: final_mean                     ! mean of final solution vector
  real(r8b)    :: R_est                          ! estimated max relative difference in solution between iterations
  integer(i4b) :: k                              ! iteration counter
  logical      :: exit_flag                      ! exit flag
  ! LAPACK Variables
  real(r8b)    :: B(1:f_obj % n,1:1)             ! right-hand side / solution vector

  ! initialize convergence flag
  f_obj % converged = .false.

  f_obj % inner = .false. ! classical iterations only
  exit_flag=.false.
  f_obj % xk(:) = f_obj % x0(:) ! initialize
  do k=0,f_obj % kmax
   if (f_obj % f_eval_flag) call f_obj % f_vec_eval(f_obj % xk) ! compute non-linear function vector (f_obj % f_vec)
   if (f_obj % J_eval_flag) call f_obj % J_eval(f_obj % xk)     ! compute Jacobian (f_obj % J)

   ! begin LAPACK operations
   B(:,1)=-f_obj % f_vec(:) ! initialize right-side vector used by LAPACK
   call linear_solve(f_obj,f_obj % J,B,f_obj % tol) ! Solve Jx=B -- x stored in B on output

   if (f_obj % refinement) then
    call f_obj % apply_refinement(.false.,f_obj % J,f_obj % xk,B(:,1),f_obj % xkp1,f_obj % f_vec) ! apply Newton step refinement to obtain next guess
   else
    f_obj % xkp1(:) = f_obj % xk(:) + B(:,1) ! update guess based on unrefined Newton step
   end if

   call check_residual_vector(f_obj,k,f_obj % xkp1,f_obj % xk,R_est,exit_flag)
   if (f_obj % out_detail) write(f_obj % unit,'(i4,3(g23.15))') k,sum(f_obj % xk)/f_obj % n,f_obj % R(0),R_est
   if (exit_flag) then ! exit loop if convergence criterion is met
    f_obj % converged = .true.
    exit
   end if

   if (f_obj % constraints) call f_obj % apply_constraints(f_obj % xk,f_obj % xkp1) ! apply constraints without interfering with the convergence criterion
   f_obj % xk(:) = f_obj % xkp1(:) ! prep for next iteration - can probably evaluate in place
  end do
  ! final output
  if (exit_flag.eqv..true.) then
   if (f_obj % out_detail) then
    final_mean=sum(f_obj % xkp1)/f_obj % n 
    write(f_obj % unit,'(i4,3(g23.15))') k+1,final_mean,f_obj % R(1)
   end if
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

  f_obj % x1(:) = f_obj % xkp1(:)
  if (f_obj % out_basic) write(f_obj % unit,*) "Convergence Error=",f_obj % R(1)

 end subroutine Newton_vector

 subroutine nested_Newton_vector(f_obj)
  ! Newton solver
  type(f_obj_type),intent(inout) :: f_obj 
  real(r8b)    :: final_mean                     ! mean value of final solution vector
  real(r8b)    :: R_est                          ! estimated max relative difference in solution between iterations
  integer(i4b) :: k,l                            ! iteration counters
  integer(i4b) :: l_total                        ! total number of inner iterations
  logical      :: exit_outer,exit_inner 
  ! LAPACK Variables
  real(r8b)    :: B(1:f_obj % n,1:1)             ! right-hand side / solution vector
  real(r8b)    :: f2mJ2xk0(1:f_obj % n)             ! right-hand side / solution vector

  ! initialize convergence flag
  f_obj % converged = .false.
 
  l_total=0
  exit_outer=.false.
  f_obj % inner=.false.        ! start with outer iterations
  f_obj % xk0(:)=f_obj % x0(:) ! initial guess
  outer: do k=0,f_obj % kmax

   if (f_obj % f2_eval_flag) call f_obj % f2_vec_eval(f_obj % xk0)
   if (f_obj % J2_eval_flag) call f_obj % J2_eval(f_obj % xk0) ! compute Jacobian
   f2mJ2xk0(:) = f_obj % f2_vec(:) - matrix_vector_product(f_obj,f_obj % J2,f_obj % xk0)
   exit_inner=.false.
   f_obj % xkp1l(:) = f_obj % xk0(:) !initial guess for inner iterations

   f_obj % inner=.true. ! inner iterations for next loop
   inner: do l=0,f_obj % lmax ! inner iterations
    if (f_obj % f1_eval_flag) call f_obj % f1_vec_eval(f_obj % xkp1l)
    if (f_obj % J1_eval_flag) call f_obj % J1_eval(f_obj % xkp1l) ! compute Jacobian
    f_obj % Jdiff(:,:) = f_obj % J1(:,:) - f_obj % J2(:,:)

    ! begin LAPACK operations
    ! initialize right-side vector used by LAPACK
    B(:,1) = f2mJ2xk0(:) - f_obj % f1_vec(:) + matrix_vector_product(f_obj,f_obj % J1,f_obj % xkp1l) 
    call linear_solve(f_obj,f_obj % Jdiff,B,f_obj % tol) ! Solve Jdiff*x=B -- x stored in B on output -- M is the # of rows/columns of A
    f_obj % xkp1lp1(:)=B(:,1) ! update guess

    if (f_obj % refinement_inner) then
     call f_obj % apply_refinement(.true.,f_obj % Jdiff,f_obj % xkp1l,B(:,1),f_obj % xkp1lp1,f_obj % f_vec)
     if (.not.f_obj % f1_eval_flag) call f_obj % f1_vec_eval(f_obj % xkp1lp1) ! ensure we have f1_vec value for refined solution
     !if (.not.f_obj % f2_eval_flag) call f_obj % f2_vec_eval(f_obj % xkp1lp1) ! assume f2_vec does not change during inner iterations
    end if

    call check_residual_vector(f_obj,l,f_obj % xkp1lp1,f_obj % xkp1l,R_est,exit_inner)
    ! print exact convergence error
    if (f_obj % out_detail) write(f_obj % unit,'(a2,i4,3(g23.15))') "  ",l,sum(f_obj % xkp1l)/f_obj % n,f_obj % R_inner(0),R_est
    if (exit_inner) exit inner
    if (f_obj % constraints_inner) then
     call f_obj % apply_constraints(f_obj % xkp1l,f_obj % xkp1lp1) ! apply constraints without interfering with the convergence criterion
    end if
    f_obj % xkp1l(:) = f_obj % xkp1lp1(:) ! set up next inner iteration
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
    if (f_obj % out_detail) write(f_obj % unit,'(a2,i4,2(g23.15))') "  ",l+1,sum(f_obj % xkp1lp1)/f_obj % n,f_obj % R_inner(1) 
    l_total=l_total+(l+1)
   else
    l_total=l_total+l
   end if

   f_obj % inner=.false.

   ! apply Newton step refinement
   if (f_obj % refinement) then
    call f_obj % apply_refinement(.true.,f_obj % Jdiff,f_obj % xk0,B(:,1),f_obj % xkp1lp1,f_obj % f_vec) 
    if (.not.f_obj % f1_eval_flag) call f_obj % f1_vec_eval(f_obj % xkp1lp1) ! ensure we have f1_vec value for refined solution
    if (.not.f_obj % f2_eval_flag) call f_obj % f2_vec_eval(f_obj % xkp1lp1) ! ensure we have f2_vec value for refined solution
   end if

   call check_residual_vector(f_obj,k,f_obj % xkp1lp1,f_obj % xk0,R_est,exit_outer)
   if (f_obj % out_detail) then ! convergence error info for iteration k
    write(f_obj % unit,'(i4,3(g23.15))') k,sum(f_obj % xk0)/f_obj % n,f_obj % R(0),R_est 
   end if
   if (exit_outer) then ! exit loop if convergence criterion is met
    f_obj % converged = .true.
    exit outer
   end if

   if (f_obj % constraints) call f_obj % apply_constraints(f_obj % xk0,f_obj % xkp1lp1) ! apply constraints without interfering with the convergence criterion
   f_obj % xk0(:) = f_obj % xkp1lp1(:)
  end do outer

  ! outer iterations counts 
  if (exit_outer) then
   ! final convergence error for outer iterations (if early loop exit occurred)
   if (f_obj % out_detail) then
    final_mean=sum(f_obj % xkp1lp1)/f_obj % n
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

  f_obj % x1(:) = f_obj % xkp1lp1(:)
  if (f_obj % out_basic) write(f_obj % unit,*) "Convergence Error=",f_obj % R(1)
  f_obj % lcount = l_total
  if (f_obj % out_detail) then
   write(f_obj % unit,*) "# of outer iterations=",f_obj % kcount
   write(f_obj % unit,*) "# of inner iterations=",l_total
  end if
 end subroutine nested_Newton_vector

 subroutine check_residual_vector(f_obj,iteration,xkp1,xk,R_est,exit_flag)
  ! *** Check residual vector for potential loop exit ***
  type(f_obj_type),intent(inout) :: f_obj 
  integer(i4b),intent(in)  :: iteration   ! interation count
  real(r8b),intent(in)     :: xkp1(1:f_obj % n)  ! current root estimate
  real(r8b),intent(in)     :: xk(1:f_obj % n)    ! previous root estimate
  logical,intent(inout)    :: exit_flag   ! exit flag
  real(r8b),intent(out)    :: R_est       ! estimated R for current iteration (computed in the previous call)
  ! local variables
  real(r8b)                :: tol         ! tolerance
  real(r8b)                :: R(-1:1)     ! maximum residual array (two previous exact values and prediction for next iteration)
  integer(i4b)             :: i                  ! index for residual vector
  real(r8b)                :: R_vec(1:f_obj % n) ! residual vector
  real(r8b)                :: b                  ! exponent used for convergence error estimation 
  character(:),allocatable :: convergence        ! convergence option string that adapts to inner and outer/classical iterations

  if (f_obj % inner) then ! inner iterations
   convergence = f_obj % convergence_inner
  else                    ! outer iterations
   convergence = f_obj % convergence
  end if

  if (convergence.eq.'custom') then ! use custom convergence criterion
   exit_flag = f_obj % custom_convergence()
   if (exit_flag)  return  ! set exit flag if criterion is satisfied
  else

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
   if (convergence.eq.'strict') then ! strict estimate
    R(1)=R(0) ! estimated residual for iteration+1
   else if (convergence.eq.'predictive') then
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
    tol = f_obj % tol_inner   ! set tolerance
    f_obj % R_inner(0) = R(0) ! store exact residual for current iteration
    f_obj % R_inner(1) = R(1) ! store estimated residual for iteration+1
   else                    ! outer/classical iterations
    tol = f_obj % tol         ! set tolerance
    f_obj % R(0) = R(0)       ! store exact residual for current iteration
    f_obj % R(1) = R(1)       ! store estimated residual for iteration+1
   end if
   if (iteration.eq.0) R_est=R(1) ! initialize R_est for iteration zero

   ! check exact error from current iteration and estimated error for next iteration
   if ((R(0).lt.tol).or.(R(1).lt.tol)) then
    exit_flag=.true.; return  ! set exit flag if criterion is satisfied
   end if

  end if
 end subroutine check_residual_vector

 function matrix_vector_product(f_obj,A,x) result(y)
  ! *** Compute matrix vector product y=A*x ***
  ! input
  type(f_obj_type),intent(in) :: f_obj        ! class object containing solver options
  real(r8b),intent(in) :: A(:,:)  ! input matrix 
  real(r8b),intent(in) :: x(1:f_obj % n)              ! input vector

  ! output
  real(r8b) :: y(1:f_obj % n)                         ! product vector

  ! local variables
  character(1),parameter :: TRANS='N'                ! option for matrix transposition
  integer(i4b) :: KL,KU                              ! # of subdiagonals and superdiagonals of A (banded storage)
  integer(i4b) :: LDA                                ! first dimension of A
  integer(i4b),parameter :: INCX=1_i4b, INCY=1_i4b   ! increment for elements of x and y vectors
  real(r8b),parameter    :: ALPHA=1._r8b,BETA=0._r8b ! scalars used in LAPACK solvers
 
  ! set LAPACK parameters for choice of matrix storage
  if (f_obj % banded) then ! banded storage
   KL=f_obj % subdiag; KU=f_obj % superdiag 
   LDA=KL + KU + 1
  else ! full matrix storage
   LDA=f_obj % n
  end if
  
  if (f_obj % banded) then ! banded storage
   call DGBMV(TRANS,f_obj % n,f_obj % n,KL,KU,ALPHA,A,LDA,x,INCX,BETA,y,INCY) ! BLAS
  else ! full matrix storage
   !y=matmul(A,x)
   call DGEMV(TRANS,f_obj % n,f_obj % n,ALPHA,A,LDA,x,INCX,BETA,y,INCY) ! BLAS
  end if
 end function matrix_vector_product

 subroutine linear_solve(f_obj,A,B,tol)
  ! *** Solve Ax=B -- x stored in B on output *** 
  type(f_obj_type),intent(inout) :: f_obj                  ! nested Newton object
  ! LAPACK Variables
  real(r8b),intent(in)    :: A(:,:)                        ! input matrix
  real(r8b),intent(inout) :: B(1:f_obj % n,1:1) ! right-hand side / solution vector
  real(r8b),intent(in)    :: tol                   ! tolerance value used by the calling routine
  ! local variables
  character(1),parameter :: FACT='E'               ! option for matrix factoring (equilibrate matrix prior to factoring)
  character(1),parameter :: TRANS='N'              ! option for matrix transposition (no transposition)
  character(1)           :: EQUED                  ! specifies equilibration type (no initial equilibration)
  integer(i4b),parameter :: NRHS = 1_i4b           ! # of right-hand-side vectors
  integer(i4b) :: INFO                             ! error code
  integer(i4b) :: IPIV(1:f_obj % n)                ! pivot index vector
  integer(i4b) :: IWORK(1:f_obj % n)               ! work integer array
  real(r8b) :: RA(1:f_obj % n),CA(1:f_obj % n)     ! row and column scale factors for A
  real(r8b) :: RCOND                               ! estimate of condition number reciprocal
  real(r8b) :: X(1:f_obj % n,1:1)                  ! solution to original (unscaled) system
  real(r8b) :: FERR(1:1),BERR(1:1)                 ! forward and backward error estimates (single right-hand side assumed)

  ! begin LAPACK operations
  if (f_obj % linear_system_solver .eq. "LAPACK_standard") then ! use standard LAPACK solver
   if (f_obj % banded) then ! banded matrix storage
    ! load banded storage matrix used by LAPACK (stores LU factors on output)
    f_obj % AF(1:f_obj % KL,:)=0._r8b; f_obj % AF(f_obj % KL+1:f_obj % LDAF,:)=A(1:f_obj % LDA,:)
    ! scale
    if (f_obj % scaling) call f_obj % custom_scaling(B) ! B will be scaled solution vector after solving
    ! solve 
    call DGBSV(f_obj % n,f_obj % KL,f_obj % KU,NRHS,f_obj % AF,f_obj % LDAF,IPIV,B,f_obj % LDB,INFO)
   else ! full matrix storage
    f_obj % AF(:,:)=A(:,:) ! load matrix used by LAPACK (stores LU factors on output) 
    ! scale
    if (f_obj % scaling) call f_obj % custom_scaling(B) ! B will be scaled solution vector after solving
    call DGESV(f_obj % n,NRHS,f_obj % AF,f_obj % LDAF,IPIV,B,f_obj % LDB,INFO) ! solve
   end if
  else if (f_obj % linear_system_solver .eq. "LAPACK_expert") then ! Use expert LAPACK solver with scaling and iterative refinement
   if (f_obj % scaling) then
     if (f_obj % out_error) then
      write(f_obj % unit,*) "LAPACK Error: expert solver does not currently support custom scaling."
     end if
     stop ! fatal error
   end if
   EQUED='N' ! note: not a parameter because LAPACK may change this value on output
   if (f_obj % banded) then ! banded matrix storage
    call DGBSVX(FACT,TRANS,f_obj % n,f_obj % KL,f_obj % KU,NRHS,A,f_obj % LDA,f_obj % AF,f_obj % LDAF,&
               &IPIV,EQUED,RA,CA,B,f_obj % LDB,X,f_obj % LDX,RCOND,FERR,BERR,f_obj % WORK,IWORK,INFO)
   else ! full matrix storage
    call DGESVX(FACT,TRANS,f_obj % n,NRHS,A,f_obj % LDA,f_obj % AF,f_obj % LDAF,IPIV,EQUED,RA,CA,B,f_obj % LDB,&
               &X,f_obj % LDX,RCOND,FERR,BERR,f_obj % WORK,IWORK,INFO)
   end if
   B(:,:)=X(:,:) ! put solution in output vector
  end if

  ! compute descaled solution if needed (not needed for Newton step refinement)
  if ((f_obj % scaling) .and. (.not.f_obj % refinement)) call f_obj % custom_descaling(B)

  ! error control
  if (f_obj % linear_system_solver .eq. "LAPACK_standard") then ! use standard LAPACK solver
   if (INFO.ne.0) then
     if (f_obj % out_error) then
      write(f_obj % unit,*) "LAPACK Error: DGESV or DGBSV exited with an error code of",info,"."
     end if
     stop ! fatal error
   end if
  else if (f_obj % linear_system_solver .eq. "LAPACK_expert") then ! Use expert LAPACK solver with scaling and iterative refinement
   if (INFO.ne.0) then
    if (INFO.eq.(f_obj % n+1_i4b)) then
     if (f_obj % out_warning) then
      write(f_obj % unit,*) "LAPACK Warning: RCOND=",RCOND,"may be too low for an accurate solution."
     end if
    else
     if (f_obj % out_error) then
      write(f_obj % unit,*) "LAPACK Error: DGESVX or DGBSVX exited with an error code of",info,"."
     end if
     stop ! fatal error
    end if
   end if
   if (f_obj % out_warning) then
    if ((FERR(1).gt.tol).or.(BERR(1).gt.tol)) then ! print error information if tolerance is not met
      write(f_obj % unit,*) "LAPACK Warning -- tolerance not met using expert solver:",RCOND,FERR,BERR 
    end if
   end if
  else
   if (f_obj % out_error) then
    write(f_obj % unit,*) "Linear system solver choice is not supported."
   end if
   stop ! fatal error
  end if

 end subroutine linear_solve

end module Newton_solvers
