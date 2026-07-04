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
    if (f_obj % dynamic) then ! for dynamic selection of classical or nested regimes
     ! note: start with classical regime and switch to nested regime if needed 
     f_obj % dynamic_classical = .true.
     call Newton_vector(f_obj) ! classical iterations using classical algorithm
     if (.not.f_obj % dynamic_classical) then ! go to nested iterations if classical iterations do not converge well

      ! need to intialize f1,f2,J1,J2 for intial condition from classical iterations
      ! note: guess vector elements from classical iterations are reused where possible
      call f_obj % f1_f2_vec_eval(f_obj % x0) ! get f1 and f2 (also initializes scaled residual and computes line search objective function)
      if (f_obj % f_error) return             ! check for function evaluation errors
      f_obj % L0 = f_obj % out_SS4HG % fNew   ! initialize line search objective function value based on computed value 
      call f_obj % J1_eval(f_obj % x0)        ! get J1
      call f_obj % J2_eval(f_obj % x0)        ! get J2

      call nested_Newton_vector(f_obj,f_obj % kmax,f_obj % lmax) ! nested iterations

     end if   
    else ! use nested regime only (original behaviour)
     call nested_Newton_vector(f_obj,f_obj % kmax,f_obj % lmax)   
    end if
   else ! classical iterations only
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
  ! *** Newton solver for vector problems ***
  use Newton_functions,only: LS_C
  type(f_obj_type),intent(inout) :: f_obj        ! nested Newton object 
  real(r8b)    :: final_mean                     ! mean of final solution vector
  real(r8b)    :: R_est                          ! estimated max relative difference in solution between iterations
  integer(i4b) :: k                              ! iteration counter
  integer(i4b) :: i                              ! loop index
  logical      :: exit_flag                      ! exit flag
  logical      :: return_flag                    ! return flag for early return from Newton solver call
  ! LAPACK Variables
  real(r8b)    :: B(1:f_obj % n,1:1)             ! right-hand side / solution vector

  ! initialize error flag
  f_obj % f_error = .false. ! error flag for the computation of f, f1, or f2

  ! initialize convergence flag
  f_obj % converged = .false.
  
  ! initialize residual values
  f_obj % R(:)      = 0._r8b ! use default residual value used as solver output for early exit

  ! determine Newton step refinement option and whether we need to evaluate the RHS B vector (can be reused from the line search)
  if (f_obj % refinement) then
   !f_obj % line_search_option = LS_C ! assignment not needed
   f_obj % evaluate_B = .false.
  else
   f_obj % evaluate_B = .true.
  end if

  f_obj % inner = .false. ! classical iterations only
  exit_flag=.false.
  f_obj % xk(:) = f_obj % x0(:) ! initialize
  do k=0,f_obj % kmax_classical
   f_obj % k = k ! store index 
   if (f_obj % f_eval_flag) then
    call f_obj % f_vec_eval(f_obj % xk); if (f_obj % f_error) return ! compute non-linear function vector (f_obj % f_vec)
   end if

   ! obtain RHS vector
   if (f_obj % evaluate_B) then ! compute (unscaled) RHS if using the 'L' scheme or not doing the line search
    do concurrent (i = 1:f_obj % n)
     B(i,1) = -f_obj % f_vec(i) ! initialize right-side vector used by LAPACK
    end do
   else ! get scaled RHS from previous line search call or initial value
    do concurrent (i = 1:f_obj % n)
     B(i,1) = -f_obj % rVecScaled(i)
    end do
   end if

   ! solve for Newton step
   if (f_obj % J_eval_flag) call f_obj % J_eval(f_obj % xk) ! compute Jacobian (f_obj % J)
   f_obj % AF(:,:) = f_obj % J(:,:) ! load matrix used for LU factors
   call linear_solve(f_obj,B,f_obj % tol) ! Solve Jx=B -- x stored in B on output
   if (f_obj % LAPACK_error) return ! check for LAPACK errors to allow recovery (if supported by the external driver)

   ! Newton step refinement and update guess
   if (f_obj % refinement) then
    call f_obj % apply_nested_line_search(LS_C,.false.,B(:,1)); if (f_obj % f_error) return
   else
    do concurrent (i = 1:f_obj % n)
     f_obj % xkp1(i) = f_obj % xk(i) + B(i,1) ! update guess based on unrefined Newton step
    end do
   end if

   call check_residual_vector(f_obj,f_obj % convergence,k,f_obj % xkp1,f_obj % xk,&
                             &R_est,exit_flag,return_flag); if (return_flag) return
   if (f_obj % out_detail) write(f_obj % unit,'(i4,3(g23.15))') f_obj % k,sum(f_obj % xk)/f_obj % n,f_obj % R(0),R_est
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
    write(f_obj % unit,'(i4,3(g23.15))') f_obj % k+1,final_mean,f_obj % R(1)
   end if
   f_obj % kcount=f_obj % k+1
  else
   f_obj % kcount=f_obj % k
  end if

  if (f_obj % k.gt.f_obj % kmax_classical) then
   if (f_obj % out_warning) then
    write(f_obj % unit,*) "Warning - classical Newton solver has reached the maximum number of iterations&
                          & - accuracy may not be sufficient."
   end if
  end if

  !f_obj % x1(:) = f_obj % xkp1(:)
  f_obj % x0(:) = f_obj % xkp1(:)
  if (f_obj % out_basic) write(f_obj % unit,*) "Convergence Error=",f_obj % R(1)

 end subroutine Newton_vector

 subroutine nested_Newton_vector(f_obj,kmax,lmax)
  ! *** Nested Newton solver for vector problems ***
  use Newton_functions,only: LS_C,LS_I,LS_O ! line search options
  use Newton_functions,only: custom!,custom_strict,strict,custom_predictive,predictive ! convergence options 
  ! arguments
  type(f_obj_type),intent(inout) :: f_obj 
  integer(i4b),intent(in) :: kmax                           ! max k index value (outer iterations)
  integer(i4b),intent(in) :: lmax                           ! max l index value (inner iterations)
  ! Newton solver variables
  real(r8b)    :: final_mean                     ! mean value of final solution vector
  real(r8b)    :: R_est                          ! estimated max relative difference in solution between iterations
  integer(i4b) :: l_total                        ! total number of inner iterations
  integer(i4b) :: i,j
  logical      :: exit_outer,exit_inner          ! exit flags for outer and inner loops
  logical      :: return_flag                    ! return flag for early return from Newton solver call
  ! LAPACK Variables
  real(r8b)    :: B(1:f_obj % n,1:1)             ! right-hand side / solution vector

  ! initialize error flag
  f_obj % f_error = .false. ! error flag for the computation of f, f1, or f2

  ! initialize convergence flag
  f_obj % converged = .false.

  ! initialize residual values
  f_obj % R(:)       = 0._r8b ! use default residual value used as solver output for early exit
  f_obj % R_inner(:) = 0._r8b ! use default residual value used as solver output for early exit
 
  l_total=0
  exit_outer=.false.
  f_obj % inner=.false.        ! start with outer iterations
  f_obj % xk0(:)=f_obj % x0(:) ! initial guess
  f_obj % k = 0_i4b ! intialize loop index
  outer: do while (f_obj % k <= kmax)

   if (f_obj % f2_eval_flag) then
    call f_obj % f2_vec_eval(f_obj % xk0); if (f_obj % f_error) return
   end if
   if (f_obj % J2_eval_flag) call f_obj % J2_eval(f_obj % xk0) ! compute Jacobian
   exit_inner=.false.
   f_obj % xkp1l(:) = f_obj % xk0(:) !initial guess for inner iterations

   f_obj % inner=.true. ! inner iterations for next loop
   f_obj % l = 0_i4b ! intialize loop index
   !f_obj % lmax_loop = f_obj % lmax ! initialize lmax value used for inner loop
   f_obj % lmax_loop = lmax ! initialize lmax value used for inner loop -- use value from input argument
   inner: do while (f_obj % l <= f_obj % lmax_loop) ! inner iterations

    if (f_obj % f1_eval_flag) then
     call f_obj % f1_vec_eval(f_obj % xkp1l); if (f_obj % f_error) return
    end if

    ! refactored for efficiency
    ! determine line search scheme and whether we need to evaluate RHS vector
    if (f_obj % refinement) then
     if (f_obj % lmax == 0_i4b) then ! classical regime -- why doesn't (f_obj % lmax_loop == 0_i4b) work here?
      f_obj % line_search_option = LS_C
      f_obj % evaluate_B = .false.
     else if (f_obj % l < f_obj % lmax_loop) then ! initial inner iterations
      f_obj % line_search_option = LS_I
      if (f_obj % k == 0_i4b) then ! first inner scheme iteration -- reuse initial value from systemSolv 
       if (f_obj % l == 0_i4b) f_obj % evaluate_B = .false.
      else if (f_obj % l > 0_i4b) then
       f_obj % evaluate_B = .false.
      else
       f_obj % evaluate_B = .true.
      end if
     else ! last inner iteration
      f_obj % line_search_option = LS_O
      f_obj % evaluate_B = .false. ! RHS computed in previous inner LS scheme
     end if
    else
     f_obj % evaluate_B = .true.
    end if

    ! obtain RHS vector
    if (f_obj % evaluate_B) then ! compute (unscaled) RHS if using the 'L' scheme or not doing the line search
      B(:,1) = f_obj % matrix_vector_product(f_obj % J2,f_obj % xkp1l - f_obj % xk0)
     do concurrent (i = 1:f_obj % n)
      B(i,1) = -(f_obj % f1_vec(i) - f_obj % f2_vec(i)) + B(i,1)
     end do
    else ! get scaled RHS from previous line search call or initial value
     do concurrent (i = 1:f_obj % n)
      B(i,1) = -f_obj % rVecScaled(i)
     end do
    end if

    if (f_obj % J1_eval_flag) call f_obj % J1_eval(f_obj % xkp1l) ! compute Jacobian
    do concurrent (i = 1:f_obj % nrow, j = 1:f_obj % n)
     f_obj % AF(i,j) = f_obj % J1(i,j) - f_obj % J2(i,j) ! difference of Jacobians (formerly Jdiff)
    end do
    call linear_solve(f_obj,B,f_obj % tol) ! Solve AF*x_step_inner=B -- inner Newton step stored in B on output
    if (f_obj % LAPACK_error) return ! check for LAPACK errors to allow recovery (if supported by the external driver)

    ! apply Newton step refinement and update guess
    if (f_obj % refinement) then
     if (f_obj % line_search_option == LS_O) then
      do concurrent (i = 1:f_obj % n) ! compute search direction for outer line search scheme
       B(i,1) = f_obj % xkp1l(i) + B(i,1) - f_obj % xk0(i)
      end do
      !call f_obj % apply_nested_line_search(f_obj % line_search_option,.true.,f_obj % xkp1l(:) + B(:,1) - f_obj % xk0(:)); if (f_obj % f_error) return
     !else
     ! call f_obj % apply_nested_line_search(f_obj % line_search_option,.true.,B(:,1)); if (f_obj % f_error) return
     end if
     call f_obj % apply_nested_line_search(f_obj % line_search_option,.true.,B(:,1)); if (f_obj % f_error) return
    else
     do concurrent (i = 1:f_obj % n)
      f_obj % xkp1lp1(i)=f_obj % xkp1l(i)+B(i,1) ! update guess if no refinement
     end do
    end if

    call check_residual_vector(f_obj,f_obj % convergence_inner,f_obj % l,f_obj % xkp1lp1,f_obj % xkp1l,&
                              &R_est,exit_inner,return_flag); if (return_flag) return
    ! print exact convergence error for iteration l
    if (f_obj % out_detail) write(f_obj % unit,'(a2,i4,3(g23.15))') "  ",f_obj % l,sum(f_obj % xkp1l)/f_obj % n,f_obj % R_inner(0),R_est
    if (exit_inner) exit inner
    if (f_obj % constraints_inner) then
     call f_obj % apply_constraints(f_obj % xkp1l,f_obj % xkp1lp1) ! apply constraints without interfering with the convergence criterion
    end if
    f_obj % xkp1l(:) = f_obj % xkp1lp1(:) ! set up next inner iteration

    f_obj % l = f_obj % l + 1_i4b
   end do inner

   if (f_obj % l.gt.f_obj % lmax_loop) then
    if (f_obj % out_detail) write(f_obj % unit,'(a2,i4,2(g23.15))') "  ",f_obj % l,sum(f_obj % xkp1lp1)/f_obj % n,f_obj % R_inner(1)
    if (f_obj % out_warning) then
     write(f_obj % unit,*) "Warning - nested Newton solver has reached the maximum number of inner iterations&
                           & - accuracy may not be sufficient."
    end if
   end if
   ! inner iteration counts 
   if (exit_inner) then
    ! final output for inner iterations
    if (f_obj % out_detail) write(f_obj % unit,'(a2,i4,2(g23.15))') "  ",f_obj % l+1,sum(f_obj % xkp1lp1)/f_obj % n,f_obj % R_inner(1) 
    l_total=l_total+(f_obj % l+1_i4b)
   else
    l_total=l_total+f_obj % l
   end if

   f_obj % inner=.false.

   call check_residual_vector(f_obj,f_obj % convergence,f_obj % k,f_obj % xkp1lp1,f_obj % xk0,&
                             &R_est,exit_outer,return_flag); if (return_flag) return
   if (f_obj % out_detail) then ! convergence error info for iteration k
    write(f_obj % unit,'(i4,3(g23.15))') f_obj % k,sum(f_obj % xk0)/f_obj % n,f_obj % R(0),R_est 
   end if
   if (exit_outer) then ! exit loop if convergence criterion is met
    f_obj % converged = .true.
    exit outer
   end if

   if (f_obj % constraints) call f_obj % apply_constraints(f_obj % xk0,f_obj % xkp1lp1) ! apply constraints without interfering with the convergence criterion
   f_obj % xk0(:) = f_obj % xkp1lp1(:)

   f_obj % k = f_obj % k + 1_i4b
  end do outer

  ! outer iterations counts 
  if (exit_outer) then
   ! final convergence error for outer iterations (if early loop exit occurred)
   if (f_obj % out_detail) then
    final_mean=sum(f_obj % xkp1lp1)/f_obj % n
    if (f_obj % convergence .ne. custom) then
     write(f_obj % unit,'(i4,2(g23.15))') f_obj % k+1,final_mean,f_obj % R(1) ! mean of final solution 
    else ! custom methods don't have f_obj % R computed
     write(f_obj % unit,'(i4,2(g23.15))') f_obj % k+1,final_mean ! mean of final solution 
    end if
   end if
   f_obj % kcount = f_obj % k+1_i4b
  else
   f_obj % kcount = f_obj % k
  end if

  if (f_obj % k.gt.kmax) then
   if (f_obj % out_warning) then
    write(f_obj % unit,*) "Warning - nested Newton solver has reached the maximum number of outer iterations&
                          & - accuracy may not be sufficient."
   end if
  end if

  !f_obj % x1(:) = f_obj % xkp1lp1(:)
  f_obj % x0(:) = f_obj % xkp1lp1(:)
  if ((f_obj % out_basic).and.(f_obj % convergence .ne. custom)) write(f_obj % unit,*) "Convergence Error=",f_obj % R(1)
  f_obj % lcount = l_total
  if (f_obj % out_detail) then
   write(f_obj % unit,*) "# of outer iterations=",f_obj % kcount
   write(f_obj % unit,*) "# of inner iterations=",l_total
  end if
 end subroutine nested_Newton_vector

 subroutine check_residual_vector(f_obj,convergence,iteration,xkp1,xk,R_est,exit_flag,return_flag)
  ! *** Check residual vector for potential loop exit ***
  use,intrinsic :: ieee_arithmetic,only: ieee_is_finite
  use Newton_functions,only: custom,custom_strict,strict,custom_predictive,predictive ! convergence options 
  type(f_obj_type),intent(inout) :: f_obj 
  integer(i4b),intent(in)  :: convergence         ! convergence option string that adapts to inner and outer/classical iterations
  integer(i4b),intent(in)  :: iteration           ! interation count
  real(r8b),intent(in)     :: xkp1(:)   ! current root estimate
  real(r8b),intent(in)     :: xk(:)     ! previous root estimate
  logical,intent(inout)    :: exit_flag           ! exit flag
  logical,intent(out)      :: return_flag         ! return flag for early return from Newton solver call
  real(r8b),intent(out)    :: R_est               ! estimated R for current iteration (computed in the previous call)
  ! local variables
  real(r8b)                :: tol                 ! tolerance
  real(r8b)                :: R(-1:1)             ! maximum residual array (two previous exact values and prediction for next iteration)
  integer(i4b)             :: i                   ! index for residual vector
  real(r8b)                :: R_vec(1:f_obj % n)  ! residual vector
  real(r8b)                :: b                   ! exponent used for convergence error estimation 
  real(r8b),parameter      :: tol_inner_LS = 1.e-4_r8b ! 10._r8b*epsilon(1._r8b) ! tolerance threshold for switching to outer line search scheme during inner iterations 

  return_flag = .false. ! initialize return flag

  if (convergence.eq.custom) then ! use custom convergence criterion

   exit_flag = f_obj % custom_convergence()
   if (exit_flag) return  ! set exit flag if criterion is satisfied

   ! if doing line search for inner iterations and inner iterate has not changed much, ensure that one more inner iteration is performed using the outer line search scheme
   ! note: this eliminates unproductive inner iterations
   if (f_obj % refinement) then
    if (f_obj % inner) then
     if (f_obj % l < f_obj % lmax_loop) then
      ! look for precise agreement within a tight tolerance
      call compute_relative_residual
      if (all(R_vec < tol_inner_LS)) then
       f_obj % lmax_loop = f_obj % l + 1_i4b; return
      end if
     end if
    end if
   end if

   ! SJT: new dynamic lmax option for switching between classical and nested regimes
   if (f_obj % nested) then
    if (f_obj % dynamic) then
     if (f_obj % dynamic_classical) then
      if (.not.f_obj % inner) then ! check classical residuals using outer iteration residuals
       call check_dynamic_mode ! classical algorithm used
       if (return_flag) return ! return if switching from classical to nested iterations
      end if
     end if
    end if
   end if

  else

   ! for hybrid of custom and built-in methods: check custom flag for possible early exit (else proceed with built-in methods)
   if ((convergence.eq.custom_strict).or.(convergence.eq.custom_predictive)) then
    exit_flag = f_obj % custom_convergence()
    if (exit_flag)  return  ! set exit flag if criterion is satisfied
   end if

   ! compute current residual
   call compute_relative_residual

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
   if ((convergence.eq.strict).or.(convergence.eq.custom_strict)) then ! strict estimate
    R(1)=R(0) ! estimated residual for iteration+1
   else if ((convergence.eq.predictive).or.(convergence.eq.custom_predictive)) then
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
    ! if doing line search for inner iterations, ensure that one more inner iteration is performed using the outer line search scheme
    if ((f_obj % inner).and.(f_obj % refinement).and.(f_obj % l < f_obj % lmax_loop)) then
     f_obj % lmax_loop = f_obj % l + 1_i4b; return
    else
     exit_flag=.true.; return  ! set exit flag if criterion is satisfied
    end if
   end if

  end if

  contains

   subroutine compute_relative_residual
    ! ** compute current residual **
    do concurrent (i=1:f_obj % n)
     if (xk(i).ne.0._r8b) then
      R_vec(i)=abs((xkp1(i)-xk(i))/xk(i))
     else if (xkp1(i).ne.0._r8b) then
      R_vec(i)=abs(xkp1(i)-xk(i)) ! avoid residuals of unity (since xk(i) equals zero)
     else
      R_vec(i)=0._r8b ! both xk and xkp1 are zero -- set the residual to zero
     end if
    end do
   end subroutine compute_relative_residual

   subroutine check_dynamic_mode
    ! ** Dynamic Newton iteration type selection mode: check convergence order of classical iterations and swith to nested if needed **
    !logical                :: accept(1:f_obj % n) ! accept classical guess as initial guess for nested iterations in dynamic mode?
    integer(i4b),parameter :: k_check=20_i4b ! k_check=2_i4b is the minimum
    
    if (f_obj % k == k_check - 2_i4b) then
     f_obj % xk_0(:) = xk(:)   ! x0
     f_obj % xk_1(:) = xkp1(:) ! x1
    else if (f_obj % k == k_check) then ! check convergence order for third classical iteration

     call check_convergence_order(f_obj % order_min,f_obj % xk_0,f_obj % xk_1,xk,xkp1,&
                                 &f_obj % accept,f_obj % dynamic_revert,f_obj % dynamic_classical)

     ! go to nested iterations if needed 
     if (.not.f_obj % dynamic_classical) then
      if (.not.f_obj % dynamic_revert) then
       f_obj % x0(:) = merge(xkp1,f_obj % x0,f_obj % accept) ! use accepted classical guess vector components for nested initial guess
      end if
      return_flag = .true.
      return 
     end if

    end if
   end subroutine check_dynamic_mode

   subroutine convergence_order_cutoff
    ! ** Dynamic Newton iteration type selection mode: check convergence order of classical iterations and swith to nested if needed **
    logical              :: accept(1:f_obj % n) ! accept classical guess as initial guess for nested iterations in dynamic mode?
    logical              :: revert    ! does solution vector need to be completely reverted to original guess before starting nested iterations?
    logical              :: success   ! is the convergence order threshold successfully met for all solution vector elements?
    real(r8b),parameter       :: order_min=0._r8b    
    integer(i4b),parameter    :: k_cutoff=20_i4b

    if (.not.f_obj % inner) then ! check classical residuals using outer iteration residuals
     if (f_obj % k == k_cutoff-2_i4b) then
      f_obj % xk_0(:) = xk(:)   ! x0
      f_obj % xk_1(:) = xkp1(:) ! x1
     else if (f_obj % k == k_cutoff) then ! check convergence order for third classical iteration

      call check_convergence_order(order_min,f_obj % xk_0,f_obj % xk_1,xk,xkp1,accept,revert,success)

      if (.not.success) return_flag = .true.

     end if
    end if
   end subroutine convergence_order_cutoff

   subroutine check_convergence_order(order_min,x0,x1,x2,x3,accept,revert,success)
    real(r8b),intent(in) :: order_min
    real(r8b),intent(in) :: x0(:),x1(:),x2(:),x3(:)
    real(r8b)            :: x3m2,x2m1,x1m0   ! absolute differences used in convergence order calculation for dynamic mode
    real(r8b)            :: num_arg,den_arg     ! arguments for numerator and denominator of convergence order formula
    real(r8b)            :: order               ! approximate convergence order
    logical,intent(out)  :: accept(:) ! accept classical guess element as initial guess element for nested iterations in dynamic mode?
    logical,intent(out)  :: revert    ! does solution vector need to be completely reverted to original guess before starting nested iterations?
    logical,intent(out)  :: success   ! is the convergence order threshold successfully met for all solution vector elements?

      success = .true.
      revert = .false. ! initialize initial condition reversion flag
      do i=1,f_obj % n

       ! initialize acceptance flag
       if (ieee_is_finite(x3(i))) then ! check for normal finite ieee value (e.g., not infinity or NaN)
        accept(i) = .true.
       else
        accept(i) = .false.
        success = .false.
        revert = .true. ! cannot use classical iteration solution as initial condition
        exit ! no need to check remaining solution values because we are reverting to the original solution
       end if

       ! cases where classical iterations have converged to machine precision
       if (x1(i) == x0(i)) cycle
       if (x2(i) == x1(i)) cycle
       if (x3(i) == x2(i)) cycle

       ! compute absolute differences
       x1m0 = x1(i) - x0(i)
       x2m1 = x2(i) - x1(i)
       x3m2 = x3(i) - x2(i)

       ! compute arguments for comvergence order formula
       num_arg = abs(x3m2/x2m1) ! numerator argument (error ratio of k=3 and k=2)
       den_arg = abs(x2m1/x1m0) ! denominator argument (error ratio of k=2 and k=1)
       if ((num_arg == 1._r8b).or.(den_arg == 1._r8b)) then ! if errors did not improve over classical iterations, used nested
        accept(i) = .false. ! do not accept solution for nested initial guess
        success = .false.
       end if

       ! compute estimate of convergence order
       order = log( num_arg ) / log( den_arg )
       if (order < order_min) then ! if convergence rate is not sufficiently high with classical go to nested
        accept(i) = .false. ! do not accept solution for nested initial guess
        success = .false.
       end if

      end do

   end subroutine check_convergence_order

 end subroutine check_residual_vector

 subroutine linear_solve(f_obj,B,tol)
  use Newton_functions,only: LAPACK_standard,LAPACK_expert
  ! *** Solve Ax=B -- x stored in B on output *** 
  type(f_obj_type),intent(inout) :: f_obj          ! nested Newton object
  ! LAPACK Variables
  real(r8b),intent(inout) :: B(:,:)                ! right-hand side / solution vector
  real(r8b),intent(in)    :: tol                   ! tolerance value used by the calling routine
  ! local variables
  character(1),parameter :: FACT='E'               ! option for matrix factoring (equilibrate matrix prior to factoring)
  character(1),parameter :: TRANS='N'              ! option for matrix transposition (no transposition)
  character(1)           :: EQUED                  ! specifies equilibration type (no initial equilibration)
  integer(i4b),parameter :: NRHS = 1_i4b           ! # of right-hand-side vectors
  integer(i4b) :: INFO                             ! error code
  !integer(i4b) :: IPIV(1:f_obj % n)                ! pivot index vector
  !integer(i4b) :: IWORK(1:f_obj % n)               ! work integer array
  !real(r8b) :: RA(1:f_obj % n),CA(1:f_obj % n)     ! row and column scale factors for A
  real(r8b) :: RCOND                               ! estimate of condition number reciprocal
  !real(r8b) :: X(1:f_obj % n,1:1)                  ! solution to original (unscaled) system
  real(r8b) :: FERR(1:1),BERR(1:1)                 ! forward and backward error estimates (single right-hand side assumed)

  ! initialize error flag (used to enable recoverable errors for external drivers)
  f_obj % LAPACK_error = .false.

  ! begin LAPACK operations
  if (f_obj % linear_system_solver .eq. LAPACK_standard) then ! use standard LAPACK solver
   if (f_obj % banded) then ! banded matrix storage
    ! load banded storage matrix used by LAPACK (stores LU factors on output)
    !f_obj % AF(:,:)=A(:,:) ! load matrix used by LAPACK (stores LU factors on output) 
    ! scale (if needed)
    if (f_obj % scaling) call f_obj % custom_scaling(B) ! B will be scaled solution vector after solving
    ! solve 
    call DGBSV(f_obj % n,f_obj % KL,f_obj % KU,NRHS,f_obj % AF,f_obj % LDAF,f_obj % IPIV,B,f_obj % LDB,INFO)
   else ! full matrix storage
    !f_obj % AF(:,:)=A(:,:) ! load matrix used by LAPACK (stores LU factors on output) 
    ! scale
    if (f_obj % scaling) call f_obj % custom_scaling(B) ! B will be scaled solution vector after solving
    call DGESV(f_obj % n,NRHS,f_obj % AF,f_obj % LDAF,f_obj % IPIV,B,f_obj % LDB,INFO) ! solve
   end if
  else if (f_obj % linear_system_solver .eq. LAPACK_expert) then ! Use expert LAPACK solver with scaling and iterative refinement
   if (f_obj % scaling) then
     if (f_obj % out_error) then
      write(f_obj % unit,*) "LAPACK Error: expert solver does not currently support custom scaling."
     end if
     stop ! fatal error
   end if
   EQUED='N' ! note: not a parameter because LAPACK may change this value on output
   if (f_obj % banded) then ! banded matrix storage ---------------- may need to update A argument in this call (tried a fix but not tested)
    call DGBSVX(FACT,TRANS,f_obj % n,f_obj % KL,f_obj % KU,NRHS,f_obj % AF(f_obj % KL+1:,:),f_obj % LDA,f_obj % AF,f_obj % LDAF,&
               &f_obj % IPIV,EQUED,f_obj % RA,f_obj % CA,B,f_obj % LDB,f_obj % X,f_obj % LDX,RCOND,FERR,BERR,f_obj % WORK,f_obj % IWORK,INFO)
   else ! full matrix storage
    call DGESVX(FACT,TRANS,f_obj % n,NRHS,f_obj % AF,f_obj % LDA,f_obj % AF,f_obj % LDAF,f_obj % IPIV,EQUED,f_obj % RA,f_obj % CA,B,f_obj % LDB,&
               &f_obj % X,f_obj % LDX,RCOND,FERR,BERR,f_obj % WORK,f_obj % IWORK,INFO)
   end if
   B(:,:)=f_obj % X(:,:) ! put solution in output vector
  end if

  ! compute descaled solution if needed
  if (f_obj % scaling) call f_obj % custom_descaling(B)

  ! error control
  if (f_obj % linear_system_solver .eq. LAPACK_standard) then ! use standard LAPACK solver
   if (INFO.ne.0) then
     if (f_obj % out_warning) then
      write(f_obj % unit,*) "LAPACK Error: DGESV or DGBSV exited with an error code of",info,"."
     end if
     f_obj % LAPACK_error = .true.
     return ! recoverable error (if supported by external driver)
   end if
  else if (f_obj % linear_system_solver .eq. LAPACK_expert) then ! Use expert LAPACK solver with scaling and iterative refinement
   if (INFO.ne.0) then
    if (INFO.eq.(f_obj % n+1_i4b)) then
     if (f_obj % out_warning) then
      write(f_obj % unit,*) "LAPACK Warning: RCOND=",RCOND,"may be too low for an accurate solution."
     end if
    else
     if (f_obj % out_warning) then
      write(f_obj % unit,*) "LAPACK Error: DGESVX or DGBSVX exited with an error code of",info,"."
     end if
     f_obj % LAPACK_error = .true.
     return ! recoverable error (if supported by external driver)
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
