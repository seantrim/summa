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
  integer(i4b) :: kmax                           ! max k index value
  integer(i4b) :: lmax                           ! max l index value
  logical(i4b),parameter :: dynamic_strict = .true. ! strict initialization of f1, f2, J1, and J2 for nested portion of dynamic mode

  if (f_obj % n .gt. 0_i4b) then
   if (f_obj % nested) then
    if (f_obj % dynamic) then ! for dynamic selection of classical or nested regimes
     ! note: start with classical regime and switch to nested regime if needed 
     f_obj % f1_vec_save(:) = f_obj % f1_vec(:); f_obj % f2_vec_save(:) = f_obj % f2_vec(:) ! save f1 and f2 for original initial condition in case of reversion
     !f_obj % L_save = f_obj % out_SS4HG % fNew; f_obj % f_vec_scaled_save(:) = f_obj % rVecScaled(:) ! save line search quantities
     f_obj % L_save = f_obj % L0; f_obj % f_vec_scaled_save(:) = f_obj % rVecScaled(:) ! save line search quantities
     f_obj % J1_save(:,:) = f_obj % J1(:,:); f_obj % J2_save(:,:) = f_obj % J2(:,:) ! save J1 and J2 for original initial condition in case of reversion
     kmax = f_obj % kmax_classical
     lmax = 0_i4b
     f_obj % dynamic_classical = .true.
     call nested_Newton_vector(f_obj,kmax,lmax) ! classical iterations
     if (.not.f_obj % dynamic_classical) then

      ! revert to original initial condition if needed
      if (f_obj % dynamic_revert) then
       f_obj % f1_vec(:) = f_obj % f1_vec_save(:); f_obj % f2_vec(:) = f_obj % f2_vec_save(:)
       !f_obj % out_SS4HG % fNew = f_obj % L_save; f_obj % rVecScaled(:) = f_obj % f_vec_scaled_save(:)
       f_obj % L0 = f_obj % L_save; f_obj % rVecScaled(:) = f_obj % f_vec_scaled_save(:)
       f_obj % J1(:,:) = f_obj % J1_save(:,:); f_obj % J2(:,:) = f_obj % J2_save(:,:)
      else if (dynamic_strict) then
       ! need to intialize f1,f2,J1,J2 for intial condition from classical iterations
       call f_obj % f1_f2_vec_eval(f_obj % x0) ! get f1 and f2 (also initializes line search objective function and scaled residual)
       call f_obj % J1_J2_eval(f_obj % x0) ! get J1 and J2
      end if
      kmax = f_obj % kmax
      lmax = f_obj % lmax
      call nested_Newton_vector(f_obj,kmax,lmax) ! nested iterations
     end if   
    else ! use nested regime only (original behaviour)
     kmax = f_obj % kmax
     lmax = f_obj % lmax
     call nested_Newton_vector(f_obj,kmax,lmax)   
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
  ! Newton solver for vector problems
  type(f_obj_type),intent(inout) :: f_obj        ! nested Newton object 
  real(r8b)    :: final_mean                     ! mean of final solution vector
  real(r8b)    :: R_est                          ! estimated max relative difference in solution between iterations
  integer(i4b) :: k                              ! iteration counter
  logical      :: exit_flag                      ! exit flag
  logical      :: return_flag                    ! return flag for early return from Newton solver call
  ! LAPACK Variables
  real(r8b)    :: B(1:f_obj % n,1:1)             ! right-hand side / solution vector

  ! initialize error flag
  f_obj % f_error = .false. ! error flag for the computation of f, f1, or f2

  ! initialize convergence flag
  f_obj % converged = .false.

  ! determine Newton step refinement option and whether we need to evaluate the RHS B vector (can be reused from the line search)
  if (f_obj % refinement) then
   f_obj % line_search_option = 'C'
   f_obj % evaluate_B = .false.
  else
   f_obj % evaluate_B = .true.
  end if

  f_obj % inner = .false. ! classical iterations only
  exit_flag=.false.
  f_obj % xk(:) = f_obj % x0(:) ! initialize
  do k=0,f_obj % kmax_classical
   f_obj % k = k ! store index 
   if (f_obj % f_eval_flag) call f_obj % f_vec_eval(f_obj % xk) ; if (f_obj % f_error) return ! compute non-linear function vector (f_obj % f_vec)
   if (f_obj % J_eval_flag) call f_obj % J_eval(f_obj % xk)     ! compute Jacobian (f_obj % J)

   ! obtain RHS vector
   if (f_obj % evaluate_B) then ! compute (unscaled) RHS if using the 'L' scheme or not doing the line search
    B(:,1) = -f_obj % f_vec(:) ! initialize right-side vector used by LAPACK
   else ! get scaled RHS from previous line search call or initial value
    B(:,1) = -f_obj % rVecScaled(:)
   end if

   ! solve for Newton step
   call linear_solve(f_obj,f_obj % J,B,f_obj % tol) ! Solve Jx=B -- x stored in B on output
   if (f_obj % LAPACK_error) return ! check for LAPACK errors to allow recovery (if supported by the external driver)
   f_obj % xkp1(:) = f_obj % xk(:) + B(:,1) ! update guess based on unrefined Newton step

   ! Newton step refinement
   if (f_obj % refinement) then
    call f_obj % apply_nested_line_search('C')
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

  f_obj % x1(:) = f_obj % xkp1(:)
  if (f_obj % out_basic) write(f_obj % unit,*) "Convergence Error=",f_obj % R(1)

 end subroutine Newton_vector

 subroutine nested_Newton_vector(f_obj,kmax,lmax)
  ! *** Nested Newton solver for vector problems ***
  ! arguments
  type(f_obj_type),intent(inout) :: f_obj 
  integer(i4b),intent(in) :: kmax                           ! max k index value (outer iterations)
  integer(i4b),intent(in) :: lmax                           ! max l index value (inner iterations)
  ! Newton solver variables
  real(r8b)    :: final_mean                     ! mean value of final solution vector
  real(r8b)    :: R_est                          ! estimated max relative difference in solution between iterations
  integer(i4b) :: l_total                        ! total number of inner iterations
  logical      :: exit_outer,exit_inner          ! exit flags for outer and inner loops
  logical      :: return_flag                    ! return flag for early return from Newton solver call
  ! LAPACK Variables
  real(r8b)    :: B(1:f_obj % n,1:1)             ! right-hand side / solution vector
  !real(r8b)    :: f2mJ2xk0(1:f_obj % n)          ! right-hand side / solution vector

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

   if (f_obj % f2_eval_flag) call f_obj % f2_vec_eval(f_obj % xk0); if (f_obj % f_error) return
   if (f_obj % J2_eval_flag) call f_obj % J2_eval(f_obj % xk0) ! compute Jacobian
   !f2mJ2xk0(:) = f_obj % f2_vec(:) - f_obj % matrix_vector_product(f_obj % J2,f_obj % xk0)
   exit_inner=.false.
   f_obj % xkp1l(:) = f_obj % xk0(:) !initial guess for inner iterations

   f_obj % inner=.true. ! inner iterations for next loop
   f_obj % l = 0_i4b ! intialize loop index
   !f_obj % lmax_loop = f_obj % lmax ! initialize lmax value used for inner loop
   f_obj % lmax_loop = lmax ! initialize lmax value used for inner loop -- use value from input argument
   inner: do while (f_obj % l <= f_obj % lmax_loop) ! inner iterations

    ! determine Newton step refinement option
    if (f_obj % refinement) then
     if (f_obj % lmax == 0_i4b) then ! classical regime
      f_obj % line_search_option = 'C'
     else if (f_obj % l < f_obj % lmax_loop) then ! initial inner iterations
      f_obj % line_search_option = 'I'
     else ! last inner iteration
      f_obj % line_search_option = 'L'
      !f_obj % line_search_option = 'C' ! testing 'C' scheme on last inner iteration
     end if
    end if

    if (f_obj % f1_eval_flag) call f_obj % f1_vec_eval(f_obj % xkp1l); if (f_obj % f_error) return
    if (f_obj % J1_eval_flag) call f_obj % J1_eval(f_obj % xkp1l) ! compute Jacobian
    f_obj % Jdiff(:,:) = f_obj % J1(:,:) - f_obj % J2(:,:)

    ! begin LAPACK operations
    ! initialize right-side vector used by LAPACK ---------- OG method (solve for updated inner iteration solution directly)
    !B(:,1) = f2mJ2xk0(:) - f_obj % f1_vec(:) + f_obj % matrix_vector_product(f_obj % J1,f_obj % xkp1l) 
    !call linear_solve(f_obj,f_obj % Jdiff,B,f_obj % tol) ! Solve Jdiff*x=B -- x stored in B on output
    !if (f_obj % LAPACK_error) return ! check for LAPACK errors to all recovery (if supported by the external driver)
    !f_obj % xkp1lp1(:)=B(:,1) ! update guess

    ! SJT: solve for inner step using LAPACK ------- testing ----------------
    ! do we need to evaluate RHS vector?
    if ((f_obj % refinement).and.(f_obj % line_search_option == 'C')) then
     f_obj % evaluate_B = .false.
    else if ((f_obj % refinement).and.(f_obj % line_search_option == 'I')) then
     if ((f_obj % k == 0_i4b).and.(f_obj % l == 0_i4b)) then ! first inner scheme iteration -- reuse initial value from systemSolv 
      f_obj % evaluate_B = .false.
     else if (f_obj % l > 0_i4b) then
      f_obj % evaluate_B = .false.
     else
      f_obj % evaluate_B = .true.
     end if
    else
     f_obj % evaluate_B = .true.
    end if

    ! obtain RHS vector
    if (f_obj % evaluate_B) then ! compute (unscaled) RHS if using the 'L' scheme or not doing the line search
     B(:,1) = -(f_obj % f1_vec(:) - f_obj % f2_vec(:)) + f_obj % matrix_vector_product(f_obj % J2,f_obj % xkp1l - f_obj % xk0)
    else ! get scaled RHS from previous line search call or initial value
     B(:,1) = -f_obj % rVecScaled(:)
    end if

    call linear_solve(f_obj,f_obj % Jdiff,B,f_obj % tol) ! Solve Jdiff*x_step_inner=B -- inner Newton step stored in B on output
    if (f_obj % LAPACK_error) return ! check for LAPACK errors to allow recovery (if supported by the external driver)
    f_obj % xkp1lp1(:)=f_obj % xkp1l(:)+B(:,1) ! update guess

    ! apply Newton step refinement
    if (f_obj % refinement) then
     call f_obj % apply_nested_line_search(f_obj % line_search_option)
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
    if (f_obj % convergence .ne. 'custom') then
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

  f_obj % x1(:) = f_obj % xkp1lp1(:)
  if ((f_obj % out_basic).and.(f_obj % convergence .ne. 'custom')) write(f_obj % unit,*) "Convergence Error=",f_obj % R(1)
  f_obj % lcount = l_total
  if (f_obj % out_detail) then
   write(f_obj % unit,*) "# of outer iterations=",f_obj % kcount
   write(f_obj % unit,*) "# of inner iterations=",l_total
  end if
 end subroutine nested_Newton_vector

 subroutine check_residual_vector(f_obj,convergence,iteration,xkp1,xk,R_est,exit_flag,return_flag)
  ! *** Check residual vector for potential loop exit ***
  use,intrinsic :: ieee_arithmetic,only: ieee_is_finite
  type(f_obj_type),intent(inout) :: f_obj 
  character(*),intent(in)  :: convergence         ! convergence option string that adapts to inner and outer/classical iterations
  integer(i4b),intent(in)  :: iteration           ! interation count
  real(r8b),intent(in)     :: xkp1(1:f_obj % n)   ! current root estimate
  real(r8b),intent(in)     :: xk(1:f_obj % n)     ! previous root estimate
  logical,intent(inout)    :: exit_flag           ! exit flag
  logical,intent(out)      :: return_flag         ! return flag for early return from Newton solver call
  real(r8b),intent(out)    :: R_est               ! estimated R for current iteration (computed in the previous call)
  ! local variables
  real(r8b)                :: tol                 ! tolerance
  real(r8b)                :: R(-1:1)             ! maximum residual array (two previous exact values and prediction for next iteration)
  integer(i4b)             :: i                   ! index for residual vector
  real(r8b)                :: R_vec(1:f_obj % n)  ! residual vector
  real(r8b)                :: b                   ! exponent used for convergence error estimation 
  real(r8b),parameter      :: tol_inner_LS = 10._r8b*epsilon(1._r8b) ! tolerance threshold for switching to outer line search scheme during inner iterations 

  return_flag = .false. ! initialize return flag

  if (convergence.eq.'custom') then ! use custom convergence criterion

   exit_flag = f_obj % custom_convergence()
   if (exit_flag) return  ! set exit flag if criterion is satisfied

   ! if doing line search for inner iterations and inner iterate has not changed, ensure that one more inner iteration is performed using the outer line search scheme
   ! note: this eliminates unproductive inner iterations
   if ((f_obj % inner).and.(f_obj % refinement).and.(f_obj % l < f_obj % lmax_loop)) then
    ! look for precise agreement within a tight tolerance (faster)
    call compute_relative_residual
    if (all(R_vec < tol_inner_LS)) then
     f_obj % lmax_loop = f_obj % l + 1_i4b; return
    end if
   end if

   ! SJT: new dynamic lmax option for switching between classical and nested regimes
   if (f_obj % nested) then
    if (f_obj % dynamic) then
     if (f_obj % dynamic_classical) then
      call check_dynamic_mode(f_obj % xkp1lp1) ! nested algorithm used with lmax=0
      !call check_dynamic_mode(f_obj % xkp1) ! classical algorithm used
     end if
    end if
   end if

  else

   ! for hybrid of custom and built-in methods: check custom flag for possible early exit (else proceed with built-in methods)
   if ((convergence.eq.'custom-strict').or.(convergence.eq.'custom-predictive')) then
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
   if ((convergence.eq.'strict').or.(convergence.eq.'custom-strict')) then ! strict estimate
    R(1)=R(0) ! estimated residual for iteration+1
   else if ((convergence.eq.'predictive').or.(convergence.eq.'custom-predictive')) then
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
    do i=1,f_obj % n
     if (xk(i).ne.0._r8b) then
      R_vec(i)=abs((xkp1(i)-xk(i))/xk(i))
     else if (xkp1(i).ne.0._r8b) then
      R_vec(i)=abs(xkp1(i)-xk(i)) ! avoid residuals of unity (since xk(i) equals zero)
     else
      R_vec(i)=0._r8b ! both xk and xkp1 are zero -- set the residual to zero
     end if
    end do
   end subroutine compute_relative_residual

   subroutine check_dynamic_mode(xkp1)
    ! ** Dynamic Newton iteration type selection mode: check convergence order of classical iterations and swith to nested if needed **
    real(r8b),intent(in) :: xkp1(:)             ! current classical iteration guess
    real(r8b)            :: xk3m2,xk2m1,xk1m0   ! absolute differences used in convergence order calculation for dynamic mode
    real(r8b)            :: num_arg,den_arg     ! arguments for numerator and denominator of convergence order formula
    real(r8b)            :: order               ! approximate convergence order
    logical              :: accept(1:f_obj % n) ! accept classical guess as initial guess for nested iterations in dynamic mode?
    
    if (.not.f_obj % inner) then ! check classical residuals using outer iteration residuals
     if (f_obj % k == 0_i4b) then
      f_obj % xk1(:) = xkp1(:) ! x1
     else if (f_obj % k == 1_i4b) then
      f_obj % xk2(:) = xkp1(:) ! x2
     else if (f_obj % k == 2_i4b) then ! check convergence rate for second classical iteration

      f_obj % dynamic_revert = .false. ! initialize initial condition reversion flag
      do i=1,f_obj % n

       ! initialize acceptance flag
       if (ieee_is_finite(xkp1(i))) then ! check for normal finite ieee value (e.g., not infinity or NaN)
        accept(i) = .true.
       else
        accept(i) = .false.
        f_obj % dynamic_classical = .false.
        f_obj % dynamic_revert = .true. ! cannot use classical iteration solution as initial condition
        exit ! no need to check remaining solution values because we are reverting to the original solution
       end if

       ! cases where classical iterations have converged to machine precision
       if (f_obj % xk1(i) == f_obj % x0(i)) cycle
       if (f_obj % xk2(i) == f_obj % xk1(i)) cycle
       if (xkp1(i) == f_obj % xk2(i)) cycle

       ! compute absolute differences
       xk1m0 = f_obj % xk1(i) - f_obj % x0(i)
       xk2m1 = f_obj % xk2(i) - f_obj % xk1(i)
       xk3m2 = xkp1(i) - f_obj % xk2(i)

       ! compute arguments for comvergence order formula
       num_arg = abs(xk3m2/xk2m1) ! numerator argument (error ratio of k=3 and k=2)
       den_arg = abs(xk2m1/xk1m0) ! denominator argument (error ratio of k=2 and k=1)
       if ((num_arg == 1._r8b).or.(den_arg == 1._r8b)) then ! if errors did not improve over classical iterations, used nested
        accept(i) = .false. ! do not accept solution for nested initial guess
        f_obj % dynamic_classical = .false.
       end if

       ! compute estimate of convergence order
       order = log( num_arg ) / log( den_arg )
       if (order < f_obj % order_min) then ! if convergence rate is not sufficiently high with classical go to nested
        accept(i) = .false. ! do not accept solution for nested initial guess
        f_obj % dynamic_classical = .false.
       end if

      end do

      ! go to nested iterations if needed 
      if (.not.f_obj % dynamic_classical) then
       if (.not.f_obj % dynamic_revert) then
        f_obj % x0(:) = merge(xkp1,f_obj % x0,accept) ! use accepted classical guess vector components for nested initial guess
       end if
       return_flag = .true.
       return 
      end if

     end if
    end if
   end subroutine check_dynamic_mode

 end subroutine check_residual_vector

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

  ! initialize error flag (used to enable recoverable errors for external drivers)
  f_obj % LAPACK_error = .false.

  ! begin LAPACK operations
  if (f_obj % linear_system_solver .eq. "LAPACK_standard") then ! use standard LAPACK solver
   if (f_obj % banded) then ! banded matrix storage
    ! load banded storage matrix used by LAPACK (stores LU factors on output)
    f_obj % AF(1:f_obj % KL,:)=0._r8b; f_obj % AF(f_obj % KL+1:f_obj % LDAF,:)=A(1:f_obj % LDA,:)
    ! scale (if needed)
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

  ! compute descaled solution if needed
  if (f_obj % scaling) call f_obj % custom_descaling(B)

  ! error control
  if (f_obj % linear_system_solver .eq. "LAPACK_standard") then ! use standard LAPACK solver
   if (INFO.ne.0) then
     if (f_obj % out_warning) then
      write(f_obj % unit,*) "LAPACK Error: DGESV or DGBSV exited with an error code of",info,"."
     end if
     f_obj % LAPACK_error = .true.
     return ! recoverable error (if supported by external driver)
   end if
  else if (f_obj % linear_system_solver .eq. "LAPACK_expert") then ! Use expert LAPACK solver with scaling and iterative refinement
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
