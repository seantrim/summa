module Richards
  ! h=pressure head
  use Richards_exact,only:S_Richards_exact ! Maple functions
  use Richards_exact,only:C_hle0_exact,dhdt_exact ! SageMath functions (temporal terms)
  use Richards_exact,only:K_hle0_exact,dKdz_hle0_exact,dhdz_exact,d2hdz2_exact ! SageMath functions (spatial terms)
  use Richards_exact,only:S_Richards_exact_SageMath ! SageMath source term function
  use kind_params,only:i4b,r8b,r16b
  implicit none
  private
  public :: Richards_obj
  real(r8b),parameter :: pi=3.1415926535897932_r8b

  type,public :: Richards_base
    logical      :: FV   ! finite volumes if true, centred FD if false
    integer(i4b) :: unit ! file unit for output
    integer(i4b) :: nz   ! # of interior points
    integer(i4b) :: i    ! layer index
    real(r8b),allocatable    :: h(:)  ! current pressure head
    real(r8b),allocatable    :: h0(:) ! previous pressure head
    real(r8b),allocatable    :: zg(:) ! array for vertical grid points (positive up)
    real(r8b)    :: t    ! current time 
    real(r8b)    :: dz   ! grid spacing
    real(r8b)    :: dt   ! time step size
    real(r8b)    :: L    ! domain length
   contains
    ! setup
    procedure :: allocate_memory => Richards_allocate_memory
    procedure :: deallocate_memory => Richards_deallocate_memory
    procedure :: create_grid =>  Richards_create_grid
    ! output
    procedure :: plot_h => Richards_plot_h
  end type Richards_base

  type,extends(Richards_base),public :: Richards_input_functions
    ! van Genuchten parameters
    real(r8b) :: alpha,theta_s,theta_r,n,K_s
    ! manufactured solution parameters
    real(r8b) :: A0,A1,q0,q1,t1,hBC
   contains
   ! ! **** For Test Richards Problem (using SageMath) ****
   ! ! pressure head procedures
   ! procedure :: initialize_h => Richards_test_initialize_h
   ! ! hydraulic conductivity
   ! procedure :: K => Richards_test_K
   ! procedure :: dKdh => Richards_test_dKdh
   ! ! water retention capacity
   ! procedure :: C => Richards_test_C
   ! procedure :: dCdh => Richards_test_dCdh
   ! ! source term
   ! procedure :: S => Richards_test_S

   ! ! **** For Celia Problem ****
   ! ! pressure head procedures
   ! procedure :: initialize_h => Richards_Celia_initialize_h
   ! procedure :: h_BCs => Richards_Celia_h_BCs
   ! ! hydraulic conductivity
   ! procedure :: K => Richards_Celia_K
   ! procedure :: dKdh => Richards_Celia_dKdh
   ! ! water retention capacity
   ! procedure :: C => Richards_Celia_C
   ! procedure :: dCdh => Richards_Celia_dCdh
   ! ! source term
   ! procedure :: S => Richards_Celia_S


   ! ! **** Casulli and Zanolli (2010) ****
   ! ! pressure head procedures
   ! procedure :: initialize_h => Richards_CZ2010_initialize_h
   ! procedure :: h_BCs => Richards_CZ2010_h_BCs
   ! procedure :: hd => h_transition ! transitional pressure head value 
   ! ! hydraulic conductivity
   ! procedure :: K => Richards_CZ2010_K
   ! procedure :: dKdh => Richards_CZ2010_dKdh
   ! ! water retention capacity
   ! procedure :: C => Richards_CZ2010_C
   ! procedure :: dCdh => Richards_CZ2010_dCdh
   ! ! source term
   ! procedure :: S => Richards_CZ2010_S

    ! **** Manufactured Richards Problem using Maple ****
    ! ** same parameterization as CZ2010 but different BCs and IC **
    ! pressure head procedures
    procedure :: h_exact => Richards_exact_h
    procedure :: initialize_h => Richards_exact_initialize_h
    procedure :: h_BCs => Richards_exact_h_BCs
    procedure :: h_error => Richards_h_error
    !procedure :: hd => h_transition ! transitional pressure head value 
    ! hydraulic conductivity
    procedure :: K => Richards_CZ2010_K
    procedure :: dKdh => Richards_CZ2010_dKdh
    ! water retention capacity
    procedure :: C => Richards_CZ2010_C
    procedure :: dCdh => Richards_CZ2010_dCdh
    ! source term
    !procedure :: S => Richards_exact_S ! Maple
    procedure :: S => Richards_exact_S_SageMath ! SageMath
    ! SageMath functions
    procedure :: S_exact => Richards_exact_S_SageMath
    procedure :: C_exact => Richards_exact_C_SageMath
    procedure :: dhdt_exact => Richards_exact_dhdt_SageMath
    procedure :: K_exact => Richards_exact_K_SageMath
    procedure :: dKdz_exact => Richards_exact_dKdz_SageMath
    procedure :: dhdz_exact => Richards_exact_dhdz_SageMath
    procedure :: d2hdz2_exact => Richards_exact_d2hdz2_SageMath

  end type Richards_input_functions

  type,extends(Richards_input_functions),public :: Richards_discrete_form
   contains
    ! objective function terms -- depends on form of Richards equation and discretization methods
    procedure :: CT => Richards_hform_BE_CT
    procedure :: dCTdh => Richards_hform_BE_dCTdh
    procedure :: KT => Richards_hform_FD2_KT
    procedure :: dKTdh => Richards_hform_FD2_dKTdh
  end type Richards_discrete_form

  type,extends(Richards_discrete_form),public :: Richards_type
   contains
   ! extra utility functions here 

  end type Richards_type

  type(Richards_type) :: Richards_obj

 contains

  ! ************************************************** General Procedures ************************************************** !

  subroutine Richards_allocate_memory(Richards_obj)
   class(Richards_base),intent(inout) :: Richards_obj
   associate(nz => Richards_obj % nz)
    allocate(Richards_obj % h(0:nz+1),Richards_obj % h0(0:nz+1),Richards_obj % zg(0:nz+1)) ! allocate pressure head arrays
   end associate
  end subroutine Richards_allocate_memory

  subroutine Richards_deallocate_memory(Richards_obj)
   class(Richards_base),intent(inout) :: Richards_obj
   associate(nz => Richards_obj % nz)
    deallocate(Richards_obj % h,Richards_obj % h0,Richards_obj % zg) ! deallocate pressure head arrays
   end associate
  end subroutine Richards_deallocate_memory

  subroutine Richards_create_grid(Richards_obj,zbot,ztop)
   ! *** Create uniform grid for 1D Richards problem ***
   class(Richards_base),intent(inout) :: Richards_obj
   real(r8b),intent(in)       :: zbot,ztop
   ! local variables
   integer(i4b)               :: i
   real(r8b)                  :: z

   associate(zg => Richards_obj % zg, nz => Richards_obj % nz, dz => Richards_obj % dz)
    if (Richards_obj % FV) then ! finite volumes

     dz=(ztop-zbot)/real(nz,r8b) ! nz=# of volumes
     zg(0)=zbot ! bottom BC

     do i=1,nz
      z=real(i-1,r8b)*dz+dz/2._r8b ! z at volume centres
      zg(i)=z
     end do

     zg(nz+1)=ztop ! top BC

    else ! centred FD

     dz=(ztop-zbot)/real(nz+1_i4b,r8b) !nz=# of interior points
     zg(0)=zbot ! bottom BC
     do i=1,nz
      z=real(i,r8b)*dz ! z at interior points (i.e., volume faces)
      zg(i)=z
     end do
     zg(nz+1)=ztop ! top BC

    end if
   end associate

   Richards_obj % L = ztop-zbot
  end subroutine Richards_create_grid

  subroutine Richards_plot_h(Richards_obj,fname)
   class(Richards_base),intent(in) :: Richards_obj
   character(11), intent(in) :: fname
   integer(i4b) :: i 
   open(unit=101,file=fname)
   do i=0,Richards_obj % nz+1
    write(101,*) Richards_obj % zg(i),Richards_obj % h(i)
   end do
   close(101)
  end subroutine Richards_plot_h

  subroutine Richards_h_error(Richards_obj)
   ! *** Compute error statistics for the pressure head relative to the exact solution ***
   ! note: assumes exact solution is available
   ! input
   class(Richards_input_functions),intent(in) :: Richards_obj
   ! output
   real(r8b)    :: error_mean ! mean error for pressure head
   real(r8b)    :: error_max  ! max error for pressure head
   integer(i4b) :: error_maxloc ! index value for max error
   ! local variables
   integer(i4b) :: i ! mesh index
   real(r8b)    :: h_exact ! exact pressure head value
   real(r8b),allocatable :: error(:)
   allocate(error(0:Richards_obj % nz+1))
   do i=0,Richards_obj % nz+1
    h_exact=Richards_obj % h_exact(i)
    if ((Richards_obj % h(i).ne.0._r8b).and.(h_exact.ne.0._r8b)) then
     error(i)=abs((Richards_obj % h(i)-h_exact)/h_exact) ! relative error
    else
     error(i)=abs(Richards_obj % h(i)-h_exact) ! absolute error
    end if
   end do
   error_mean=sum(error)/real(Richards_obj % nz+2,r8b) ! mean error (proportional to 1-norm)
   error_max =maxval(error) ! max relative error (infinity norm)
   error_maxloc=maxloc(error,DIM=1)-1_i4b
   write(Richards_obj % unit,'(a39,2(g25.16),i6)') "Mean Error, Max Error, Max Error Index=",error_mean,error_max,error_maxloc
  end subroutine Richards_h_error

  ! ****** General Discretization Procedures ******

  real(r8b) function Richards_hform_BE_CT(Richards_obj) result(CT)
   ! water retention capacity term
   ! BE time integration
   class(Richards_discrete_form),intent(in) :: Richards_obj
   associate(dt => Richards_obj % dt, h0 => Richards_obj % h0, h => Richards_obj % h, i => Richards_obj % i)
    CT=Richards_obj % C(h(i))*(h(i)-h0(i))/dt
   end associate
  end function Richards_hform_BE_CT

  real(r8b) function Richards_hform_BE_dCTdh(Richards_obj,j) result(dCTdh)
   ! water retention capacity term derivative WRT h(i)
   class(Richards_discrete_form),intent(in) :: Richards_obj
   integer(i4b),intent(in)         :: j ! derivative WRT h(j)
   associate(dt => Richards_obj % dt, h0 => Richards_obj % h0, h => Richards_obj % h, i => Richards_obj % i)
    if (j.eq.i) then ! WRT h(i)
     dCTdh=Richards_obj % dCdh(h(i))*(h(i)-h0(i))/dt + Richards_obj % C(h(i))/dt 
    else ! if h(j) is outside of the stencil
     dCTdh=0._r8b
    end if
   end associate
  end function Richards_hform_BE_dCTdh

  real(r8b) function Richards_hform_FD2_KT(Richards_obj) result(KT)
   ! hydraulic conductivity term
   ! second-order centred FD
   class(Richards_discrete_form),intent(in) :: Richards_obj
   associate(dz => Richards_obj % dz, h => Richards_obj % h, i => Richards_obj % i)
    KT=(Richards_obj % K(h(i+1))-Richards_obj % K(h(i-1)))/(2._r8b*dz)*((h(i+1)-h(i-1))/(2._r8b*dz)+1._r8b)&
      & + Richards_obj % K(h(i))*(h(i+1)-2_i4b*h(i)+h(i-1))/(dz**2_i4b)
   end associate
  end function Richards_hform_FD2_KT

  real(r8b) function Richards_hform_FD2_dKTdh(Richards_obj,j) result(dKTdh)
   ! hydraulic conductivity term derivative WRT h(i)
   class(Richards_discrete_form),intent(in) :: Richards_obj
   integer(i4b),intent(in)         :: j ! derivative WRT h(j)
   associate(dz => Richards_obj % dz, h => Richards_obj % h, i => Richards_obj % i)
    if (j.eq.i) then ! WRT h(i)
     dKTdh=Richards_obj % dKdh(h(i))*(h(i+1)-2_i4b*h(i)+h(i-1))/(dz**2_i4b)&
          & + Richards_obj % K(h(i))*(-2_i4b)/(dz**2_i4b)
    elseif (j.eq.(i-1)) then ! WRT h(i-1)
     dKTdh=-Richards_obj % dKdh(h(i-1))/(2._r8b*dz)*((h(i+1)-h(i-1))/(2._r8b*dz)+1._r8b)&
          &+ (Richards_obj % K(h(i+1))-Richards_obj % K(h(i-1)))/(2._r8b*dz)*(-1._r8b/(2._r8b*dz))&
          &+ Richards_obj % K(h(i))*(1._r8b/(dz**2_i4b))
    elseif (j.eq.(i+1)) then ! WRT h(i+1)
     dKTdh=Richards_obj % dKdh(h(i+1))/(2._r8b*dz)*((h(i+1)-h(i-1))/(2._r8b*dz)+1._r8b)&
          &+ (Richards_obj % K(h(i+1))-Richards_obj % K(h(i-1)))/(2._r8b*dz)*(1._r8b/(2._r8b*dz))&
          &+ Richards_obj % K(h(i))*(1._r8b/(dz**2_i4b))
    else ! if h(j) is outside of the stencil
     dKTdh=0._r8b
    end if
   end associate
  end function Richards_hform_FD2_dKTdh

  ! ****** For non-Maple Richards test problem (does not use the full van Genuchten model) ******

  real(r8b) function Richards_test_h_exact(Richards_obj,i) result(h_exact)
   ! exact formula for h (from the manufactured solution of the discrete system)
   class(Richards_base),intent(in) :: Richards_obj
   integer(i4b) :: i ! grid index
   real(r8b)                       :: z
   associate(t => Richards_obj % t)
    z=Richards_obj % zg(i)
    h_exact=-sin(pi*z)*exp(-t/1.e5_r8b)
   end associate
  end function Richards_test_h_exact

  real(r8b) function Richards_test_K(Richards_obj,hval) result(K)
   ! hydraulic conductivity
   class(Richards_input_functions),intent(in) :: Richards_obj
   real(r8b),intent(in) :: hval ! pressure head value
   !K=(sqrt(-1._r8b/(h**2 + 1._r8b) + 1._r8b) - 1._r8b)**2_i4b/sqrt(sqrt(h**2_i4b + 1._r8b)) ! van Genuchten model
   K=(8._r8b*hval+9._r8b)/10._r8b ! simplified but similar to van Genuchten model
  end function Richards_test_K

  real(r8b) function Richards_test_dKdh(Richards_obj,hval) result(dKdh)
   ! hydraulic conductivity derivative WRT h
   class(Richards_input_functions),intent(in) :: Richards_obj
   real(r8b) :: hval ! pressure head value  
   dKdH=8._r8b/10._r8b ! simplified but similar to van Genuchten model
  end function Richards_test_dKdh

  real(r8b) function Richards_test_C(Richards_obj,hval) result(C)
   ! water retention capacity
   class(Richards_input_functions),intent(in) :: Richards_obj
   real(r8b) :: hval ! pressure head value  
   C=-hval*(sqrt(1._r8b+hval**2_i4b))**(-3_i4b)
  end function Richards_test_C

  real(r8b) function Richards_test_dCdh(Richards_obj,hval) result(dCdh)
   ! water retention capacity derivative WRT h
   class(Richards_input_functions),intent(in) :: Richards_obj
   real(r8b) :: hval ! pressure head value  
   dCdh=3._r8b*hval**2_i4b/(hval**2_i4b + 1._r8b)**(5._r8b/2._r8b) - 1._r8b/(hval**2_i4b + 1._r8b)**(3._r8b/2._r8b)
  end function Richards_test_dCdh

  real(r8b) function Richards_test_S(Richards_obj) result(S)
   ! discrete source term (from SageMath script -- factored form used)
   class(Richards_input_functions),intent(in) :: Richards_obj
   real(r8b)                       :: z   
   associate(dz => Richards_obj % dz, dt => Richards_obj % dt, t => Richards_obj % t, i => Richards_obj % i)
    z=Richards_obj % zg(i)
    S=1._r8b/5._r8b*(exp(-1._r8b/100000._r8b*t)*sin(pi*dz + pi*z)&
     & - exp(-1._r8b/100000._r8b*t)*sin(-pi*dz + pi*z))*((exp(-1._r8b/100000._r8b*t)*sin(pi*dz + pi*z)&
     & - exp(-1._r8b/100000._r8b*t)*sin(-pi*dz + pi*z))/dz - 2._r8b)/dz - (exp(1._r8b/100000._r8b*dt - 1._r8b/100000._r8b*t)*sin(pi*z)&
     & - exp(-1._r8b/100000._r8b*t)*sin(pi*z))*exp(-1._r8b/100000._r8b*t)&
     &*sin(pi*z)/((exp(-1._r8b/50000._r8b*t)*sin(pi*z)**2_i4b + 1._r8b)**(3._r8b/2._r8b)*dt)&
     & + 1._r8b/10._r8b*(exp(-1._r8b/100000._r8b*t)*sin(pi*dz + pi*z) + exp(-1._r8b/100000._r8b*t)*sin(-pi*dz + pi*z)&
     & - 2._r8b*exp(-1._r8b/100000._r8b*t)*sin(pi*z))*(8._r8b*exp(-1._r8b/100000._r8b*t)*sin(pi*z) - 9._r8b)/dz**2_i4b
   end associate
  end function Richards_test_S

  subroutine Richards_test_initialize_h(Richards_obj)
   ! *** Initialize the pressure head arrays with values for the previous time step ***
   class(Richards_input_functions),intent(inout) :: Richards_obj
   real(r8b)                  :: t0 ! time from previous time step
   real(r8b)                  :: t_save ! original time value (i.e., current time)
   integer(i4b)               :: i

   associate(nz => Richards_obj % nz, dz => Richards_obj % dz, dt => Richards_obj % dt, h0 => Richards_obj % h0, h => Richards_obj % h,&
            & t => Richards_obj %  t)
    t_save=t ! save current time value
    t=t-dt   ! use time from previous time step for h_exact() function
    h0(0)=0._r8b ! bottom BC
    do i=1,nz
     h0(i)=Richards_test_h_exact(Richards_obj,i)!Richards_obj % h_exact(i)
    end do
    h0(nz+1)=0._r8b ! top BC
    h=h0 ! use h0 as a default initial guess for h
    t=t_save ! restore current time value
   end associate
  end subroutine Richards_test_initialize_h

  function Richards_test_expected_h(Richards_obj) result(h)
   ! *** Initialize the pressure head arrays with values for the previous time step ***
   class(Richards_type),intent(inout) :: Richards_obj
   real(r8b)                  :: h(0:Richards_obj % nz+1)
   real(r8b)                  :: t0 ! time from previous time step
   integer(i4b)               :: i

   associate(nz => Richards_obj % nz, dz => Richards_obj % dz, dt => Richards_obj % dt,t => Richards_obj %  t)
    h(0)=0._r8b ! bottom BC
    do i=1,nz
     !Richards_obj % i =i
     h(i)=Richards_test_h_exact(Richards_obj,i)!Richards_obj % h_exact(i)
    end do
    h(nz+1)=0._r8b ! top BC
   end associate
  end function Richards_test_expected_h

  ! ************************** Class Procedures for the Celia Problem ************************** !

  real(r8b) function Richards_Celia_C(Richards_obj,hval) result(C)
   ! hydraulic conductivity (assumes h<0)
   class(Richards_input_functions),intent(in) :: Richards_obj
   real(r8b),intent(in) :: hval ! pressure head value
   real(r8b),parameter :: alpha=1.611e6_r8b,theta_s=0.287_r8b,theta_r=0.075_r8b,beta=3.96_r8b 
   C=alpha*(theta_s-theta_r)*beta*((-hval)**(beta-1._r8b))/(alpha+(-hval)**beta)**2_i4b ! Celia Eq. 10a
  end function Richards_Celia_C

  real(r8b) function Richards_Celia_dCdh(Richards_obj,hval) result(dCdh)
   ! hydraulic conductivity (assumes h<0)
   class(Richards_input_functions),intent(in) :: Richards_obj
   real(r8b),intent(in) :: hval ! pressure head value
   real(r8b),parameter :: alpha=1.611e6_r8b,theta_s=0.287_r8b,theta_r=0.075_r8b,beta=3.96_r8b 
   dCdh=(2._r8b*alpha*(theta_s-theta_r)*beta**2_i4b)*((-hval)**(2._r8b*beta-2._r8b))/(alpha+(-hval)**beta)**3_i4b&
       &-alpha*(theta_s-theta_r)*beta*(beta-1._r8b)*((-hval)**(beta-2._r8b))/(alpha+(-hval)**beta)**2_i4b ! Celia Eq. 10a
  end function Richards_Celia_dCdh

  real(r8b) function Richards_Celia_K(Richards_obj,hval) result(K)
   ! hydraulic conductivity
   class(Richards_input_functions),intent(in) :: Richards_obj
   real(r8b),intent(in) :: hval ! pressure head value
   real(r8b),parameter :: Ks=0.00944_r8b,A=1.175e6_r8b,gam=4.74_r8b
   K=Ks*A/(A+abs(hval)**gam) ! Celia Eq. 10b
  end function Richards_Celia_K

  real(r8b) function Richards_Celia_dKdh(Richards_obj,hval) result(dKdh)
   ! hydraulic conductivity derivative WRT h (h<0 assumed)
   class(Richards_input_functions),intent(in) :: Richards_obj
   real(r8b) :: hval ! pressure head value  
   real(r8b),parameter :: Ks=0.00944_r8b,A=1.175e6_r8b,gam=4.74_r8b
   dKdH=Ks*A*gam*((-hval)**(gam-1._r8b))/(A+(-hval)**gam)**2_i4b 
  end function Richards_Celia_dKdh

  real(r8b) function Richards_Celia_S(Richards_obj) result(S)
   ! discrete source term (from SageMath script -- factored form used)
   class(Richards_input_functions),intent(in) :: Richards_obj
    S=0._r8b
  end function Richards_Celia_S

  subroutine Richards_Celia_initialize_h(Richards_obj)
   ! *** Initialize the pressure head arrays with values for the previous time step ***
   class(Richards_input_functions),intent(inout) :: Richards_obj
   integer(i4b)               :: i
   real(r8b),parameter        :: t=0._r8b ! initial time         

   call Richards_obj % h_BCs(t) ! enforce BCs
   associate(nz => Richards_obj % nz,  h0 => Richards_obj % h0, h => Richards_obj % h)            
    h0(0)=h(0) ! bottom BC
    do i=1,nz
     h0(i)=h0(0)
    end do
    h0(nz+1)=h(nz+1) ! top BC
    h=h0 ! use h0 as a default initial guess for h
   end associate
  end subroutine Richards_Celia_initialize_h

  subroutine Richards_Celia_h_BCs(Richards_obj,t)
   class(Richards_input_functions),intent(inout) :: Richards_obj
   real(r8b),intent(in) :: t
   associate(nz => Richards_obj % nz)
    Richards_obj % h(0)    = -61.5_r8b ! [cm]
    Richards_obj % h(nz+1) = -20.7_r8b ! [cm]
   end associate
  end subroutine Richards_Celia_h_BCs

  ! ******** Casulli and Zanolli (2010) ******** !

  subroutine Richards_CZ2010_h_BCs(Richards_obj,t)
   class(Richards_input_functions),intent(inout) :: Richards_obj
   real(r8b),intent(in) :: t
   associate(nz => Richards_obj % nz)
    Richards_obj % h(0)    = 0._r8b  ! [m]
    if (t.le.1.e5_r8b) then
     Richards_obj % h(nz+1) = -0.05_r8b + 0.03_r8b*sin(2._r8b*pi*t/1.e5_r8b)
    elseif (t.le.1.8e5_r8b) then
     Richards_obj % h(nz+1) = 0.1_r8b ! [m]
    else
     Richards_obj % h(nz+1) = -0.05_r8b + 2952.45_r8b*exp(-t/18204.8_r8b)
    end if
   end associate
  end subroutine Richards_CZ2010_h_BCs

  subroutine Richards_CZ2010_initialize_h(Richards_obj)
   ! *** Initialize the pressure head arrays with values for the previous time step ***
   class(Richards_input_functions),intent(inout) :: Richards_obj
   integer(i4b)               :: i
   real(r8b),parameter        :: t=0._r8b ! initial time         

   associate(nz => Richards_obj % nz,  h0 => Richards_obj % h0, h => Richards_obj % h)            
    do i=1,nz
     !h0(i)=Richards_obj % zg(i)-2._r8b
     h(i)=-Richards_obj % zg(i)
    end do
    call Richards_obj % h_BCs(t) ! enforce BCs
    h0=h ! use h0 as a default initial guess for h
   end associate
  end subroutine Richards_CZ2010_initialize_h

  real(r8b) function Richards_CZ2010_K(Richards_obj,h) result(K)
   ! *** Hydraulic conductivity for CZ2010 test ***
   class(Richards_input_functions),intent(in) :: Richards_obj
   real(r8b) :: h ! pressure head value
   real(r8b) :: m
   associate(n => Richards_obj % n, K_s => Richards_obj % K_s, alpha => Richards_obj % alpha)
    m=1._r8b-1._r8b/n
    if (h.le.0._r8b) then
     K=K_s*((-(1._r8b/(((-alpha*h)**n + 1._r8b)**m))**(1._r8b/m) + 1._r8b)**m - 1._r8b)**2_i4b/sqrt(((-alpha*h)**n + 1._r8b)**m)
    else
     K=K_s
    end if
   end associate
  end function Richards_CZ2010_K

  real(r8b) function Richards_CZ2010_dKdh(Richards_obj,h) result(dKdh)
   ! *** Derivative of  Hydraulic conductivity for CZ2010 test ***
   class(Richards_input_functions),intent(in) :: Richards_obj
   real(r8b) :: h ! pressure head value
   real(r8b) :: m
   associate(n => Richards_obj % n, K_s => Richards_obj % K_s, alpha => Richards_obj % alpha)
    !if (h.le.0._r8b) then ! OG
    if (h.lt.0._r8b) then ! SageMath gives dKdh(h=0)=0 after simplifying
     m=1._r8b-1._r8b/n

     ! OG -- appears to be from the factored form from SageMath
     !dKdh=-1._r8b/2._r8b*(5._r8b*(-(1._r8b/(((-alpha*h)**n + 1._r8b)**m))**(1._r8b/m) + 1._r8b)**m&
     !    &*(1._r8b/(((-alpha*h)**n + 1._r8b)**m))**(1._r8b/m) - (-(1._r8b/(((-alpha*h)**n + 1._r8b)**m))**(1._r8b/m)&
     !    & + 1._r8b)**m - (1._r8b/(((-alpha*h)**n + 1._r8b)**m))**(1._r8b/m) + 1._r8b)*(-alpha*h)**n&
     !    &*K_s*m*n*((-(1._r8b/(((-alpha*h)**n + 1._r8b)**m))**(1._r8b/m) + 1._r8b)**m - 1._r8b)/&
     !    &(((-alpha*h)**n + 1._r8b)*sqrt(((-alpha*h)**n + 1._r8b)**m)*h*((1._r8b/(((-alpha*h)**n + 1._r8b)**m))**(1._r8b/m) - 1._r8b))
     ! NaN result at h=0: division by zero when h=0 due to last line (factor of h near the centre)
     ! ** may also be a negative value raised to a fractional power

     ! unfactored form from SageMath
     dKdh=-2._r8b*(-alpha*h)**(n - 1._r8b)*((-alpha*h)**n + 1._r8b)**(m - 1._r8b)*(1._r8b/(((-alpha*h)**n + 1._r8b)**m))&
         &**(1._r8b/m - 1._r8b)*K_s*alpha*m*n*(-(1._r8b/(((-alpha*h)**n + 1._r8b)**m))**(1._r8b/m) + 1._r8b)**(m - 1._r8b)*&
         &((-(1._r8b/(((-alpha*h)**n + 1._r8b)**m))**(1._r8b/m) + 1._r8b)**m - 1._r8b)/(((-alpha*h)**n + 1._r8b)**m)&
         &**(5._r8b/2._r8b) + 1._r8b/2._r8b*(-alpha*h)**(n - 1._r8b)*((-alpha*h)**n + 1._r8b)**(m - 1._r8b)*K_s*alpha*m*n*&
         &((-(1._r8b/(((-alpha*h)**n + 1._r8b)**m))**(1._r8b/m) + 1._r8b)**m - 1._r8b)**2/(((-alpha*h)**n + 1._r8b)**m)&
         &**(3._r8b/2._r8b)
    else
     dKdh=0._r8b
    end if
   end associate
  end function Richards_CZ2010_dKdh

  real(r8b) function Richards_CZ2010_C(Richards_obj,h) result(C)
   ! *** water retention capacity for CZ2010 test ***
   class(Richards_input_functions),intent(in) :: Richards_obj
   real(r8b) :: h ! pressure head value
   real(r8b) :: m
   associate(n => Richards_obj % n, alpha => Richards_obj % alpha, theta_s => Richards_obj % theta_s, theta_r => Richards_obj % theta_r)
    m=1._r8b-1._r8b/n
    if (h.le.0._r8b) then
     !C=(-alpha*h)**n*((-alpha*h)**n + 1._r8b)**(-m - 1._r8b)*m*n*(theta_r - theta_s)/h !! OG from SageMath
     C=-(alpha)**n*(-h)**(n-1._r8b)*((-alpha*h)**n + 1._r8b)**(-m - 1._r8b)*m*n*(theta_r - theta_s) ! simplified to avoid division by zero
    else
     C=0._r8b
    end if
   end associate
  end function Richards_CZ2010_C

  real(r8b) function Richards_CZ2010_dCdh(Richards_obj,h) result(dCdh)
   ! *** water retention capacity for CZ2010 test ***
   class(Richards_input_functions),intent(in) :: Richards_obj
   real(r8b) :: h ! pressure head value
   real(r8b) :: m
   associate(n => Richards_obj % n, alpha => Richards_obj % alpha, theta_s => Richards_obj % theta_s, theta_r => Richards_obj % theta_r)
    m=1._r8b-1._r8b/n
    if (h.le.0._r8b) then ! rearranged to avoid division by zero when h=0
     !dCdh=-((-alpha*h)**n*m*n + (-alpha*h)**n - n + 1._r8b)*(-alpha*h)**n*((-alpha*h)**n + 1._r8b)**(-m - 2._r8b)*m*n*(theta_r - theta_s)/h**2_i4b
     dCdh=-((-alpha*h)**n*m*n + (-alpha*h)**n - n + 1._r8b)*(alpha)**n*(-h)**(n-2._r8b)*((-alpha*h)**n + 1._r8b)**(-m - 2._r8b)*m*n*(theta_r - theta_s)
    else
     dCdh=0._r8b
    end if
   end associate
  end function Richards_CZ2010_dCdh

  real(r8b) function Richards_CZ2010_S(Richards_obj) result(S)
   ! *** Source term for CZ2010 test ***
   class(Richards_input_functions),intent(in) :: Richards_obj
   S=0._r8b
  end function Richards_CZ2010_S

  real(r8b) function h_transition(Richards_obj)
   ! *** Transitional pressure head value used for Jordan decomposition in CZ2010 ***
   class(Richards_input_functions),intent(in) :: Richards_obj
   real(r8b) :: m
   m=1._r8b-1._r8b/Richards_obj % n
   h_transition=-((Richards_obj % n-1.0_r8b)/(Richards_obj % n*m+1.0_r8b))**(1.0_r8b/ Richards_obj % n)/ Richards_obj % alpha
  end function h_transition

  ! ******** Manufactured Richards Problem using Maple ******** !
 
  subroutine Richards_exact_h_BCs(Richards_obj,t)
   class(Richards_input_functions),intent(inout) :: Richards_obj
   real(r8b),intent(in) :: t
   associate(nz => Richards_obj % nz)
    Richards_obj % h(0)    = Richards_obj % hBC ! [m]
    Richards_obj % h(nz+1) = Richards_obj % hBC ! [m] 
   end associate
  end subroutine Richards_exact_h_BCs

  function Richards_exact_h(Richards_obj,i) result(h)
   ! *** Exact pressure head solution ***
   class(Richards_input_functions),intent(in) :: Richards_obj
   integer(i4b) :: i ! grid index

   real(r8b) :: z ! vertical position [m]
   real(r8b) :: h ! pressure head [m]
   ! amplitude
   real(r8b) :: A
   ! exponent
   real(r8b) :: q

   associate(t => Richards_obj % t, L => Richards_obj % L,&
            &A0 => Richards_obj % A0,A1 => Richards_obj % A1,q0 => Richards_obj % q0,q1 => Richards_obj % q1,&
            &t1 => Richards_obj % t1,hBC => Richards_obj % hBC)
    z=Richards_obj % zg(i)
    A = A0*exp(-log(A0/A1)*(t/t1)**2_i4b)
    q = q0 + (q1-q0)*(t/t1)**3_i4b
    h = -A*sin(pi*((L-z)/L)**q)+hBC
   end associate
  end function Richards_exact_h

  subroutine Richards_exact_initialize_h(Richards_obj) 
   ! *** Initialize the pressure head arrays with values for the previous time step ***
   class(Richards_input_functions),intent(inout) :: Richards_obj
   integer(i4b)               :: i
   real(r8b)                  :: t_save

   associate(nz => Richards_obj % nz, dt => Richards_obj % dt, h0 => Richards_obj % h0, h => Richards_obj % h,&
            & t => Richards_obj %  t)
    t_save=t ! save current time value
    t=t-dt   ! use time from previous time step for h_exact() function
    do i=1,nz
     h(i)=Richards_obj % h_exact(i)
    end do
    call Richards_obj % h_BCs(t) ! enforce BCs
    h0=h ! use h0 as a default initial guess for h
    t=t_save ! restore current time value
   end associate
  end subroutine Richards_exact_initialize_h

  real(r8b) function Richards_exact_S(Richards_obj) result(S)
   ! *** Source term for manufactured Richards problem using Maple ***
   class(Richards_input_functions),intent(in) :: Richards_obj
   real(r8b) :: z ! vertical position [m] (positive up)
   associate(alpha => Richards_obj % alpha, theta_s => Richards_obj % theta_s, theta_r => Richards_obj % theta_r,&
            &n => Richards_obj % n, Ks => Richards_obj % K_s, t => Richards_obj % t, L => Richards_obj % L,&
            & i => Richards_obj % i,&
            &A0 => Richards_obj % A0,A1 => Richards_obj % A1,q0 => Richards_obj % q0,q1 => Richards_obj % q1,&
            &t1 => Richards_obj % t1,hBC => Richards_obj % hBC)
    z=Richards_obj % zg(i)
    S=S_Richards_exact(z,t,L,alpha,n,Ks,theta_r,theta_s,A0,A1,q0,q1,t1)
   end associate
  end function Richards_exact_S

  ! *** SageMath Functions for Manufactured Richards Problem

  real(r8b) function Richards_exact_C_SageMath(Richards_obj) result(C)
   ! *** C for manufactured Richards problem using SageMath ***
   class(Richards_input_functions),intent(in) :: Richards_obj
   real(r8b) :: z ! vertical position [m] (positive up)
   associate(alpha => Richards_obj % alpha, theta_s => Richards_obj % theta_s, theta_r => Richards_obj % theta_r,&
            &n => Richards_obj % n, t => Richards_obj % t, L => Richards_obj % L,&
            & i => Richards_obj % i,&
            &A0 => Richards_obj % A0,A1 => Richards_obj % A1,q0 => Richards_obj % q0,q1 => Richards_obj % q1,&
            &t1 => Richards_obj % t1,hBC => Richards_obj % hBC)
    z=Richards_obj % zg(i)
    C=C_hle0_exact(z,t,L,alpha,n,theta_r,theta_s,A0,A1,q0,q1,t1,hBC)
   end associate
  end function Richards_exact_C_SageMath

  real(r8b) function Richards_exact_dhdt_SageMath(Richards_obj) result(dhdt)
   ! *** dhdt for manufactured Richards problem using SageMath ***
   class(Richards_input_functions),intent(in) :: Richards_obj
   real(r8b) :: z ! vertical position [m] (positive up)
   associate(t => Richards_obj % t, L => Richards_obj % L,i => Richards_obj % i,&
            &A0 => Richards_obj % A0,A1 => Richards_obj % A1,q0 => Richards_obj % q0,q1 => Richards_obj % q1,&
            &t1 => Richards_obj % t1,hBC => Richards_obj % hBC)
    z=Richards_obj % zg(i)
    dhdt=dhdt_exact(z,t,L,A0,A1,q0,q1,t1)
   end associate
  end function Richards_exact_dhdt_SageMath

  real(r8b) function Richards_exact_K_SageMath(Richards_obj) result(K)
   ! *** K for manufactured Richards problem using SageMath ***
   class(Richards_input_functions),intent(in) :: Richards_obj
   real(r8b) :: z ! vertical position [m] (positive up)
   associate(alpha => Richards_obj % alpha, n => Richards_obj % n, K_s => Richards_obj % K_s,&
            &t => Richards_obj % t, L => Richards_obj % L, i => Richards_obj % i,&
            &A0 => Richards_obj % A0,A1 => Richards_obj % A1,q0 => Richards_obj % q0,q1 => Richards_obj % q1,&
            &t1 => Richards_obj % t1,hBC => Richards_obj % hBC)
    z=Richards_obj % zg(i)
    K=K_hle0_exact(z,t,L,alpha,n,K_s,A0,A1,q0,q1,t1,hBC)
   end associate
  end function Richards_exact_K_SageMath

  real(r8b) function Richards_exact_dKdz_SageMath(Richards_obj) result(dKdz)
   ! *** K for manufactured Richards problem using SageMath ***
   class(Richards_input_functions),intent(in) :: Richards_obj
   real(r8b) :: z ! vertical position [m] (positive up)
   associate(alpha => Richards_obj % alpha, n => Richards_obj % n, K_s => Richards_obj % K_s,&
            &t => Richards_obj % t, L => Richards_obj % L, i => Richards_obj % i,&
            &A0 => Richards_obj % A0,A1 => Richards_obj % A1,q0 => Richards_obj % q0,q1 => Richards_obj % q1,&
            &t1 => Richards_obj % t1,hBC => Richards_obj % hBC)
    z=Richards_obj % zg(i)
    dKdz=dKdz_hle0_exact(z,t,L,alpha,n,K_s,A0,A1,q0,q1,t1,hBC)
   end associate
  end function Richards_exact_dKdz_SageMath

  real(r8b) function Richards_exact_dhdz_SageMath(Richards_obj) result(dhdz)
   ! *** dhdt for manufactured Richards problem using SageMath ***
   class(Richards_input_functions),intent(in) :: Richards_obj
   real(r8b) :: z ! vertical position [m] (positive up)
   associate(t => Richards_obj % t, L => Richards_obj % L,i => Richards_obj % i,&
            &A0 => Richards_obj % A0,A1 => Richards_obj % A1,q0 => Richards_obj % q0,q1 => Richards_obj % q1,&
            &t1 => Richards_obj % t1,hBC => Richards_obj % hBC)
    z=Richards_obj % zg(i)
    dhdz=dhdz_exact(z,t,L,A0,A1,q0,q1,t1)
   end associate
  end function Richards_exact_dhdz_SageMath

  real(r8b) function Richards_exact_d2hdz2_SageMath(Richards_obj) result(d2hdz2)
   ! *** dhdt for manufactured Richards problem using SageMath ***
   class(Richards_input_functions),intent(in) :: Richards_obj
   real(r8b) :: z ! vertical position [m] (positive up)
   associate(t => Richards_obj % t, L => Richards_obj % L,i => Richards_obj % i,&
            &A0 => Richards_obj % A0,A1 => Richards_obj % A1,q0 => Richards_obj % q0,q1 => Richards_obj % q1,&
            &t1 => Richards_obj % t1,hBC => Richards_obj % hBC)
    z=Richards_obj % zg(i)
    d2hdz2=d2hdz2_exact(z,t,L,A0,A1,q0,q1,t1)
   end associate
  end function Richards_exact_d2hdz2_SageMath

  real(r8b) function Richards_exact_S_SageMath(Richards_obj) result(S)
   ! *** S for manufactured Richards problem using SageMath ***
   class(Richards_input_functions),intent(in) :: Richards_obj
   real(r8b) :: z ! vertical position [m] (positive up)
   associate(alpha => Richards_obj % alpha, theta_s => Richards_obj % theta_s, theta_r => Richards_obj % theta_r,&
            &n => Richards_obj % n, Ks => Richards_obj % K_s, t => Richards_obj % t, L => Richards_obj % L,&
            &i => Richards_obj % i,&
            &A0 => Richards_obj % A0,A1 => Richards_obj % A1,q0 => Richards_obj % q0,q1 => Richards_obj % q1,&
            &t1 => Richards_obj % t1,hBC => Richards_obj % hBC)
    z=Richards_obj % zg(i)
    S=S_Richards_exact_SageMath(z,t,L,alpha,n,Ks,theta_r,theta_s,A0,A1,q0,q1,t1,hBC)
   end associate
  end function Richards_exact_S_SageMath

end module Richards
