module Richards_exact
  use iso_fortran_env, only: int32,int64,real32,real64,real128 !!integer and real kind parameters
  implicit none
  private
  ! Maple functions
  public :: S_Richards_exact
  ! SageMath functions
  public :: C_hle0_exact,dhdt_exact
  public :: K_hle0_exact,dKdz_hle0_exact,dhdz_exact,d2hdz2_exact
  public :: S_Richards_exact_SageMath

  ! kind parameters
  integer, parameter :: i4b=int32
  integer, parameter :: i8b=int64
  integer, parameter :: r4b=real32
  integer, parameter :: r8b=real64
  integer, parameter :: r16b=real128

  ! other parameters
  real(r8b), parameter :: pi=3.1415926535897932_r8b
 contains

  ! ************************************* SageMath Functions ************************************* !
  function C_hle0_exact(z,t,L,alpha,n,theta_r,theta_s,A0,A1,q0,q1,t1,h0) result(C)
   ! *** Water retention capacity for the manufactured Richards problem (from SageMath) *** 
   ! input
   real(r8b),intent(in) :: z,t ! vertical position [m] (positive up) and time [s]
   real(r8b),intent(in) :: L   ! domain length
   real(r8b),intent(in) :: alpha,n,theta_r,theta_s ! van Genuchten parameters
   real(r8b),intent(in) :: A0,A1,q0,q1,t1,h0 ! manufactured solution parameters

   ! output 
   real(r8b) :: C

   ! local variables
   real(r8b) :: m

   m=1._r8b-1._r8b/n

   C=-(-alpha*(h0 - A0*sin(pi*(-z/L + 1)**q0*(-z/L + 1)**(q1*t**3/t1**3)/(-z/L + 1)**(q0*t**3/t1**3))/(A0/A1)**(t**2/t1**2)))**(n - 1)*((-alpha*(h0 - A0*sin(pi*(-z/L + 1)**q0*(-z/L + 1)**(q1*t**3/t1**3)/(-z/L + 1)**(q0*t**3/t1**3))/(A0/A1)**(t**2/t1**2)))**n + 1)**(m - 1)*alpha*m*n*(theta_r - theta_s)/((-alpha*(h0 - A0*sin(pi*(-z/L + 1)**q0*(-z/L + 1)**(q1*t**3/t1**3)/(-z/L + 1)**(q0*t**3/t1**3))/(A0/A1)**(t**2/t1**2)))**n + 1)**(2*m)
  end function C_hle0_exact

  function dhdt_exact(z,t,L,A0,A1,q0,q1,t1) result(dhdt)
   ! *** Time derivative of pressure head for the manufactured Richards problem (from SageMath) ***
   ! input
   real(r8b),intent(in) :: z,t ! vertical position [m] (positive up) and time [s]
   real(r8b),intent(in) :: L   ! domain length
   real(r8b),intent(in) :: A0,A1,q0,q1,t1 ! manufactured solution parameters

   ! output
   real(r8b) :: dhdt

   dhdt=(3*(pi*A0*q0 - pi*A0*q1)*t**2*((L - z)/L)**(q0 + q1*t**3/t1**3)*cos(pi*(-z/L + 1)**q0*(-z/L + 1)**(q1*t**3/t1**3)/(-z/L + 1)**(q0*t**3/t1**3))*log((L - z)/L) + 2*A0*t*t1*((L - z)/L)**(q0*t**3/t1**3)*log(A0/A1)*sin(pi*(-z/L + 1)**q0*(-z/L + 1)**(q1*t**3/t1**3)/(-z/L + 1)**(q0*t**3/t1**3)))/(t1**3*(A0/A1)**(t**2/t1**2)*((L - z)/L)**(q0*t**3/t1**3))
  end function dhdt_exact

  function K_hle0_exact(z,t,L,alpha,n,K_s,A0,A1,q0,q1,t1,h0) result(K)
   ! *** Hydraulic conductivity for the manufactured Richards problem (from SageMath) ***
   ! input
   real(r8b),intent(in) :: z,t ! vertical position [m] (positive up) and time [s]
   real(r8b),intent(in) :: L   ! domain length
   real(r8b),intent(in) :: alpha,n,K_s ! van Genuchten parameters
   real(r8b),intent(in) :: A0,A1,q0,q1,t1,h0 ! manufactured solution parameters

   ! output
   real(r8b) :: K

   ! local variables
   real(r8b) :: m

   m=1._r8b-1._r8b/n

   K=K_s*((-(1/(((-alpha*(h0 - A0*sin(pi*(-z/L + 1)**q0*(-z/L + 1)**(q1*t**3/t1**3)/(-z/L + 1)**(q0*t**3/t1**3))/(A0/A1)**(t**2/t1**2)))**n + 1)**m))**(1/m) + 1)**m - 1)**2/sqrt(((-alpha*(h0 - A0*sin(pi*(-z/L + 1)**q0*(-z/L + 1)**(q1*t**3/t1**3)/(-z/L + 1)**(q0*t**3/t1**3))/(A0/A1)**(t**2/t1**2)))**n + 1)**m)
  end function K_hle0_exact

  function dKdz_hle0_exact(z,t,L,alpha,n,K_s,A0,A1,q0,q1,t1,h0) result(dKdz)
   ! *** Derivative of hydraulic conductivity with respect to z for the manufactured Richards problem (from SageMath) ***
   ! input
   real(r8b),intent(in) :: z,t ! vertical position [m] (positive up) and time [s]
   real(r8b),intent(in) :: L   ! domain length
   real(r8b),intent(in) :: alpha,n,K_s ! van Genuchten parameters
   real(r8b),intent(in) :: A0,A1,q0,q1,t1,h0 ! manufactured solution parameters

   ! output
   real(r8b) :: dKdz

   ! local variables
   real(r8b) :: m

   m=1._r8b-1._r8b/n

   dKdz=1._r8b/2._r8b*pi*(q0*t**3 - q1*t**3 - q0*t1**3)*(5*(-(1/(((-(h0*(A0/A1)**(t**2/t1**2) - A0*sin(pi*(-z/L + 1)**q0*(-z/L + 1)**(q1*t**3/t1**3)/(-z/L + 1)**(q0*t**3/t1**3)))*alpha/(A0/A1)**(t**2/t1**2))**n + 1)**m))**(1/m) + 1)**m*(1/(((-(h0*(A0/A1)**(t**2/t1**2) - A0*sin(pi*(-z/L + 1)**q0*(-z/L + 1)**(q1*t**3/t1**3)/(-z/L + 1)**(q0*t**3/t1**3)))*alpha/(A0/A1)**(t**2/t1**2))**n + 1)**m))**(1/m) - (-(1/(((-(h0*(A0/A1)**(t**2/t1**2) - A0*sin(pi*(-z/L + 1)**q0*(-z/L + 1)**(q1*t**3/t1**3)/(-z/L + 1)**(q0*t**3/t1**3)))*alpha/(A0/A1)**(t**2/t1**2))**n + 1)**m))**(1/m) + 1)**m - (1/(((-(h0*(A0/A1)**(t**2/t1**2) - A0*sin(pi*(-z/L + 1)**q0*(-z/L + 1)**(q1*t**3/t1**3)/(-z/L + 1)**(q0*t**3/t1**3)))*alpha/(A0/A1)**(t**2/t1**2))**n + 1)**m))**(1/m) + 1)*A0*K_s*m*n*(-(h0*(A0/A1)**(t**2/t1**2) - A0*sin(pi*(-z/L + 1)**q0*(-z/L + 1)**(q1*t**3/t1**3)/(-z/L + 1)**(q0*t**3/t1**3)))*alpha/(A0/A1)**(t**2/t1**2))**n*((-(1/(((-(h0*(A0/A1)**(t**2/t1**2) - A0*sin(pi*(-z/L + 1)**q0*(-z/L + 1)**(q1*t**3/t1**3)/(-z/L + 1)**(q0*t**3/t1**3)))*alpha/(A0/A1)**(t**2/t1**2))**n + 1)**m))**(1/m) + 1)**m - 1)*((L - z)/L)**(q0 - q0*t**3/t1**3 + q1*t**3/t1**3)*cos(pi*(-z/L + 1)**q0*(-z/L + 1)**(q1*t**3/t1**3)/(-z/L + 1)**(q0*t**3/t1**3))/((h0*(A0/A1)**(t**2/t1**2) - A0*sin(pi*(-z/L + 1)**q0*(-z/L + 1)**(q1*t**3/t1**3)/(-z/L + 1)**(q0*t**3/t1**3)))*(L - z)*t1**3*sqrt(((-(h0*(A0/A1)**(t**2/t1**2) - A0*sin(pi*(-z/L + 1)**q0*(-z/L + 1)**(q1*t**3/t1**3)/(-z/L + 1)**(q0*t**3/t1**3)))*alpha/(A0/A1)**(t**2/t1**2))**n + 1)**m)*((-(h0*(A0/A1)**(t**2/t1**2) - A0*sin(pi*(-z/L + 1)**q0*(-z/L + 1)**(q1*t**3/t1**3)/(-z/L + 1)**(q0*t**3/t1**3)))*alpha/(A0/A1)**(t**2/t1**2))**n + 1)*((1/(((-(h0*(A0/A1)**(t**2/t1**2) - A0*sin(pi*(-z/L + 1)**q0*(-z/L + 1)**(q1*t**3/t1**3)/(-z/L + 1)**(q0*t**3/t1**3)))*alpha/(A0/A1)**(t**2/t1**2))**n + 1)**m))**(1/m) - 1))

  end function dKdz_hle0_exact

  function dhdz_exact(z,t,L,A0,A1,q0,q1,t1) result(dhdz)
   ! *** First derivative of pressure head with respect to z for the manufactured Richards problem (from SageMath) ***
   ! input
   real(r8b),intent(in) :: z,t ! vertical position [m] (positive up) and time [s]
   real(r8b),intent(in) :: L   ! domain length
   real(r8b),intent(in) :: A0,A1,q0,q1,t1 ! manufactured solution parameters

   ! output
   real(r8b) :: dhdz
 
   dhdz=-pi*(q0*t**3 - q1*t**3 - q0*t1**3)*A0*((L - z)/L)**(q0 - q0*t**3/t1**3 + q1*t**3/t1**3)*cos(pi*(-z/L + 1)**q0*(-z/L + 1)**(q1*t**3/t1**3)/(-z/L + 1)**(q0*t**3/t1**3))/((L - z)*t1**3*(A0/A1)**(t**2/t1**2))
  end function dhdz_exact
 
  function d2hdz2_exact(z,t,L,A0,A1,q0,q1,t1) result(d2hdz2)
   ! *** Second derivative of pressure head with respect to z for the manufactured Richards problem (from SageMath) ***
   ! input
   real(r8b),intent(in) :: z,t ! vertical position [m] (positive up) and time [s]
   real(r8b),intent(in) :: L   ! domain length
   real(r8b),intent(in) :: A0,A1,q0,q1,t1 ! manufactured solution parameters

   ! output
   real(r8b) :: d2hdz2
   
   ! local variables
   real(r8b) :: term1,term2,term3,term4,term5,term6
   real(r8b) :: denominator

   ! OG from SageMath
   !d2hdz2=-(((pi*A0*q0**2 - 2*pi*A0*q0*q1 + pi*A0*q1**2)*t**6 - (2*pi*A0*q0**2 - pi*A0*q0 - (2*pi*A0*q0 - pi*A0)*q1)*t**3*t1**3 + (pi*A0*q0**2 - pi*A0*q0)*t1**6)*((L - z)/L)**(q0 + q0*t**3/t1**3 + q1*t**3/t1**3)*cos(pi*(-z/L + 1)**q0*(-z/L + 1)**(q1*t**3/t1**3)/(-z/L + 1)**(q0*t**3/t1**3)) - (pi**2*A0*q0**2*t1**6 + (pi**2*A0*q0**2 - 2*pi**2*A0*q0*q1 + pi**2*A0*q1**2)*t**6 - 2*(pi**2*A0*q0**2 - pi**2*A0*q0*q1)*t**3*t1**3)*((L - z)/L)**(2*q0 + 2*q1*t**3/t1**3)*sin(pi*(-z/L + 1)**q0*(-z/L + 1)**(q1*t**3/t1**3)/(-z/L + 1)**(q0*t**3/t1**3)))/((L**2*t1**6*(A0/A1)**(t**2/t1**2) - 2*L*t1**6*z*(A0/A1)**(t**2/t1**2) + t1**6*z**2*(A0/A1)**(t**2/t1**2))*((L - z)/L)**(2*q0*t**3/t1**3))

   ! manual optimization for repeated terms
   !term1=-z/L + 1
   !term2=(L - z)/L 
   !d2hdz2=-(((pi*A0*q0**2 - 2*pi*A0*q0*q1 + pi*A0*q1**2)*t**6 - (2*pi*A0*q0**2 - pi*A0*q0 - (2*pi*A0*q0 - pi*A0)*q1)*t**3*t1**3 + (pi*A0*q0**2 - pi*A0*q0)*t1**6)*(term2)**(q0 + q0*t**3/t1**3 + q1*t**3/t1**3)*cos(pi*(term1)**q0*(term1)**(q1*t**3/t1**3)/(term1)**(q0*t**3/t1**3)) - (pi**2*A0*q0**2*t1**6 + (pi**2*A0*q0**2 - 2*pi**2*A0*q0*q1 + pi**2*A0*q1**2)*t**6 - 2*(pi**2*A0*q0**2 - pi**2*A0*q0*q1)*t**3*t1**3)*(term2)**(2*q0 + 2*q1*t**3/t1**3)*sin(pi*(term1)**q0*(term1)**(q1*t**3/t1**3)/(term1)**(q0*t**3/t1**3)))/((L**2*t1**6*(A0/A1)**(t**2/t1**2) - 2*L*t1**6*z*(A0/A1)**(t**2/t1**2) + t1**6*z**2*(A0/A1)**(t**2/t1**2))*(term2)**(2*q0*t**3/t1**3))

   ! further manual optimization for repeated terms
   term1=-z/L + 1  !OG
   term2=(L - z)/L !OG
   !term3=(A0/A1)**(t**2/t1**2) !OG
   term3=(A0/A1)**((t/t1)**2) !simplified
   term4=t1**6 !OG
   !term5=t**3/t1**3 !OG
   term5=(t/t1)**3 !simplified
   term6=term3*term4 ! OG

   ! expression in denominator
   ! minimal simplification
   !denominator=((L**2*t1**6*(A0/A1)**(t**2/t1**2) - 2*L*t1**6*z*(A0/A1)**(t**2/t1**2) + t1**6*z**2*(A0/A1)**(t**2/t1**2))*(term2)**(2*q0*t**3/t1**3)) 
   ! quadratic polynomial in z in denominator
   !denominator=((L**2 - 2*L*z + z**2)*term6*(term2)**(2*q0*term5)) 
   ! factored form for quadratic polynomial in z in denominator
   ! note: using this form allows consistent results between compilers without using quad precision
   denominator=((z-L)**2*term6*(term2)**(2*q0*term5))
   ! derivative expression 
   d2hdz2=-(((pi*A0*q0**2 - 2*pi*A0*q0*q1 + pi*A0*q1**2)*t**6 - (2*pi*A0*q0**2 - pi*A0*q0 - (2*pi*A0*q0 - pi*A0)*q1)*t**3*t1**3 + (pi*A0*q0**2 - pi*A0*q0)*term4)*(term2)**(q0 + q0*term5 + q1*term5)*cos(pi*(term1)**q0*(term1)**(q1*term5)/(term1)**(q0*term5)) - (pi**2*A0*q0**2*term4 + (pi**2*A0*q0**2 - 2*pi**2*A0*q0*q1 + pi**2*A0*q1**2)*t**6 - 2*(pi**2*A0*q0**2 - pi**2*A0*q0*q1)*t**3*t1**3)*(term2)**(2*q0 + 2*q1*term5)*sin(pi*(term1)**q0*(term1)**(q1*term5)/(term1)**(q0*term5)))/denominator

  end function d2hdz2_exact


  function S_Richards_exact_SageMath(z,t,L,alpha,n,K_s,theta_r,theta_s,A0,A1,q0,q1,t1,h0) result(S)
   ! *** Source term for manufactured Richards solution using functions generated by SageMath ***

   ! input
   real(r8b),intent(in) :: z,t ! vertical position [m] (positive up) and time [s]
   real(r8b),intent(in) :: L   ! domain length
   real(r8b),intent(in) :: alpha,n,K_s,theta_r,theta_s ! van Genuchten parameters
   real(r8b),intent(in) :: A0,A1,q0,q1,t1,h0 ! manufactured solution parameters

   ! output 
   real(r8b) :: S

   ! local variables
   real(r8b) :: C,dhdt,dKdz,dhdz,K,d2hdz2

   C=C_hle0_exact(z,t,L,alpha,n,theta_r,theta_s,A0,A1,q0,q1,t1,h0)

   K=K_hle0_exact(z,t,L,alpha,n,K_s,A0,A1,q0,q1,t1,h0)
   dKdz=dKdz_hle0_exact(z,t,L,alpha,n,K_s,A0,A1,q0,q1,t1,h0)

   dhdt=dhdt_exact(z,t,L,A0,A1,q0,q1,t1)
   dhdz=dhdz_exact(z,t,L,A0,A1,q0,q1,t1)
   d2hdz2=d2hdz2_exact(z,t,L,A0,A1,q0,q1,t1)

   S=dKdz*(dhdz+1._r8b)+K*d2hdz2-C*dhdt

  end function S_Richards_exact_SageMath

  ! ************************************* Maple Functions ************************************* !
  function S_Richards_exact(z,t,L,alpha,n,Ks,theta_r,theta_s,A0,A1,q0,q1,t1) result(S)
   ! *** Source term for manufactured Richards solution generated by Maple ***

   ! input
   real(r8b),intent(in) :: z,t ! vertical position [m] (positive up) and time [s]
   real(r8b),intent(in) :: L   ! domain length
   real(r8b),intent(in) :: alpha,n,Ks,theta_r,theta_s ! van Genuchten parameters
   real(r8b),intent(in) :: A0,A1,q0,q1,t1 ! manufactured solution parameters

   ! output 
   real(r8b) :: S

   ! local variables generated by Maple
   real(r8b) :: t4
   real(r8b) :: t5
   real(r8b) :: t7
   real(r8b) :: t8
   real(r8b) :: t10 
   real(r8b) :: t11 
   real(r8b) :: t13 
   real(r8b) :: t14 
   real(r8b) :: t18 
   real(r8b) :: t20 
   real(r8b) :: t21 
   real(r8b) :: t22 
   real(r8b) :: t23 
   real(r8b) :: t26
   real(r8b) :: t27 
   real(r8b) :: t29 
   real(r8b) :: t30 
   real(r8b) :: t31 
   real(r8b) :: t32 
   real(r8b) :: t35 
   real(r8b) :: t37 
   real(r8b) :: t39 
   real(r8b) :: t40 
   real(r8b) :: t41 
   real(r8b) :: t42 
   real(r8b) :: t44 
   real(r8b) :: t45 
   real(r8b) :: t46 
   real(r8b) :: t47 
   real(r8b) :: t53 
   real(r8b) :: t55 
   real(r8b) :: t82

   real(r8b) :: cg,cg1,cg3,cg5,cg7,cg9,cg11,cg13 

   ! interface between user variables and Maple variables
   cg=A0; cg1=A1
   cg3=Ks
   cg5=q0; cg7=q1
   cg9=t1
   cg11=theta_r; cg13=theta_s   

   t4 = log(cg / cg1)
   t5 = (t ** 2)
   t7 = (cg9 ** 2)
   t8 = 0.1D1 / t7
   t10 = exp(-t8 * t5 * t4)
   t11 = (L - z)
   t13 = 0.1D1 / L * t11
   t14 = (cg7 - cg5)
   t18 = 0.1D1 / cg9 / t7
   t20 = t18 * t * t5 * t14 + cg5
   t21 = t13 ** t20
   t22 = t21 * 0.31415926535897932D1
   t23 = sin(t22)
   t26 = (t23 * t10 * alpha * cg) ** n
   t27 = 1 + t26
   t29 = (0.1D1 - 0.1D1 / n)
   t30 = t27 ** t29
   t31 = 0.1D1 / t30
   t32 = sqrt(t31)
   t35 = t31 ** (0.1D1 / t29)
   t37 = (1 - t35) ** t29
   t39 = (1 - t37) ** 2
   t40 = t10 * cg
   t41 = 0.31415926535897932D1 * t40
   t42 = t20 ** 2
   t44 = t11 ** 2
   t45 = 0.1D1 / t44
   t46 = cos(t22)
   t47 = t46 * t45
   t53 = (0.31415926535897932D1 ** 2)
   t55 = t21 ** 2
   t82 = log(t13)
   S = (t23 * t45 * t42 * t55 * t53 * t40 + t47 * t20 * t21 * t41 - t& 
     &47 * t42 * t21 * t41) * t39 * t32 * cg3 - (-3 * t46 * t82 * t18 *& 
     &t5 * t14 * t22 * t40 + 2 * t23 * t10 * t8 * t * t4 * cg) / t27 / t& 
     &23 / t10 / cg * n * t26 * t29 * t31 * (cg13 - cg11) 
  end function S_Richards_exact

end module Richards_exact
