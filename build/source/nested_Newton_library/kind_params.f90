module kind_params
 use, intrinsic :: iso_fortran_env, only: int32,int64,real32,real64,real128 !!integer and real kind parameters
 implicit none
 integer, parameter :: i4b=int32
 integer, parameter :: i8b=int64
 integer, parameter :: r4b=real32
 integer, parameter :: r8b=real64
 integer, parameter :: r16b=real128
end module kind_params

