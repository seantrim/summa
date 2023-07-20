! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module bigAquifer_module
! -----------------------------------------------------------------------------------------------------------
! numerical recipes data types
USE nrtype

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! physical constants
USE multiconst,only:&
                    LH_vap,  &      ! latent heat of vaporization   (J kg-1)
                    iden_water      ! intrinsic density of water    (kg m-3)
! -----------------------------------------------------------------------------------------------------------
implicit none
private
public :: bigAquifer
contains

 ! ************************************************************************************************************************
 ! public subroutine bigAquifer: compute aquifer transpiration, recharge, baseflow, and derivative in the net aquifer flux  
 ! ************************************************************************************************************************
 subroutine bigAquifer(&
                       ! input: state variables and fluxes
                       in_data,                      & ! intent(in):  state variables and fluxes
                       ! input: diagnostic variables and parameters
                       mpar_data,                    & ! intent(in):  model parameter structure
                       diag_data,                    & ! intent(in):  diagnostic variable structure
                       ! output: fluxes and error control
                       out_data)                       ! intent(out): fluxes and error control
 ! named variables
 USE var_lookup,only:iLookDIAG                         ! named variables for structure elements
 USE var_lookup,only:iLookPARAM                        ! named variables for structure elements
 ! data types
 USE data_types,only:var_dlength                       ! x%var(:)%dat   (dp)
 USE data_types,only:data_bin                          ! x%bin(:)%{lgt(:),i4b(:),rkind(:),rmatrix(:,:),string}, x%err [i4b], x%msg [character]
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 implicit none
 ! input: state variables, fluxes, and parameters
 type(data_bin),intent(in)           :: in_data                      ! state variables, fluxes, and parameters
 ! input: diagnostic variables and parameters
 type(var_dlength),intent(in)        :: mpar_data                    ! model parameters
 type(var_dlength),intent(in)        :: diag_data                    ! diagnostic variables for a local HRU
 ! output: fluxes and error control
 type(data_bin),intent(out)          :: out_data                     ! fluxes and error control
 ! -----------------------------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 real(rkind)                         :: aquiferTranspireFrac         ! fraction of total transpiration that comes from the aquifer (-)
 real(rkind)                         :: xTemp                        ! temporary variable (-)
 ! -------------------------------------------------------------------------------------------------------------------------------------------------

 ! make association between local variables and the information in the data structures
 associate(&
 ! input: state variables, fluxes, and parameters
 scalarAquiferStorageTrial   => in_data%bin(1)%rkind(1),                           & ! intent(in): [dp] trial value of aquifer storage (m)
 scalarCanopyTranspiration   => in_data%bin(1)%rkind(2),                           & ! intent(in): [dp] canopy transpiration (kg m-2 s-1)
 scalarSoilDrainage          => in_data%bin(1)%rkind(3),                           & ! intent(in): [dp] soil drainage (m s-1)
 ! model diagnostic variables: contribution of the aquifer to transpiration
 scalarTranspireLim     => diag_data%var(iLookDIAG%scalarTranspireLim)%dat(1),     & ! intent(in): [dp] weighted average of the transpiration limiting factor (-)
 scalarAquiferRootFrac  => diag_data%var(iLookDIAG%scalarAquiferRootFrac)%dat(1),  & ! intent(in): [dp] fraction of roots below the lowest soil layer (-)
 scalarTranspireLimAqfr => diag_data%var(iLookDIAG%scalarTranspireLimAqfr)%dat(1), & ! intent(in): [dp] transpiration limiting factor for the aquifer (-)
 ! model parameters: baseflow flux
 aquiferBaseflowRate    => mpar_data%var(iLookPARAM%aquiferBaseflowRate)%dat(1),   & ! intent(in): [dp] tbaseflow rate when aquiferStorage = aquiferScaleFactor (m s-1)
 aquiferScaleFactor     => mpar_data%var(iLookPARAM%aquiferScaleFactor)%dat(1),    & ! intent(in): [dp] scaling factor for aquifer storage in the big bucket (m)
 aquiferBaseflowExp     => mpar_data%var(iLookPARAM%aquiferBaseflowExp)%dat(1),    & ! intent(in): [dp] baseflow exponent (-)
 ! output: fluxes
 scalarAquiferTranspire => out_data%bin(1)%rkind(1),                               & ! intent(out): [dp] transpiration loss from the aquifer (m s-1)
 scalarAquiferRecharge  => out_data%bin(1)%rkind(2),                               & ! intent(out): [dp] recharge to the aquifer (m s-1)
 scalarAquiferBaseflow  => out_data%bin(1)%rkind(3),                               & ! intent(out): [dp] total baseflow from the aquifer (m s-1)
 dBaseflow_dAquifer     => out_data%bin(1)%rkind(4),                               & ! intent(out): [dp] change in baseflow flux w.r.t. aquifer storage (s-1)
 ! output: error control
 err                    => out_data%err,                                           & ! intent(out): [i4b] error code
 message                => out_data%msg                                            & ! intent(out): [character] error message
 )  ! end associating local variables with the information in the data structures
 err=0; message='bigAquifer/' ! initialize error control

 ! compute aquifer transpiration (m s-1)
 aquiferTranspireFrac   = scalarAquiferRootFrac*scalarTranspireLimAqfr/scalarTranspireLim   ! fraction of total transpiration that comes from the aquifer (-)
 scalarAquiferTranspire = aquiferTranspireFrac*scalarCanopyTranspiration/iden_water         ! aquifer transpiration (kg m-2 s-1 --> m s-1)

 ! compute aquifer recharge (transfer variables -- included for generality for basin-wide aquifer)
 scalarAquiferRecharge = scalarSoilDrainage ! m s-1

 ! compute the aquifer baseflow (m s-1)
 xTemp                 = scalarAquiferStorageTrial/aquiferScaleFactor
 scalarAquiferBaseflow = aquiferBaseflowRate*(xTemp**aquiferBaseflowExp)

 ! compute the derivative in the net aquifer flux
 dBaseflow_dAquifer    = -(aquiferBaseflowExp*aquiferBaseflowRate*(xTemp**(aquiferBaseflowExp - 1._rkind)))/aquiferScaleFactor

 end associate ! end association to data in structures

 end subroutine bigAquifer

end module bigAquifer_module
