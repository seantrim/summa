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

module ssdNrgFlux_module

! data types
USE nrtype

! data types
USE data_types,only:var_d           ! x%var(:)       (dp)
USE data_types,only:var_dlength     ! x%var(:)%dat   (dp)
USE data_types,only:var_ilength     ! x%var(:)%dat   (i4b)
USE data_types,only:data_bin        ! x%bin(:)%{lgt(:),i4b(:),rkind(:),rmatrix(:,:),string}, x%err [i4b], x%msg [character]

! physical constants
USE multiconst,only:&
                    sb,          &  ! Stefan Boltzman constant      (W m-2 K-4)
                    Em_Sno,      &  ! emissivity of snow            (-)
                    Cp_air,      &  ! specific heat of air          (J kg-1 K-1)
                    Cp_water,    &  ! specifric heat of water       (J kg-1 K-1)
                    LH_fus,      &  ! latent heat of fusion         (J kg-1)
                    LH_vap,      &  ! latent heat of vaporization   (J kg-1)
                    LH_sub,      &  ! latent heat of sublimation    (J kg-1)
                    gravity,     &  ! gravitational acceleteration  (m s-2)
                    Tfreeze,     &  ! freezing point of pure water  (K)
                    iden_air,    &  ! intrinsic density of air      (kg m-3)
                    iden_ice,    &  ! intrinsic density of ice      (kg m-3)
                    iden_water      ! intrinsic density of water    (kg m-3)

! missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! named variables for snow and soil
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! named variables
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements

! model decisions
USE globalData,only:model_decisions ! model decision structure
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure

! provide access to look-up values for model decisions
USE mDecisions_module,only:      &
 ! look-up values for method used to compute derivative
 numerical,                      &  ! numerical solution
 analytical,                     &  ! analytical solution
 ! look-up values for choice of boundary conditions for thermodynamics
 prescribedTemp,                 &  ! prescribed temperature
 energyFlux,                     &  ! energy flux
 zeroFlux,                       &  ! zero flux
 ! look-up values for choice of boundary conditions for soil hydrology
 prescribedHead                     ! prescribed head

! -------------------------------------------------------------------------------------------------
implicit none
private
public :: ssdNrgFlux
! global parameters
real(rkind),parameter            :: valueMissing=-9999._rkind   ! missing value parameter
contains

 ! ************************************************************************************************
 ! public subroutine ssdNrgFlux: compute energy fluxes and derivatives at layer interfaces
 ! ************************************************************************************************
 subroutine ssdNrgFlux(&
                       ! input: model control, fluxes, derivatives and trial model state variables
                       in_data,                            & ! intent(in):    model control, fluxes, derivatives, and trial state variables
                       ! input-output: data structures
                       mpar_data,                          & ! intent(in):    model parameters
                       indx_data,                          & ! intent(in):    model indices
                       prog_data,                          & ! intent(in):    model prognostic variables for a local HRU
                       diag_data,                          & ! intent(in):    model diagnostic variables for a local HRU
                       flux_data,                          & ! intent(inout): model fluxes for a local HRU
                       ! output: fluxes/derivatives at layer interfaces, and error control
                       out_data)                             ! intent(out):   fluxes/derivatives at layer interfaces, and error control
 implicit none
 ! input: model control, fluxes, derivatives, and trial state variables
 type(data_bin),intent(in)          :: in_data                    ! input data structure (model options, fluxes, derivatives, and trial state variables)
 ! input-output: data structures
 type(var_dlength),intent(in)       :: mpar_data                  ! model parameters
 type(var_ilength),intent(in)       :: indx_data                  ! state vector geometry
 type(var_dlength),intent(in)       :: prog_data                  ! prognostic variables for a local HRU
 type(var_dlength),intent(in)       :: diag_data                  ! diagnostic variables for a local HRU
 type(var_dlength),intent(inout)    :: flux_data                  ! model fluxes for a local HRU
 ! output: fluxes and derivatives at layer interfaces, and error control
 type(data_bin),intent(out)         :: out_data                   ! output data structure (fluxes/derivatives at layer interfaces, and error control)
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                       :: nLayers                    ! total number of layers
 integer(i4b)                       :: iLayer                     ! index of model layers
 integer(i4b)                       :: ixLayerDesired(1)          ! layer desired (scalar solution)
 integer(i4b)                       :: ixTop                      ! top layer in subroutine call
 integer(i4b)                       :: ixBot                      ! bottom layer in subroutine call
 real(rkind)                        :: qFlux                      ! liquid flux at layer interfaces (m s-1)
 real(rkind)                        :: dz                         ! height difference (m)
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! allocate data structure for storing intent(out) arguments
 nLayers=indx_data%var(iLookINDEX%nLayers)%dat(1); allocate(out_data%bin(1:3))
 allocate(out_data%bin(1)%rkind(0:nLayers),out_data%bin(2)%rkind(0:nLayers),out_data%bin(3)%rkind(0:nLayers)) 
 ! make association of local variables with information in the data structures
 associate(&
  ! input: model control, fluxes, derivatives, and trial state variables
  scalarSolution             => in_data%bin(1)%lgt(1),                          & ! intent(in): flag to denote if implementing the scalar solution
  groundNetFlux              => in_data%bin(1)%rkind(1),                        & ! intent(in): net energy flux for the ground surface (W m-2)
  dGroundNetFlux_dGroundTemp => in_data%bin(1)%rkind(2),                        & ! intent(in): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
  iLayerLiqFluxSnow          => in_data%bin(2)%rkind,                           & ! intent(in): liquid water flux at the interface of each snow layer (m s-1) 
  iLayerLiqFluxSoil          => in_data%bin(3)%rkind,                           & ! intent(in): liquid water flux at the interface of each soil layer (m s-1)
  mLayerTempTrial            => in_data%bin(4)%rkind,                           & ! intent(in): trial temperature of each snow/soil layer at the current iteration (K)
  ! input: model decisions for calculation of flux derivatives and boundary conditions
  ix_fDerivMeth        => model_decisions(iLookDECISIONS%fDerivMeth)%iDecision, & ! intent(in): method used to calculate flux derivatives
  ix_bcLowrTdyn        => model_decisions(iLookDECISIONS%bcLowrTdyn)%iDecision, & ! intent(in): method used to calculate the lower boundary condition for thermodynamics
  ! input: model coordinates
  nSnow                => indx_data%var(iLookINDEX%nSnow)%dat(1),               & ! intent(in): number of snow layers
  layerType            => indx_data%var(iLookINDEX%layerType)%dat,              & ! intent(in): layer type (iname_soil or iname_snow)
  ixLayerState         => indx_data%var(iLookINDEX%ixLayerState)%dat,           & ! intent(in): list of indices for all model layers
  ixSnowSoilNrg        => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat,          & ! intent(in): index in the state subset for energy state variables in the snow+soil domain
  ! input: thermal properties
  mLayerDepth          => prog_data%var(iLookPROG%mLayerDepth)%dat,             & ! intent(in): depth of each layer (m)
  mLayerHeight         => prog_data%var(iLookPROG%mLayerHeight)%dat,            & ! intent(in): height at the mid-point of each layer (m)
  iLayerThermalC       => diag_data%var(iLookDIAG%iLayerThermalC)%dat,          & ! intent(in): thermal conductivity at the interface of each layer (W m-1 K-1)
  lowerBoundTemp       => mpar_data%var(iLookPARAM%lowerBoundTemp)%dat(1),      & ! intent(in): temperature of the lower boundary (K)
  ! output: diagnostic fluxes
  iLayerConductiveFlux => flux_data%var(iLookFLUX%iLayerConductiveFlux)%dat,    & ! intent(out): conductive energy flux at layer interfaces at end of time step (W m-2)
  iLayerAdvectiveFlux  => flux_data%var(iLookFLUX%iLayerAdvectiveFlux)%dat,     & ! intent(out): advective energy flux at layer interfaces at end of time step (W m-2)
  ! output: fluxes and derivatives at all layer interfaces 
  iLayerNrgFlux        => out_data%bin(1)%rkind,                                & ! intent(out): energy flux at the layer interfaces (W m-2)
  dFlux_dTempAbove     => out_data%bin(2)%rkind,                                & ! intent(out): flux derivatives w.r.t. temperature in the layer above (J m-2 s-1 K-1)
  dFlux_dTempBelow     => out_data%bin(3)%rkind,                                & ! intent(out): flux derivatives w.r.t. temperature in the layer below (J m-2 s-1 K-1)
  ! output: error control
  err                  => out_data%err,                                         & ! intent(out): error code
  message              => out_data%msg                                          & ! intent(out): error message
 )  ! end associate statement of local variables with information in the data structures
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 err=0; message='ssdNrgFlux/' ! initialize error control

 ! set conductive and advective fluxes to missing in the upper boundary
 ! NOTE: advective flux at the upper boundary is included in the ground heat flux
 iLayerConductiveFlux(0) = valueMissing
 iLayerAdvectiveFlux(0)  = valueMissing

 ! get the indices for the snow+soil layers
 if (scalarSolution) then
  ixLayerDesired = pack(ixLayerState, ixSnowSoilNrg/=integerMissing)
  ixTop = ixLayerDesired(1)
  ixBot = ixLayerDesired(1)
 else
  ixTop = 1
  ixBot = nLayers
 end if

 ! -------------------------------------------------------------------------------------------------------------------------
 ! ***** compute the conductive fluxes at layer interfaces *****
 ! -------------------------------------------------------------------------------------------------------------------------
 do iLayer=ixTop,ixBot ! loop through model layers
  if (iLayer==nLayers) then ! compute fluxes at the lower boundary -- positive downwards
   ! flux depends on the type of lower boundary condition
   select case(ix_bcLowrTdyn) ! identify the lower boundary condition for thermodynamics
    case(prescribedTemp); iLayerConductiveFlux(nLayers) = -iLayerThermalC(iLayer)*(lowerBoundTemp - mLayerTempTrial(iLayer))/(mLayerDepth(iLayer)*0.5_rkind)
    case(zeroFlux);       iLayerConductiveFlux(nLayers) = 0._rkind
    case default;         err=20; message=trim(message)//'unable to identify lower boundary condition for thermodynamics'; return
   end select  ! end identifying the lower boundary condition for thermodynamics
  else ! compute fluxes within the domain -- positive downwards
    iLayerConductiveFlux(iLayer)  = -iLayerThermalC(iLayer)*(mLayerTempTrial(iLayer+1) - mLayerTempTrial(iLayer)) / &
                                    (mLayerHeight(iLayer+1) - mLayerHeight(iLayer))
  end if ! end if at lower boundary
 end do

 ! -------------------------------------------------------------------------------------------------------------------------
 ! ***** compute the advective fluxes at layer interfaces *****
 ! -------------------------------------------------------------------------------------------------------------------------
 do iLayer=ixTop,ixBot
  ! get the liquid flux at layer interfaces
  select case(layerType(iLayer))
   case(iname_snow); qFlux = iLayerLiqFluxSnow(iLayer)
   case(iname_soil); qFlux = iLayerLiqFluxSoil(iLayer-nSnow)
   case default; err=20; message=trim(message)//'unable to identify layer type'; return
  end select
  ! compute fluxes at the lower boundary -- positive downwards
  if(iLayer==nLayers)then
   iLayerAdvectiveFlux(iLayer) = -Cp_water*iden_water*qFlux*(lowerBoundTemp - mLayerTempTrial(iLayer))
  ! compute fluxes within the domain -- positive downwards
  else
   iLayerAdvectiveFlux(iLayer) = -Cp_water*iden_water*qFlux*(mLayerTempTrial(iLayer+1) - mLayerTempTrial(iLayer))
  end if
 end do  ! looping through layers

 ! -------------------------------------------------------------------------------------------------------------------------
 ! ***** compute the total fluxes at layer interfaces *****
 ! -------------------------------------------------------------------------------------------------------------------------
 ! NOTE: ignore advective fluxes for now
 iLayerNrgFlux(0)           = groundNetFlux
 iLayerNrgFlux(ixTop:ixBot) = iLayerConductiveFlux(ixTop:ixBot)

 ! -------------------------------------------------------------------------------------------------------------------------
 ! ***** compute the derivative in fluxes at layer interfaces w.r.t temperature in the layer above and the layer below *****
 ! -------------------------------------------------------------------------------------------------------------------------

 ! initialize un-used elements
 dFlux_dTempBelow(nLayers) = -huge(lowerBoundTemp)  ! don't expect this to be used, so deliberately set to a ridiculous value to cause problems

 ! ***** the upper boundary
 dFlux_dTempBelow(0) = dGroundNetFlux_dGroundTemp

 ! ***** validate method selected for derivative computations (only analytical allowed -- numerical no longer used)
 if (ix_fDerivMeth/=analytical) then; err=20; message=trim(message)//'only analytical derivatives are allowed -- numerical derivatives no longer used'; return; end if

 ! loop through INTERFACES...
 do iLayer=ixTop,ixBot
  ! ***** the lower boundary
  if (iLayer==nLayers) then  ! if at lower boundary
   select case(ix_bcLowrTdyn) ! identify the lower boundary condition
    case(prescribedTemp) ! * prescribed temperature at the lower boundary
     dz = mLayerDepth(iLayer)*0.5_rkind
     dFlux_dTempAbove(iLayer) = iLayerThermalC(iLayer)/dz ! ** analytical derivatives
    case(zeroFlux) ! * zero flux at the lower boundary
     dFlux_dTempAbove(iLayer) = 0._rkind
    case default; err=20; message=trim(message)//'unable to identify lower boundary condition for thermodynamics'; return
   end select  ! end identifying the lower boundary condition for thermodynamics
  else ! ***** internal layers
   dz = (mLayerHeight(iLayer+1) - mLayerHeight(iLayer))
   dFlux_dTempAbove(iLayer) =  iLayerThermalC(iLayer)/dz ! ** analytical derivatives
   dFlux_dTempBelow(iLayer) = -iLayerThermalC(iLayer)/dz
  end if  ! end if for type of layer (upper, internal, or lower)
 end do  ! end looping through layers

 end associate ! end association of local variables with information in the data structures

 end subroutine ssdNrgFlux

end module ssdNrgFlux_module

