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

module computResid_module

! data types
USE nr_type

! derived types to define the data structures
USE data_types,only:&
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength     ! data vector with variable length dimension (rkind)

! named variables
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookDERIV      ! named variables for structure elements

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! define access to state variables to print
USE globalData,only: iJac1          ! first layer of the Jacobian to print
USE globalData,only: iJac2          ! last layer of the Jacobian to print

! domain types
USE globalData,only:iname_veg       ! named variables for vegetation
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers

! constants
USE multiconst,only:&
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)
! privacy
implicit none
private
public::computResid
contains

! **********************************************************************************************************
! public subroutine computResid: compute the residual vector
! **********************************************************************************************************
subroutine computResid(&
                      ! input: model control
                      dt,                        & ! intent(in):  length of the time step (seconds)
                      nSnow,                     & ! intent(in):  number of snow layers
                      nSoil,                     & ! intent(in):  number of soil layers
                      nLayers,                   & ! intent(in):  total number of layers
                      mixdformNrg,               & ! intent(in):  flag to use enthalpy formulation
                      mass_flag,                 & ! intent(in):  flag to compute mass terms
                      energy_flag,               & ! intent(in):  flag to compute energy terms
                      f_flag,f1_mass,f1_energy,f2_mass,f2_energy, & ! flags to compute f, f1, and f2 for nested Newton
                      ! input: flux vectors
                      sMul,                      & ! intent(in):  state vector multiplier (used in the residual calculations)
                      fVec,                      & ! intent(in):  flux vector
                      ! input: state variables (already disaggregated into scalars and vectors)
                      scalarCanairTempTrial,     & ! intent(in):  trial value for the temperature of the canopy air space (K)
                      scalarCanopyTempTrial,     & ! intent(in):  trial value for the temperature of the vegetation canopy (K)
                      scalarCanopyWatTrial,      & ! intent(in):  trial value for the water on the vegetation canopy (kg m-2)
                      mLayerTempTrial,           & ! intent(in):  trial value for the temperature of each snow and soil layer (K)
                      scalarAquiferStorageTrial, & ! intent(in):  trial value of storage of water in the aquifer (m)
                      ! input: diagnostic variables defining the liquid water and ice content (function of state variables)
                      scalarCanopyIceTrial,      & ! intent(in):  trial value for the ice on the vegetation canopy (kg m-2)
                      scalarCanopyLiqTrial,      & ! intent(in):  trial value for the liq on the vegetation canopy (kg m-2)
                      mLayerVolFracIceTrial,     & ! intent(in):  trial value for the volumetric ice in each snow and soil layer (-)
                      mLayerVolFracWatTrial,     & ! intent(in):  trial value for the volumetric water in each snow and soil layer (-)
                      mLayerVolFracLiqTrial,     & ! intent(in):  trial value for the volumetric liq in each snow and soil layer (-)
                      ! input: enthalpy terms
                      scalarCanopyCmTrial,       & ! intent(in):  Cm for vegetation canopy (J kg-1)
                      mLayerCmTrial,             & ! intent(in):  Cm for each snow+soil layer (J m-3)
                      scalarCanairEnthalpyTrial, & ! intent(in):  trial value for  enthalpy of the canopy air space (J m-3)
                      scalarCanopyEnthTempTrial, & ! intent(in):  trial value for temperature component of enthalpy of the vegetation canopy (J m-3)
                      mLayerEnthTempTrial,       & ! intent(in):  trial vector of temperature component of enthalpy of each snow+soil layer (J m-3)  
                      ! input: data structures
                      prog_data,                 & ! intent(in):  model prognostic variables for a local HRU
                      diag_data,                 & ! intent(in):  model diagnostic variables for a local HRU
                      flux_data,                 & ! intent(in):  model fluxes for a local HRU
                      indx_data,                 & ! intent(in):  index data
                      deriv_data,                & ! intent(in):  derivatives in model fluxes w.r.t. relevant state variables
                      ! output
                      f,f1,f2,                   & ! intent(out): f, f1, and f2 vectors for nested Newton objects
                      fRHS,                      & ! intent(out): right-hand-side function for ARKODE
                      rAdd,                      & ! intent(out): additional (sink) terms on the RHS of the state equation
                      rVec,                      & ! intent(out): residual vector
                      err,message)                 ! intent(out): error control
  ! --------------------------------------------------------------------------------------------------------------------------------
  implicit none
  ! input: model control
  real(rkind),intent(in)             :: dt                        ! length of the time step (seconds)
  integer(i4b),intent(in)            :: nSnow                     ! number of snow layers
  integer(i4b),intent(in)            :: nSoil                     ! number of soil layers
  integer(i4b),intent(in)            :: nLayers                   ! total number of layers in the snow+soil domain
  logical(lgt),intent(in)            :: mixdformNrg               ! flag to use enthalpy formulation
  logical(lgt),intent(in)            :: mass_flag,energy_flag     ! flags to compute mass and energy terms
  logical(lgt),intent(in)            :: f_flag,f1_mass,f1_energy,f2_mass,f2_energy ! flags to compute f, f1, and f2 for nested Newton
  ! input: flux vectors
  real(qp),intent(in)                :: sMul(:)   ! NOTE: qp      ! state vector multiplier (used in the residual calculations)
  real(rkind),intent(in)             :: fVec(:)                   ! flux vector
  ! input: state variables (already disaggregated into scalars and vectors)
  real(rkind),intent(in)             :: scalarCanairTempTrial     ! trial value for temperature of the canopy air space (K)
  real(rkind),intent(in)             :: scalarCanopyTempTrial     ! trial value for temperature of the vegetation canopy (K)
  real(rkind),intent(in)             :: scalarCanopyWatTrial      ! trial value for canopy total water content (kg m-2)
  real(rkind),intent(in)             :: mLayerTempTrial(:)        ! trial value for temperature of each snow/soil layer (K)
  real(rkind),intent(in)             :: scalarAquiferStorageTrial ! trial value of aquifer storage (m)
  ! input: diagnostic variables defining the liquid water and ice content (function of state variables)
  real(rkind),intent(in)             :: scalarCanopyIceTrial      ! trial value for mass of ice on the vegetation canopy (kg m-2)
  real(rkind),intent(in)             :: scalarCanopyLiqTrial      ! trial value for the liq on the vegetation canopy (kg m-2)
  real(rkind),intent(in)             :: mLayerVolFracIceTrial(:)  ! trial value for volumetric fraction of ice (-)
  real(rkind),intent(in)             :: mLayerVolFracWatTrial(:)  ! trial value for the volumetric water in each snow and soil layer (-)
  real(rkind),intent(in)             :: mLayerVolFracLiqTrial(:)  ! trial value for the volumetric water in each snow and soil layer (-)
  ! input: enthalpy terms
  real(rkind),intent(in)             :: scalarCanopyCmTrial       ! Cm for vegetation canopy (J kg-1)
  real(rkind),intent(in)             :: mLayerCmTrial(:)          ! Cm for each snow+soil layer (J m-3)
  real(rkind),intent(in)             :: scalarCanairEnthalpyTrial ! trial value for enthalpy of the canopy air space (J m-3)
  real(rkind),intent(in)             :: scalarCanopyEnthTempTrial ! trial value for temperature component of enthalpy of the vegetation canopy (J m-3)
  real(rkind),intent(in)             :: mLayerEnthTempTrial(:)    ! trial vector of temperature component of enthalpy of each snow+soil layer (J m-3)
  ! input: data structures
  type(var_dlength),intent(in)       :: prog_data                 ! prognostic variables for a local HRU
  type(var_dlength),intent(in)       :: diag_data                 ! diagnostic variables for a local HRU
  type(var_dlength),intent(in)       :: flux_data                 ! model fluxes for a local HRU
  type(var_ilength),intent(in)       :: indx_data                 ! indices defining model states and layers
  type(var_dlength),intent(in)       :: deriv_data                ! derivatives in model fluxes w.r.t. relevant state variables
  ! output
  real(rkind),intent(inout)          :: f(:),f1(:),f2(:)          ! f, f1, and f2 vectors for nested Newton objects
  real(rkind),intent(out)            :: fRHS(:)                   ! right-hand-side function for ARKODE
  real(rkind),intent(out)            :: rAdd(:)                   ! additional (sink) terms on the RHS of the state equation
  real(qp),intent(out)               :: rVec(:)   ! NOTE: qp      ! residual vector
  integer(i4b),intent(out)           :: err                       ! error code
  character(*),intent(out)           :: message                   ! error message
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  ! --------------------------------------------------------------------------------------------------------------------------------
  integer(i4b)                       :: iLayer                    ! index of layer within the snow+soil domain
  integer(i4b),parameter             :: ixVegVolume=1             ! index of the desired vegetation control volumne (currently only one veg layer)
  real(rkind)                        :: scalarCanopyHydTrial      ! trial value of canopy water content (kg m-2), either liquid water content or total water content
  real(rkind)                        :: scalarCanopyHyd           ! canopy water content (kg m-2), either liquid water content or total water content
  real(rkind),dimension(nLayers)     :: mLayerVolFracHydTrial     ! trial vector of volumetric water content (-), either liquid water content or total water content
  real(rkind),dimension(nLayers)     :: mLayerVolFracHyd          ! vector of volumetric water content (-), either liquid water content or total water content
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! link to the necessary variables for the residual computations
  associate(&
    ! model state variables (vegetation canopy)
    scalarCanairTemp        => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)       ,& ! intent(in): [dp]     temperature of the canopy air space (K)
    scalarCanopyTemp        => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,& ! intent(in): [dp]     temperature of the vegetation canopy (K)
    scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)        ,& ! intent(in): [dp]     mass of ice on the vegetation canopy (kg m-2)
    scalarCanopyLiq         => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)        ,& ! intent(in): [dp]     mass of liquid water on the vegetation canopy (kg m-2)
    scalarCanopyWat         => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)        ,& ! intent(in): [dp]     mass of total water on the vegetation canopy (kg m-2)
    ! model state variables (snow and soil domains)
    mLayerTemp              => prog_data%var(iLookPROG%mLayerTemp)%dat                ,& ! intent(in): [dp(:)]  temperature of each snow/soil layer (K)
    mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat          ,& ! intent(in): [dp(:)]  volumetric fraction of ice (-)
    mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,& ! intent(in): [dp(:)]  volumetric fraction of liquid water (-)
    mLayerVolFracWat        => prog_data%var(iLookPROG%mLayerVolFracWat)%dat          ,& ! intent(in): [dp(:)]  volumetric fraction of total water (-)
    ! enthalpy terms
    scalarCanairEnthalpy    => prog_data%var(iLookPROG%scalarCanairEnthalpy)%dat(1)   ,& ! intent(in): [dp]     enthalpy of the canopy air space (J m-3)
    scalarCanopyEnthTemp    => diag_data%var(iLookDIAG%scalarCanopyEnthTemp)%dat(1)   ,& ! intent(in): [dp]     temperature component of enthalpy of the vegetation canopy (J m-3)
    mLayerEnthTemp          => diag_data%var(iLookDIAG%mLayerEnthTemp)%dat            ,& ! intent(in): [dp(:)]  temperature component of enthalpy of the snow+soil layers (J m-3)
    ! model state variables (aquifer)
    scalarAquiferStorage    => prog_data%var(iLookPROG%scalarAquiferStorage)%dat(1)   ,& ! intent(in): [dp]     storage of water in the aquifer (m)
    ! canopy and layer depth
    canopyDepth             => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)      ,& ! intent(in): [dp]     canopy depth (m)
    mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,& ! intent(in): [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
    ! model fluxes (sink terms in the soil domain)
    mLayerTranspire         => flux_data%var(iLookFLUX%mLayerTranspire)%dat           ,& ! intent(in): [dp(:)]  transpiration loss from each soil layer (m s-1)
    mLayerBaseflow          => flux_data%var(iLookFLUX%mLayerBaseflow)%dat            ,& ! intent(in): [dp(:)]  baseflow from each soil layer (m s-1)
    mLayerCompress          => diag_data%var(iLookDIAG%mLayerCompress)%dat            ,& ! intent(in): [dp(:)]  change in storage associated with compression of the soil matrix (-)
    ! number of state variables of a specific type
    nSnowSoilNrg            => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in): [i4b]    number of energy state variables in the snow+soil domain
    nSnowSoilHyd            => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in): [i4b]    number of hydrology variables in the snow+soil domain
    nSoilOnlyHyd            => indx_data%var(iLookINDEX%nSoilOnlyHyd )%dat(1)         ,& ! intent(in): [i4b]    number of hydrology variables in the soil domain
    ! model indices
    ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in): [i4b]    index of canopy air space energy state variable
    ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in): [i4b]    index of canopy energy state variable
    ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,& ! intent(in): [i4b]    index of canopy hydrology state variable (mass)
    ixAqWat                 => indx_data%var(iLookINDEX%ixAqWat)%dat(1)               ,& ! intent(in): [i4b]    index of water storage in the aquifer
    ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in): [i4b(:)] indices for energy states in the snow+soil subdomain
    ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in): [i4b(:)] indices for hydrology states in the snow+soil subdomain
    ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat            ,& ! intent(in): [i4b(:)] indices for hydrology states in the soil subdomain
    ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat              ,& ! intent(in): [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
    ixHydCanopy             => indx_data%var(iLookINDEX%ixHydCanopy)%dat              ,& ! intent(in): [i4b(:)] index of the hydrology states in the canopy domain
    ixHydType               => indx_data%var(iLookINDEX%ixHydType)%dat                ,& ! intent(in): [i4b(:)] named variables defining the type of hydrology states in snow+soil domain
    layerType               => indx_data%var(iLookINDEX%layerType)%dat                 & ! intent(in): [i4b(:)] named variables defining the type of layer in snow+soil domain
    ) ! association to necessary variables for the residual computations
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message="computResid/"

    ! initialize rVec
    rVec(:) = 0._qp
    if (f_flag) f(:) = 0._rkind
    if (f1_mass.or.f1_energy) f1(:) = 0._rkind
    if (f2_mass.or.f2_energy) f2(:) = 0._rkind

    ! ---
    ! * compute sink terms...
    ! -----------------------

    ! intialize additional terms on the RHS as zero
    rAdd(:) = 0._rkind

    if (energy_flag) then
      ! compute energy associated with melt freeze for the vegetation canopy (J m-3)
      if (ixVegNrg/=integerMissing) rAdd(ixVegNrg) = rAdd(ixVegNrg) + LH_fus*( scalarCanopyIceTrial - scalarCanopyIce )/canopyDepth
 
      ! compute energy associated with melt/freeze for snow
      ! NOTE: allow expansion of ice during melt-freeze for snow; deny expansion of ice during melt-freeze for soil
      if (nSnowSoilNrg>0) then
        ! loop through non-missing energy state variables in the snow+soil domain
        do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)           
          select case( layerType(iLayer) )
            case(iname_snow); rAdd( ixSnowSoilNrg(iLayer) ) = rAdd( ixSnowSoilNrg(iLayer) )&
                                                          & + LH_fus*iden_ice  *( mLayerVolFracIceTrial(iLayer) - mLayerVolFracIce(iLayer) )
            case(iname_soil); rAdd( ixSnowSoilNrg(iLayer) ) = rAdd( ixSnowSoilNrg(iLayer) )&
                                                          & + LH_fus*iden_water*( mLayerVolFracIceTrial(iLayer) - mLayerVolFracIce(iLayer) )
          end select
        end do 
      end if
    end if

    if (mass_flag) then
      ! sink terms soil hydrology (-)
      ! NOTE 1: state variable is volumetric water content, so melt-freeze is not included
      ! NOTE 2: ground evaporation was already included in the flux at the upper boundary
      ! NOTE 3: rAdd(ixSnowOnlyWat)=0, and is defined in the initialization above
      ! NOTE 4: same sink terms for matric head and liquid matric potential
      if (nSoilOnlyHyd>0) then ! loop through non-missing hydrology state variables in the snow+soil domain
        do concurrent (iLayer=1:nSoil,ixSoilOnlyHyd(iLayer)/=integerMissing) 
          rAdd( ixSoilOnlyHyd(iLayer) ) = rAdd( ixSoilOnlyHyd(iLayer) )&
                                      & + ( ( mLayerTranspire(iLayer) - mLayerBaseflow(iLayer) )/mLayerDepth(iLayer+nSnow) - mLayerCompress(iLayer) )*dt
        end do
      end if
    end if

    ! ---
    ! * compute the residual vector...
    ! --------------------------------

    ! compute the residual vector for the vegetation canopy
    ! NOTE: sMul(ixVegHyd) = 1, but include as it converts all variables to quadruple precision
    ! --> energy balance
    if (energy_flag) then
      if (mixdformNrg) then
        if (ixCasNrg/=integerMissing) then
          fRHS(ixCasNrg) = ( fVec(ixCasNrg) + rAdd(ixCasNrg)/dt )
          rVec(ixCasNrg) = ( scalarCanairEnthalpyTrial - scalarCanairEnthalpy )&
                       & - ( fVec(ixCasNrg)*dt + rAdd(ixCasNrg) )
          if (f_flag) f(ixCasNrg) = real(rVec(ixCasNrg),rkind) ! nested Newton solver
          if (f1_energy) f1(ixCasNrg) = f(ixCasNrg) ! nested Newton solver
          if (f2_energy) f2(ixCasNrg) = -f(ixCasNrg) ! nested Newton solver
        end if
        if (ixVegNrg/=integerMissing) then
          fRHS(ixVegNrg) = ( fVec(ixVegNrg) + rAdd(ixVegNrg)/dt )
          rVec(ixVegNrg) = ( scalarCanopyEnthTempTrial - scalarCanopyEnthTemp )&
                         & - ( fVec(ixVegNrg)*dt + rAdd(ixVegNrg) )
          if (f_flag) f(ixVegNrg) = real(rVec(ixVegNrg),rkind) ! nested Newton solver
          if (f1_energy) f1(ixVegNrg) = f(ixVegNrg) ! nested Newton solver
          if (f2_energy) f2(ixVegNrg) = -f(ixVegNrg) ! nested Newton solver
        end if
      else
        if (ixCasNrg/=integerMissing) then
          fRHS(ixCasNrg) = ( fVec(ixCasNrg) + rAdd(ixCasNrg)/dt )/real(sMul(ixCasNrg),rkind)
          rVec(ixCasNrg) = sMul(ixCasNrg)*( scalarCanairTempTrial - scalarCanairTemp )&
                       & - ( fVec(ixCasNrg)*dt + rAdd(ixCasNrg) )
          if (f_flag) f(ixCasNrg) = real(rVec(ixCasNrg),rkind) ! nested Newton solver
          if (f1_energy) f1(ixCasNrg) = f(ixCasNrg) ! nested Newton solver
          if (f2_energy) f2(ixCasNrg) = -f(ixCasNrg) ! nested Newton solver
        end if
        if (ixVegNrg/=integerMissing) then
          fRHS(ixVegNrg) = 0._rkind ! not clear how to isolate the RHS function for ARKODE due to multiple time derivatives 
          rVec(ixVegNrg) = sMul(ixVegNrg)*( scalarCanopyTempTrial - scalarCanopyTemp )&
                       & + scalarCanopyCmTrial*( scalarCanopyWatTrial - scalarCanopyWat )/canopyDepth &
                       & - ( fVec(ixVegNrg)*dt + rAdd(ixVegNrg) )
          if (f_flag) f(ixVegNrg) = real(rVec(ixVegNrg),rkind) ! nested Newton solver
          if (f1_energy) f1(ixVegNrg) = f(ixVegNrg) ! nested Newton solver
          if (f2_energy) f2(ixVegNrg) = -f(ixVegNrg) ! nested Newton solver
        end if
      end if
    end if
    ! --> mass balance
    if (mass_flag) then
      if (ixVegHyd/=integerMissing) then
        scalarCanopyHydTrial = merge(scalarCanopyWatTrial, scalarCanopyLiqTrial, (ixStateType( ixHydCanopy(ixVegVolume) )==iname_watCanopy) )
        scalarCanopyHyd      = merge(scalarCanopyWat,      scalarCanopyLiq,      (ixStateType( ixHydCanopy(ixVegVolume) )==iname_watCanopy) )
        fRHS(ixVegHyd) = ( fVec(ixVegHyd) + rAdd(ixVegHyd)/dt )/real(sMul(ixVegHyd),rkind)
        rVec(ixVegHyd) = sMul(ixVegHyd)*scalarCanopyHydTrial - ( sMul(ixVegHyd)*scalarCanopyHyd + fVec(ixVegHyd)*dt + rAdd(ixVegHyd) )
        !rVec(ixVegHyd) = sMul(ixVegHyd)*scalarCanopyHydTrial - sMul(ixVegHyd)*scalarCanopyHyd&
        !             & - real(fVec(ixVegHyd)*dt + rAdd(ixVegHyd),qp) ! may need all terms to be qp before doing the sum to match original 
        if (f_flag) f(ixVegHyd) = real(rVec(ixVegHyd),rkind) ! nested Newton solver
        if (f1_mass) f1(ixVegHyd) = f(ixVegHyd) ! nested Newton solver
        if (f2_mass) f2(ixVegHyd) = -f(ixVegHyd) ! nested Newton solver
      end if
    end if

    ! compute the residual vector for the snow and soil sub-domains for energy
    if (energy_flag) then
      if (nSnowSoilNrg>0) then
        ! loop through non-missing energy state variables in the snow+soil domain
        do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   
          if (mixdformNrg) then
            fRHS( ixSnowSoilNrg(iLayer) ) = ( fVec( ixSnowSoilNrg(iLayer) ) + rAdd( ixSnowSoilNrg(iLayer) )/dt )
            rVec( ixSnowSoilNrg(iLayer) ) = ( mLayerEnthTempTrial(iLayer) - mLayerEnthTemp(iLayer) )&
                                        & - ( fVec( ixSnowSoilNrg(iLayer) )*dt + rAdd( ixSnowSoilNrg(iLayer) ) )
            if (f_flag) f(ixSnowSoilNrg(iLayer)) = real(rVec(ixSnowSoilNrg(iLayer)),rkind) ! nested Newton solver
            if (f1_energy) f1(ixSnowSoilNrg(iLayer)) = f(ixSnowSoilNrg(iLayer)) ! nested Newton solver
            if (f2_energy) f2(ixSnowSoilNrg(iLayer)) = -f(ixSnowSoilNrg(iLayer)) ! nested Newton solver
          else
            ! not clear how to isolate the RHS function for ARKODE due to multiple time derivatives
            ! so temperature formulation will not be used for ARKODE
            fRHS( ixSnowSoilNrg(iLayer) ) = 0._rkind ! set unused value to zero
            rVec( ixSnowSoilNrg(iLayer) ) = sMul( ixSnowSoilNrg(iLayer) )*( mLayerTempTrial(iLayer) - mLayerTemp(iLayer) )&
                                        & + mLayerCmTrial(iLayer)*( mLayerVolFracWatTrial(iLayer) - mLayerVolFracWat(iLayer) )&
                                          - ( fVec( ixSnowSoilNrg(iLayer) )*dt + rAdd( ixSnowSoilNrg(iLayer) ) )
            if (f_flag) f(ixSnowSoilNrg(iLayer)) = real(rVec(ixSnowSoilNrg(iLayer)),rkind) ! nested Newton solver
            if (f1_energy) f1(ixSnowSoilNrg(iLayer)) = f(ixSnowSoilNrg(iLayer)) ! nested Newton solver
            if (f2_energy) f2(ixSnowSoilNrg(iLayer)) = -f(ixSnowSoilNrg(iLayer)) ! nested Newton solver
          end if
        end do 
      end if
    end if

    ! compute the residual vector for the snow and soil sub-domains for hydrology
    ! NOTE: residual depends on choice of state variable
    if (mass_flag) then
      if (nSnowSoilHyd>0) then
        ! loop through non-missing hydrology state variables in the snow+soil domain
        do concurrent (iLayer=1:nLayers,ixSnowSoilHyd(iLayer)/=integerMissing) 
          ! get the correct state variable
          mLayerVolFracHydTrial(iLayer) = merge(mLayerVolFracWatTrial(iLayer), mLayerVolFracLiqTrial(iLayer) ,&
                                        & (ixHydType(iLayer)==iname_watLayer .or. ixHydType(iLayer)==iname_matLayer) )
          mLayerVolFracHyd(iLayer)      = merge(mLayerVolFracWat(iLayer),      mLayerVolFracLiq(iLayer),&
                                        & (ixHydType(iLayer)==iname_watLayer .or. ixHydType(iLayer)==iname_matLayer) )
          fRHS( ixSnowSoilHyd(iLayer) ) = 0._rkind !! SJT: temporary value that needs to be fixed
          rVec( ixSnowSoilHyd(iLayer) ) = ( mLayerVolFracHydTrial(iLayer) -  mLayerVolFracHyd(iLayer) )&
                                      & - ( fVec( ixSnowSoilHyd(iLayer) )*dt + rAdd( ixSnowSoilHyd(iLayer) ) )
          if (f_flag) f(ixSnowSoilHyd(iLayer)) = real(rVec(ixSnowSoilHyd(iLayer)),rkind) ! nested Newton solver
          if (f1_mass) f1(ixSnowSoilHyd(iLayer)) = f(ixSnowSoilHyd(iLayer)) ! nested Newton solver
          if (f2_mass) f2(ixSnowSoilHyd(iLayer)) = -f(ixSnowSoilHyd(iLayer)) ! nested Newton solver
        end do 
      end if
    end if

    ! compute the residual vector for the aquifer
    if (mass_flag) then ! aquifer residual included in mass portion of residual vector
      if (ixAqWat/=integerMissing) then
        fRHS(ixAqWat) = ( fVec(ixAqWat) + rAdd(ixAqWat)/dt )
        rVec(ixAqWat) = sMul(ixAqWat)*( scalarAquiferStorageTrial - scalarAquiferStorage )&
                    & - ( fVec(ixAqWat)*dt + rAdd(ixAqWat) )
        if (f_flag) f(ixAqWat) = real(rVec(ixAqWat),rkind) ! nested Newton solver
        if (f1_mass) f1(ixAqWat) = f(ixAqWat) ! nested Newton solver
        if (f2_mass) f2(ixAqWat) = -f(ixAqWat) ! nested Newton solver
      end if
    end if

    ! print the state variables if requested
    if(globalPrintFlag)then
      write(*,'(a)') 'In computResid:'
      write(*,'(a,i4)') '  nSnow = ', nSnow
      write(*,'(a,i4)') '  nSoil = ', nSoil
      write(*,'(a,i4)') '  nLayers = ', nLayers
      write(*,'(a,f12.5)') '  dt = ', dt
      write(*,'(a,e12.5)') '  scalarCanairTempTrial = ', scalarCanairTempTrial
      write(*,'(a,e12.5)') '  scalarCanopyTempTrial = ', scalarCanopyTempTrial
      write(*,'(a,e12.5)') '  scalarCanopyWatTrial = ', scalarCanopyWatTrial
      write(*,'(a,1x,100(e12.5,1x))') '  mLayerTempTrial = ', mLayerTempTrial(min(iJac1,size(mLayerTempTrial)):min(iJac2,size(mLayerTempTrial)))
      write(*,'(a,e12.5)') '  scalarAquiferStorageTrial = ', scalarAquiferStorageTrial
      write(*,'(a,e12.5)') '  scalarCanopyIceTrial = ', scalarCanopyIceTrial
      write(*,'(a,e12.5)') '  scalarCanopyLiqTrial = ', scalarCanopyLiqTrial
      write(*,'(a,1x,100(e12.5,1x))') '  mLayerVolFracIceTrial = ', mLayerVolFracIceTrial(min(iJac1,size(mLayerVolFracIceTrial)):min(iJac2,size(mLayerVolFracIceTrial)))
      write(*,'(a,1x,100(e12.5,1x))') '  mLayerVolFracWatTrial = ', mLayerVolFracWatTrial(min(iJac1,size(mLayerVolFracWatTrial)):min(iJac2,size(mLayerVolFracWatTrial)))
      write(*,'(a,1x,100(e12.5,1x))') '  mLayerVolFracLiqTrial = ', mLayerVolFracLiqTrial(min(iJac1,size(mLayerVolFracLiqTrial)):min(iJac2,size(mLayerVolFracLiqTrial)))
      write(*,'(a,e12.5)') '  scalarCanopyCmTrial = ', scalarCanopyCmTrial
      write(*,'(a,1x,100(e12.5,1x))') '  mLayerCmTrial = ', mLayerCmTrial(min(iJac1,size(mLayerCmTrial)):min(iJac2,size(mLayerCmTrial)))
      write(*,'(a,e12.5)') '  scalarCanairEnthalpyTrial = ', scalarCanairEnthalpyTrial 
      write(*,'(a,e12.5)') '  scalarCanopyEnthTempTrial = ', scalarCanopyEnthTempTrial
      write(*,'(a,1x,100(e12.5,1x))') '  mLayerEnthTempTrial = ', mLayerEnthTempTrial(min(iJac1,size(mLayerEnthTempTrial)):min(iJac2,size(mLayerEnthTempTrial)))
      write(*,'(a,1x,100(e12.5,1x))') 'sMul = ', sMul(min(iJac1,size(sMul)):min(iJac2,size(sMul)))
    endif

    ! print result
    if(globalPrintFlag .or. any(isNan(rVec)))then
      write(*,'(a,1x,100(e12.5,1x))') 'fVec = ', fVec(min(iJac1,size(rVec)):min(iJac2,size(rVec)))
      write(*,'(a,1x,100(e12.5,1x))') 'rVec = ', rVec(min(iJac1,size(rVec)):min(iJac2,size(rVec)))
    endif
    if(any(isNan(rVec)))then; message=trim(message)//'NaN in residuals'; err=20; return; endif

  end associate

end subroutine computResid

end module computResid_module
