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

module updateVars_module

! data types
USE nrtype

! missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! access the global print flag
USE globalData,only:globalPrintFlag

! domain types
USE globalData,only:iname_cas       ! named variables for canopy air space
USE globalData,only:iname_veg       ! named variables for vegetation canopy
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil
USE globalData,only:iname_aquifer   ! named variables for the aquifer

! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of total water on the vegetation canopy
USE globalData,only:iname_liqCanopy ! named variable defining the mass of liquid water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers

! metadata for information in the data structures
USE globalData,only:indx_meta       ! metadata for the variables in the index structure

! constants
USE multiconst,only:&
                    Tfreeze,      & ! temperature at freezing              (K)
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (rkind)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    zLookup,      & ! data vector with variable length dimension (rkind)
                    model_options,& ! defines the model decisions
                    var_dlength     ! data vector with variable length dimension (rkind)

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookDIAG             ! named variables for structure elements
USE var_lookup,only:iLookPROG             ! named variables for structure elements
USE var_lookup,only:iLookDERIV            ! named variables for structure elements
USE var_lookup,only:iLookPARAM            ! named variables for structure elements
USE var_lookup,only:iLookINDEX            ! named variables for structure elements
USE var_lookup,only:iLookFLUX             ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS        ! named variables for elements of the decision structure

! look-up values for the numerical method
USE mDecisions_module,only:         &
                    homegrown      ,& ! homegrown backward Euler solution using concepts from numerical recipes
                    kinsol         ,& ! SUNDIALS backward Euler solution using Kinsol
                    ida               ! SUNDIALS solution using IDA

! constants
USE multiconst,only:&
                    Tfreeze,        & ! freezing temperature                 (K)
                    LH_fus,         & ! latent heat of fusion                (J kg-1)
                    LH_vap,         & ! latent heat of vaporization          (J kg-1)
                    iden_ice,       & ! intrinsic density of ice             (kg m-3)
                    iden_water        ! intrinsic density of liquid water    (kg m-3)

! provide access to routines to update states
USE updatState_module,only:updateSnow     ! update snow states
USE updatState_module,only:updateSoil     ! update soil states

! provide access to functions for the constitutive functions and derivatives
USE snow_utils_module,only:fracliquid          ! compute the fraction of liquid water (snow)
USE snow_utils_module,only:dFracLiq_dTk        ! differentiate the freezing curve w.r.t. temperature (snow)
USE soil_utils_module,only:dTheta_dTk          ! differentiate the freezing curve w.r.t. temperature (soil)
USE soil_utils_module,only:dTheta_dPsi         ! derivative in the soil water characteristic (soil)
USE soil_utils_module,only:matricHead          ! compute the matric head based on volumetric water content
USE soil_utils_module,only:volFracLiq          ! compute volumetric fraction of liquid water
USE soil_utils_module,only:crit_soilT          ! compute critical temperature below which ice exists
USE soil_utils_module,only:liquidHead          ! compute the liquid water matric potential
USE enthalpyTemp_module,only:T2enthTemp_cas    ! convert temperature to enthalpy for canopy air space
USE enthalpyTemp_module,only:T2enthTemp_veg   ! convert temperature to enthalpy for vegetation
USE enthalpyTemp_module,only:T2enthTemp_snow   ! convert temperature to enthalpy for snow
USE enthalpyTemp_module,only:T2enthTemp_soil   ! convert temperature to enthalpy for soil 

! IEEE check
USE, intrinsic :: ieee_arithmetic         ! check values (NaN, etc.)

implicit none
private
public::updateVars
public::updateProg

contains

! **********************************************************************************************************
! public subroutine updateVars: compute diagnostic variables and derivatives
! **********************************************************************************************************
subroutine updateVars(&
                      ! input
                      computeEnthTemp,                           & ! intent(in):    flag if computing temperature compoment of enthalpy
                      use_lookup,                                & ! intent(in):    flag to use the lookup table for soil enthalpy 
                      do_adjustTemp,                             & ! intent(in):    logical flag to adjust temperature to account for the energy used in melt+freeze
                      mpar_data,                                 & ! intent(in):    model parameters for a local HRU
                      indx_data,                                 & ! intent(in):    indices defining model states and layers
                      prog_data,                                 & ! intent(in):    model prognostic variables for a local HRU
                      diag_data,                                 & ! intent(inout): model diagnostic variables for a local HRU
                      deriv_data,                                & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                      lookup_data,                               & ! intent(in):    lookup table data structure
                      scalarCanairTempTrial,                     & ! intent(in):    trial value of canopy air space temperature (K)
                      ! output: variables for the vegetation canopy
                      scalarCanopyTempTrial,                     & ! intent(inout): trial value of canopy temperature (K)
                      scalarCanopyWatTrial,                      & ! intent(inout): trial value of canopy total water (kg m-2)
                      scalarCanopyLiqTrial,                      & ! intent(inout): trial value of canopy liquid water (kg m-2)
                      scalarCanopyIceTrial,                      & ! intent(inout): trial value of canopy ice content (kg m-2)
                      ! output: variables for the snow-soil domain
                      mLayerTempTrial,                           & ! intent(inout): trial vector of layer temperature (K)
                      mLayerVolFracWatTrial,                     & ! intent(inout): trial vector of volumetric total water content (-)
                      mLayerVolFracLiqTrial,                     & ! intent(inout): trial vector of volumetric liquid water content (-)
                      mLayerVolFracIceTrial,                     & ! intent(inout): trial vector of volumetric ice water content (-)
                      mLayerMatricHeadTrial,                     & ! intent(inout): trial vector of total water matric potential (m)
                      mLayerMatricHeadLiqTrial,                  & ! intent(inout): trial vector of liquid water matric potential (m)
                      ! output: enthalpy variables  
                      scalarCanairEnthalpyTrial,                 & ! intent(inout): trial value for enthalpy of the canopy air space (J m-3)
                      scalarCanopyEnthTempTrial,                 & ! intent(inout): trial value for temperature component of enthalpy of the vegetation canopy (J m-3)
                      mLayerEnthTempTrial,                       & ! intent(inout): trial vector of temperature component of enthalpy of each snow+soil layer (J m-3)                          
                      ! output: error control
                      err,message)                                 ! intent(out):   error control
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------------------------------------
  implicit none
  ! input
  logical(lgt)     ,intent(in)       :: computeEnthTemp                 ! flag if computing temperature compoment of enthalpy
  logical(lgt)     ,intent(in)       :: use_lookup                      ! flag to use the lookup table for soil enthalpy
  logical(lgt)     ,intent(in)       :: do_adjustTemp                   ! flag to adjust temperature to account for the energy used in melt+freeze
  type(var_dlength),intent(in)       :: mpar_data                       ! model parameters for a local HRU
  type(var_ilength),intent(in)       :: indx_data                       ! indices defining model states and layers
  type(var_dlength),intent(in)       :: prog_data                       ! prognostic variables for a local HRU
  type(var_dlength),intent(inout)    :: diag_data                       ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout)    :: deriv_data                      ! derivatives in model fluxes w.r.t. relevant state variables
  type(zLookup)    ,intent(in)       :: lookup_data                     ! lookup tables
  real(rkind),intent(in)             :: scalarCanairTempTrial           ! trial value of canopy air space temperature (K)
  ! output: variables for the vegetation canopy
  real(rkind),intent(inout)          :: scalarCanopyTempTrial           ! trial value of canopy temperature (K)
  real(rkind),intent(inout)          :: scalarCanopyWatTrial            ! trial value of canopy total water (kg m-2)
  real(rkind),intent(inout)          :: scalarCanopyLiqTrial            ! trial value of canopy liquid water (kg m-2)
  real(rkind),intent(inout)          :: scalarCanopyIceTrial            ! trial value of canopy ice content (kg m-2)
  ! output: variables for the snow-soil domain
  real(rkind),intent(inout)          :: mLayerTempTrial(:)              ! trial vector of layer temperature (K)
  real(rkind),intent(inout)          :: mLayerVolFracWatTrial(:)        ! trial vector of volumetric total water content (-)
  real(rkind),intent(inout)          :: mLayerVolFracLiqTrial(:)        ! trial vector of volumetric liquid water content (-)
  real(rkind),intent(inout)          :: mLayerVolFracIceTrial(:)        ! trial vector of volumetric ice water content (-)
  real(rkind),intent(inout)          :: mLayerMatricHeadTrial(:)        ! trial vector of total water matric potential (m)
  real(rkind),intent(inout)          :: mLayerMatricHeadLiqTrial(:)     ! trial vector of liquid water matric potential (m)
  ! output: enthalpy variables
  real(rkind),intent(inout)          :: scalarCanairEnthalpyTrial       ! trial value for enthalpy of the canopy air space (J m-3)
  real(rkind),intent(inout)          :: scalarCanopyEnthTempTrial       ! trial value for temperature component of enthalpy of the vegetation canopy (J m-3)
  real(rkind),intent(inout)          :: mLayerEnthTempTrial(:)          ! trial vector of temperature component of enthalpy of each snow+soil layer (J m-3)
  ! output: error control
  integer(i4b),intent(out)           :: err                             ! error code
  character(*),intent(out)           :: message                         ! error message
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! general local variables
  integer(i4b)                       :: iState                          ! index of model state variable
  integer(i4b)                       :: iLayer                          ! index of layer within the snow+soil domain
  integer(i4b)                       :: ixFullVector                    ! index within full state vector
  integer(i4b)                       :: ixDomainType                    ! name of a given model domain
  integer(i4b)                       :: ixControlIndex                  ! index within a given model domain
  integer(i4b)                       :: ixOther,ixOtherLocal            ! index of the coupled state variable within the (full, local) vector
  logical(lgt)                       :: isCoupled                       ! .true. if a given variable shared another state variable in the same control volume
  logical(lgt)                       :: isNrgState                      ! .true. if a given variable is an energy state
  logical(lgt),allocatable           :: computedCoupling(:)             ! .true. if computed the coupling for a given state variable
  real(rkind)                        :: scalarVolFracLiq                ! volumetric fraction of liquid water (-)
  real(rkind)                        :: scalarVolFracIce                ! volumetric fraction of ice (-)
  real(rkind)                        :: Tcrit                           ! critical soil temperature below which ice exists (K)
  real(rkind)                        :: xTemp                           ! temporary temperature (K)
  real(rkind)                        :: effSat                          ! effective saturation (-)
  real(rkind)                        :: avPore                          ! available pore space (-)
  character(len=256)                 :: cMessage                        ! error message of downwind routine
  logical(lgt),parameter             :: printFlag=.false.               ! flag to turn on printing
  ! iterative solution for temperature
  real(rkind)                        :: meltNrg                         ! energy for melt+freeze (J m-3)
  real(rkind)                        :: residual                        ! residual in the energy equation (J m-3)
  real(rkind)                        :: derivative                      ! derivative in the energy equation (J m-3 K-1)
  real(rkind)                        :: tempInc                         ! iteration increment (K)
  integer(i4b)                       :: iter                            ! iteration index
  integer(i4b)                       :: niter                           ! number of iterations
  integer(i4b),parameter             :: maxiter=100                     ! maximum number of iterations
  real(rkind),parameter              :: nrgConvTol=1.e-4_rkind          ! convergence tolerance for energy (J m-3)
  real(rkind),parameter              :: tempConvTol=1.e-6_rkind         ! convergence tolerance for temperature (K)
  real(rkind)                        :: critDiff                        ! temperature difference from critical (K)
  real(rkind)                        :: tempMin                         ! minimum bracket for temperature (K)
  real(rkind)                        :: tempMax                         ! maximum bracket for temperature (K)
  logical(lgt)                       :: bFlag                           ! flag to denote that iteration increment was constrained using bi-section
  real(rkind),parameter              :: epsT=1.e-7_rkind                ! small interval above/below critical temperature (K)
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! make association with variables in the data structures
  associate(&
    ! number of model layers, and layer type
    nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in):  [i4b]    total number of snow layers
    nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in):  [i4b]    total number of soil layers
    nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)               ,& ! intent(in):  [i4b]    total number of snow and soil layers
    mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,& ! intent(in):  [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
    ! indices defining model states and layers
    ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in):  [i4b]    index of canopy energy state variable
    ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,& ! intent(in):  [i4b]    index of canopy hydrology state variable (mass)
    ! indices in the full vector for specific domains
    ixNrgCanair             => indx_data%var(iLookINDEX%ixNrgCanair)%dat              ,& ! intent(in):  [i4b(:)] indices IN THE FULL VECTOR for energy states in canopy air space domain
    ixNrgCanopy             => indx_data%var(iLookINDEX%ixNrgCanopy)%dat              ,& ! intent(in):  [i4b(:)] indices IN THE FULL VECTOR for energy states in the canopy domain
    ixHydCanopy             => indx_data%var(iLookINDEX%ixHydCanopy)%dat              ,& ! intent(in):  [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the canopy domain
    ixNrgLayer              => indx_data%var(iLookINDEX%ixNrgLayer)%dat               ,& ! intent(in):  [i4b(:)] indices IN THE FULL VECTOR for energy states in the snow+soil domain
    ixHydLayer              => indx_data%var(iLookINDEX%ixHydLayer)%dat               ,& ! intent(in):  [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the snow+soil domain
    ! mapping between the full state vector and the state subset
    ixMapFull2Subset        => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat         ,& ! intent(in):  [i4b(:)] list of indices in the state subset for each state in the full state vector
    ixMapSubset2Full        => indx_data%var(iLookINDEX%ixMapSubset2Full)%dat         ,& ! intent(in):  [i4b(:)] [state subset] list of indices of the full state vector in the state subset
    ! type of domain, type of state variable, and index of control volume within domain
    ixDomainType_subset     => indx_data%var(iLookINDEX%ixDomainType_subset)%dat      ,& ! intent(in):  [i4b(:)] [state subset] id of domain for desired model state variables
    ixControlVolume         => indx_data%var(iLookINDEX%ixControlVolume)%dat          ,& ! intent(in):  [i4b(:)] index of the control volume for different domains (veg, snow, soil)
    ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat              ,& ! intent(in):  [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
    ! snow parameters
    snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)         ,& ! intent(in):  [dp]     scaling parameter for the snow freezing curve (K-1)
    ! depth-varying model parameters (heat capacity, enthalpy)
    specificHeatVeg         => mpar_data%var(iLookPARAM%specificHeatVeg)%dat(1)       ,& ! intent(in):  [dp   ]  specific heat of vegetation (J kg-1 K-1)
    maxMassVegetation       => mpar_data%var(iLookPARAM%maxMassVegetation)%dat(1)     ,& ! intent(in):  [dp   ]  maximum mass of vegetation (kg m-2)
    soil_dens_intr          => mpar_data%var(iLookPARAM%soil_dens_intr)%dat           ,& ! intent(in):  [dp(:)]  intrinsic soil density (kg m-3)
    vGn_m                   => diag_data%var(iLookDIAG%scalarVGn_m)%dat               ,& ! intent(in):  [dp(:)]  van Genutchen "m" parameter (-)
    vGn_n                   => mpar_data%var(iLookPARAM%vGn_n)%dat                    ,& ! intent(in):  [dp(:)]  van Genutchen "n" parameter (-)
    vGn_alpha               => mpar_data%var(iLookPARAM%vGn_alpha)%dat                ,& ! intent(in):  [dp(:)]  van Genutchen "alpha" parameter (m-1)
    theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat                ,& ! intent(in):  [dp(:)]  soil porosity (-)
    theta_res               => mpar_data%var(iLookPARAM%theta_res)%dat                ,& ! intent(in):  [dp(:)]  soil residual volumetric water content (-)
    ! model diagnostic variables (heat capacity)
    canopyDepth             => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)      ,& ! intent(in):  [dp   ]  canopy depth (m)
    scalarBulkVolHeatCapVeg => diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1),& ! intent(in):  [dp   ]  volumetric heat capacity of the vegetation (J m-3 K-1)
    mLayerVolHtCapBulk      => diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat        ,& ! intent(in):  [dp(:)]  volumetric heat capacity in each layer (J m-3 K-1)
    ! model diagnostic variables (fraction of liquid water)
    scalarFracLiqVeg        => diag_data%var(iLookDIAG%scalarFracLiqVeg)%dat(1)       ,& ! intent(out): [dp]     fraction of liquid water on vegetation (-)
    mLayerFracLiqSnow       => diag_data%var(iLookDIAG%mLayerFracLiqSnow)%dat         ,& ! intent(out): [dp(:)]  fraction of liquid water in each snow layer (-)
    ! model states from a previous solution
    scalarCanopyTemp        => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,& ! intent(in):  [dp]     temperature of the vegetation canopy (K)
    mLayerTemp              => prog_data%var(iLookPROG%mLayerTemp)%dat                ,& ! intent(in):  [dp(:)]  temperature of each snow/soil layer (K)
    scalarCanopyWat         => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)        ,& ! intent(in):  [dp]     mass of total water on the vegetation canopy (kg m-2)
    mLayerVolFracWat        => prog_data%var(iLookPROG%mLayerVolFracWat)%dat          ,& ! intent(in):  [dp(:)]  volumetric fraction of total water (-)
    ! model diagnostic variables from a previous solution
    scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)        ,& ! intent(in):  [dp(:)]  mass of ice on the vegetation canopy (kg m-2)
    mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat          ,& ! intent(in):  [dp(:)]  volumetric fraction of ice (-)
    ! derivatives
    dVolTot_dPsi0           => deriv_data%var(iLookDERIV%dVolTot_dPsi0   )%dat        ,& ! intent(out): [dp(:)]  derivative in total water content w.r.t. total water matric potential
    dPsiLiq_dPsi0           => deriv_data%var(iLookDERIV%dPsiLiq_dPsi0   )%dat        ,& ! intent(out): [dp(:)]  derivative in liquid water matric pot w.r.t. the total water matric pot (-)
    dPsiLiq_dTemp           => deriv_data%var(iLookDERIV%dPsiLiq_dTemp   )%dat        ,& ! intent(out): [dp(:)]  derivative in the liquid water matric potential w.r.t. temperature
    mLayerdTheta_dTk        => deriv_data%var(iLookDERIV%mLayerdTheta_dTk)%dat        ,& ! intent(out): [dp(:)]  derivative of volumetric liquid water content w.r.t. temperature
    dTheta_dTkCanopy        => deriv_data%var(iLookDERIV%dTheta_dTkCanopy)%dat(1)     ,& ! intent(out): [dp]     derivative of volumetric liquid water content w.r.t. temperature
    dFracLiqWat_dTk        => deriv_data%var(iLookDERIV%dFracLiqWat_dTk)%dat          ,& ! intent(out): [dp(:)]  derivative in fraction of liquid water w.r.t. temperature
    dFracLiqVeg_dTkCanopy   => deriv_data%var(iLookDERIV%dFracLiqVeg_dTkCanopy)%dat(1),& ! intent(out): [dp   ]  derivative in fraction of (throughfall + drainage) w.r.t. temperature
    ! derivatives inside solver for Jacobian only
    mLayerdTemp_dt          => deriv_data%var(iLookDERIV%mLayerdTemp_dt )%dat         ,& ! intent(out): [dp(:)]  timestep change in layer temperature
    scalarCanopydTemp_dt    => deriv_data%var(iLookDERIV%scalarCanopydTemp_dt)%dat(1) ,& ! intent(out): [dp   ]  timestep change in canopy temperature
    mLayerdWat_dt           => deriv_data%var(iLookDERIV%mLayerdWat_dt)%dat           ,& ! intent(out): [dp(:)]  timestep change in layer volumetric fraction of total water
    scalarCanopydWat_dt     => deriv_data%var(iLookDERIV%scalarCanopydWat_dt)%dat(1)   & ! intent(out): [dp   ]  timestep change in canopy total water
    ) ! association with variables in the data structures

    ! --------------------------------------------------------------------------------------------------------------------------------
    ! --------------------------------------------------------------------------------------------------------------------------------

    ! initialize error control
    err=0; message='updateVars/'

    ! allocate space and assign values to the flag vector
    allocate(computedCoupling(size(ixMapSubset2Full)),stat=err)        ! .true. if computed the coupling for a given state variable
    if(err/=0)then; message=trim(message)//'problem allocating computedCoupling'; return; endif
    computedCoupling(:)=.false.

    ! loop through model state variables
    do iState=1,size(ixMapSubset2Full)

      ! check the need for the computations
      if(computedCoupling(iState)) cycle

      ! -----
      ! - compute indices...
      ! --------------------

      ! get domain type, and index of the control volume within the domain
      ixFullVector   = ixMapSubset2Full(iState)       ! index within full state vector
      ixDomainType   = ixDomainType_subset(iState)    ! named variables defining the domain (iname_cas, iname_veg, etc.)
      ixControlIndex = ixControlVolume(ixFullVector)  ! index within a given domain

      ! get the layer index
      select case(ixDomainType)
        case(iname_cas);     iLayer = 0
        case(iname_veg);     iLayer = 0
        case(iname_snow);    iLayer = ixControlIndex
        case(iname_soil);    iLayer = ixControlIndex + nSnow
        case(iname_aquifer); cycle ! aquifer: do nothing
        case default; err=20; message=trim(message)//'expect case to be iname_cas, iname_veg, iname_snow, iname_soil, iname_aquifer'; return
      end select

      ! get the index of the other (energy or mass) state variable within the full state vector
      select case(ixDomainType)
        case(iname_cas)             ; ixOther = integerMissing
        case(iname_veg)             ; ixOther = merge(ixHydCanopy(1),    ixNrgCanopy(1),    ixStateType(ixFullVector)==iname_nrgCanopy)
        case(iname_snow, iname_soil); ixOther = merge(ixHydLayer(iLayer),ixNrgLayer(iLayer),ixStateType(ixFullVector)==iname_nrgLayer)
        case default; err=20; message=trim(message)//'expect case to be iname_veg, iname_snow, iname_soil'; return
      end select

      ! get the index in the local state vector
      if(ixDomainType==iname_cas)then
        ixOtherLocal = integerMissing
      else
        ixOtherLocal = ixMapFull2Subset(ixOther)  ! ixOtherLocal could equal integerMissing
      endif
      if(ixOtherLocal/=integerMissing) computedCoupling(ixOtherLocal)=.true.

      ! check if we have a coupled solution
      isCoupled    = (ixOtherLocal/=integerMissing)

      ! check if we are an energy state
      isNrgState   = (ixStateType(ixFullVector)==iname_nrgCanopy .or. ixStateType(ixFullVector)==iname_nrgLayer)

      if(printFlag)then
        print*, 'iState         = ', iState, size(ixMapSubset2Full)
        print*, 'ixFullVector   = ', ixFullVector
        print*, 'ixDomainType   = ', ixDomainType
        print*, 'ixControlIndex = ', ixControlIndex
        print*, 'ixOther        = ', ixOther
        print*, 'ixOtherLocal   = ', ixOtherLocal
        print*, 'do_adjustTemp  = ', do_adjustTemp
        print*, 'isCoupled      = ', isCoupled
        print*, 'isNrgState     = ', isNrgState
      endif

      ! calculate temperature component of enthalpy for canopy air space
      if(ixDomainType==iname_cas)then
        if(computeEnthTemp)then
          call T2enthTemp_cas(&
                      scalarCanairTempTrial,       & ! intent(in): canopy air temperature (K)
                      scalarCanairEnthalpyTrial)     ! intent(out): enthalpy of the canopy air space (J m-3)
        else
          scalarCanairEnthalpyTrial = realMissing
        endif
        cycle ! no more to do on canopy air space
      endif

      ! update hydrology state variables for the uncoupled solution
      if(.not.isNrgState .and. .not.isCoupled)then

        ! update the total water from volumetric liquid water
        if(ixStateType(ixFullVector)==iname_liqCanopy .or. ixStateType(ixFullVector)==iname_liqLayer)then
          select case(ixDomainType)
            case(iname_veg);    scalarCanopyWatTrial          = scalarCanopyLiqTrial          + scalarCanopyIceTrial
            case(iname_snow);   mLayerVolFracWatTrial(iLayer) = mLayerVolFracLiqTrial(iLayer) + mLayerVolFracIceTrial(iLayer)*iden_ice/iden_water
            case(iname_soil);   mLayerVolFracWatTrial(iLayer) = mLayerVolFracLiqTrial(iLayer) + mLayerVolFracIceTrial(iLayer) ! no volume expansion
            case default; err=20; message=trim(message)//'expect case to be iname_veg, iname_snow, or iname_soil'; return
          end select
        endif

        ! update the total water and the total water matric potential
        if(ixDomainType==iname_soil)then
          select case( ixStateType(ixFullVector) )
            ! --> update the total water from the liquid water matric potential
            case(iname_lmpLayer)
              effSat = volFracLiq(mLayerMatricHeadLiqTrial(ixControlIndex),vGn_alpha(ixControlIndex),0._rkind,1._rkind,vGn_n(ixControlIndex),vGn_m(ixControlIndex))  ! effective saturation
              avPore = theta_sat(ixControlIndex) - mLayerVolFracIceTrial(iLayer) - theta_res(ixControlIndex)  ! available pore space
              mLayerVolFracLiqTrial(iLayer) = effSat*avPore + theta_res(ixControlIndex)
              mLayerVolFracWatTrial(iLayer) = mLayerVolFracLiqTrial(iLayer) + mLayerVolFracIceTrial(iLayer) ! no volume expansion
              mLayerMatricHeadTrial(ixControlIndex) = matricHead(mLayerVolFracWatTrial(iLayer),vGn_alpha(ixControlIndex),theta_res(ixControlIndex),theta_sat(ixControlIndex),vGn_n(ixControlIndex),vGn_m(ixControlIndex))
            ! --> update the total water from the total water matric potential
            case(iname_matLayer)
              mLayerVolFracWatTrial(iLayer) = volFracLiq(mLayerMatricHeadTrial(ixControlIndex),vGn_alpha(ixControlIndex),theta_res(ixControlIndex),theta_sat(ixControlIndex),vGn_n(ixControlIndex),vGn_m(ixControlIndex))
            ! --> update the total water matric potential (assume already have mLayerVolFracWatTrial given block above)
            case(iname_liqLayer, iname_watLayer)
              mLayerMatricHeadTrial(ixControlIndex) = matricHead(mLayerVolFracWatTrial(iLayer),vGn_alpha(ixControlIndex),theta_res(ixControlIndex),theta_sat(ixControlIndex),vGn_n(ixControlIndex),vGn_m(ixControlIndex))
            case default; err=20; message=trim(message)//'expect iname_lmpLayer, iname_matLayer, iname_liqLayer, or iname_watLayer'; return
          end select
        endif  ! if in the soil domain

      endif  ! if hydrology state variable or uncoupled solution

      ! compute the critical soil temperature below which ice exists
      select case(ixDomainType)
        case(iname_veg, iname_snow); Tcrit = Tfreeze
        case(iname_soil);            Tcrit = crit_soilT( mLayerMatricHeadTrial(ixControlIndex) )
        case default; err=20; message=trim(message)//'expect case to be iname_veg, iname_snow, iname_soil'; return
      end select

      ! initialize temperature
      select case(ixDomainType)
        case(iname_veg);              xTemp = scalarCanopyTempTrial
        case(iname_snow, iname_soil); xTemp = mLayerTempTrial(iLayer)
        case default; err=20; message=trim(message)//'expect case to be iname_veg, iname_snow, iname_soil'; return
      end select

      ! define brackets for the root
      ! NOTE: start with an enormous range; updated quickly in the iterations
      tempMin = xTemp - 10._rkind
      tempMax = xTemp + 10._rkind

      ! get iterations (set to maximum iterations if adjusting the temperature)
      niter = merge(maxiter, 1, do_adjustTemp)

      ! iterate
      iterations: do iter=1,niter

        ! restrict temperature
        if(xTemp <= tempMin .or. xTemp >= tempMax)then
          xTemp = 0.5_rkind*(tempMin + tempMax)  ! new value
          bFlag = .true.
        else
          bFlag = .false.
        endif

        ! -----
        ! - compute derivatives...
        ! ------------------------

        ! compute temperature time derivatives
        select case(ixDomainType)
          case(iname_veg); scalarCanopydTemp_dt = xTemp - scalarCanopyTemp
          case(iname_snow, iname_soil); mLayerdTemp_dt(iLayer) = xTemp - mLayerTemp(iLayer)
        end select

        ! compute the derivative in total water content w.r.t. total water matric potential (m-1)
        ! NOTE 1: valid for frozen and unfrozen conditions
        ! NOTE 2: for case "iname_lmpLayer", dVolTot_dPsi0 = dVolLiq_dPsi
        if(ixDomainType==iname_soil)then
          select case( ixStateType(ixFullVector) )
            case(iname_lmpLayer); dVolTot_dPsi0(ixControlIndex) = dTheta_dPsi(mLayerMatricHeadLiqTrial(ixControlIndex),vGn_alpha(ixControlIndex),0._rkind,1._rkind,vGn_n(ixControlIndex),vGn_m(ixControlIndex))*avPore
            case default;         dVolTot_dPsi0(ixControlIndex) = dTheta_dPsi(mLayerMatricHeadTrial(ixControlIndex),vGn_alpha(ixControlIndex),theta_res(ixControlIndex),theta_sat(ixControlIndex),vGn_n(ixControlIndex),vGn_m(ixControlIndex))
          end select
        endif

        ! compute the derivative in liquid water content w.r.t. temperature
        ! --> partially frozen: dependence of liquid water on temperature
        if(xTemp<Tcrit)then
          select case(ixDomainType)
            case(iname_veg)
              dFracLiqVeg_dTkCanopy = dFracLiq_dTk(xTemp,snowfrz_scale)
              dTheta_dTkCanopy = dFracLiqVeg_dTkCanopy * scalarCanopyWatTrial/(iden_water*canopyDepth)
            case(iname_snow)
              dFracLiqWat_dTk(iLayer) = dFracLiq_dTk(xTemp,snowfrz_scale)
              mLayerdTheta_dTk(iLayer) = dFracLiqWat_dTk(iLayer) * mLayerVolFracWatTrial(iLayer)
            case(iname_soil)
              dFracLiqWat_dTk(iLayer) = 0._rkind !dTheta_dTk(xTemp,theta_res(ixControlIndex),theta_sat(ixControlIndex),vGn_alpha(ixControlIndex),vGn_n(ixControlIndex),vGn_m(ixControlIndex))/ mLayerVolFracWatTrial(iLayer)
              mLayerdTheta_dTk(iLayer) = dTheta_dTk(xTemp,theta_res(ixControlIndex),theta_sat(ixControlIndex),vGn_alpha(ixControlIndex),vGn_n(ixControlIndex),vGn_m(ixControlIndex))
            case default; err=20; message=trim(message)//'expect case to be iname_veg, iname_snow, iname_soil'; return
          end select  ! domain type

        ! --> unfrozen: no dependence of liquid water on temperature
        else
          select case(ixDomainType)
            case(iname_veg);              dTheta_dTkCanopy         = 0._rkind; dFracLiqVeg_dTkCanopy   = 0._rkind
            case(iname_snow, iname_soil); mLayerdTheta_dTk(iLayer) = 0._rkind; dFracLiqWat_dTk(iLayer) = 0._rkind
            case default; err=20; message=trim(message)//'expect case to be iname_veg, iname_snow, iname_soil'; return
          end select  ! domain type
        endif


        ! -----
        ! - update volumetric fraction of liquid water and ice...
        !    => case of hydrology state uncoupled with energy (and when not adjusting the temperature)...
        ! -----------------------------------------------------------------------------------------------

        ! case of hydrology state uncoupled with energy (and when not adjusting the temperature)
        if(.not.do_adjustTemp .and. .not.isNrgState .and. .not.isCoupled)then

          ! compute the fraction of snow
          select case(ixDomainType)
            case(iname_veg);  scalarFracLiqVeg          = fracliquid(xTemp,snowfrz_scale)
            case(iname_snow); mLayerFracLiqSnow(iLayer) = fracliquid(xTemp,snowfrz_scale)
            case(iname_soil)  ! do nothing
            case default; err=20; message=trim(message)//'expect case to be iname_veg, iname_snow, iname_soil'; return
          end select  ! domain type

          ! -----
          ! - update volumetric fraction of liquid water and ice...
          !    => case of energy state or coupled solution (or adjusting the temperature)...
          ! --------------------------------------------------------------------------------

          ! case of energy state OR coupled solution (or adjusting the temperature)
        elseif(do_adjustTemp .or. ( (isNrgState .or. isCoupled) ) )then

          ! identify domain type
          select case(ixDomainType)

            ! *** vegetation canopy
            case(iname_veg)

              ! compute volumetric fraction of liquid water and ice
              call updateSnow(xTemp,                                        & ! intent(in):  temperature (K)
                              scalarCanopyWatTrial/(iden_water*canopyDepth),& ! intent(in):  volumetric fraction of total water (-)
                              snowfrz_scale,                                & ! intent(in):  scaling parameter for the snow freezing curve (K-1)
                              scalarVolFracLiq,                             & ! intent(out): trial volumetric fraction of liquid water (-)
                              scalarVolFracIce,                             & ! intent(out): trial volumetric fraction if ice (-)
                              scalarFracLiqVeg,                             & ! intent(out): fraction of liquid water (-)
                              err,cmessage)                                   ! intent(out): error control
              if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

              ! compute mass of water on the canopy
              ! NOTE: possibilities for speed-up here
              scalarCanopyLiqTrial =             scalarFracLiqVeg *scalarCanopyWatTrial !(kg m-2), scalarVolFracLiq*iden_water*canopyDepth
              scalarCanopyIceTrial = (1._rkind - scalarFracLiqVeg)*scalarCanopyWatTrial !(kg m-2), scalarVolFracIce* iden_ice *canopyDepth

            ! *** snow layers
            case(iname_snow)

              ! compute volumetric fraction of liquid water and ice
              call updateSnow(xTemp,                          & ! intent(in):  temperature (K)
                              mLayerVolFracWatTrial(iLayer),  & ! intent(in):  mass state variable = trial volumetric fraction of water (-)
                              snowfrz_scale,                  & ! intent(in):  scaling parameter for the snow freezing curve (K-1)
                              mLayerVolFracLiqTrial(iLayer),  & ! intent(out): trial volumetric fraction of liquid water (-)
                              mLayerVolFracIceTrial(iLayer),  & ! intent(out): trial volumetric fraction if ice (-)
                              mLayerFracLiqSnow(iLayer),      & ! intent(out): fraction of liquid water (-)
                              err,cmessage)                     ! intent(out): error control
              if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

              ! *** soil layers
            case(iname_soil)

              ! compute volumetric fraction of liquid water and ice
              call updateSoil(xTemp,                                  & ! intent(in):  temperature (K)
                              mLayerMatricHeadTrial(ixControlIndex),  & ! intent(in):  total water matric potential (m)
                              vGn_alpha(ixControlIndex),vGn_n(ixControlIndex),theta_sat(ixControlIndex),theta_res(ixControlIndex),vGn_m(ixControlIndex), & ! intent(in): soil parameters
                              mLayerVolFracWatTrial(iLayer),          & ! intent(in):  mass state variable = trial volumetric fraction of water (-)
                              mLayerVolFracLiqTrial(iLayer),          & ! intent(out): trial volumetric fraction of liquid water (-)
                              mLayerVolFracIceTrial(iLayer),          & ! intent(out): trial volumetric fraction if ice (-)
                              err,cmessage)                             ! intent(out): error control
              if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

              ! check
              case default; err=20; message=trim(message)//'expect case to be iname_veg, iname_snow, iname_soil'; return

          end select  ! domain type

        ! final check
        else

          ! do nothing (input = output) -- and check that we got here correctly
          if( (isNrgState .or. isCoupled) )then
            scalarVolFracLiq = realMissing
            scalarVolFracIce = realMissing
          else
            message=trim(message)//'unexpected else branch'
            err=20; return
          endif

        endif  ! if energy state or solution is coupled

        ! compute water time derivatives
        select case(ixDomainType)
          case(iname_veg); scalarCanopydWat_dt = scalarCanopyWatTrial - scalarCanopyWat
          case(iname_snow, iname_soil); mLayerdWat_dt(iLayer) = mLayerVolFracWatTrial(iLayer) - mLayerVolFracWat(iLayer)
        end select

        ! -----
        ! - update temperatures...
        ! ------------------------

        ! check the need to adjust temperature
        if(do_adjustTemp)then

          ! get the melt energy
          meltNrg = merge(LH_fus*iden_ice, LH_fus*iden_water, ixDomainType==iname_snow)

          ! compute the residual and the derivative
          select case(ixDomainType)

            ! * vegetation
            case(iname_veg)
              call xTempSolve(&
                              ! constant over iterations
                              meltNrg         = meltNrg                                 ,&  ! intent(in):    energy for melt+freeze (J m-3)
                              heatCap         = scalarBulkVolHeatCapVeg                 ,&  ! intent(in):    volumetric heat capacity (J m-3 K-1)
                              tempInit        = scalarCanopyTemp                        ,&  ! intent(in):    initial temperature (K)
                              volFracIceInit  = scalarCanopyIce/(iden_water*canopyDepth),&  ! intent(in):    initial volumetric fraction of ice (-)
                              ! trial values
                              xTemp           = xTemp                                   ,&  ! intent(inout): trial value of temperature
                              dLiq_dT         = dTheta_dTkCanopy                        ,&  ! intent(in):    derivative in liquid water content w.r.t. temperature (K-1)
                              volFracIceTrial = scalarVolFracIce                        ,&  ! intent(in):    trial value for volumetric fraction of ice
                              ! residual and derivative
                              residual        = residual                                ,&  ! intent(out):   residual (J m-3)
                              derivative      = derivative                               )  ! intent(out):   derivative (J m-3 K-1)

                    ! * snow and soil
            case(iname_snow, iname_soil)
              call xTempSolve(&
                              ! constant over iterations
                              meltNrg         = meltNrg                        ,&  ! intent(in):    energy for melt+freeze (J m-3)
                              heatCap         = mLayerVolHtCapBulk(iLayer)     ,&  ! intent(in):    volumetric heat capacity (J m-3 K-1)
                              tempInit        = mLayerTemp(iLayer)             ,&  ! intent(in):    initial temperature (K)
                              volFracIceInit  = mLayerVolFracIce(iLayer)       ,&  ! intent(in):    initial volumetric fraction of ice (-)
                              ! trial values
                              xTemp           = xTemp                          ,&  ! intent(inout): trial value of temperature
                              dLiq_dT         = mLayerdTheta_dTk(iLayer)       ,&  ! intent(in):    derivative in liquid water content w.r.t. temperature (K-1)
                              volFracIceTrial = mLayerVolFracIceTrial(iLayer)  ,&  ! intent(in):    trial value for volumetric fraction of ice
                              ! residual and derivative
                              residual        = residual                       ,&  ! intent(out):   residual (J m-3)
                              derivative      = derivative                      )  ! intent(out):   derivative (J m-3 K-1)

            ! * check
            case default; err=20; message=trim(message)//'expect case to be iname_veg, iname_snow, iname_soil'; return

          end select  ! domain type

          ! check validity of residual
          if( ieee_is_nan(residual) )then
            message=trim(message)//'residual is not valid'
            err=20; return
          endif

          ! update bracket
          if(residual < 0._rkind)then
            tempMax = min(xTemp,tempMax)
          else
            tempMin = max(tempMin,xTemp)
          end if

          ! compute iteration increment
          tempInc    = residual/derivative  ! K

          ! check
          if(globalPrintFlag)&
          write(*,'(i4,1x,e20.10,1x,5(f20.10,1x),L1)') iter, residual, xTemp-Tcrit, tempInc, Tcrit, tempMin, tempMax, bFlag

          ! check convergence
          if(abs(residual) < nrgConvTol .or. abs(tempInc) < tempConvTol) exit iterations

          ! add constraints for snow temperature
          if(ixDomainType==iname_veg .or. ixDomainType==iname_snow)then
            if(tempInc > Tcrit - xTemp) tempInc=(Tcrit - xTemp)*0.5_rkind  ! simple bi-section method
          endif  ! if the domain is vegetation or snow

          ! deal with the discontinuity between partially frozen and unfrozen soil
          if(ixDomainType==iname_soil)then
            ! difference from the temperature below which ice exists
            critDiff = Tcrit - xTemp
            ! --> initially frozen (T < Tcrit)
            if(critDiff > 0._rkind)then
              if(tempInc > critDiff) tempInc = critDiff + epsT  ! set iteration increment to slightly above critical temperature
            ! --> initially unfrozen (T > Tcrit)
            else
              if(tempInc < critDiff) tempInc = critDiff - epsT  ! set iteration increment to slightly below critical temperature
            endif
          endif  ! if the domain is soil

          ! update the temperature trial
          xTemp = xTemp + tempInc

          ! check failed convergence
          if(iter==maxiter)then
            message=trim(message)//'failed to converge'
            err=-20; return ! negative error code = try to recover
          endif

        endif   ! if adjusting the temperature

      end do iterations ! iterating

      ! save temperature
      select case(ixDomainType)
        case(iname_veg);              scalarCanopyTempTrial   = xTemp
        case(iname_snow, iname_soil); mLayerTempTrial(iLayer) = xTemp
      end select

      ! calculate temperature component of enthalpy for remaining domains
      if(ixDomainType==iname_veg)then
        if(computeEnthTemp)then
          call T2enthTemp_veg(&
                      canopyDepth,                 & ! intent(in): canopy depth (m)
                      specificHeatVeg,             & ! intent(in): specific heat of vegetation (J kg-1 K-1)
                      maxMassVegetation,           & ! intent(in): maximum mass of vegetation (kg m-2)
                      snowfrz_scale,               & ! intent(in): scaling parameter for the snow freezing curve  (K-1)
                      scalarCanopyTempTrial,       & ! intent(in): canopy temperature (K)
                      scalarCanopyWatTrial,        & ! intent(in): canopy water content (kg m-2)
                      scalarCanopyEnthTempTrial)     ! intent(out): temperature component of enthalpy of the vegetation canopy (J m-3)
        else
          scalarCanopyEnthTempTrial = realMissing
        endif
      elseif(ixDomainType==iname_snow)then
        if(computeEnthTemp)then
          call T2enthTemp_snow(&
                      snowfrz_scale,                   & ! intent(in):  scaling parameter for the snow freezing curve  (K-1)
                      mLayerTempTrial(iLayer),         & ! intent(in):  layer temperature (K)
                      mLayerVolFracWatTrial(iLayer),   & ! intent(in):  volumetric total water content (-)
                      mLayerEnthTempTrial(iLayer))       ! intent(out): temperature component of enthalpy of each snow layer (J m-3)
        else
          mLayerEnthTempTrial(iLayer) = realMissing
        endif
      elseif(ixDomainType==iname_soil)then
        if(computeEnthTemp)then
          call T2enthTemp_soil(&
                      use_lookup,                            & ! intent(in):  flag to use the lookup table for soil enthalpy
                      soil_dens_intr(ixControlIndex),        & ! intent(in):  intrinsic soil density (kg m-3)
                      vGn_alpha(ixControlIndex),vGn_n(ixControlIndex),theta_sat(ixControlIndex),theta_res(ixControlIndex),vGn_m(ixControlIndex), & ! intent(in): soil parameters
                      ixControlIndex,                        & ! intent(in):  index of the control volume within the domain
                      lookup_data,                           & ! intent(in):  lookup table data structure
                      realMissing,                           & ! intent(in):  lower value of integral (not computed)
                      mLayerTempTrial(iLayer),               & ! intent(in):  layer temperature (K)
                      mLayerMatricHeadTrial(ixControlIndex), & ! intent(in):  matric head (m)
                      mLayerEnthTempTrial(iLayer),           & ! intent(out): temperature component of enthalpy soil layer (J m-3)
                      err,cmessage)                            ! intent(out): error control
          if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
        else
          mLayerEnthTempTrial(iLayer) = realMissing
        endif
      endif

      ! =======================================================================================================================================
      ! =======================================================================================================================================

      ! -----
      ! - compute the liquid water matric potential (and necessary derivatives)...
      ! -------------------------------------------------------------------------

      ! only for soil
      if(ixDomainType==iname_soil)then

        ! check liquid water (include tolerance)
        if(mLayerVolFracLiqTrial(iLayer) > theta_sat(ixControlIndex)+epsT )then
          message=trim(message)//'liquid water greater than porosity'
          print*,'---------------'
          print*,'porosity(theta_sat)=', theta_sat(ixControlIndex)
          print*,'liq water =',mLayerVolFracLiqTrial(iLayer)
          print*,'layer =',iLayer
          print*,'---------------'
          err=20; return
        endif

        ! case of hydrology state uncoupled with energy
        if(.not.isNrgState .and. .not.isCoupled)then

          ! derivatives relating liquid water matric potential to total water matric potential and temperature
          dPsiLiq_dPsi0(ixControlIndex) = 1._rkind  ! exact correspondence (psiLiq=psi0)
          dPsiLiq_dTemp(ixControlIndex) = 0._rkind  ! no relationship between liquid water matric potential and temperature

        ! case of energy state or coupled solution
        else

          ! compute the liquid matric potential (and the derivatives w.r.t. total matric potential and temperature)
          call liquidHead(&
                        ! input
                        mLayerMatricHeadTrial(ixControlIndex)     ,& ! intent(in):  total water matric potential (m)
                        mLayerVolFracLiqTrial(iLayer)             ,& ! intent(in):  volumetric fraction of liquid water (-)
                        mLayerVolFracIceTrial(iLayer)             ,& ! intent(in):  volumetric fraction of ice (-)
                        vGn_alpha(ixControlIndex),vGn_n(ixControlIndex),theta_sat(ixControlIndex),theta_res(ixControlIndex),vGn_m(ixControlIndex), & ! intent(in): soil parameters
                        dVolTot_dPsi0(ixControlIndex)             ,& ! intent(in):  derivative in the soil water characteristic (m-1)
                        mLayerdTheta_dTk(iLayer)                  ,& ! intent(in):  derivative in volumetric total water w.r.t. temperature (K-1)
                        ! output
                        mLayerMatricHeadLiqTrial(ixControlIndex)  ,& ! intent(out): liquid water matric potential (m)
                        dPsiLiq_dPsi0(ixControlIndex)             ,& ! intent(out): derivative in the liquid water matric potential w.r.t. the total water matric potential (-)
                        dPsiLiq_dTemp(ixControlIndex)             ,& ! intent(out): derivative in the liquid water matric potential w.r.t. temperature (m K-1)
                        err,cmessage)                                ! intent(out): error control
          if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

        endif  ! switch between hydrology and energy state

      endif  ! if domain is soil

    end do ! looping through state variables

    deallocate(computedCoupling,stat=err)        ! .true. if computed the coupling for a given state variable
    if(err/=0)then; message=trim(message)//'problem deallocating computedCoupling'; return; endif

  end associate

 end subroutine updateVars


! **********************************************************************************************************
! private subroutine xTempSolve: compute residual and derivative for temperature
! **********************************************************************************************************
subroutine xTempSolve(&
                      ! input: constant over iterations
                      meltNrg          ,&  ! intent(in):    energy for melt+freeze (J m-3)
                      heatCap          ,&  ! intent(in):    volumetric heat capacity (J m-3 K-1)
                      tempInit         ,&  ! intent(in):    initial temperature (K)
                      volFracIceInit   ,&  ! intent(in):    initial volumetric fraction of ice (-)
                      ! input-output: trial values
                      xTemp            ,&  ! intent(inout): trial value of temperature
                      dLiq_dT          ,&  ! intent(in):    derivative in liquid water content w.r.t. temperature (K-1)
                      volFracIceTrial  ,&  ! intent(in):    trial value for volumetric fraction of ice
                      ! output: residual and derivative
                      residual         ,&  ! intent(out):   residual (J m-3)
                      derivative        )  ! intent(out):   derivative (J m-3 K-1)
  implicit none
  ! input: constant over iterations
  real(rkind),intent(in)             :: meltNrg          ! energy for melt+freeze (J m-3)
  real(rkind),intent(in)             :: heatCap          ! volumetric heat capacity (J m-3 K-1)
  real(rkind),intent(in)             :: tempInit         ! initial temperature (K)
  real(rkind),intent(in)             :: volFracIceInit   ! initial volumetric fraction of ice (-)
  ! input-output: trial values
  real(rkind),intent(inout)          :: xTemp            ! trial value for temperature
  real(rkind),intent(in)             :: dLiq_dT          ! derivative in liquid water content w.r.t. temperature (K-1)
  real(rkind),intent(in)             :: volFracIceTrial  ! trial value for the volumetric fraction of ice (-)
  ! output: residual and derivative
  real(rkind),intent(out)            :: residual         ! residual (J m-3)
  real(rkind),intent(out)            :: derivative       ! derivative (J m-3 K-1)
  ! subroutine starts here
  residual   = -heatCap*(xTemp - tempInit) + meltNrg*(volFracIceTrial - volFracIceInit)  ! J m-3
  derivative = heatCap + LH_fus*iden_water*dLiq_dT  ! J m-3 K-1
end subroutine xTempSolve

! **********************************************************************************************************
! public subroutine updateProg: update prognostic variables
! **********************************************************************************************************
subroutine updateProg(dt,nSnow,nSoil,nLayers,untappedMelt,stateVecTrial,stateVecPrime,                                           & ! input: states
                      doAdjustTemp,computeVegFlux,computMassBalance,computNrgBalance,computeEnthTemp,enthalpyStateVec,use_lookup,& ! input: model control
                      model_decisions,lookup_data,mpar_data,indx_data,flux_data,prog_data,diag_data,deriv_data,                  & ! input-output: data structures
                      fluxVec,resVec,balance,waterBalanceError,nrgFluxModified,err,message)                                        ! input-output: balances, flags, and error control
USE getVectorz_module,only:varExtract                              ! extract variables from the state vector
#ifdef SUNDIALS_ACTIVE
  USE updateVarsWithPrime_module,only:updateVarsWithPrime          ! update prognostic variables
#endif
  !USE updateVars_module,only:updateVars                            ! update prognostic variables
  USE enthalpyTemp_module,only:enthTemp_or_enthalpy                ! add phase change terms to delta temperature component of enthalpy
  implicit none
  ! model control
  real(rkind)      ,intent(in)    :: dt                            ! time step (s)
  integer(i4b)     ,intent(in)    :: nSnow                         ! number of snow layers
  integer(i4b)     ,intent(in)    :: nSoil                         ! number of soil layers
  integer(i4b)     ,intent(in)    :: nLayers                       ! total number of layers
  logical(lgt)     ,intent(in)    :: doAdjustTemp                  ! flag to indicate if we adjust the temperature
  logical(lgt)     ,intent(in)    :: computeVegFlux                ! flag to compute the vegetation flux
  real(rkind)      ,intent(in)    :: untappedMelt(:)               ! un-tapped melt energy (J m-3 s-1)
  real(rkind)      ,intent(in)    :: stateVecTrial(:)              ! trial state vector (mixed units)
  real(rkind)      ,intent(in)    :: stateVecPrime(:)              ! trial state vector (mixed units)
  logical(lgt)     ,intent(in)    :: computMassBalance             ! flag to check the mass balance
  logical(lgt)     ,intent(in)    :: computNrgBalance              ! flag to check the energy balance
  logical(lgt)     ,intent(in)    :: computeEnthTemp               ! flag to compute enthalpy
  logical(lgt)     ,intent(in)    :: enthalpyStateVec              ! flag if enthalpy is a state variable (ida)
  logical(lgt)     ,intent(in)    :: use_lookup                    ! flag to use the lookup table for soil enthalpy, otherwise use analytical solution
  ! data structures
  type(model_options),intent(in)  :: model_decisions(:)            ! model decisions
  type(zLookup),intent(in)        :: lookup_data                   ! lookup tables
  type(var_dlength),intent(in)    :: mpar_data                     ! model parameters
  type(var_ilength),intent(in)    :: indx_data                     ! indices for a local HRU
  type(var_dlength),intent(inout) :: flux_data                     ! model fluxes for a local HRU
  type(var_dlength),intent(inout) :: prog_data                     ! prognostic variables for a local HRU
  type(var_dlength),intent(inout) :: diag_data                     ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout) :: deriv_data                    ! derivatives in model fluxes w.r.t. relevant state variables
  ! balances, flags, and error control
  real(rkind)      ,intent(in)    :: fluxVec(:)                    ! flux vector (mixed units)
  real(qp)         ,intent(in)    :: resVec(:)    ! NOTE: qp       ! residual vector
  real(rkind)      ,intent(inout) :: balance(:)                    ! balance of energy per domain per second
  logical(lgt)     ,intent(out)   :: waterBalanceError             ! flag to denote that there is a water balance error
  logical(lgt)     ,intent(out)   :: nrgFluxModified               ! flag to denote that the energy fluxes were modified
  integer(i4b)     ,intent(out)   :: err                           ! error code
  character(*)     ,intent(out)   :: message                       ! error message
  ! ==================================================================================================================
  ! general
  integer(i4b)                    :: i                             ! indices
  integer(i4b)                    :: iState                        ! index of model state variable
  integer(i4b)                    :: ixSubset                      ! index within the state subset
  integer(i4b)                    :: ixFullVector                  ! index within full state vector
  integer(i4b)                    :: ixControlIndex                ! index within a given domain
  real(rkind)                     :: volMelt                       ! volumetric melt (kg m-3)
  real(rkind),parameter           :: eps=epsilon(1._rkind)         ! a very small number (deal with precision issues)
  real(rkind)                     :: eps_veg                       ! precision needs to vary based on set canopy water tolerance for IDA
  real(rkind)                     :: eps_snow                      ! precision needs to vary based on set snow water tolerance for IDA
  ! mass balance
  real(rkind)                     :: canopyBalance0,canopyBalance1 ! canopy storage at start/end of time step
  real(rkind)                     :: soilBalance0,soilBalance1     ! soil storage at start/end of time step
  real(rkind)                     :: vertFlux                      ! change in storage due to vertical fluxes
  real(rkind)                     :: tranSink,baseSink,compSink    ! change in storage due to sink terms
  real(rkind)                     :: liqError                      ! water balance error
  real(rkind)                     :: fluxNet                       ! net water fluxes (kg m-2 s-1)
  real(rkind)                     :: superflousWat                 ! superflous water used for evaporation (kg m-2 s-1)
  real(rkind)                     :: superflousNrg                 ! superflous energy that cannot be used for evaporation (W m-2 [J m-2 s-1])
  character(LEN=256)              :: cmessage                      ! error message of downwind routine
  logical(lgt),parameter          :: printFlag=.false.             ! flag to print water balance error information
  ! trial state variables
  real(rkind)                     :: scalarCanairTempTrial         ! trial value for temperature of the canopy air space (K)
  real(rkind)                     :: scalarCanopyTempTrial         ! trial value for temperature of the vegetation canopy (K)
  real(rkind)                     :: scalarCanopyWatTrial          ! trial value for liquid water storage in the canopy (kg m-2)
  real(rkind),dimension(nLayers)  :: mLayerTempTrial               ! trial vector of temperature of layers in the snow and soil domains (K)
  real(rkind),dimension(nLayers)  :: mLayerVolFracWatTrial         ! trial vector of volumetric fraction of total water (-)
  real(rkind),dimension(nSoil)    :: mLayerMatricHeadTrial         ! trial vector of total water matric potential (m)
  real(rkind),dimension(nSoil)    :: mLayerMatricHeadLiqTrial      ! trial vector of liquid water matric potential (m)
  real(rkind)                     :: scalarAquiferStorageTrial     ! trial value for storage of water in the aquifer (m)
  real(rkind)                     :: scalarCanairEnthalpyTrial     ! trial value for enthalpy of the canopy air space (J m-3)
  real(rkind)                     :: scalarCanopyEnthTempTrial     ! trial value for temperature component of enthalpy of the vegetation canopy (J m-3)
  real(rkind),dimension(nLayers)  :: mLayerEnthTempTrial           ! trial vector of temperature component of enthalpy of snow + soil (J m-3)
  real(rkind)                     :: scalarCanopyEnthalpyTrial     ! trial value for enthalpy of the vegetation canopy (J m-3)
  real(rkind),dimension(nLayers)  :: mLayerEnthalpyTrial           ! trial vector of enthalpy of each snow and soil layer (J m-3)
  ! diagnostic variables
  real(rkind)                     :: scalarCanopyLiqTrial          ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
  real(rkind)                     :: scalarCanopyIceTrial          ! trial value for mass of ice on the vegetation canopy (kg m-2)
  real(rkind),dimension(nLayers)  :: mLayerVolFracLiqTrial         ! trial vector of volumetric fraction of liquid water (-)
  real(rkind),dimension(nLayers)  :: mLayerVolFracIceTrial         ! trial vector of volumetric fraction of ice (-)
  ! prime state variables
  real(rkind)                     :: scalarCanairTempPrime         ! trial value for temperature of the canopy air space (K)
  real(rkind)                     :: scalarCanopyTempPrime         ! trial value for temperature of the vegetation canopy (K)
  real(rkind)                     :: scalarCanopyWatPrime          ! trial value for liquid water storage in the canopy (kg m-2)
  real(rkind),dimension(nLayers)  :: mLayerTempPrime               ! trial vector of temperature of layers in the snow and soil domains (K)
  real(rkind),dimension(nLayers)  :: mLayerVolFracWatPrime         ! trial vector of volumetric fraction of total water (-)
  real(rkind),dimension(nSoil)    :: mLayerMatricHeadPrime         ! trial vector of total water matric potential (m)
  real(rkind),dimension(nSoil)    :: mLayerMatricHeadLiqPrime      ! trial vector of liquid water matric potential (m)
  real(rkind)                     :: scalarAquiferStoragePrime     ! trial value for storage of water in the aquifer (m)
  ! diagnostic prime or delta variables
  real(rkind)                     :: scalarCanopyLiqPrime          ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
  real(rkind)                     :: scalarCanopyIcePrime          ! trial value for mass of ice on the vegetation canopy (kg m-2)
  real(rkind)                     :: scalarCanopyIceDelta          ! delta value for mass of ice on the vegetation canopy (kg m-2)
  real(rkind)                     :: scalarCanopyHDelta            ! delta value for enthalpy of the vegetation canopy (J m-3)
  real(rkind),dimension(nLayers)  :: mLayerVolFracLiqPrime         ! trial vector of volumetric fraction of liquid water (-)
  real(rkind),dimension(nLayers)  :: mLayerVolFracIcePrime         ! trial vector of volumetric fraction of ice (-)
  real(rkind),dimension(nLayers)  :: mLayerVolFracIceDelta         ! delta vector volumetric fraction of ice of snow + soil (-)
  real(rkind),dimension(nLayers)  :: mLayerHDelta                  ! delta vector of enthalpy of snow+soil (J m-3)
  ! dummy state variables
  real(rkind)                     :: scalarCanairNrgTrial        ! trial value for energy of the canopy air space
  real(rkind)                     :: scalarCanopyNrgTrial        ! trial value for energy of the vegetation canopy
  real(rkind),dimension(nLayers)  :: mLayerNrgTrial              ! trial vector of energy of each snow and soil layer
  real(rkind)                     :: scalarCanairNrgPrime        ! prime value for energy of the canopy air space
  real(rkind)                     :: scalarCanopyNrgPrime        ! prime value for energy of the vegetation canopy
  real(rkind),dimension(nLayers)  :: mLayerNrgPrime              ! prime vector of energy of each snow and soil layer
  ! -------------------------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------------------------
  ! point to flux variables in the data structure
  associate(&
    ! model decisions
    ixNumericalMethod         => model_decisions(iLookDECISIONS%num_method)%iDecision       ,& ! intent(in):  [i4b] choice of numerical solver
    ! get indices for balances
    ixCasNrg                  => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)                  ,& ! intent(in)   : [i4b]    index of canopy air space energy state variable
    ixVegNrg                  => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)                  ,& ! intent(in)   : [i4b]    index of canopy energy state variable
    ixVegHyd                  => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)                  ,& ! intent(in)   : [i4b]    index of canopy hydrology state variable (mass)
    !ixTopNrg                  => indx_data%var(iLookINDEX%ixTopNrg)%dat(1)                  ,& ! intent(in)   : [i4b]    index of upper-most energy state in the snow+soil subdomain
    !ixTopHyd                  => indx_data%var(iLookINDEX%ixTopHyd)%dat(1)                  ,& ! intent(in)   : [i4b]    index of upper-most hydrology state in the snow+soil subdomain
    ixAqWat                   => indx_data%var(iLookINDEX%ixAqWat)%dat(1)                   ,& ! intent(in)   : [i4b]    index of water storage in the aquifer
    ixSoilOnlyHyd             => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat                ,& ! intent(in)   : [i4b(:)] index in the state subset for hydrology state variables in the soil domain
    ixSnowSoilNrg             => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat                ,& ! intent(in)   : [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
    ixSnowSoilHyd             => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat                ,& ! intent(in)   : [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
    nSnowSoilNrg              => indx_data%var(iLookINDEX%nSnowSoilNrg)%dat(1)              ,& ! intent(in)   : [i4b]    number of energy state variables in the snow+soil domain
    nSnowSoilHyd              => indx_data%var(iLookINDEX%nSnowSoilHyd)%dat(1)              ,& ! intent(in)   : [i4b]    number of hydrology state variables in the snow+soil domain
    ! get indices for the un-tapped melt
    ixNrgOnly                 => indx_data%var(iLookINDEX%ixNrgOnly)%dat                    ,& ! intent(in)   : [i4b(:)] list of indices for all energy states
    ixDomainType              => indx_data%var(iLookINDEX%ixDomainType)%dat                 ,& ! intent(in)   : [i4b(:)] indices defining the domain of the state (iname_veg, iname_snow, iname_soil)
    ixControlVolume           => indx_data%var(iLookINDEX%ixControlVolume)%dat              ,& ! intent(in)   : [i4b(:)] index of the control volume for different domains (veg, snow, soil)
    ixMapSubset2Full          => indx_data%var(iLookINDEX%ixMapSubset2Full)%dat             ,& ! intent(in)   : [i4b(:)] [state subset] list of indices of the full state vector in the state subset
    ! water fluxes
    scalarRainfall            => flux_data%var(iLookFLUX%scalarRainfall)%dat(1)             ,& ! intent(in)   : [dp]     rainfall rate (kg m-2 s-1)
    scalarThroughfallRain     => flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1)      ,& ! intent(in)   : [dp]     rain reaches ground without touching the canopy (kg m-2 s-1)
    scalarCanopyEvaporation   => flux_data%var(iLookFLUX%scalarCanopyEvaporation)%dat(1)    ,& ! intent(in)   : [dp]     canopy evaporation/condensation (kg m-2 s-1)
    scalarCanopyLiqDrainage   => flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1)    ,& ! intent(in)   : [dp]     drainage liquid water from vegetation canopy (kg m-2 s-1)
    iLayerLiqFluxSoil         => flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat             ,& ! intent(in)   : [dp(0:)] vertical liquid water flux at soil layer interfaces (-)
    iLayerNrgFlux             => flux_data%var(iLookFLUX%iLayerNrgFlux)%dat                 ,& ! intent(in)   :
    mLayerNrgFlux             => flux_data%var(iLookFLUX%mLayerNrgFlux)%dat                 ,& ! intent(out)  : [dp]     net energy flux for each layer within the snow+soil domain (J m-3 s-1)
    mLayerTranspire           => flux_data%var(iLookFLUX%mLayerTranspire)%dat               ,& ! intent(in)   : [dp(:)]  transpiration loss from each soil layer (m s-1)
    mLayerBaseflow            => flux_data%var(iLookFLUX%mLayerBaseflow)%dat                ,& ! intent(in)   : [dp(:)]  baseflow from each soil layer (m s-1)
    mLayerCompress            => diag_data%var(iLookDIAG%mLayerCompress)%dat                ,& ! intent(in)   : [dp(:)]  change in storage associated with compression of the soil matrix (-)
    ! energy fluxes
    scalarLatHeatCanopyEvap   => flux_data%var(iLookFLUX%scalarLatHeatCanopyEvap)%dat(1)    ,& ! intent(in)   : [dp]     latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
    scalarSenHeatCanopy       => flux_data%var(iLookFLUX%scalarSenHeatCanopy)%dat(1)        ,& ! intent(in)   : [dp]     sensible heat flux from the canopy to the canopy air space (W m-2)
    ! domain depth
    canopyDepth               => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)          ,& ! intent(in)   : [dp   ]  canopy depth (m)
    mLayerDepth               => prog_data%var(iLookPROG%mLayerDepth)%dat                   ,& ! intent(in)   : [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
    ! model state variables (vegetation canopy)
    scalarCanairTemp          => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)           ,& ! intent(inout): [dp]     temperature of the canopy air space (K)
    scalarCanopyTemp          => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)           ,& ! intent(inout): [dp]     temperature of the vegetation canopy (K)
    scalarCanopyIce           => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)            ,& ! intent(inout): [dp]     mass of ice on the vegetation canopy (kg m-2)
    scalarCanopyLiq           => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)            ,& ! intent(inout): [dp]     mass of liquid water on the vegetation canopy (kg m-2)
    scalarCanopyWat           => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)            ,& ! intent(inout): [dp]     mass of total water on the vegetation canopy (kg m-2)
    ! model state variables (snow and soil domains)
    mLayerTemp                => prog_data%var(iLookPROG%mLayerTemp)%dat                    ,& ! intent(inout): [dp(:)]  temperature of each snow/soil layer (K)
    mLayerVolFracIce          => prog_data%var(iLookPROG%mLayerVolFracIce)%dat              ,& ! intent(inout): [dp(:)]  volumetric fraction of ice (-)
    mLayerVolFracLiq          => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat              ,& ! intent(inout): [dp(:)]  volumetric fraction of liquid water (-)
    mLayerVolFracWat          => prog_data%var(iLookPROG%mLayerVolFracWat)%dat              ,& ! intent(inout): [dp(:)]  volumetric fraction of total water (-)
    mLayerMatricHead          => prog_data%var(iLookPROG%mLayerMatricHead)%dat              ,& ! intent(inout): [dp(:)]  matric head (m)
    mLayerMatricHeadLiq       => diag_data%var(iLookDIAG%mLayerMatricHeadLiq)%dat           ,& ! intent(inout): [dp(:)]  matric potential of liquid water (m)
    ! enthalpy
    scalarCanairEnthalpy      => prog_data%var(iLookPROG%scalarCanairEnthalpy)%dat(1)       ,& ! intent(inout): [dp]     enthalpy of the canopy air space (J m-3)
    scalarCanopyEnthalpy      => prog_data%var(iLookPROG%scalarCanopyEnthalpy)%dat(1)       ,& ! intent(inout): [dp]     enthalpy of the vegetation canopy (J m-3)
    scalarCanopyEnthTemp      => diag_data%var(iLookDIAG%scalarCanopyEnthTemp)%dat(1)       ,& ! intent(inout): [dp]     temperature component of enthalpy of the vegetation canopy (J m-3)
    mLayerEnthalpy            => prog_data%var(iLookPROG%mLayerEnthalpy)%dat                ,& ! intent(inout): [dp(:)]  enthalpy of the snow+soil layers (J m-3)
    mLayerEnthTemp            => diag_data%var(iLookDIAG%mLayerEnthTemp)%dat                ,& ! intent(inout): [dp(:)]  temperature component of enthalpy of the snow+soil layers (J m-3)
    ! model state variables (aquifer)
    scalarAquiferStorage      => prog_data%var(iLookPROG%scalarAquiferStorage)%dat(1)       ,& ! intent(inout): [dp(:)]  storage of water in the aquifer (m)
    ! error tolerance
    absConvTol_liquid         => mpar_data%var(iLookPARAM%absConvTol_liquid)%dat(1)          & ! intent(in)   : [dp]     absolute convergence tolerance for vol frac liq water (-)
    ) ! associating flux variables in the data structure
    ! -------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message='updateProg/'

    ! initialize flags for water balance error and energy flux modification
    waterBalanceError=.false.
    nrgFluxModified = .false.

    ! get storage at the start of the step
    canopyBalance0 = merge(scalarCanopyLiq + scalarCanopyIce, realMissing, computeVegFlux)
    soilBalance0   = sum( (mLayerVolFracLiq(nSnow+1:nLayers) + mLayerVolFracIce(nSnow+1:nLayers)  )*mLayerDepth(nSnow+1:nLayers) )

    ! -----
    ! * update states...
    ! ------------------

    ! initialize to state variable from the last update
    scalarCanairTempTrial     = scalarCanairTemp
    scalarCanairEnthalpyTrial = scalarCanairEnthalpy
    scalarCanopyTempTrial     = scalarCanopyTemp
    scalarCanopyEnthalpyTrial = scalarCanopyEnthalpy
    scalarCanopyEnthTempTrial = scalarCanopyEnthTemp
    scalarCanopyWatTrial      = scalarCanopyWat
    scalarCanopyLiqTrial      = scalarCanopyLiq
    scalarCanopyIceTrial      = scalarCanopyIce
    mLayerTempTrial           = mLayerTemp
    mLayerEnthalpyTrial       = mLayerEnthalpy
    mLayerEnthTempTrial       = mLayerEnthTemp
    mLayerVolFracWatTrial     = mLayerVolFracWat
    mLayerVolFracLiqTrial     = mLayerVolFracLiq
    mLayerVolFracIceTrial     = mLayerVolFracIce
    mLayerMatricHeadTrial     = mLayerMatricHead
    mLayerMatricHeadLiqTrial  = mLayerMatricHeadLiq
    scalarAquiferStorageTrial = scalarAquiferStorage

    if(enthalpyStateVec)then ! use state variable as enthalpy
      scalarCanairNrgTrial = scalarCanairEnthalpy
      scalarCanopyNrgTrial = realMissing ! currently not splitting in ida so no need to update
      mLayerNrgTrial       = realMissing ! currently not splitting in ida so no need to update
    else
      scalarCanairNrgTrial = scalarCanairTemp
      scalarCanopyNrgTrial = scalarCanopyTemp
      mLayerNrgTrial       = mLayerTemp
    endif
      
    ! extract states from the state vector
    call varExtract(&
                    ! input
                    stateVecTrial,             & ! intent(in):    model state vector (mixed units)
                    indx_data,                 & ! intent(in):    indices defining model states and layers
                    ! output: variables for the vegetation canopy
                    scalarCanairNrgTrial,      & ! intent(inout): trial value of energy of the canopy air space, temperature (K) or enthalpy (J m-3)
                    scalarCanopyNrgTrial,      & ! intent(inout): trial value of energy of the vegetation canopy, temperature (K) or enthalpy (J m-3)
                    scalarCanopyWatTrial,      & ! intent(inout): trial value of canopy total water (kg m-2)
                    scalarCanopyLiqTrial,      & ! intent(inout): trial value of canopy liquid water (kg m-2)
                    ! output: variables for the snow-soil domain
                    mLayerNrgTrial,            & ! intent(inout): trial vector of energy, temperature (K) or enthalpy (J m-3)
                    mLayerVolFracWatTrial,     & ! intent(inout): trial vector of volumetric total water content (-)
                    mLayerVolFracLiqTrial,     & ! intent(inout): trial vector of volumetric liquid water content (-)
                    mLayerMatricHeadTrial,     & ! intent(inout): trial vector of total water matric potential (m)
                    mLayerMatricHeadLiqTrial,  & ! intent(inout): trial vector of liquid water matric potential (m)
                    ! output: variables for the aquifer
                    scalarAquiferStorageTrial, & ! intent(inout): trial value of storage of water in the aquifer (m)
                    ! output: error control
                    err,cmessage)               ! intent(out):   error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)
  
    if(enthalpyStateVec)then ! use state variable as enthalpy
      scalarCanairEnthalpyTrial = scalarCanairNrgTrial
      scalarCanopyEnthalpyTrial = scalarCanopyNrgTrial
      mLayerEnthalpyTrial       = mLayerNrgTrial
    else
      scalarCanairTempTrial = scalarCanairNrgTrial
      scalarCanopyTempTrial = scalarCanopyNrgTrial
      mLayerTempTrial       = mLayerNrgTrial
    endif

    ! Placeholder: if we decide to use splitting, we need to pass all the previous values of the state variables
    scalarCanairNrgPrime      = realMissing
    scalarCanopyNrgPrime      = realMissing
    scalarCanopyWatPrime      = realMissing
    scalarCanopyLiqPrime      = realMissing
    scalarCanopyIcePrime      = realMissing
    mLayerNrgPrime            = realMissing
    mLayerVolFracWatPrime     = realMissing
    mLayerVolFracLiqPrime     = realMissing
    mLayerVolFracIcePrime     = realMissing
    mLayerMatricHeadPrime     = realMissing
    mLayerMatricHeadLiqPrime  = realMissing
    scalarAquiferStoragePrime = realMissing

    ! set the default precision
    eps_veg  = eps*2._rkind
    eps_snow = eps*2._rkind

    select case(ixNumericalMethod)
      case(ida)
#ifdef SUNDIALS_ACTIVE
        ! IDA precision needs to vary based on set tolerances
        eps_veg = mpar_data%var(iLookPARAM%absTolWatVeg)%dat(1)*2._rkind
        eps_snow = mpar_data%var(iLookPARAM%absTolWatSnow)%dat(1)*2._rkind

        ! extract the derivatives from the state vector
        call varExtract(&
                  ! input
                  stateVecPrime,             & ! intent(in):    derivative of model state vector (mixed units)
                  indx_data,                 & ! intent(in):    indices defining model states and layers
                  ! output: variables for the vegetation canopy
                  scalarCanairNrgPrime,      & ! intent(inout): derivative of energy of the canopy air space, temperature (K s-1) or enthalpy (W m-3)
                  scalarCanopyNrgPrime,      & ! intent(inout): derivative of energy of the vegetation canopy, temperature (K s-1) or enthalpy (W m-3)
                  scalarCanopyWatPrime,      & ! intent(inout): derivative of canopy total water (kg m-2 s-1)
                  scalarCanopyLiqPrime,      & ! intent(inout): derivative of canopy liquid water (kg m-2 s-1)
                  ! output: variables for the snow-soil domain
                  mLayerNrgPrime,            & ! intent(inout): derivative of energy of each snow and soil layer, temperature (K s-1) or enthalpy (W m-3)
                  mLayerVolFracWatPrime,     & ! intent(inout):   derivative of volumetric total water content (-)
                  mLayerVolFracLiqPrime,     & ! intent(inout):   derivative of volumetric liquid water content (-)
                  mLayerMatricHeadPrime,     & ! intent(inout):   derivative of total water matric potential (m)
                  mLayerMatricHeadLiqPrime,  & ! intent(inout):   derivative of liquid water matric potential (m)
                  ! output: variables for the aquifer
                  scalarAquiferStoragePrime, & ! intent(inout):   derivative of storage of water in the aquifer (m)
                  ! output: error control
                  err,cmessage)               ! intent(out):   error control
        if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

        if(enthalpyStateVec)then ! use state variable as enthalpy, need to compute temperature
          ! do not use these variables
          scalarCanairTempPrime = realMissing
          scalarCanopyTempPrime = realMissing
          mLayerTempPrime       = realMissing
        else ! use state variable as temperature
          scalarCanairTempPrime = scalarCanairNrgPrime
          scalarCanopyTempPrime = scalarCanopyNrgPrime
          mLayerTempPrime       = mLayerNrgPrime   
        endif !(choice of how conservation of energy is implemented)
    
        ! update diagnostic variables
        call updateVarsWithPrime(&
                    ! input
                    enthalpyStateVec,                 & ! intent(in):    flag if enthalpy is used as state variable
                    use_lookup,                       & ! intent(in):    flag to use the lookup table for soil enthalpy
                    .false.,                          & ! intent(in):    logical flag if computing for Jacobian update
                    doAdjustTemp,                     & ! intent(in):    logical flag to adjust temperature to account for the energy used in melt+freeze
                    mpar_data,                        & ! intent(in):    model parameters for a local HRU
                    indx_data,                        & ! intent(in):    indices defining model states and layers
                    prog_data,                        & ! intent(in):    model prognostic variables for a local HRU
                    diag_data,                        & ! intent(inout): model diagnostic variables for a local HRU
                    deriv_data,                       & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                    lookup_data,                      & ! intent(in):    lookup table data structure
                    ! input: enthalpy state variables  
                    scalarCanairEnthalpyTrial,        & ! intent(in):    trial value for enthalpy of the canopy air space (J m-3)
                    scalarCanopyEnthalpyTrial,        & ! intent(in):    trial value for enthalpy of the vegetation canopy (J m-3)
                    mLayerEnthalpyTrial,              & ! intent(in):    trial vector of enthalpy of each snow+soil layer (J m-3)                      
                    ! output: variables for the vegetation canopy
                    scalarCanairTempTrial,            & ! intent(inout): trial value of canopy air space temperature (K)
                    scalarCanopyTempTrial,            & ! intent(inout): trial value of canopy temperature (K)
                    scalarCanopyWatTrial,             & ! intent(inout): trial value of canopy total water (kg m-2)
                    scalarCanopyLiqTrial,             & ! intent(inout): trial value of canopy liquid water (kg m-2)
                    scalarCanopyIceTrial,             & ! intent(inout): trial value of canopy ice content (kg m-2)
                    scalarCanopyTempPrime,            & ! intent(inout): trial value of canopy temperature (K)
                    scalarCanopyWatPrime,             & ! intent(inout): trial value of canopy total water (kg m-2)
                    scalarCanopyLiqPrime,             & ! intent(inout): trial value of canopy liquid water (kg m-2)
                    scalarCanopyIcePrime,             & ! intent(inout): trial value of canopy ice content (kg m-2)
                    ! output: variables for the snow-soil domain
                    mLayerTempTrial,                  & ! intent(inout): trial vector of layer temperature (K)
                    mLayerVolFracWatTrial,            & ! intent(inout): trial vector of volumetric total water content (-)
                    mLayerVolFracLiqTrial,            & ! intent(inout): trial vector of volumetric liquid water content (-)
                    mLayerVolFracIceTrial,            & ! intent(inout): trial vector of volumetric ice water content (-)
                    mLayerMatricHeadTrial,            & ! intent(inout): trial vector of total water matric potential (m)
                    mLayerMatricHeadLiqTrial,         & ! intent(inout): trial vector of liquid water matric potential (m)
                    mLayerTempPrime,                  & ! intent(inout): Prime vector of layer temperature (K)
                    mLayerVolFracWatPrime,            & ! intent(inout): Prime vector of volumetric total water content (-)
                    mLayerVolFracLiqPrime,            & ! intent(inout): Prime vector of volumetric liquid water content (-)
                    mLayerVolFracIcePrime,            & ! intent(inout): Prime vector of volumetric ice water content (-)
                    mLayerMatricHeadPrime,            & ! intent(inout): Prime vector of total water matric potential (m)
                    mLayerMatricHeadLiqPrime,         & ! intent(inout): Prime vector of liquid water matric potential (m)
                    ! output: error control
                    err,cmessage)                       ! intent(out):   error control
#endif
      case(kinsol, homegrown)
        ! update diagnostic variables
        call updateVars(&
                 ! input
                 computeEnthTemp,           & ! intent(in):    flag if computing temperature component of enthalpy
                 use_lookup,                & ! intent(in):    flag to use the lookup table for soil enthalpy
                 doAdjustTemp,              & ! intent(in):    logical flag to adjust temperature to account for the energy used in melt+freeze
                 mpar_data,                 & ! intent(in):    model parameters for a local HRU
                 indx_data,                 & ! intent(in):    indices defining model states and layers
                 prog_data,                 & ! intent(in):    model prognostic variables for a local HRU
                 diag_data,                 & ! intent(inout): model diagnostic variables for a local HRU
                 deriv_data,                & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                 lookup_data,               & ! intent(in):    lookup table data structure
                 scalarCanairTempTrial,     & ! intent(in):    trial value of canopy air space temperature (K)
                 ! output: variables for the vegetation canopy
                 scalarCanopyTempTrial,     & ! intent(inout): trial value of canopy temperature (K)
                 scalarCanopyWatTrial,      & ! intent(inout): trial value of canopy total water (kg m-2)
                 scalarCanopyLiqTrial,      & ! intent(inout): trial value of canopy liquid water (kg m-2)
                 scalarCanopyIceTrial,      & ! intent(inout): trial value of canopy ice content (kg m-2)
                 ! output: variables for the snow-soil domain
                 mLayerTempTrial,           & ! intent(inout): trial vector of layer temperature (K)
                 mLayerVolFracWatTrial,     & ! intent(inout): trial vector of volumetric total water content (-)
                 mLayerVolFracLiqTrial,     & ! intent(inout): trial vector of volumetric liquid water content (-)
                 mLayerVolFracIceTrial,     & ! intent(inout): trial vector of volumetric ice water content (-)
                 mLayerMatricHeadTrial,     & ! intent(inout): trial vector of total water matric potential (m)
                 mLayerMatricHeadLiqTrial,  & ! intent(inout): trial vector of liquid water matric potential (m)
                 ! output: enthalpy state variables  
                 scalarCanairEnthalpyTrial, & ! intent(inout): trial value for enthalpy of the canopy air space (J m-3)
                 scalarCanopyEnthTempTrial, & ! intent(inout): trial value for temperature component of enthalpy of the vegetation canopy (J m-3)
                 mLayerEnthTempTrial,       & ! intent(inout): trial vector of temperature component of enthalpy of each snow+soil layer (J m-3)                     
                 ! output: error control
                 err,cmessage)                ! intent(out):   error control

    end select
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    if(computNrgBalance)then
      ! compute energy balance if didn't do inside solver substeps
      select case(ixNumericalMethod)
        case(ida); ! do nothing, already computed
        case(kinsol, homegrown)
          ! calculate delta ice
          scalarCanopyIceDelta  = scalarCanopyIceTrial - scalarCanopyIce
          mLayerVolFracIceDelta = mLayerVolFracIceTrial - mLayerVolFracIce(1:nLayers)

          ! initialize delta enthalpy (HDelta) to delta temperature component of enthalpy, no difference in canopy air space
          scalarCanopyHDelta = scalarCanopyEnthTempTrial - scalarCanopyEnthTemp
          mLayerHDelta       = mLayerEnthTempTrial - mLayerEnthTemp(1:nLayers)
          
          ! compute mixture enthalpy for current values, do on delta value so only have to do once
          call enthTemp_or_enthalpy(&
                            ! input: data structures
                            .true.,                & ! intent(in):    flag to convert enthTemp to enthalpy
                            diag_data,             & ! intent(in):    model diagnostic variables for a local HRU
                            indx_data,             & ! intent(in):    model indices
                            ! input: ice content change
                            scalarCanopyIceDelta,  & ! intent(in):    delta value for canopy ice content (kg m-2)
                            mLayerVolFracIceDelta, & ! intent(in):    delta vector of volumetric ice water content (-)
                            ! input/output: enthalpy
                            scalarCanopyHDelta,    & ! intent(inout): delta value for enthalpy of the vegetation canopy (J m-3)
                            mLayerHDelta,          & ! intent(inout): delta vector of enthalpy of each snow+soil layer (J m-3)
                            ! output: error control    
                            err,cmessage)             ! intent(out): error control
          if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

          ! compute energy balance, maybe should use to check for step reduction
          if(ixCasNrg/=integerMissing) balance(ixCasNrg) = (scalarCanairEnthalpyTrial - scalarCanairEnthalpy)/dt - fluxVec(ixCasNrg)
          if(ixVegNrg/=integerMissing) balance(ixVegNrg) = scalarCanopyHDelta/dt - fluxVec(ixVegNrg)
          if(nSnowSoilNrg>0)then
            do concurrent (i=1:nLayers,ixSnowSoilNrg(i)/=integerMissing)
              balance(ixSnowSoilNrg(i)) = mLayerHDelta(i)/dt - fluxVec(ixSnowSoilNrg(i))
            enddo
          endif
          ! This is equivalent to above if, and only if, ixNrgConserv.ne.closedForm
          !!if(ixCasNrg/=integerMissing) balance(ixCasNrg) = resVec(ixCasNrg)/dt
          !if(ixVegNrg/=integerMissing) balance(ixVegNrg) = resVec(ixVegNrg)/dt
          !if(nSnowSoilNrg>0)then
          !  do concurrent (i=1:nLayers,ixSnowSoilNrg(i)/=integerMissing)
          !    balance(ixSnowSoilNrg(i)) = resVec(ixSnowSoilNrg(i))/dt
          !  enddo
          !endif

      end select
    else ! if not checking energy balance set balance to missing
      if(ixCasNrg/=integerMissing) balance(ixCasNrg) = realMissing
      if(ixVegNrg/=integerMissing) balance(ixVegNrg) = realMissing
      if(nSnowSoilNrg>0)then
        do concurrent (i=1:nLayers,ixSnowSoilNrg(i)/=integerMissing)
          balance(ixSnowSoilNrg(i)) = realMissing
        enddo
      endif
    endif  ! if checking energy balance

    ! -----
    ! * check mass balance...
    ! -----------------------

    ! NOTE: currently this will only fail with kinsol solver, since mass balance is checked in the homegrown solver and not checked for ida solver
    !   Negative error code will mean step will be failed and retried with smaller step size
    if(computMassBalance)then

      if(ixVegHyd/=integerMissing)then ! check for complete drainage

        ! handle cases where fluxes empty the canopy
        fluxNet = scalarRainfall + scalarCanopyEvaporation - scalarThroughfallRain - scalarCanopyLiqDrainage
        if(-fluxNet*dt > canopyBalance0)then

          ! --> first add water
          canopyBalance1 = canopyBalance0 + (scalarRainfall - scalarThroughfallRain)*dt

          ! --> next, remove canopy evaporation -- put the unsatisfied evap into sensible heat
          canopyBalance1 = canopyBalance1 + scalarCanopyEvaporation*dt
          if(canopyBalance1 < 0._rkind)then
            ! * get superfluous water and energy
            superflousWat = -canopyBalance1/dt     ! kg m-2 s-1
            superflousNrg = superflousWat*LH_vap   ! W m-2 (J m-2 s-1)
            ! * update fluxes and states
            canopyBalance1          = 0._rkind
            scalarCanopyEvaporation = scalarCanopyEvaporation + superflousWat
            scalarLatHeatCanopyEvap = scalarLatHeatCanopyEvap + superflousNrg
            scalarSenHeatCanopy     = scalarSenHeatCanopy - superflousNrg
          endif

          ! --> next, remove canopy drainage
          canopyBalance1 = canopyBalance1 -scalarCanopyLiqDrainage*dt
          if(canopyBalance1 < 0._rkind)then
            superflousWat           = -canopyBalance1/dt     ! kg m-2 s-1
            canopyBalance1          = 0._rkind
            scalarCanopyLiqDrainage = scalarCanopyLiqDrainage + superflousWat
          endif

          ! update the trial state
          scalarCanopyWatTrial = canopyBalance1

          ! set the modification flag
          nrgFluxModified = .true.

        else
          canopyBalance1  = canopyBalance0 + fluxNet*dt
          nrgFluxModified = .false.
        endif  ! cases where fluxes empty the canopy
      
      endif ! check for complete drainage

      ! compute mass balance if didn't do inside solver substeps
      select case(ixNumericalMethod)
        case(ida); ! do nothing
        case(kinsol, homegrown)
          ! old mass balance checks
          if(ixVegHyd/=integerMissing)then
            ! check the mass balance for the canopy for step reduction (ida and kinsol should have done this already unless modified canopy water above)
            fluxNet  = scalarRainfall + scalarCanopyEvaporation - scalarThroughfallRain - scalarCanopyLiqDrainage
            liqError = (canopyBalance0 + fluxNet*dt) - scalarCanopyWatTrial
            if(abs(liqError) > absConvTol_liquid*10._rkind)then  ! *10 because of precision issues
              if(printFlag)then
                write(*,'(a,1x,f20.10)') 'dt = ', dt
                write(*,'(a,1x,f20.10)') 'scalarCanopyWatTrial         = ', scalarCanopyWatTrial
                write(*,'(a,1x,f20.10)') 'canopyBalance0               = ', canopyBalance0
                write(*,'(a,1x,f20.10)') 'canopyBalance1               = ', canopyBalance1
                write(*,'(a,1x,f20.10)') 'scalarRainfall*dt            = ', scalarRainfall*dt
                write(*,'(a,1x,f20.10)') 'scalarCanopyLiqDrainage*dt   = ', scalarCanopyLiqDrainage*dt
                write(*,'(a,1x,f20.10)') 'scalarCanopyEvaporation*dt   = ', scalarCanopyEvaporation*dt
                write(*,'(a,1x,f20.10)') 'scalarThroughfallRain*dt     = ', scalarThroughfallRain*dt
                write(*,'(a,1x,f20.10)') 'liqError                     = ', liqError
              endif
              waterBalanceError = .true.
              return
            endif  ! if there is a water balance error
          endif  ! if veg canopy

          ! check mass balance for soil domain for step reduction (ida and kinsol should have done this already 
          if(count(ixSoilOnlyHyd/=integerMissing)==nSoil)then
            soilBalance1 = sum( (mLayerVolFracLiqTrial(nSnow+1:nLayers) + mLayerVolFracIceTrial(nSnow+1:nLayers) )*mLayerDepth(nSnow+1:nLayers) )
            vertFlux     = -(iLayerLiqFluxSoil(nSoil) - iLayerLiqFluxSoil(0))*dt           ! m s-1 --> m
            tranSink     = sum(mLayerTranspire)*dt                                         ! m s-1 --> m
            baseSink     = sum(mLayerBaseflow)*dt                                          ! m s-1 --> m
            compSink     = sum(mLayerCompress(1:nSoil) * mLayerDepth(nSnow+1:nLayers) )*dt ! m s-1 --> m
            liqError     = soilBalance1 - (soilBalance0 + vertFlux + tranSink - baseSink - compSink)
            if(abs(liqError) > absConvTol_liquid*10._rkind)then   ! *10 because of precision issues
              if(printFlag)then
                write(*,'(a,1x,f20.10)') 'dt = ', dt
                write(*,'(a,1x,f20.10)') 'soilBalance0      = ', soilBalance0
                write(*,'(a,1x,f20.10)') 'soilBalance1      = ', soilBalance1
                write(*,'(a,1x,f20.10)') 'vertFlux          = ', vertFlux
                write(*,'(a,1x,f20.10)') 'tranSink          = ', tranSink
                write(*,'(a,1x,f20.10)') 'baseSink          = ', baseSink
                write(*,'(a,1x,f20.10)') 'compSink          = ', compSink
                write(*,'(a,1x,f20.10)') 'liqError          = ', liqError
              endif
              waterBalanceError = .true.
              return
            endif  ! if there is a water balance error
          endif  ! if hydrology states exist in the soil domain

          ! compute mass balance, maybe should use to check for step reduction
          ! resVec is the residual vector from the solver over dt
          if(ixVegHyd/=integerMissing) balance(ixVegHyd) = resVec(ixVegHyd)/dt
          if(nSnowSoilHyd>0)then
            do concurrent (i=1:nLayers,ixSnowSoilHyd(i)/=integerMissing)
              balance(ixSnowSoilHyd(i)) = resVec(ixSnowSoilHyd(i))/dt
            end do
          endif
          if(ixAqWat/=integerMissing) balance(ixAqWat) = resVec(ixAqWat)/dt

      end select
    else ! if not checking mass balance set balance to missing
      if(ixVegHyd/=integerMissing) balance(ixVegHyd) = realMissing
      if(nSnowSoilHyd>0)then
        do concurrent (i=1:nLayers,ixSnowSoilHyd(i)/=integerMissing)
          balance(ixSnowSoilHyd(i)) = realMissing
        end do
      endif
      if(ixAqWat/=integerMissing) balance(ixAqWat) = realMissing
    endif  ! if checking the mass balance

    ! -----
    ! * remove untapped melt energy... always 0 at the moment but if use should be in solved as affects state
    ! --------------------------------

    ! only work with energy state variables
    if(size(ixNrgOnly)>0)then  ! energy state variables exist

      ! loop through energy state variables
      do iState=1,size(ixNrgOnly)

        ! get index of the control volume within the domain
        ixSubset       = ixNrgOnly(iState)             ! index within the state subset
        ixFullVector   = ixMapSubset2Full(ixSubset)    ! index within full state vector
        ixControlIndex = ixControlVolume(ixFullVector) ! index within a given domain

        ! compute volumetric melt (kg m-3)
        volMelt = dt*untappedMelt(ixSubset)/LH_fus  ! (kg m-3)

        ! update ice content
        select case( ixDomainType(ixFullVector) )
          case(iname_cas);  cycle ! do nothing, since there is no snow stored in the canopy air space
          case(iname_veg);  scalarCanopyIceTrial                        = scalarCanopyIceTrial                        - volMelt*canopyDepth  ! (kg m-2)
          case(iname_snow); mLayerVolFracIceTrial(ixControlIndex)       = mLayerVolFracIceTrial(ixControlIndex)       - volMelt/iden_ice     ! (-)
          case(iname_soil); mLayerVolFracIceTrial(ixControlIndex+nSnow) = mLayerVolFracIceTrial(ixControlIndex+nSnow) - volMelt/iden_water   ! (-)
          case default; err=20; message=trim(message)//'unable to identify domain type [remove untapped melt energy]'; return
        end select

        ! update liquid water content
        select case( ixDomainType(ixFullVector) )
          case(iname_cas);  cycle ! do nothing, since there is no snow stored in the canopy air space
          case(iname_veg);  scalarCanopyLiqTrial                        = scalarCanopyLiqTrial                        + volMelt*canopyDepth  ! (kg m-2)
          case(iname_snow); mLayerVolFracLiqTrial(ixControlIndex)       = mLayerVolFracLiqTrial(ixControlIndex)       + volMelt/iden_water   ! (-)
          case(iname_soil); mLayerVolFracLiqTrial(ixControlIndex+nSnow) = mLayerVolFracLiqTrial(ixControlIndex+nSnow) + volMelt/iden_water   ! (-)
          case default; err=20; message=trim(message)//'unable to identify domain type [remove untapped melt energy]'; return
        end select

      end do  ! looping through energy variables

      ! ========================================================================================================

      ! *** ice

      ! --> check if we removed too much water
      if(scalarCanopyIceTrial < 0._rkind  .or. any(mLayerVolFracIceTrial < 0._rkind) )then

        ! **
        ! canopy within numerical precision
        if(scalarCanopyIceTrial < 0._rkind)then

          if(scalarCanopyIceTrial > -eps_veg)then
            scalarCanopyLiqTrial = scalarCanopyLiqTrial - scalarCanopyIceTrial
            scalarCanopyIceTrial = 0._rkind

          ! encountered an inconsistency: spit the dummy
          else
            print*, 'dt = ', dt
            print*, 'untappedMelt          = ', untappedMelt
            print*, 'untappedMelt*dt       = ', untappedMelt*dt
            print*, 'scalarCanopyiceTrial  = ', scalarCanopyIceTrial
            message=trim(message)//'melted more than the available water'
            err=20; return
          endif  ! (inconsistency)

        endif  ! if checking the canopy
        ! **
        ! snow+soil within numerical precision
        do iState=1,size(mLayerVolFracIceTrial)

          ! snow layer within numerical precision
          if(mLayerVolFracIceTrial(iState) < 0._rkind)then

            if(mLayerVolFracIceTrial(iState) > -eps_snow)then
              mLayerVolFracLiqTrial(iState) = mLayerVolFracLiqTrial(iState) - mLayerVolFracIceTrial(iState)
              mLayerVolFracIceTrial(iState) = 0._rkind

            ! encountered an inconsistency: spit the dummy
            else
              print*, 'dt = ', dt
              print*, 'untappedMelt          = ', untappedMelt
              print*, 'untappedMelt*dt       = ', untappedMelt*dt
              print*, 'mLayerVolFracIceTrial = ', mLayerVolFracIceTrial
              message=trim(message)//'melted more than the available water'
              err=20; return
            endif  ! (inconsistency)

          endif  ! if checking a snow layer

        end do ! (looping through state variables)

      endif  ! (if we removed too much water)

      ! ========================================================================================================

      ! *** liquid water

      ! --> check if we removed too much water
      if(scalarCanopyLiqTrial < 0._rkind  .or. any(mLayerVolFracLiqTrial < 0._rkind) )then

        ! **
        ! canopy within numerical precision
        if(scalarCanopyLiqTrial < 0._rkind)then

          if(scalarCanopyLiqTrial > -eps_veg)then
            scalarCanopyIceTrial = scalarCanopyIceTrial - scalarCanopyLiqTrial
            scalarCanopyLiqTrial = 0._rkind

          ! encountered an inconsistency: spit the dummy
          else
            print*, 'dt = ', dt
            print*, 'untappedMelt          = ', untappedMelt
            print*, 'untappedMelt*dt       = ', untappedMelt*dt
            print*, 'scalarCanopyLiqTrial  = ', scalarCanopyLiqTrial
            message=trim(message)//'frozen more than the available water'
            err=20; return
          endif  ! (inconsistency)
        endif  ! checking the canopy

        ! **
        ! snow+soil within numerical precision
        do iState=1,size(mLayerVolFracLiqTrial)

          ! snow layer within numerical precision
          if(mLayerVolFracLiqTrial(iState) < 0._rkind)then

            if(mLayerVolFracLiqTrial(iState) > -eps_snow)then
              mLayerVolFracIceTrial(iState) = mLayerVolFracIceTrial(iState) - mLayerVolFracLiqTrial(iState)
              mLayerVolFracLiqTrial(iState) = 0._rkind

            ! encountered an inconsistency: spit the dummy
            else
              print*, 'dt = ', dt
              print*, 'untappedMelt          = ', untappedMelt
              print*, 'untappedMelt*dt       = ', untappedMelt*dt
              print*, 'mLayerVolFracLiqTrial = ', mLayerVolFracLiqTrial
              message=trim(message)//'frozen more than the available water'
              err=20; return
            endif  ! (inconsistency)

          endif  ! checking a snow layer

        end do ! (looping through state variables)

      endif  ! (if we removed too much water)

    endif  ! (if energy state variables exist)

    ! -----
    ! * update enthalpy as a diagnostic variable... 
    !   if computeEnthTemp then enthTemp will change, if enthalpyStateVec then enthalpy will change
    ! --------------------------------
    scalarCanairEnthalpy = scalarCanairEnthalpyTrial ! equivalent to scalarCanairEnthTemp
    scalarCanopyEnthTemp = scalarCanopyEnthTempTrial
    scalarCanopyEnthalpy = scalarCanopyEnthalpyTrial
    mLayerEnthTemp       = mLayerEnthTempTrial
    mLayerEnthalpy       = mLayerEnthalpyTrial

    ! -----
    ! * update prognostic variables...
    ! --------------------------------
    ! update state variables for the vegetation canopy
    scalarCanairTemp    = scalarCanairTempTrial    ! trial value of canopy air temperature (K)
    scalarCanopyTemp    = scalarCanopyTempTrial    ! trial value of canopy temperature (K)
    scalarCanopyWat     = scalarCanopyWatTrial     ! trial value of canopy total water (kg m-2)
    scalarCanopyLiq     = scalarCanopyLiqTrial     ! trial value of canopy liquid water (kg m-2)
    scalarCanopyIce     = scalarCanopyIceTrial     ! trial value of canopy ice content (kg m-2)

    ! update state variables for the snow+soil domain
    mLayerTemp          = mLayerTempTrial          ! trial vector of layer temperature (K)
    mLayerVolFracWat    = mLayerVolFracWatTrial    ! trial vector of volumetric total water content (-)
    mLayerVolFracLiq    = mLayerVolFracLiqTrial    ! trial vector of volumetric liquid water content (-)
    mLayerVolFracIce    = mLayerVolFracIceTrial    ! trial vector of volumetric ice water content (-)
    mLayerMatricHead    = mLayerMatricHeadTrial    ! trial vector of matric head (m)
    mLayerMatricHeadLiq = mLayerMatricHeadLiqTrial ! trial vector of matric head (m)

    ! update state variables for the aquifer
    scalarAquiferStorage = scalarAquiferStorageTrial

    ! end associations to info in the data structures
  end associate

end subroutine updateProg

end module updateVars_module
