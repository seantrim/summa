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

module enthalpyTemp_module

! constants
USE multiconst, only: gravity, &                          ! gravitational acceleration (m s-1)
                      Tfreeze, &                          ! freezing point of water (K)
                      Cp_soil,Cp_water,Cp_ice,Cp_air,&    ! specific heat of soil, water and ice (J kg-1 K-1)
                      iden_water,iden_ice,iden_air,&      ! intrinsic density of water and ice (kg m-3)
                      LH_fus                              ! latent heat of fusion (J kg-1)

! data types
USE nrtype
USE data_types,only:var_iLength                    ! var(:)%dat(:)
USE data_types,only:var_dLength                    ! var(:)%dat(:)
USE data_types,only:zLookup                        ! z(:)%var(:)%lookup(:)

! indices within parameter structure
USE var_lookup,only:iLookPARAM                     ! named variables to define structure element
USE var_lookup,only:iLookINDEX                     ! named variables to define structure element
USE var_lookup,only:iLookLOOKUP                    ! named variables to define structure element
USE var_lookup,only:iLookDIAG                      ! named variables for structure elements

! data dimensions
USE var_lookup,only:maxvarLookup                   ! maximum number of variables in the lookup tables

! domain types
USE globalData,only:iname_cas                      ! named variables for canopy air space
USE globalData,only:iname_veg                      ! named variables for vegetation canopy
USE globalData,only:iname_snow                     ! named variables for snow
USE globalData,only:iname_soil                     ! named variables for soil
USE globalData,only:iname_aquifer                  ! named variables for the aquifer

! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair                ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy                ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_nrgLayer                 ! named variable defining the energy state variable for snow+soil layers

! missing values
USE globalData,only:integerMissing                 ! missing integer
USE globalData,only:realMissing                    ! missing real number

! privacy

implicit none
public::T2H_lookup_snow
public::T2L_lookup_soil
public::H2T_snow
public::T2H_snow
public::T2enthTemp
public::enthTemp2H
private::hyp_2F1_real

! define the snow look-up table used to compute temperature based on enthalpy
integer(i4b),parameter               :: nlook=10001       ! number of elements in the lookup table
real(rkind),dimension(nlook),public  :: H_lookup          ! enthalpy values (J kg-1)
real(rkind),dimension(nlook),public  :: T_lookup          ! temperature values (K)
contains


! ************************************************************************************************************************
! public subroutine T2H_lookup_snow:: define a look-up table to mixture enthalpy based on temperature
!                                    appropriate when no dry mass, as in snow
! ************************************************************************************************************************
subroutine T2H_lookup_snow(mpar_data,                     &  ! intent(in):    parameter data structure
                           err,message)
  ! -------------------------------------------------------------------------------------------------------------------------
  ! downwind routines 
  USE nr_utility_module,only:arth                       ! use to build vectors with regular increments
  USE spline_int_module,only:spline,splint              ! use for cubic spline interpolation
  implicit none
  ! -------------------------------------------------------------------------------------------------------------------------
  ! declare dummy variables
  type(var_dlength),intent(in)  :: mpar_data            ! model parameters
  integer(i4b),intent(out)      :: err                  ! error code
  character(*),intent(out)      :: message              ! error message
  ! declare local variables
  character(len=128)            :: cmessage             ! error message in downwind routine
  real(rkind),parameter         :: T_start=260.0_rkind  ! start temperature value where all liquid water is assumed frozen (K)
  real(rkind)                   :: T_incr,H_incr        ! temperature/enthalpy increments
  real(rkind),dimension(nlook)  :: Tk                   ! initial temperature vector
  real(rkind),dimension(nlook)  :: Hy                   ! initial enthalpy vector
  real(rkind),parameter         :: waterWght=1._rkind   ! weight applied to total water (kg m-3) --- cancels out
  real(rkind),dimension(nlook)  :: H2                   ! 2nd derivatives of the interpolating function at tabulated points
  real(rkind)                   :: dT                   ! derivative of temperature with enthalpy at H_lookup
  integer(i4b)                  :: ilook                ! loop through lookup table
  ! -------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message="T2H_lookup_snow/"

  ! associate
  associate( snowfrz_scale => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1) )

    ! define initial temperature vector
    T_incr = (Tfreeze - T_start) / real(nlook-1, kind(rkind))  ! temperature increment
    Tk     = arth(T_start,T_incr,nlook)
    ! ***** compute specific enthalpy (NOTE: J m-3 --> J kg-1) *****

    do ilook=1,nlook
      Hy(ilook) = T2H_snow(Tk(ilook),waterWght,snowfrz_scale)/waterWght  ! (J m-3 --> J kg-1)
    end do

    ! define the final enthalpy vector
    H_incr   = (-Hy(1)) / real(nlook-1, kind(rkind))  ! enthalpy increment
    H_lookup = arth(Hy(1),H_incr,nlook)

    ! use cubic spline interpolation to obtain temperature values at the desired values of enthalpy
    call spline(Hy,Tk,1.e30_rkind,1.e30_rkind,H2,err,cmessage)  ! get the second derivatives
    if(err/=0) then; message=trim(message)//trim(cmessage); return; end if

    do ilook=1,nlook
      call splint(Hy,Tk,H2,H_lookup(ilook),T_lookup(ilook),dT,err,cmessage)
      if(err/=0) then; message=trim(message)//trim(cmessage); return; end if
    end do

  end associate

 end subroutine T2H_lookup_snow

! ************************************************************************************************************************
! public subroutine T2L_lookup_soil: define a look-up table to compute integral of soil Clapeyron equation liquid water
!                                    matric potential from temperature
! ************************************************************************************************************************
subroutine T2L_lookup_soil(nSoil,                         &  ! intent(in):    number of soil layers
                           mpar_data,                     &  ! intent(in):    parameter data structure
                           lookup_data,                   &  ! intent(inout): lookup table data structure
                           err,message)
  ! -------------------------------------------------------------------------------------------------------------------------
  ! downwind routines                    
  USE nr_utility_module,only:arth                       ! use to build vectors with regular increments
  USE spline_int_module,only:spline,splint              ! use for cubic spline interpolation
  USE soil_utils_module,only:volFracLiq                 ! use to compute the volumetric fraction of liquid water
  implicit none
  ! -------------------------------------------------------------------------------------------------------------------------
  ! declare dummy variables
  integer(i4b),intent(in)       :: nSoil
  type(var_dlength),intent(in)  :: mpar_data            ! model parameters
  type(zLookup),intent(inout)   :: lookup_data          ! lookup tables
  integer(i4b),intent(out)      :: err                  ! error code
  character(*),intent(out)      :: message              ! error message
  ! declare local variables
  character(len=128)            :: cmessage             ! error message in downwind routine
  integer(i4b),parameter        :: nLook=500            ! number of elements in the lookup table
  integer(i4b),parameter        :: nIntegr8=10000       ! number of points used in the numerical integration
  real(rkind),parameter         :: T_lower=260.0_rkind  ! lowest temperature value where all liquid water is assumed frozen (K)
  real(rkind),dimension(nLook)  :: xTemp                ! temporary vector
  real(rkind)                   :: xIncr                ! temporary increment
  real(rkind)                   :: T_incr               ! temperature increment
  real(rkind)                   :: dL                   ! derivative of integral with temperature at T_test
  integer(i4b)                  :: iVar                 ! loop through variables
  integer(i4b)                  :: iSoil                ! loop through soil layers
  integer(i4b)                  :: iLook                ! loop through lookup table
  integer(i4b)                  :: jIntegr8             ! index for numerical integration
  logical(lgt)                  :: check                ! flag to check allocation
  real(rkind)                   :: vGn_m                ! van Genuchten "m" parameter (-)
  real(rkind)                   :: vFracLiq             ! volumetric fraction of liquid water (-)
  real(rkind)                   :: matricHead           ! matric head (m)
  ! -------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message="T2L_lookup_soil/"

  ! get the values of temperature for the lookup table
  xIncr = 1._rkind/real(nLook-1, kind(rkind))
  xTemp = T_lower + (Tfreeze - T_lower)*sqrt(sqrt(arth(0._rkind,xIncr,nLook))) ! use sqrt(sqrt()) to give more values near freezing

  ! -----
  ! * allocate space for the lookup table...
  ! ----------------------------------------

  ! initialize checks
  check=.false.

  ! allocate space for soil layers
  if(allocated(lookup_data%z))then; check=.true.; else; allocate(lookup_data%z(nSoil), stat=err); endif
  if(check) then; err=20; message=trim(message)//'lookup table z dimension was unexpectedly allocated already'; return; end if
  if(err/=0)then; err=20; message=trim(message)//'problem allocating lookup table z dimension dimension'; return; end if

  ! allocate space for the variables in the lookup table
  do iSoil=1,nSoil
    if(allocated(lookup_data%z(iSoil)%var))then; check=.true.; else; allocate(lookup_data%z(iSoil)%var(maxvarLookup), stat=err); endif
    if(check) then; err=20; message=trim(message)//'lookup table var dimension was unexpectedly allocated already'; return; end if
    if(err/=0)then; err=20; message=trim(message)//'problem allocating lookup table var dimension dimension'; return; end if

    ! allocate space for the values in the lookup table
    do iVar=1,maxvarLookup
      if(allocated(lookup_data%z(iSoil)%var(iVar)%lookup))then; check=.true.; else; allocate(lookup_data%z(iSoil)%var(iVar)%lookup(nLook), stat=err); endif
      if(check) then; err=20; message=trim(message)//'lookup table value dimension was unexpectedly allocated already'; return; end if
      if(err/=0)then; err=20; message=trim(message)//'problem allocating lookup table vaule dimension dimension'; return; end if

    end do ! (looping through variables)
  end do ! (looping through soil layers)

  ! loop through soil layers
  do iSoil=1,nSoil

    ! -----
    ! * make association to variables in the data structures...
    ! ---------------------------------------------------------

    associate(&

      ! associate model parameters
      snowfrz_scale  => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)           , & ! scaling parameter for freezing     (K-1)
      soil_dens_intr => mpar_data%var(iLookPARAM%soil_dens_intr)%dat(iSoil)      , & ! intrinsic soil density             (kg m-3)
      theta_sat      => mpar_data%var(iLookPARAM%theta_sat)%dat(iSoil)           , & ! soil porosity                      (-)
      theta_res      => mpar_data%var(iLookPARAM%theta_res)%dat(iSoil)           , & ! volumetric residual water content  (-)
      vGn_alpha      => mpar_data%var(iLookPARAM%vGn_alpha)%dat(iSoil)           , & ! van Genuchten "alpha" parameter    (m-1)
      vGn_n          => mpar_data%var(iLookPARAM%vGn_n)%dat(iSoil)               , & ! van Genuchten "n" parameter        (-)

      ! associate values in the lookup table
      Tk            => lookup_data%z(iSoil)%var(iLookLOOKUP%temperature)%lookup  , & ! temperature (K)
      Ly            => lookup_data%z(iSoil)%var(iLookLOOKUP%psiLiq_int)%lookup   , & ! integral of mLayerPsiLiq from Tfreeze to Tk (K)
      L2            => lookup_data%z(iSoil)%var(iLookLOOKUP%deriv2)%lookup         & ! second derivative of the interpolating function

      ) ! end associate statement

      ! compute vGn_m
      vGn_m = 1._rkind - 1._rkind/vGn_n

      ! -----
      ! * populate the lookup table...
      ! ------------------------------

      ! initialize temperature and integral
      Tk(nLook) = Tfreeze
      Ly(nLook) = 0._rkind

      ! loop through lookup table
      do iLook=(nLook-1),1,-1

        ! update temperature and integral
        Tk(iLook) = Tk(iLook+1)
        Ly(iLook) = Ly(iLook+1)

        ! get the temperature increment for the numerical integration
        T_incr = (xTemp(iLook)-xTemp(iLook+1))/real(nIntegr8, kind(rkind))

        ! numerical integration between different values of the lookup table
        do jIntegr8=1,nIntegr8

          ! update temperature
          Tk(iLook)  = Tk(iLook) + T_incr

          ! compute the volumetric liquid water and ice content at the mid point of the temperature increment
          matricHead = (LH_fus/gravity)*(Tk(iLook) - Tfreeze - T_incr/2._rkind)/Tfreeze
          vFracLiq   = volFracLiq(matricHead,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)

          ! compute integral
          Ly(iLook)  = Ly(iLook) + vFracLiq*T_incr
  
        end do  ! numerical integration

      end do  ! loop through lookup table

      ! use cubic spline interpolation to obtain integral values at the desired values of temperature
      call spline(Tk,Ly,1.e30_rkind,1.e30_rkind,L2,err,cmessage)  ! get the second derivatives
      if(err/=0) then; message=trim(message)//trim(cmessage); return; end if

    ! end asssociation to variables in the data structures
    end associate

  end do  ! (looping through soil layers)
end subroutine T2L_lookup_soil


! ************************************************************************************************************************
! public subroutine H2T_snow: compute temperature based on specific temperature component of enthalpy 
!                             appropriate when no dry mass, as in snow
! ************************************************************************************************************************
subroutine H2T_snow(Hy,BulkDenWater,fc_param,Tk,err,message)
  ! -------------------------------------------------------------------------------------------------------------------------
  implicit none
  ! -------------------------------------------------------------------------------------------------------------------------
  ! declare dummy variables
  real(rkind),intent(in)      :: Hy            ! total temperature component of enthalpy (J m-3)
  real(rkind),intent(in)      :: BulkDenWater  ! bulk density of water (kg m-3)
  real(rkind),intent(in)      :: fc_param      ! freezing curve parameter (K-1)
  real(rkind),intent(out)     :: Tk            ! initial temperature guess / final temperature value (K)
  integer(i4b),intent(out)    :: err           ! error code
  character(*),intent(out)    :: message       ! error message
  ! declare local variables
  real(rkind),parameter       :: dx=1.d-8      ! finite difference increment (J kg-1)
  real(rkind),parameter       :: atol=1.d-12   ! convergence criteria (J kg-1)
  real(rkind)                 :: H_spec        ! specific enthalpy (J kg-1)
  real(rkind)                 :: H_incr        ! enthalpy increment
  integer(i4b)                :: niter=15      ! maximum number of iterations
  integer(i4b)                :: iter          ! iteration index
  integer(i4b)                :: i0            ! position in lookup table
  real(rkind)                 :: Tg0,Tg1       ! trial temperatures (K)
  real(rkind)                 :: Ht0,Ht1       ! specific enthalpy, based on the trial temperatures (J kg-1)
  real(rkind)                 :: f0,f1         ! function evaluations (difference between enthalpy guesses)
  real(rkind)                 :: dh            ! enthalpy derivative
  real(rkind)                 :: dT            ! temperature increment
  ! -------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message="H2T_snow/"
  ! convert input of total enthalpy (J m-3) to total specific enthalpy (J kg-1)
  H_spec = Hy/BulkDenWater ! (NOTE: no soil)
 
  ! ***** get initial guess and derivative assuming all water is frozen
  if(H_spec<H_lookup(1))then ! process cases below the limit of the look-up table
    ! get temperature guess
    Tg0 = (H_spec - H_lookup(1))/Cp_ice + T_lookup(1)
    Tg1 = Tg0+dx
    ! compute enthalpy
    Ht0 = T2H_snow(Tg0,1._rkind,fc_param)
    Ht1 = T2H_snow(Tg1,1._rkind,fc_param)
    ! compute function evaluations
    f0  = Ht0 - H_spec
    f1  = Ht1 - H_spec

  ! ***** get initial guess and derivative from the look-up table
  else
    ! get enthalpy increment
    H_incr = H_lookup(2) - H_lookup(1)
    ! get position in lookup table
    i0 = ceiling( (H_spec - H_lookup(1)) / H_incr, kind(i4b) )
    ! check found the appropriate value in the look-up table
    if(H_spec < H_lookup(i0) .or. H_spec > H_lookup(i0+1) .or. &
       i0 < 1 .or. i0+1 > nlook)then
     err=10; message=trim(message)//'problem finding appropriate value in lookup table'; return
    end if
    ! get temperature guess
    Tg0 = T_lookup(i0)
    Tg1 = T_lookup(i0+1)
    ! compute function evaluations
    f0  = H_lookup(i0) - H_spec
    f1  = H_lookup(i0+1) - H_spec
  end if

  ! compute initial derivative
  dh  = (f1 - f0) / (Tg1 - Tg0)
  ! compute initial change in T
  dT  = -f0/dh
  ! exit if already close enough
  if(abs(dT)<atol)then
    Tk = Tg0+dT
    return
  end if

  ! **** iterate a little
  do iter=1,niter
    ! comute new value of Tg
    Tg1 = Tg0+dT
    ! get new function evaluation
    Ht1 = T2H_snow(Tg1,1._rkind,fc_param)
    f1  = Ht1 - H_spec
    ! compute derivative if dT
    dh  = (f1 - f0)/dT
    ! compute change in T
    dT  = -f1/dh
    ! exit if converged
    if(abs(dT)<atol)then
      Tk = Tg1+dT
      return
    end if
    ! get ready for next iteration -- save old function evaluation and temperature
    f0  = f1
    Tg0 = Tg1
    ! and check for convergence
    if(iter==niter)then; err=20; message=trim(message)//"failedToConverge"; return; end if
  end do  ! (iteration loop)
end subroutine H2T_snow


! ************************************************************************************************************************
! public function T2H_snow: compute liquid and ice mixture enthalpy based on temperature and mass (J m-3) for a
!                           layer only where the layer has no dry mass, as in snow
!                           NOTE: enthalpy is a relative value, defined as zero at Tfreeze where all water is liquid
! ************************************************************************************************************************
function T2H_snow(Tk,BulkDenWater,fc_param)
  ! -------------------------------------------------------------------------------------------------------------------------
  implicit none
  ! declare dummy variables
  real(rkind),intent(in)  :: Tk              ! layer temperature (K)
  real(rkind),intent(in)  :: BulkDenWater    ! bulk density of water (kg m-3)
  real(rkind),intent(in)  :: fc_param        ! freezing curve parameter (K-1)
  real(rkind)             :: T2H_snow        ! return value of the function, total specific enthalpy (J m-3)
  ! declare local variables
  real(rkind)             :: frac_liq        ! fraction of liquid water
  real(rkind)             :: enthTempWater   ! temperature component of specific enthalpy for total water (liquid and ice) (J kg-1)
  real(rkind)             :: enthMass        ! mass component of specific enthalpy (J kg-1)
  ! -------------------------------------------------------------------------------------------------------------------------
  ! compute the fraction of liquid water in the given layer
  frac_liq     = 1._rkind / ( 1._rkind + ( fc_param*( Tfreeze - min(Tk,Tfreeze) ) )**2_i4b )

  ! compute the temperature component of enthalpy for total water (J kg-1)
  ! NOTE: negative enthalpy means require energy to bring to Tfreeze
  if(Tk< Tfreeze) enthTempWater = Cp_ice*(Tk - Tfreeze) - (Cp_water - Cp_ice)*(atan(fc_param*(Tfreeze - Tk))/fc_param)
  if(Tk>=Tfreeze) enthTempWater = Cp_water*(Tk - Tfreeze)

  ! compute the mass component of enthalpy -- energy required to melt ice (J kg-1)
  ! NOTE: negative enthalpy means require energy to bring to Tfreeze
  enthMass = -LH_fus*(1._rkind - frac_liq)

  ! finally, compute the total enthalpy (J m-3)
  T2H_snow = BulkDenWater*(enthTempWater + enthMass) !+ BulkDenSoil*enthTempSoil
end function T2H_snow


! ************************************************************************************************************************
! public subroutine T2enthTemp: compute temperature component of enthalpy from temperature and total water content
! ************************************************************************************************************************
subroutine T2enthTemp(&
                      use_lookup,                        & ! intent(in):  flag to use the lookup table for soil enthalpy
                      ! input: data structures
                      diag_data,                         & ! intent(in):  model diagnostic variables for a local HRU
                      mpar_data,                         & ! intent(in):  parameter data structure
                      indx_data,                         & ! intent(in):  model indices
                      lookup_data,                       & ! intent(in):  lookup table data structure
                      ! input: state variables for the vegetation canopy
                      scalarCanairTempTrial,             & ! intent(in):  trial value of canopy air temperature (K)
                      scalarCanopyTempTrial,             & ! intent(in):  trial value of canopy temperature (K)
                      scalarCanopyWatTrial,              & ! intent(in):  trial value of canopy total water (kg m-2)
                      ! input: variables for the snow-soil domain
                      mLayerTempTrial,                   & ! intent(in):  trial vector of layer temperature (K)
                      mLayerVolFracWatTrial,             & ! intent(in):  trial vector of volumetric total water content (-)
                      mLayerMatricHeadTrial,             & ! intent(in):  trial vector of total water matric potential (m)
                      ! output: enthalpy
                      scalarCanairEnthalpy,              & ! intent(out): enthalpy of the canopy air space (J m-3)
                      scalarCanopyEnthTemp,              & ! intent(out): temperature component of enthalpy of the vegetation canopy (J m-3)
                      mLayerEnthTemp,                    & ! intent(out): temperature component of enthalpy of each snow+soil layer (J m-3)
                      ! output: error control
                      err,message)                         ! intent(out): error control
  ! -------------------------------------------------------------------------------------------------------------------------
  ! downwind routines
  USE soil_utils_module,only:crit_soilT     ! compute critical temperature below which ice exists
  USE soil_utils_module,only:volFracLiq     ! compute volumetric fraction of liquid water
  USE spline_int_module,only:splint         ! use for cubic spline interpolation
  implicit none
  ! delare dummy variables
  ! -------------------------------------------------------------------------------------------------------------------------
  logical(lgt),intent(in)          :: use_lookup                ! flag to use the lookup table for soil enthalpy, otherwise use hypergeometric function
  ! input: data structures
  type(var_dlength),intent(in)     :: diag_data                 ! diagnostic variables for a local HRU
  type(var_dlength),intent(in)     :: mpar_data                 ! model parameters
  type(var_ilength),intent(in)     :: indx_data                 ! model indices
  type(zLookup),intent(in)         :: lookup_data               ! lookup tables
  ! input: state variables for the vegetation canopy
  real(rkind),intent(in)           :: scalarCanairTempTrial     ! trial value of canopy air temperature (K)
  real(rkind),intent(in)           :: scalarCanopyTempTrial     ! trial value of canopy temperature (K)
  real(rkind),intent(in)           :: scalarCanopyWatTrial      ! trial value of canopy total water (kg m-2)
  ! input: variables for the snow-soil domain
  real(rkind),intent(in)           :: mLayerTempTrial(:)        ! trial vector of layer temperature (K)
  real(rkind),intent(in)           :: mLayerVolFracWatTrial(:)  ! trial vector of volumetric total water content (-)
  real(rkind),intent(in)           :: mLayerMatricHeadTrial(:)  ! trial vector of total water matric potential (m)
  ! output: enthalpy
  real(rkind),intent(out)          :: scalarCanairEnthalpy      ! enthalpy of the canopy air space (J m-3)
  real(rkind),intent(out)          :: scalarCanopyEnthTemp      ! temperature component of enthalpy of the vegetation canopy (J m-3)
  real(rkind),intent(out)          :: mLayerEnthTemp(:)         ! temperature component of enthalpy of each snow+soil layer (J m-3)
  ! output: error control
  integer(i4b),intent(out)         :: err                       ! error code
  character(*),intent(out)         :: message                   ! error message
  ! -------------------------------------------------------------------------------------------------------------------------
  ! declare local variables
  character(len=128)               :: cmessage                  ! error message in downwind routine
  integer(i4b)                     :: iState                    ! index of model state variable
  integer(i4b)                     :: iLayer                    ! index of model layer
  integer(i4b)                     :: ixFullVector              ! index within full state vector
  integer(i4b)                     :: ixDomainType              ! name of a given model domain
  integer(i4b)                     :: ixControlIndex            ! index within a given model domain
  real(rkind)                      :: vGn_m                     ! van Genuchten "m" parameter (-)
  real(rkind)                      :: Tcrit                     ! temperature where all water is unfrozen (K)
  real(rkind)                      :: volFracWat                ! volumetric fraction of total water, liquid+ice (-)
  real(rkind)                      :: diff0                     ! temperature difference of Tcrit from Tfreeze
  real(rkind)                      :: diffT                     ! temperature difference of temp soil from Tfreeze
  real(rkind)                      :: integral                  ! integral of snow freezing curve
  real(rkind)                      :: dTcrit_dPsi0              ! derivative of temperature where all water is unfrozen (K) with matric head
  real(rkind)                      :: dL                        ! derivative of enthalpy with temperature at layer temperature
  real(rkind)                      :: arg                       ! argument of hypergeometric function
  real(rkind)                      :: gauss_hg_T                ! hypergeometric function result
  real(rkind)                      :: integral_unf              ! integral of unfrozen soil water content (from Tfreeze to Tcrit)
  real(rkind)                      :: integral_frz_low          ! lower limit of integral of frozen soil water content (from Tfreeze to Tcrit)
  real(rkind)                      :: integral_frz_upp          ! upper limit of integral of frozen soil water content (from Tfreeze to soil temperature)
  real(rkind)                      :: xConst                    ! constant in the freezing curve function (m K-1)
  real(rkind)                      :: mLayerPsiLiq              ! liquid water matric potential (m)
  ! enthalpy
  real(rkind)                      :: enthVeg                   ! enthalpy of the vegetation (J m-3)
  real(rkind)                      :: enthSoil                  ! enthalpy of soil particles (J m-3)
  real(rkind)                      :: enthLiq                   ! enthalpy of the liquid region (J m-3)
  real(rkind)                      :: enthIce                   ! enthalpy of the ice region (J m-3)
  real(rkind)                      :: enthAir                   ! enthalpy of air (J m-3)
  real(rkind)                      :: enthPhase                 ! enthalpy associated with phase change (J m-3)
  real(rkind)                      :: enthWater                 ! enthalpy of total water (J m-3)
  logical(lgt),parameter           :: doTest=.false.            ! flag to run unit test
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! make association with variables in the data structures
  generalVars: associate(&
    ! number of model layers, and layer type
    nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in):  [i4b]    total number of snow layers
    nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in):  [i4b]    total number of soil layers
    nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)               ,& ! intent(in):  [i4b]    total number of snow and soil layers
    ! mapping between the full state vector and the state subset
    ixMapFull2Subset        => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat         ,& ! intent(in):  [i4b(:)] list of indices in the state subset for each state in the full state vector
    ixMapSubset2Full        => indx_data%var(iLookINDEX%ixMapSubset2Full)%dat         ,& ! intent(in):  [i4b(:)] [state subset] list of indices of the full state vector in the state subset
    ! type of domain, type of state variable, and index of control volume within domain
    ixDomainType_subset     => indx_data%var(iLookINDEX%ixDomainType_subset)%dat      ,& ! intent(in):  [i4b(:)] [state subset] id of domain for desired model state variables
    ixControlVolume         => indx_data%var(iLookINDEX%ixControlVolume)%dat          ,& ! intent(in):  [i4b(:)] index of the control volume for different domains (veg, snow, soil)
    ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat               & ! intent(in):  [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
    ) ! end associate statement
    ! ------------------------------------------------------------------------------------------------------------------------------

    ! initialize error control
    err=0; message="T2enthTemp/"

    ! loop through model state variables
    do iState=1,size(ixMapSubset2Full)

      ! -----
      ! - compute indices...
      ! --------------------

      ! get domain type, and index of the control volume within the domain
      ixFullVector   = ixMapSubset2Full(iState)       ! index within full state vector
      ixDomainType   = ixDomainType_subset(iState)    ! named variables defining the domain (iname_cas, iname_veg, etc.)
      ixControlIndex = ixControlVolume(ixFullVector)  ! index within a given domain

      ! check an energy state
      if(ixStateType(ixFullVector)==iname_nrgCanair .or. ixStateType(ixFullVector)==iname_nrgCanopy .or. ixStateType(ixFullVector)==iname_nrgLayer)then

        ! get the layer index
        select case(ixDomainType)
          case(iname_cas);     iLayer = integerMissing
          case(iname_veg);     iLayer = integerMissing
          case(iname_snow);    iLayer = ixControlIndex
          case(iname_soil);    iLayer = ixControlIndex + nSnow
          case(iname_aquifer); cycle ! aquifer: do nothing (no thermodynamics in the aquifer)
          case default; err=20; message=trim(message)//'expect case to be iname_cas, iname_veg, iname_snow, iname_soil, iname_aquifer'; return
        end select

        ! identify domain
        select case(ixDomainType)
          case(iname_cas)
            scalarCanairEnthalpy = Cp_air * iden_air * (scalarCanairTempTrial - Tfreeze)

          case(iname_veg)
            ! association to necessary variables for vegetation
            vegVars: associate(&
              canopyDepth       => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1),   & ! canopy depth                                   (m)
              specificHeatVeg   => mpar_data%var(iLookPARAM%specificHeatVeg)%dat(1),    & ! specific heat of vegetation                    (J kg-1 K-1)
              maxMassVegetation => mpar_data%var(iLookPARAM%maxMassVegetation)%dat(1),  & ! maximum mass of vegetation                     (kg m-2)
              snowfrz_scale     => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)       & ! scaling parameter for the snow freezing curve  (K-1)
              )

              diffT = scalarCanopyTempTrial - Tfreeze
              enthVeg = specificHeatVeg * maxMassVegetation * diffT / canopyDepth

              if(diffT>=0._rkind)then
                enthLiq = Cp_water * scalarCanopyWatTrial * diffT / canopyDepth
                enthIce = 0._rkind
              else
                integral = (1._rkind/snowfrz_scale) * atan(snowfrz_scale * diffT)
                enthLiq = Cp_water * scalarCanopyWatTrial * integral / canopyDepth
                enthIce = Cp_ice * scalarCanopyWatTrial * ( diffT - integral ) / canopyDepth
              endif

              scalarCanopyEnthTemp = enthVeg + enthLiq + enthIce

            end associate vegVars

          case(iname_snow)

            ! association to necessary variables for snow
            snowVars: associate(&
              snowfrz_scale => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)   & ! intent(in):  [dp] scaling parameter for the snow freezing curve (K-1)
              )

              diffT = mLayerTempTrial(iLayer) - Tfreeze  ! diffT<0._rkind because snow is frozen
              integral = (1._rkind/snowfrz_scale) * atan(snowfrz_scale * diffT)
              enthLiq = iden_water * Cp_water * mLayerVolFracWatTrial(iLayer) * integral
              enthIce = iden_water * Cp_ice * mLayerVolFracWatTrial(iLayer) * ( diffT - integral )
              enthAir = iden_air * Cp_air * ( diffT - mLayerVolFracWatTrial(iLayer) * ( (iden_water/iden_ice)*(diffT-integral) + integral ) )
 
              mLayerEnthTemp(iLayer) = enthLiq + enthIce + enthAir

            end associate snowVars

          case(iname_soil)
            ! make association to variables for soil
            soilVars: associate(&
              soil_dens_intr => mpar_data%var(iLookPARAM%soil_dens_intr)%dat(ixControlIndex)      , & ! intrinsic soil density             (kg m-3)
              theta_sat      => mpar_data%var(iLookPARAM%theta_sat)%dat(ixControlIndex)           , & ! soil porosity                      (-)
              theta_res      => mpar_data%var(iLookPARAM%theta_res)%dat(ixControlIndex)           , & ! volumetric residual water content  (-)
              vGn_alpha      => mpar_data%var(iLookPARAM%vGn_alpha)%dat(ixControlIndex)           , & ! van Genuchten "alpha" parameter    (m-1)
              vGn_n          => mpar_data%var(iLookPARAM%vGn_n)%dat(ixControlIndex)                 & ! van Genuchten "n" parameter        (-)
              ) ! end associate statement
              
              ! diagnostic variables
              vGn_m    = 1._rkind - 1._rkind/vGn_n
              Tcrit    = crit_soilT( mLayerMatricHeadTrial(ixControlIndex) )
              volFracWat = volFracLiq(mLayerMatricHeadTrial(ixControlIndex),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
              diffT = mLayerTempTrial(iLayer) - Tfreeze
              diff0 = Tcrit - Tfreeze
              ! *** compute enthalpy of water for unfrozen conditions
              if(mlayerTempTrial(iLayer)>=Tcrit)then
                enthWater = iden_water * Cp_water * volFracWat * diffT ! valid for temperatures below freezing also

              ! *** compute enthalpy of water for frozen conditions
              else
                ! *** compute integral of mLayerPsiLiq from Tfreeze to layer temperature
                ! get the unfrozen water content
                integral_unf = ( Tcrit - Tfreeze ) * volFracWat

                ! get the frozen  water content
                ! initialize for case Tcrit=Tfreeze, i.e. mLayerMatricHeadTrial(ixControlIndex)>0
                integral_frz_low = 0._rkind 

                if(use_lookup)then ! cubic spline interpolation for integral of mLayerPsiLiq from Tfreeze to layer temperature
                  ! make associate to the the lookup table
                  lookVars: associate(&
                    Tk => lookup_data%z(ixControlIndex)%var(iLookLOOKUP%temperature)%lookup  , & ! temperature (K)
                    Ly => lookup_data%z(ixControlIndex)%var(iLookLOOKUP%psiLiq_int)%lookup   , & ! integral of mLayerPsiLiq from Tfreeze to Tk (K)
                    L2 => lookup_data%z(ixControlIndex)%var(iLookLOOKUP%deriv2)%lookup         & ! second derivative of the interpolating function
                    ) ! end associate statement

                    ! get the lower limit of the integral
                    if(diff0<0._rkind)then
                      call splint(Tk,Ly,L2,Tcrit,integral_frz_low,dL,err,cmessage)
                      if(err/=0) then; message=trim(message)//trim(cmessage); return; end if
                    end if
                    ! get the upper limit of the integral
                    call splint(Tk,Ly,L2,mlayerTempTrial(iLayer),integral_frz_upp,dL,err,cmessage)
                    if(err/=0) then; message=trim(message)//trim(cmessage); return; end if

                  end associate lookVars

                else ! hypergeometric function for integral of mLayerPsiLiq from Tfreeze to layer temperature
                  ! get the lower limit of the integral
                  if(diff0<0._rkind)then
                    arg = (vGn_alpha * mLayerMatricHeadTrial(ixControlIndex))**vGn_n
                    gauss_hg_T = hyp_2F1_real(vGn_m,1._rkind/vGn_n,1._rkind + 1._rkind/vGn_n,-arg)
                    integral_frz_low = diff0 * ( (theta_sat - theta_res)*gauss_hg_T + theta_res )
                  end if
                  ! get the upper limit of the integral
                  xConst        = LH_fus/(gravity*Tfreeze)        ! m K-1 (NOTE: J = kg m2 s-2)
                  mLayerPsiLiq  = xConst*diffT   ! liquid water matric potential from the Clapeyron eqution, DIFFERENT from the liquid water matric potential used in the flux calculations
                  arg = (vGn_alpha * mLayerPsiLiq)**vGn_n
                  gauss_hg_T = hyp_2F1_real(vGn_m,1._rkind/vGn_n,1._rkind + 1._rkind/vGn_n,-arg)
                  integral_frz_upp = diffT * ( (theta_sat - theta_res)*gauss_hg_T + theta_res )
                endif

                enthLiq = iden_water * Cp_water * (integral_unf + integral_frz_upp - integral_frz_low)
                enthIce = iden_ice * Cp_ice * ( volFracWat * diffT - (integral_unf + integral_frz_upp - integral_frz_low) )
                enthWater = enthIce + enthLiq

                if(doTest)then
          
                  ! write values
                  print*, 'doTest    = ', doTest
                  print*, 'T,Tcrit   = ', mlayerTempTrial(iLayer),Tcrit    ! temperature (K)
                  print*, 'integral unf,frz_upp,frz_low = ', integral_unf, integral_frz_upp, integral_frz_low    ! integral (K)
                  print*, 'theta_sat = ', theta_sat ! soil porosity                      (-)
                  print*, 'theta_res = ', theta_res ! volumetric residual water content  (-)
                  print*, 'vGn_alpha = ', vGn_alpha ! van Genuchten "alpha" parameter    (m-1)
                  print*, 'vGn_n     = ', vGn_n     ! van Genuchten "n" parameter        (-)
                  print*, trim(message)//'PAUSE: Set doTest=.false. to complete simulations'
                  read(*,*)
          
                endif  ! if testing

              endif ! (if frozen conditions)

              enthSoil = soil_dens_intr*Cp_soil*(1._rkind - theta_sat)*diffT
              enthAir = iden_air*Cp_air*(1._rkind - theta_sat - volFracWat)*diffT

              mLayerEnthTemp(iLayer) = enthWater + enthSoil + enthAir

            end associate soilVars

          ! -----
          ! - checks...
          ! -----------
          case(iname_aquifer); cycle ! aquifer: do nothing
          case default; err=20; message=trim(message)//'expect case to be iname_cas, iname_veg, iname_snow, iname_soil, iname_aquifer'; return
        end select

      end if  ! if an energy layer
    end do  ! looping through state variables

  end associate generalVars

end subroutine T2enthTemp


! ************************************************************************************************************************
! public subroutine enthTemp2H: add energy associated with thaw/freeze to temperature component of enthalpy to get total enthalpy, H
! ************************************************************************************************************************
subroutine enthTemp2H(&
                      ! input: data structures
                      diag_data,               & ! intent(in):    model diagnostic variables for a local HRU
                      indx_data,               & ! intent(in):    model indices
                      ! input: ice content change
                      scalarCanopyIce,         & ! intent(in):    value of canopy ice content (kg m-2) or prime ice content (kg m-2 s-1)
                      mLayerVolFracIce,        & ! intent(in):    vector of volumetric fraction of ice (-) or prime volumetric fraction of ice (s-1)
                      ! input/output: enthalpy
                      scalarCanopyH,           & ! intent(inout): enthalpy of the vegetation canopy (J m-3)
                      mLayerH,                 & ! intent(inout): enthalpy of each snow+soil layer (J m-3)
                      ! output: error control
                      err,message)               ! intent(out): error control
  ! -------------------------------------------------------------------------------------------------------------------------
  implicit none
  ! delare dummy variables
  ! -------------------------------------------------------------------------------------------------------------------------
  ! input: data structures
  type(var_dlength),intent(in)     :: diag_data                  ! diagnostic variables for a local HRU
  type(var_ilength),intent(in)     :: indx_data                  ! model indices
  ! input: ice content change
  real(rkind),intent(in)           :: scalarCanopyIce            ! value for canopy ice content (kg m-2) or prime ice content (kg m-2 s-1)
  real(rkind),intent(in)           :: mLayerVolFracIce(:)        ! vector of volumetric fraction of ice (-) or prime volumetric fraction of ice (s-1)
  ! input output: enthalpy
  real(rkind),intent(inout)        :: scalarCanopyH              ! value for enthalpy of the vegetation canopy (J m-3 s-1)
  real(rkind),intent(inout)        :: mLayerH(:)                 ! vector of enthalpy of each snow+soil layer (J m-3 s-1)
  ! output: error control
  integer(i4b),intent(out)         :: err                        ! error code
  character(*),intent(out)         :: message                    ! error message
  ! -------------------------------------------------------------------------------------------------------------------------
  ! declare local variables
  character(len=128)               :: cmessage                   ! error message in downwind routine
  integer(i4b)                     :: iState                     ! index of model state variable
  integer(i4b)                     :: iLayer                     ! index of model layer
  integer(i4b)                     :: ixFullVector               ! index within full state vector
  integer(i4b)                     :: ixDomainType               ! name of a given model domain
  integer(i4b)                     :: ixControlIndex             ! index within a given model domain
   ! ------------------------------------------------------------------------------------------------------------------------
  ! make association with variables in the data structures
  generalVars: associate(&
    ! number of model layers, and layer type
    nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)            ,& ! intent(in): [i4b]    total number of snow layers
    nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)            ,& ! intent(in): [i4b]    total number of soil layers
    nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)          ,& ! intent(in): [i4b]    total number of snow and soil layers
    ! mapping between the full state vector and the state subset
    ixMapFull2Subset        => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat    ,& ! intent(in): [i4b(:)] list of indices in the state subset for each state in the full state vector
    ixMapSubset2Full        => indx_data%var(iLookINDEX%ixMapSubset2Full)%dat    ,& ! intent(in): [i4b(:)] [state subset] list of indices of the full state vector in the state subset
    ! type of domain, type of state variable, and index of control volume within domain
    ixDomainType_subset     => indx_data%var(iLookINDEX%ixDomainType_subset)%dat ,& ! intent(in): [i4b(:)] [state subset] id of domain for desired model state variables
    ixControlVolume         => indx_data%var(iLookINDEX%ixControlVolume)%dat     ,& ! intent(in): [i4b(:)] index of the control volume for different domains (veg, snow, soil)
    ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat         ,& ! intent(in): [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
   ! canopy depth
    canopyDepth             => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)  & ! intent(in): [dp]     canopy depth (m)
    ) ! end associate statement
    ! -----------------------------------------------------------------------------------------------------------------------

    ! initialize error control
    err=0; message="enthTemp2H/"

    ! loop through model state variables
    do iState=1,size(ixMapSubset2Full)

      ! -----
      ! - compute indices...
      ! --------------------

      ! get domain type, and index of the control volume within the domain
      ixFullVector   = ixMapSubset2Full(iState)       ! index within full state vector
      ixDomainType   = ixDomainType_subset(iState)    ! named variables defining the domain (iname_cas, iname_veg, etc.)
      ixControlIndex = ixControlVolume(ixFullVector)  ! index within a given domain

      ! check an energy state
      if(ixStateType(ixFullVector)==iname_nrgCanair .or. ixStateType(ixFullVector)==iname_nrgCanopy .or. ixStateType(ixFullVector)==iname_nrgLayer)then

        ! get the layer index
        select case(ixDomainType)
          case(iname_cas);     cycle ! canopy air space: do nothing (no water stored in canopy air space)
          case(iname_veg)
            scalarCanopyH= scalarCanopyH - LH_fus * scalarCanopyIce/ canopyDepth
          case(iname_snow)
            iLayer = ixControlIndex
            mLayerH(iLayer) = mLayerH(iLayer) - iden_ice   * LH_fus * mLayerVolFracIce(iLayer)
          case(iname_soil)
            iLayer = ixControlIndex + nSnow
            mLayerH(iLayer) = mLayerH(iLayer) - iden_water * LH_fus * mLayerVolFracIce(iLayer)
          case(iname_aquifer); cycle ! aquifer: do nothing (no thermodynamics in the aquifer)
          case default; err=20; message=trim(message)//'expect case to be iname_cas, iname_veg, iname_snow, iname_soil, iname_aquifer'; return
        end select

      end if  ! if an energy layer
    end do  ! looping through state variables

  end associate generalVars

end subroutine enthTemp2H

!----------------------------------------------------------------------
! private function: compute hypergeometric function with real arguments into real result
!----------------------------------------------------------------------
 function hyp_2F1_real(a_real, b_real, c_real, z_real)
  !--------------------------------------------------------------------
  USE hyp_2F1_module,only:HYP_2F1 ! use for hypergeometric function
  implicit none
  real(rkind),intent(in) :: a_real, b_real, c_real, z_real
  complex(rkind)         :: a_complex, b_complex, c_complex, z_complex, result
  real(rkind)            :: hyp_2F1_real
  
  a_complex = CMPLX(a_real, 0._rkind, rkind)
  b_complex = CMPLX(b_real, 0._rkind, rkind)
  c_complex = CMPLX(c_real, 0._rkind, rkind)
  z_complex = CMPLX(z_real, 0._rkind, rkind)
  result = HYP_2F1(a_complex, b_complex, c_complex, z_complex)
  hyp_2F1_real = REAL(result, rkind)
   
end function hyp_2F1_real

end module enthalpyTemp_module