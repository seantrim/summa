module computPrime_module

! data types
USE nrtype

! derived types to define the data structures
USE data_types,only:&
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength     ! data vector with variable length dimension (rkind)

! named variables
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookPARAM            ! named variables for structure elements

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
USE globalData,only:iname_cas       ! named variables for canopy air space
USE globalData,only:iname_aquifer   ! named variables for the aquifer

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

USE soil_utils_module,only:dTheta_dPsi    ! derivative in the soil water characteristic (soil)
USE soil_utils_module,only:dPsi_dTheta    ! derivative in the soil water characteristic (soil)


public::computePrime

contains

  
  subroutine computePrime(&
       dt, &
       fVec, &
       PrimeVec, &
       sMul, &
       indx_data, &
       diag_data, &
       nLayers, &
       scalarCanopyCmTrial, &
       mLayerCmTrial, &
       mLayerMatricHeadTrial, &
       mpar_data, &
       mLayerVolFracWatTrial, scalarCanopyIce, scalarCanopyIceTrial, mLayerVolFracIce, mLayerVolFracIceTrial, mLayerDepth, flux_data)
    implicit none

    real(rkind), intent(in) :: dt
    real(rkind),intent(in) :: fVec(:)
    real(rkind), intent(out) :: PrimeVec(:)
    real(qp), intent(in) :: sMul(:)
    type(var_ilength), intent(in) :: indx_data
    type(var_dlength), intent(in) :: diag_data
    
    integer(i4b),intent(in)         :: nLayers                   ! total number of layers in the snow+soil domain
      real(qp),intent(in)             :: scalarCanopyCmTrial       ! Cm of vegetation canopy (-)
      real(qp),intent(in)             :: mLayerCmTrial(:)          ! Cm of each snow and soil layer (-)
        real(rkind),intent(inout)          :: mLayerMatricHeadTrial(:)        ! trial vector of total water matric potential (m)
  type(var_dlength),intent(in)       :: mpar_data                       ! model parameters for a local HRU
  real(rkind), intent(in) :: mLayerVolFracWatTrial(:)
    real(rkind),intent(in)          :: scalarCanopyIce         ! previous value for mass of ice on the vegetation canopy (kg m-2)
    real(rkind),intent(in)          :: scalarCanopyIceTrial         ! previous value for mass of ice on the vegetation canopy (kg m-2)
    real(rkind), intent(in) :: mLayerVolFracIce(:)
    real(rkind), intent(in) :: mLayerVolFracIceTrial(:)
    real(rkind),intent(in) :: mLayerDepth(:)
    type(var_dlength),intent(in) :: flux_data




    real(rkind) :: scalarCanairTempPrime
    real(rkind) :: scalarCanopyHydPrime
    real(rkind) :: scalarCanopyTempPrime
    real(rkind) :: scalarAquiferStoragePrime

    real(rkind), allocatable :: mLayerVolFracHydPrime (:)
    real(rkind), allocatable :: mLayerMatricHeadPrime (:)
    real(rkind), allocatable :: mLayerTempPrime (:)
      logical(lgt),allocatable           :: computedCoupling(:)             ! .true. if computed the coupling for a given state variable


      integer(i4b)                     :: iLayer                    ! index of layer within the snow+soil domain
  logical(lgt)                       :: isCoupled                       ! .true. if a given variable shared another state variable in the same control volume
  logical(lgt)                       :: isNrgState                      ! .true. if a given variable is an energy state
  integer(i4b)                       :: iState                          ! index of model state variable
  integer(i4b)                       :: ixControlIndex                  ! index within a given model domain
  integer(i4b)                       :: ixDomainType                    ! name of a given model domain
  integer(i4b)                       :: ixFullVector                    ! index within full state vector
  integer(i4b)                       :: ixOther,ixOtherLocal            ! index of the coupled state variable within the (full, local) vector


  associate(&
       nSoil => indx_data%var(iLookINDEX%nSoil)%dat(1), &
           nSoilOnlyHyd            => indx_data%var(iLookINDEX%nSoilOnlyHyd )%dat(1)         ,& ! intent(in): [i4b]    number of hydrology variables in the soil domain

           nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in):  [i4b]    total number of snow layers
    ixMapSubset2Full        => indx_data%var(iLookINDEX%ixMapSubset2Full)%dat         ,& ! intent(in):  [i4b(:)] [state subset] list of indices of the full state vector in the state subset
    ixDomainType_subset     => indx_data%var(iLookINDEX%ixDomainType_subset)%dat      ,& ! intent(in):  [i4b(:)] [state subset] id of domain for desired model state variables
    ixControlVolume         => indx_data%var(iLookINDEX%ixControlVolume)%dat          ,& ! intent(in):  [i4b(:)] index of the control volume for different domains (veg, snow, soil)
    ixHydCanopy             => indx_data%var(iLookINDEX%ixHydCanopy)%dat              ,& ! intent(in):  [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the canopy domain
    ixNrgCanopy             => indx_data%var(iLookINDEX%ixNrgCanopy)%dat              ,& ! intent(in):  [i4b(:)] indices IN THE FULL VECTOR for energy states in the canopy domain
    ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat              ,& ! intent(in):  [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
    ixHydLayer              => indx_data%var(iLookINDEX%ixHydLayer)%dat               ,& ! intent(in):  [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the snow+soil domain
    ixNrgLayer              => indx_data%var(iLookINDEX%ixNrgLayer)%dat               ,& ! intent(in):  [i4b(:)] indices IN THE FULL VECTOR for energy states in the snow+soil domain
    ixMapFull2Subset        => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat         ,& ! intent(in):  [i4b(:)] list of indices in the state subset for each state in the full state vector
    vGn_alpha               => mpar_data%var(iLookPARAM%vGn_alpha)%dat                ,& ! intent(in):  [dp(:)] van Genutchen "alpha" parameter (m-1)
    theta_res               => mpar_data%var(iLookPARAM%theta_res)%dat                ,& ! intent(in):  [dp(:)] soil residual volumetric water content (-)
    theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat                ,& ! intent(in):  [dp(:)] soil porosity (-)
    vGn_n                   => mpar_data%var(iLookPARAM%vGn_n)%dat                    ,& ! intent(in):  [dp(:)] van Genutchen "n" parameter (-)
    vGn_m                   => diag_data%var(iLookDIAG%scalarVGn_m)%dat               ,& ! intent(in):  [dp(:)] van Genutchen "m" parameter (-)
    mLayerBaseflow          => flux_data%var(iLookFLUX%mLayerBaseflow)%dat            ,& ! intent(in): [dp(:)]  baseflow from each soil layer (m s-1)

               ixAqWat                 => indx_data%var(iLookINDEX%ixAqWat)%dat(1)               ,& ! intent(in): [i4b]    index of water storage in the aquifer
    ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in): [i4b]    index of canopy air space energy state variable
    ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,& ! intent(in): [i4b]    index of canopy hydrology state variable (mass)
    ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in): [i4b]    index of canopy energy state variable
    nSnowSoilNrg            => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in): [i4b]    number of energy state variables in the snow+soil domain
        ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in): [i4b(:)] indices for hydrology states in the snow+soil subdomain
    ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in): [i4b(:)] indices for energy states in the snow+soil subdomain
    ixHydType           => indx_data%var(iLookINDEX%ixHydType)%dat                ,& ! intent(in) : [i4b(:)] index of the type of hydrology states in snow+soil domain


    canopyDepth             => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1) ,     &
    nSnowSoilHyd            => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1),         & ! intent(in): [i4b]    number of hydrology variables in the snow+soil domain
        ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat            ,& ! intent(in): [i4b(:)] indices for hydrology states in the soil subdomain
    mLayerTranspire         => flux_data%var(iLookFLUX%mLayerTranspire)%dat           ,& ! intent(in): [dp]     transpiration loss from each soil layer (m s-1)

    mLayerCompress          => diag_data%var(iLookDIAG%mLayerCompress)%dat            ,& ! intent(in): [dp(:)]  change in storage associated with compression of the soil matrix (-)

    layerType               => indx_data%var(iLookINDEX%layerType)%dat                 & ! intent(in): [i4b(:)] named variables defining the type of layer in snow+soil domain

         )

    !if (ixVegNrg/=integerMissing) fVec(ixVegNrg) = fVec(ixVegNrg) + LH_fus* (scalarCanopyIceTrial - scalarCanopyIce) / canopyDepth / dt

    ! if (nSnowSoilNrg > 0) then
    !          do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
    !     select case( layerType(iLayer) )
    !       case(iname_snow); fVec( ixSnowSoilNrg(iLayer) ) = fVec( ixSnowSoilNrg(iLayer) ) + (LH_fus*iden_ice  *( mLayerVolFracIceTrial(iLayer) - mLayerVolFracIce(iLayer) )) / dt
    !       case(iname_soil); fVec( ixSnowSoilNrg(iLayer) ) = fVec( ixSnowSoilNrg(iLayer) ) + (LH_fus*iden_water*( mLayerVolFracIceTrial(iLayer) - mLayerVolFracIce(iLayer) )) / dt
    !     end select
    !   end do  ! looping through non-missing energy state variables in the snow+soil domain
    ! endif


    ! if(nSoilOnlyHyd>0)then
    !   do concurrent (iLayer=1:nSoil,ixSoilOnlyHyd(iLayer)/=integerMissing)   ! (loop through non-missing hydrology state variables in the snow+soil domain)
    !     fVec( ixSoilOnlyHyd(iLayer) ) = fVec( ixSoilOnlyHyd(iLayer) ) + ( ( mLayerTranspire(iLayer) - mLayerBaseflow(iLayer) )/mLayerDepth(iLayer+nSnow) - mLayerCompress(iLayer) )
    !   end do  ! looping through non-missing energy state variables in the snow+soil domain
    ! endif


    allocate(mLayerVolFracHydPrime(nLayers))
    allocate(mLayerTempPrime (nLayers))
    allocate(mLayerMatricHeadPrime (nLayers))

      if (ixCasNrg/=integerMissing) scalarCanairTempPrime = fVec(ixCasNrg) / sMul(ixCasNrg)
      if(ixVegHyd/=integerMissing)then
      scalarCanopyHydPrime =  fVec(ixVegHyd) / sMul(ixVegHyd)
   endif
   if(ixVegNrg/=integerMissing) scalarCanopyTempPrime = (- scalarCanopyCmTrial * scalarCanopyHydPrime/canopyDepth + fVec(ixVegNrg) ) / sMul(ixVegNrg)

   if(nSnowSoilHyd>0)then
      do concurrent (iLayer=1:nLayers,ixSnowSoilHyd(iLayer)/=integerMissing)   ! (loop through non-missing hydrology state variables in the snow+soil domain)
        ! (get the correct state variable)
        mLayerVolFracHydPrime(iLayer) = (Max(epsilon(fVec(ixSnowSoilHyd(iLayer))), 1._qp * fVec( ixSnowSoilHyd(iLayer) ))) / sMul(ixSnowSoilHyd(iLayer))!fVec( ixSnowSoilHyd(iLayer) )
        ! (compute the residual)
      end do  ! looping through non-missing energy state variables in the snow+soil domain
   endif

   if(nSnowSoilNrg>0)then
      do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
       mLayerTempPrime(iLayer) = -1/sMul( ixSnowSoilNrg(iLayer) ) * (mLayerCmTrial(iLayer) * mLayerVolFracHydPrime(iLayer) -  fVec( ixSnowSoilNrg(iLayer) ))
      end do  ! looping through non-missing energy state variables in the snow+soil domain
   endif


    if(ixAqWat/=integerMissing)  scalarAquiferStoragePrime  = 1/sMul(ixAqWat) *   fVec(ixAqWat)

        allocate(computedCoupling(size(ixMapSubset2Full)))        ! .true. if computed the coupling for a given state variable
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
        case(iname_cas);     cycle ! canopy air space: do nothing
        case(iname_veg);     iLayer = 0
        case(iname_snow);    iLayer = ixControlIndex
        case(iname_soil);    iLayer = ixControlIndex + nSnow
        case(iname_aquifer); cycle ! aquifer: do nothing
        case default;  return
      end select

      ! get the index of the other (energy or mass) state variable within the full state vector
      select case(ixDomainType)
        case(iname_veg)             ; ixOther = merge(ixHydCanopy(1),    ixNrgCanopy(1),    ixStateType(ixFullVector)==iname_nrgCanopy)
        case(iname_snow, iname_soil); ixOther = merge(ixHydLayer(iLayer),ixNrgLayer(iLayer),ixStateType(ixFullVector)==iname_nrgLayer)
        case default;  return
      end select

      ! get the index in the local state vector
      ixOtherLocal = ixMapFull2Subset(ixOther)  ! ixOtherLocal could equal integerMissing
      !if(ixOtherLocal/=integerMissing) computedCoupling(ixOtherLocal)=.true.

      ! check if we have a coupled solution
      isCoupled    = .false.!(ixOtherLocal/=integerMissing)

      ! check if we are an energy state
      isNrgState   = (ixStateType(ixFullVector)==iname_nrgCanopy .or. ixStateType(ixFullVector)==iname_nrgLayer)

      !print*, isCoupled, isNrgState, ixControlIndex, ixStateType(ixFullVector), ixOtherLocal

      ! update hydrology state variables for the uncoupled solution
      if(.not.isNrgState .and. .not.isCoupled)then

        ! update the total water and the total water matric potential
        if(ixDomainType==iname_soil)then
          select case( ixStateType(ixFullVector) )
          case(iname_lmpLayer)
             print*, 'lmpLayer'
             cycle
              ! --> update the total water from the total water matric potential
          case(iname_matLayer)
             
             !print*, iLayer, mLayerVolFracHydPrime(iLayer), dTheta_dPsi(mLayerMatricHeadTrial(ixControlIndex), vGn_alpha(ixControlIndex), theta_res(ixControlIndex), theta_sat(ixControlIndex), vGn_n(ixControlIndex), vGn_m(ixControlIndex))
             !mLayerMatricHeadPrime(ixControlIndex) = mLayerVolFracHydPrime(iLayer) / dTheta_dPsi(mLayerMatricHeadTrial(ixControlIndex),vGn_alpha(ixControlIndex),theta_res(ixControlIndex),theta_sat(ixControlIndex),vGn_n(ixControlIndex),vGn_m(ixControlIndex))
             !if (mLayerMatricHeadTrial(iLayer) <= 0._rkind) then
                mLayerMatricHeadPrime(ixControlIndex) = mLayerVolFracHydPrime(iLayer) / dTheta_dPsi(mLayerMatricHeadTrial(ixControlIndex),vGn_alpha(ixControlIndex),theta_res(ixControlIndex),theta_sat(ixControlIndex),vGn_n(ixControlIndex),vGn_m(ixControlIndex))
             !else
             !   mLayerMatricHeadPrime(ixControlIndex) = dPsi_dTheta(mLayerVolFracWatTrial(iLayer),vGn_alpha(ixControlIndex),theta_res(ixControlIndex),theta_sat(ixControlIndex),vGn_n(ixControlIndex),vGn_m(ixControlIndex)) * mLayerVolFracHydPrime(iLayer)
             !endif
              ! --> update the total water matric potential (assume already have mLayerVolFracWatTrial given block above)
           case(iname_liqLayer, iname_watLayer)
                           !print*, iLayer, mLayerVolFracHydPrime(iLayer), dPsi_dTheta(mLayerMatricHeadTrial(ixControlIndex), vGn_alpha(ixControlIndex), theta_res(ixControlIndex), theta_sat(ixControlIndex), vGn_n(ixControlIndex), vGn_m(ixControlIndex))

              mLayerMatricHeadPrime(ixControlIndex) = dPsi_dTheta(mLayerVolFracWatTrial(iLayer),vGn_alpha(ixControlIndex),theta_res(ixControlIndex),theta_sat(ixControlIndex),vGn_n(ixControlIndex),vGn_m(ixControlIndex)) * mLayerVolFracHydPrime(iLayer)
            case default;  print*, 'default'; return
          end select
        endif  ! if in the soil domain

      endif  ! if hydrology state variable or uncoupled solution

   end do ! looping through state variables

    ! deallocate space
    deallocate(computedCoupling)        ! .true. if computed the coupling for a given state variable

    !PrimeVec(:) = 0


    if(ixCasNrg/=integerMissing .or. ixVegNrg/=integerMissing .or. ixVegHyd/=integerMissing)then

       
      ! extract temperature of the canopy air space
      if(ixCasNrg/=integerMissing) PrimeVec(ixCasNrg) = scalarCanairTempPrime

      ! extract canopy temperature
      if(ixVegNrg/=integerMissing) PrimeVec(ixVegNrg) = scalarCanopyTempPrime

      ! extract intercepted water
      if(ixVegHyd/=integerMissing)then
         PrimeVec(ixVegHyd) = scalarCanopyHydPrime
      endif

    endif  ! not computing the vegetation flux

    ! *** extract state variables from the snow+soil sub-domain


    ! overwrite with the energy values from the state vector
    if(nSnowSoilNrg>0)then
      do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
         PrimeVec(ixSnowSoilNrg(iLayer)) = mLayerTempPrime(iLayer)
      end do  ! looping through non-missing energy state variables in the snow+soil domain
    endif


    ! overwrite with the hydrology values from the state vector
    if(nSnowSoilHyd>0)then
      do concurrent (iLayer=1:nLayers,ixSnowSoilHyd(iLayer)/=integerMissing)   ! (loop through non-missing hydrology state variables in the snow+soil domain)
        select case( ixHydType(iLayer) )
          case(iname_watLayer); PrimeVec(ixSnowSoilHyd(iLayer)) = mLayerVolFracHydPrime(iLayer)
          case(iname_liqLayer); PrimeVec(ixSnowSoilHyd(iLayer)) = mLayerVolFracHydPrime(iLayer)
          case(iname_matLayer); PrimeVec(ixSnowSoilHyd(iLayer)) = mLayerMatricHeadPrime(iLayer-nSnow)
          case(iname_lmpLayer)!; mLayerMatricHeadLiqTrial(iLayer-nSnow) = stateVec( ixSnowSoilHyd(iLayer) ) ! liquid matric potential state variable for soil layers
          case default ! do nothing
        end select
      end do  ! looping through non-missing energy state variables in the snow+soil domain
    endif


    ! extract temperature of the canopy air space
    if(ixAqWat/=integerMissing) PrimeVec(ixAqWat) = scalarAquiferStoragePrime


    
    
  end associate

  
end subroutine computePrime

end module computPrime_module
