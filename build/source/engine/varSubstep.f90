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

module varSubstep_module

! data types
USE nrtype
USE globalData,only: verySmall ! a very small number used as an additive constant to check if substantial difference among real numbers

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number
USE globalData,only:quadMissing     ! missing quadruple precision number

! access the global print flag
USE globalData,only:globalPrintFlag

! domain types
USE globalData,only:iname_cas       ! named variables for the canopy air space
USE globalData,only:iname_veg       ! named variables for vegetation
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! global metadata
USE globalData,only:flux_meta       ! metadata on the model fluxes

! derived types to define the data structures
USE data_types,only:&
                    var_i,              & ! data vector (i4b)
                    var_d,              & ! data vector (rkind)
                    var_flagVec,        & ! data vector with variable length dimension (i4b)
                    var_ilength,        & ! data vector with variable length dimension (i4b)
                    var_dlength,        & ! data vector with variable length dimension (rkind)
                    zLookup,            & ! lookup tables
                    model_options,      & ! defines the model decisions
                    in_type_varSubstep, & ! class for intent(in) arguments
                    io_type_varSubstep, & ! class for intent(inout) arguments
                    out_type_varSubstep   ! class for intent(out) arguments

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookDERIV      ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure

! look up structure for variable types
USE var_lookup,only:iLookVarType

! constants
USE multiconst,only:&
                    Tfreeze,        & ! freezing temperature                 (K)
                    LH_fus,         & ! latent heat of fusion                (J kg-1)
                    LH_vap,         & ! latent heat of vaporization          (J kg-1)
                    iden_ice,       & ! intrinsic density of ice             (kg m-3)
                    iden_water        ! intrinsic density of liquid water    (kg m-3)

! look-up values for the numerical method
USE mDecisions_module,only:         &
                    homegrown      ,& ! homegrown backward Euler solution using concepts from numerical recipes
                    kinsol         ,& ! SUNDIALS backward Euler solution using Kinsol
                    ida               ! SUNDIALS solution using IDA

! look-up values for the choice of variable in energy equations (BE residual or IDA state variable)
USE mDecisions_module,only:         &
                    closedForm,     & ! use temperature with closed form heat capacity
                    enthalpyFormLU, & ! use enthalpy with soil temperature-enthalpy lookup tables
                    enthalpyForm      ! use enthalpy with soil temperature-enthalpy analytical solution

! note: updateProg was relocated to updateVars_module for interoperability with the nested Newton solver
USE updateVars_module,only: updateProg ! update prognostic variables

! safety: set private unless specified otherwise
implicit none
private
public::varSubstep
contains


! **********************************************************************************************************
! public subroutine varSubstep: run the model for a collection of substeps for a given state subset
! **********************************************************************************************************
subroutine varSubstep(&
                      ! input: model control
                      in_varSubstep,     & ! intent(in)    : model control
                      io_varSubstep,     & ! intent(inout) : model control
                      ! input/output: data structures
                      model_decisions,   & ! intent(in)    : model decisions
                      lookup_data,       & ! intent(in)    : lookup tables
                      type_data,         & ! intent(in)    : type of vegetation and soil
                      attr_data,         & ! intent(in)    : spatial attributes
                      forc_data,         & ! intent(in)    : model forcing data
                      mpar_data,         & ! intent(in)    : model parameters
                      indx_data,         & ! intent(inout) : index data
                      prog_data,         & ! intent(inout) : model prognostic variables for a local HRU
                      diag_data,         & ! intent(inout) : model diagnostic variables for a local HRU
                      flux_data,         & ! intent(inout) : model fluxes for a local HRU
                      flux_mean,         & ! intent(inout) : mean model fluxes for a local HRU
                      deriv_data,        & ! intent(inout) : derivatives in model fluxes w.r.t. relevant state variables
                      bvar_data,         & ! intent(in)    : model variables for the local basin
                      ! output: model control
                      out_varSubstep)      ! intent(out)   : model control
  ! ---------------------------------------------------------------------------------------
  ! structure allocations
  USE allocspace_module,only:allocLocal                ! allocate local data structures
  ! simulation of fluxes and residuals given a trial state vector
  USE getVectorz_module,only:popStateVec                ! populate the state vector
  USE getVectorz_module,only:varExtract                 ! extract variables from the state vector
  USE systemSolv_module,only:systemSolv                 ! solve the system of equations for one time step
  ! identify name of variable type (for error message)
  USE get_ixName_module,only:get_varTypeName           ! to access type strings for error messages
  implicit none
  ! ---------------------------------------------------------------------------------------
  ! * dummy variables
  ! ---------------------------------------------------------------------------------------
  ! input: model control
  type(in_type_varSubstep),intent(in)    :: in_varSubstep             ! model control
  type(io_type_varSubstep),intent(inout) :: io_varSubstep             ! model control
  ! input/output: data structures
  type(model_options),intent(in)         :: model_decisions(:)        ! model decisions
  type(zLookup),intent(in)               :: lookup_data               ! lookup tables
  type(var_i),intent(in)                 :: type_data                 ! type of vegetation and soil
  type(var_d),intent(in)                 :: attr_data                 ! spatial attributes
  type(var_d),intent(in)                 :: forc_data                 ! model forcing data
  type(var_dlength),intent(in)           :: mpar_data                 ! model parameters
  type(var_ilength),intent(inout)        :: indx_data                 ! indices for a local HRU
  type(var_dlength),intent(inout)        :: prog_data                 ! prognostic variables for a local HRU
  type(var_dlength),intent(inout)        :: diag_data                 ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout)        :: flux_data                 ! model fluxes for a local HRU
  type(var_dlength),intent(inout)        :: flux_mean                 ! mean model fluxes for a local HRU
  type(var_dlength),intent(inout)        :: deriv_data                ! derivatives in model fluxes w.r.t. relevant state variables
  type(var_dlength),intent(in)           :: bvar_data                 ! model variables for the local basin
  ! output: model control
  type(out_type_varSubstep),intent(out)  :: out_varSubstep            ! model control
  ! ---------------------------------------------------------------------------------------
  ! * general local variables
  ! ---------------------------------------------------------------------------------------
  ! error control
  character(LEN=256)                 :: cmessage                               ! error message of downwind routine
  ! general local variables
  integer(i4b)                       :: nLayers_in                             ! number of layers for input to systemSolv
  integer(i4b)                       :: iVar                                   ! index of variables in data structures
  integer(i4b)                       :: iSoil                                  ! index of soil layers
  integer(i4b)                       :: ixLayer                                ! index in a given domain
  integer(i4b),dimension(1)          :: ixMin,ixMax                            ! bounds of a given flux vector
  ! time stepping
  real(rkind)                        :: dtSum                                  ! sum of time from successful steps (seconds)
  real(rkind)                        :: dt_wght                                ! weight given to a given flux calculation
  real(rkind)                        :: dtSubstep                              ! length of a substep (s)
  real(rkind)                        :: maxstep                                ! maximum time step length (seconds)
  integer(i4b)                       :: nSteps                                 ! number of time steps taken in solver
  ! adaptive sub-stepping for the solution
  logical(lgt)                       :: failedSubstep                          ! flag to denote success of substepping for a given split
  integer(i4b)                       :: niter                                  ! number of iterations taken
  integer(i4b),parameter             :: n_inc=5                                ! minimum number of iterations to increase time step
  integer(i4b),parameter             :: n_dec=15                               ! maximum number of iterations to decrease time step
  real(rkind),parameter              :: F_inc = 1.25_rkind                     ! factor used to increase time step
  real(rkind),parameter              :: F_dec = 0.90_rkind                     ! factor used to decrease time step
  ! state and flux vectors (Note: nstate = in_varSubstep % nSubset)
  real(rkind)                        :: untappedMelt(in_varSubstep % nSubset)  ! un-tapped melt energy (J m-3 s-1)
  real(rkind)                        :: stateVecInit(in_varSubstep % nSubset)  ! initial state vector (mixed units)
  real(rkind)                        :: stateVecTrial(in_varSubstep % nSubset) ! trial state vector (mixed units)
  real(rkind)                        :: stateVecPrime(in_varSubstep % nSubset) ! trial state vector (mixed units)
  type(var_dlength)                  :: flux_temp                              ! temporary model fluxes
  ! flags
  logical(lgt)                       :: firstSplitOper                         ! flag to indicate if we are processing the first flux call in a splitting operation
  logical(lgt)                       :: waterBalanceError                      ! flag to denote that there is a water balance error
  logical(lgt)                       :: nrgFluxModified                        ! flag to denote that the energy fluxes were modified
  ! energy fluxes
  real(rkind)                        :: sumCanopyEvaporation                   ! sum of canopy evaporation/condensation (kg m-2 s-1)
  real(rkind)                        :: sumLatHeatCanopyEvap                   ! sum of latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
  real(rkind)                        :: sumSenHeatCanopy                       ! sum of sensible heat flux from the canopy to the canopy air space (W m-2)
  real(rkind)                        :: sumSoilCompress                        ! sum of total soil compression
  real(rkind),allocatable            :: sumLayerCompress(:)                    ! sum of soil compression by layer
  ! balances and residual vectors
  real(rkind)                        :: fluxVec(in_varSubstep % nSubset)       ! substep flux vector (mixed units)
  real(rkind)                        :: resSink(in_varSubstep % nSubset)       ! substep sink terms on the RHS of the state equation
  real(qp)                           :: resVec(in_varSubstep % nSubset)        ! substep residual vector
  real(rkind)                        :: balance(in_varSubstep % nSubset)       ! substep balance per second
  real(rkind)                        :: sumBalance(in_varSubstep % nSubset)    ! sum of substeps balance
  logical(lgt),parameter             :: computMassBalance = .true.             ! flag to compute the mass balance, will affect step length, default true
  logical(lgt),parameter             :: computNrgBalance = .true.              ! flag to compute the energy balance, will not effect solution but will not compute energy balance if false (saves expense)
  logical(lgt)                       :: computeEnthTemp                        ! flag to compute enthalpy regardless of the model decision
  logical(lgt)                       :: enthalpyStateVec                       ! flag if enthalpy is a state variable (ida)
  logical(lgt)                       :: use_lookup                             ! flag to use the lookup table for soil enthalpy, otherwise use analytical solution
  ! test variables for nested Newton -- SJT: to be removed or retained (if needed) in a future update
  logical(lgt),parameter :: nested_Newton_test=.false. ! test output

  ! ---------------------------------------------------------------------------------------
  ! initialize error control
  out_varSubstep % err=0; out_varSubstep % cmessage='varSubstep/'
  ! ---------------------------------------------------------------------------------------
  ! point to variables in the data structures
  ! ---------------------------------------------------------------------------------------
  globalVars: associate(&
    ! input: model control
    dt             => in_varSubstep % dt,             & ! intent(in): time step (seconds)
    dtInit         => in_varSubstep % dtInit,         & ! intent(in): initial time step (seconds)
    dt_min         => in_varSubstep % dt_min,         & ! intent(in): minimum time step (seconds)
    whole_step     => in_varSubstep % whole_step,     & ! intent(in): length of whole step for surface drainage and average flux
    nState         => in_varSubstep % nSubset,        & ! intent(in): total number of state variables
    doAdjustTemp   => in_varSubstep % doAdjustTemp,   & ! intent(in): flag to indicate if we adjust the temperature
    firstSubStep   => in_varSubstep % firstSubStep,   & ! intent(in): flag to indicate if processing the first sub-step
    computeVegFlux => in_varSubstep % computeVegFlux, & ! intent(in): flag to indicate if computing fluxes over vegetation (.false. means veg is buried with snow)
    scalarSolution => in_varSubstep % scalarSolution, & ! intent(in): flag to denote implementing the scalar solution
    iStateSplit    => in_varSubstep % iStateSplit,    & ! intent(in): index of the state in the splitting operation
    fluxMask       => in_varSubstep % fluxMask,       & ! intent(in): flags to denote if the flux is calculated in the given state subset
    firstFluxCall  => io_varSubstep % firstFluxCall,  & ! intent(inout): flag to define the first flux call
    fluxCount      => io_varSubstep % fluxCount,      & ! intent(inout): number of times that the flux is updated (should equal nSubsteps)
    ixSaturation   => io_varSubstep % ixSaturation,   & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
    ! model decisions
    ixNumericalMethod       => model_decisions(iLookDECISIONS%num_method)%iDecision   ,& ! intent(in):    [i4b]    choice of numerical solver
    ixNrgConserv            => model_decisions(iLookDECISIONS%nrgConserv)%iDecision   ,& ! intent(in):    [i4b]    choice of variable in either energy backward Euler residual or IDA state variable
    ! number of layers
    mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,& ! intent(in):    [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
    ! model state variables (vegetation canopy)
    scalarCanairTemp        => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)       ,& ! intent(inout): [dp]     temperature of the canopy air space (K)
    scalarCanopyTemp        => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,& ! intent(inout): [dp]     temperature of the vegetation canopy (K)
    scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)        ,& ! intent(inout): [dp]     mass of ice on the vegetation canopy (kg m-2)
    scalarCanopyLiq         => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)        ,& ! intent(inout): [dp]     mass of liquid water on the vegetation canopy (kg m-2)
    scalarCanopyWat         => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)        ,& ! intent(inout): [dp]     mass of total water on the vegetation canopy (kg m-2)
    ! model state variables (snow and soil domains)
    mLayerTemp              => prog_data%var(iLookPROG%mLayerTemp)%dat                ,& ! intent(inout): [dp(:)]  temperature of each snow/soil layer (K)
    mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat          ,& ! intent(inout): [dp(:)]  volumetric fraction of ice (-)
    mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,& ! intent(inout): [dp(:)]  volumetric fraction of liquid water (-)
    mLayerVolFracWat        => prog_data%var(iLookPROG%mLayerVolFracWat)%dat          ,& ! intent(inout): [dp(:)]  volumetric fraction of total water (-)
    mLayerMatricHead        => prog_data%var(iLookPROG%mLayerMatricHead)%dat          ,& ! intent(inout): [dp(:)]  matric head (m)
    mLayerMatricHeadLiq     => diag_data%var(iLookDIAG%mLayerMatricHeadLiq)%dat       ,& ! intent(inout): [dp(:)]  matric potential of liquid water (m)
    ! model control
    dtMultiplier      => out_varSubstep % dtMultiplier             ,& ! intent(out): substep multiplier (-)
    nSubsteps         => out_varSubstep % nSubsteps                ,& ! intent(out): number of substeps taken for a given split
    failedMinimumStep => out_varSubstep % failedMinimumStep        ,& ! intent(out): flag to denote success of substepping for a given split
    reduceCoupledStep => out_varSubstep % reduceCoupledStep        ,& ! intent(out): flag to denote need to reduce the length of the coupled step
    tooMuchMelt       => out_varSubstep % tooMuchMelt              ,& ! intent(out): flag to denote that ice is insufficient to support melt
    err               => out_varSubstep % err                      ,& ! intent(out): error code
    message           => out_varSubstep % cmessage                  & ! intent(out): error message
    )  ! end association with variables in the data structures
    ! *********************************************************************************************************************************************************


    ! initialize flag for the success of the substepping
    failedMinimumStep=.false.

    ! set the flag to compute enthalpy, may want to have this true always if want to output enthalpy
    computeEnthTemp  = .false.
    enthalpyStateVec = .false.
    use_lookup       = .false.
    if((ixNrgConserv .ne. closedForm .or. computNrgBalance) .and. ixNumericalMethod .ne. ida) computeEnthTemp = .true. ! use enthTemp to conserve energy or compute energy balance
    if(ixNrgConserv .ne. closedForm .and. ixNumericalMethod==ida) enthalpyStateVec = .true. ! enthalpy as state variable
    if(ixNrgConserv==enthalpyFormLU) use_lookup = .true. ! use lookup tables for soil enthalpy instead of analytical solution

    ! initialize the length of the substep
    dtSubstep = dtInit

    ! change maxstep with hard code here to make only the newton step loop in systemSolv* happen more frequently
    !   NOTE: this may just be amplifying the splitting error if maxstep is smaller than the full possible step
    maxstep = mpar_data%var(iLookPARAM%maxstep)%dat(1)  ! maximum time step (s).

    ! associate block for indx_data components
    associate(&
     nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in):    [i4b]    number of snow layers
     nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)                  & ! intent(in):    [i4b]    number of soil layers
    &)
      ! allocate space for the temporary model flux structure
      call allocLocal(flux_meta(:),flux_temp,nSnow,nSoil,err,cmessage)
      if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

      ! initialize the model fluxes (some model fluxes are not computed in the iterations)
      do iVar=1,size(flux_data%var)
        flux_temp%var(iVar)%dat(:) = flux_data%var(iVar)%dat(:)
      end do

      ! initialize the total energy fluxes (modified in updateProg)
      sumCanopyEvaporation = 0._rkind  ! canopy evaporation/condensation (kg m-2 s-1)
      sumLatHeatCanopyEvap = 0._rkind  ! latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
      sumSenHeatCanopy     = 0._rkind  ! sensible heat flux from the canopy to the canopy air space (W m-2)
      sumSoilCompress      = 0._rkind  ! total soil compression
      allocate(sumLayerCompress(nSoil)); sumLayerCompress = 0._rkind ! soil compression by layer
    end associate

    ! initialize balances
    sumBalance = 0._rkind

    ! define the first flux call in a splitting operation
    firstSplitOper = (.not.scalarSolution .or. iStateSplit==1)

    ! initialize subStep
    dtSum     = 0._rkind  ! keep track of the portion of the time step that is completed
    nSubsteps = 0


    ! loop through substeps
    ! NOTE: continuous do statement with exit clause
    substeps: do
      dtSubstep = min(dtSubstep,maxstep)

      ! -----
      ! * populate state vectors...
      ! ---------------------------

      ! initialize state vectors
      call popStateVec(&
                      ! input
                      nState,           & ! intent(in):  number of desired state variables
                      enthalpyStateVec, & ! intent(in):  flag to use enthalpy as a state variable
                      prog_data,        & ! intent(in):  model prognostic variables for a local HRU
                      diag_data,        & ! intent(in):  model diagnostic variables for a local HRU
                      indx_data,        & ! intent(in):  indices defining model states and layers
                      ! output
                      stateVecInit,     & ! intent(out): initial model state vector (mixed units)
                      err,cmessage)       ! intent(out): error control
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

      ! # of layers to be used as input to systemSolv
      ! note: using a separate variable rather than an association to resolve conflicting intent attributes for indx_data
      nLayers_in=indx_data%var(iLookINDEX%nLayers)%dat(1)
      ! -----
      ! * iterative solution...
      ! -----------------------
      ! solve the system of equations for a given state subset
      call systemSolv(&
                      ! input: model control
                      dtSubstep,         & ! intent(in):    time step (s)
                      whole_step,        & ! intent(in):    entire time step (s)
                      nState,            & ! intent(in):    total number of state variables
                      nLayers_in,        & ! intent(in):    total number of layers
                      firstSubStep,      & ! intent(in):    flag to denote first sub-step
                      firstFluxCall,     & ! intent(inout): flag to indicate if we are processing the first flux call
                      firstSplitOper,    & ! intent(in):    flag to indicate if we are processing the first flux call in a splitting operation
                      computeVegFlux,    & ! intent(in):    flag to denote if computing energy flux over vegetation
                      scalarSolution,    & ! intent(in):    flag to denote if implementing the scalar solution
                      computMassBalance, & ! intent(in):    flag to compute mass balance
                      computNrgBalance,  & ! intent(in):    flag to compute energy balance
                      ! input/output: data structures
                      lookup_data,       & ! intent(in):    lookup tables
                      type_data,         & ! intent(in):    type of vegetation and soil
                      attr_data,         & ! intent(in):    spatial attributes
                      forc_data,         & ! intent(in):    model forcing data
                      mpar_data,         & ! intent(in):    model parameters
                      indx_data,         & ! intent(inout): index data
                      prog_data,         & ! intent(inout): model prognostic variables for a local HRU
                      diag_data,         & ! intent(inout): model diagnostic variables for a local HRU
                      flux_temp,         & ! intent(inout): model fluxes for a local HRU
                      bvar_data,         & ! intent(in):    model variables for the local basin
                      model_decisions,   & ! intent(in):    model decisions
                      stateVecInit,      & ! intent(in):    initial state vector
                      ! output: model control
                      deriv_data,        & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                      ixSaturation,      & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                      stateVecTrial,     & ! intent(out):   updated state vector
                      stateVecPrime,     & ! intent(out):   updated state vector if need the prime space (ida)
                      fluxVec,           & ! intent(out):   model flux vector
                      resSink,           & ! intent(out):   additional (sink) terms on the RHS of the state equation
                      resVec,            & ! intent(out):   residual vector
                      untappedMelt,      & ! intent(out):   un-tapped melt energy (J m-3 s-1)
                      ! output: balances (only computed at this level for ida)
                      balance,           & ! intent(out):   balance per state variable
                      ! output  model control
                      niter,             & ! intent(out):   number of iterations taken (homegrown solver)
                      nSteps,            & ! intent(out):   number of time steps taken in solver
                      reduceCoupledStep, & ! intent(out):   flag to reduce the length of the coupled step
                      tooMuchMelt,       & ! intent(out):   flag to denote that ice is insufficient to support melt
                      err,cmessage)        ! intent(out):   error code and error message
      if (nested_Newton_test) then
       print *, "varSubstep Test A: after systemSolv call"
       print *, "err=",err
       print *, "reduceCoupledStep=",reduceCoupledStep
       print *, "tooMuchMelt=",tooMuchMelt
       print *, "indx_data nSoil=",indx_data%var(iLookINDEX%nSoil)%dat(1)
       print *, "indx_data nLayers=",indx_data%var(iLookINDEX%nLayers)%dat(1)
       print *, ""
      end if

      if(err/=0)then ! (check for errors, but do not fail yet)
        message=trim(message)//trim(cmessage)
        if(err>0) return
      endif
 
      ! if too much melt or need to reduce length of the coupled step then return
      ! NOTE: need to go all the way back to coupled_em and merge snow layers, as all splitting operations need to occur with the same layer geometry
      if(tooMuchMelt .or. reduceCoupledStep)then 
        deallocate(sumLayerCompress)
        return
      endif

      ! identify failure
      failedSubstep = (err<0)

      ! check
      if(globalPrintFlag)then
        print*, 'niter, failedSubstep, dtSubstep = ', niter, failedSubstep, dtSubstep
        print*, trim(cmessage)
      endif

      ! reduce step based on failure
      if(failedSubstep)then
        err=0; message='varSubstep/'  ! recover from failed convergence
        dtMultiplier  = 0.5_rkind        ! system failure: step halving
      else
        ! ** implicit Euler: adjust step length based on iteration count
        if(niter<n_inc)then
          dtMultiplier = F_inc
        elseif(niter>n_dec)then
          dtMultiplier = F_dec
        else
          dtMultiplier = 1._rkind
        endif
      endif  ! switch between failure and success

      ! check if we failed the substep
      if(failedSubstep)then

        ! check that the substep is greater than the minimum step
        if(dtSubstep*dtMultiplier<dt_min)then
          ! --> exit, and either (1) try another solution method; or (2) reduce coupled step
          failedMinimumStep=.true.
          exit subSteps

        else ! step is still OK
          dtSubstep = dtSubstep*dtMultiplier
          cycle subSteps
        endif  ! if step is less than the minimum

      endif  ! if failed the substep

      ! -----
      ! * update model fluxes...
      ! ------------------------

      ! NOTE: if we get to here then we are accepting the step of dtSubstep
      if(err/=0)then
        message=trim(message)//'expect err=0 if updating fluxes'
        return
      endif

      ! associate block for indx_data components
      associate(&
       nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in):    [i4b]    number of snow layers
       nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in):    [i4b]    number of soil layers
       nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)                & ! intent(in):    [i4b]    total number of layers
      &)
        if (nested_Newton_test) then
         print *, "varSubstep Test B: before updateProg call"
         print *, "nSnow=",nSnow
         print *, "nSoil=",nSoil
         print *, "indx_data nSoil=",indx_data%var(iLookINDEX%nSoil)%dat(1)
         print *, "nLayers=",nLayers
         print *, ""
        end if
        ! update prognostic variables, update balances, and check them for possible step reduction if homegrown or kinsol solver
        call updateProg(dtSubstep,nSnow,nSoil,nLayers,untappedMelt,stateVecTrial,stateVecPrime,                                    & ! input: states
                        doAdjustTemp,computeVegFlux,computMassBalance,computNrgBalance,computeEnthTemp,enthalpyStateVec,use_lookup,& ! input: model control
                        model_decisions,lookup_data,mpar_data,indx_data,flux_temp,prog_data,diag_data,deriv_data,                  & ! input-output: data structures
                        fluxVec,resVec,balance,waterBalanceError,nrgFluxModified,err,message)                                        ! input-output: balances, flags, and error control
        if (nested_Newton_test) then
         print *, "varSubstep Test C: after updateProg call"
         print *, "err=",err
         print *, "waterBalanceError=",waterBalanceError
         print *, ""
        end if
      end associate

      ! check for errors -- return if non-recoverable
      if(err/=0)then
        message=trim(message)//trim(cmessage)
        if(err>0) return
      endif

      ! if water balance error then reduce the length of the coupled step
      if(waterBalanceError)then
        message=trim(message)//'water balance error'
        reduceCoupledStep=.true.
        deallocate(sumLayerCompress)
        err=-20; return
      endif

      if(globalPrintFlag)&
      print*, trim(cmessage)//': dt = ', dtSubstep

      ! recover from errors in prognostic update
      if(err<0)then

        ! modify step
        err=0  ! error recovery
        dtSubstep = dtSubstep/2._rkind 

        ! check minimum: fail minimum step if there is an error in the update
        if(dtSubstep<dt_min)then
          failedMinimumStep=.true.
          exit subSteps
        ! minimum OK -- try again
        else
          cycle substeps
        endif

      endif  ! if errors in prognostic update

      ! associate block for indx_data components
      associate(&
       nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in):    [i4b]    number of soil layers
       nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)               ,& ! intent(in):    [i4b]    total number of layers
       ! get indices for balances
       ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy air space energy state variable
       ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy energy state variable
       ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy hydrology state variable (mass)
       ixAqWat                 => indx_data%var(iLookINDEX%ixAqWat)%dat(1)               ,& ! intent(in):    [i4b]    index of water storage in the aquifer
       ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat            ,& ! intent(in):    [i4b(:)] index in the state subset for hydrology state variables in the soil domain
       ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in):    [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
       ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in):    [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
       nSnowSoilNrg            => indx_data%var(iLookINDEX%nSnowSoilNrg)%dat(1)          ,& ! intent(in):    [i4b]    number of energy state variables in the snow+soil domain
       nSnowSoilHyd            => indx_data%var(iLookINDEX%nSnowSoilHyd)%dat(1)          ,& ! intent(in):    [i4b]    number of hydrology state variables in the snow+soil domain
       ! mapping between state vectors and control volumes
       ixLayerActive           => indx_data%var(iLookINDEX%ixLayerActive)%dat             & ! intent(in):    [i4b(:)] list of indices for all active layers (inactive=integerMissing)
      &)
      
        ! add balances to the total balances
        if(ixCasNrg/=integerMissing) sumBalance(ixCasNrg) = sumBalance(ixCasNrg) + dtSubstep*balance(ixCasNrg)
        if(ixVegNrg/=integerMissing) sumBalance(ixVegNrg) = sumBalance(ixVegNrg) + dtSubstep*balance(ixVegNrg)
        if(nSnowSoilNrg>0) then
          do concurrent (ixLayer=1:nLayers,ixSnowSoilNrg(ixLayer)/=integerMissing)
            if(ixSnowSoilNrg(ixLayer)/=integerMissing) sumBalance(ixSnowSoilNrg(ixLayer)) = sumBalance(ixSnowSoilNrg(ixLayer)) + dtSubstep*balance(ixSnowSoilNrg(ixLayer))
          end do
        endif
        if(ixVegHyd/=integerMissing) sumBalance(ixVegHyd) = sumBalance(ixVegHyd) + dtSubstep*balance(ixVegHyd)
        if(nSnowSoilHyd>0) then
          do concurrent (ixLayer=1:nLayers,ixSnowSoilHyd(ixLayer)/=integerMissing)
            if(ixSnowSoilHyd(ixLayer)/=integerMissing) sumBalance(ixSnowSoilHyd(ixLayer)) = sumBalance(ixSnowSoilHyd(ixLayer)) + dtSubstep*balance(ixSnowSoilHyd(ixLayer))
          end do
        endif
        if(ixAqWat/=integerMissing) sumBalance(ixAqWat) = sumBalance(ixAqWat) + dtSubstep*balance(ixAqWat)

        ! get the total energy fluxes (modified in updateProg), have to do differently
        if(nrgFluxModified .or. ixVegNrg/=integerMissing)then
          sumCanopyEvaporation  = sumCanopyEvaporation  + dtSubstep*flux_temp%var(iLookFLUX%scalarCanopyEvaporation)%dat(1)  ! canopy evaporation/condensation (kg m-2 s-1)
          sumLatHeatCanopyEvap  = sumLatHeatCanopyEvap  + dtSubstep*flux_temp%var(iLookFLUX%scalarLatHeatCanopyEvap)%dat(1)  ! latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
          sumSenHeatCanopy      = sumSenHeatCanopy      + dtSubstep*flux_temp%var(iLookFLUX%scalarSenHeatCanopy)%dat(1)      ! sensible heat flux from the canopy to the canopy air space (W m-2)
        else
          sumCanopyEvaporation  = sumCanopyEvaporation  + dtSubstep*flux_data%var(iLookFLUX%scalarCanopyEvaporation)%dat(1)  ! canopy evaporation/condensation (kg m-2 s-1)
          sumLatHeatCanopyEvap  = sumLatHeatCanopyEvap  + dtSubstep*flux_data%var(iLookFLUX%scalarLatHeatCanopyEvap)%dat(1)  ! latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
          sumSenHeatCanopy      = sumSenHeatCanopy      + dtSubstep*flux_data%var(iLookFLUX%scalarSenHeatCanopy)%dat(1)      ! sensible heat flux from the canopy to the canopy air space (W m-2)
        endif  ! if energy fluxes were modified

        ! get the total soil compression
        if (count(ixSoilOnlyHyd/=integerMissing)>0) then
          ! scalar compression
          if(.not.scalarSolution .or. iStateSplit==nSoil)&
          sumSoilCompress = sumSoilCompress + dtSubstep*diag_data%var(iLookDIAG%scalarSoilCompress)%dat(1) ! total soil compression
          ! vector compression
          do iSoil=1,nSoil
            if(ixSoilOnlyHyd(iSoil)/=integerMissing)&
            sumLayerCompress(iSoil) = sumLayerCompress(iSoil) + dtSubstep*diag_data%var(iLookDIAG%mLayerCompress)%dat(iSoil) ! soil compression in layers
          end do
        endif

        ! print progress
        if(globalPrintFlag)&
        write(*,'(a,1x,3(f13.2,1x))') 'updating: dtSubstep, dtSum, dt = ', dtSubstep, dtSum, dt

       ! increment fluxes
        dt_wght = dtSubstep/dt ! define weight applied to each sub-step
        do iVar=1,size(flux_meta)
          if(count(fluxMask%var(iVar)%dat)>0) then

            ! ** no domain splitting
            if(count(ixLayerActive/=integerMissing)==nLayers)then
              flux_mean%var(iVar)%dat(:) = flux_mean%var(iVar)%dat(:) + flux_temp%var(iVar)%dat(:)*dt_wght
              fluxCount%var(iVar)%dat(:) = fluxCount%var(iVar)%dat(:) + 1

            ! ** domain splitting
            else
              ixMin=lbound(flux_data%var(iVar)%dat)
              ixMax=ubound(flux_data%var(iVar)%dat)
              do ixLayer=ixMin(1),ixMax(1)
                if(fluxMask%var(iVar)%dat(ixLayer)) then
                  ! special case of the transpiration sink from soil layers: only computed for the top soil layer
                  if(iVar==iLookFLUX%mLayerTranspire)then
                    if(ixLayer==1) flux_mean%var(iVar)%dat(:) = flux_mean%var(iVar)%dat(:) + flux_temp%var(iVar)%dat(:)*dt_wght
                  ! standard case
                  else
                    flux_mean%var(iVar)%dat(ixLayer) = flux_mean%var(iVar)%dat(ixLayer) + flux_temp%var(iVar)%dat(ixLayer)*dt_wght
                  endif
                  fluxCount%var(iVar)%dat(ixLayer) = fluxCount%var(iVar)%dat(ixLayer) + 1
                endif
              end do
            endif  ! (domain splitting)

          endif   ! (if the flux is desired)
        end do  ! (loop through fluxes)

      end associate

      ! increment the number of substeps
      nSubsteps = nSubsteps + nSteps

      ! increment the sub-step legth
      dtSum = dtSum + dtSubstep

      ! check that we have completed the sub-step
      if(dtSum >= dt-verySmall)then
        failedMinimumStep=.false.
        exit subSteps
      endif

      ! adjust length of the sub-step (make sure that we don't exceed the step)
      dtSubstep = min(dt - dtSum, max(dtSubstep*dtMultiplier, dt_min) )

    end do substeps  ! time steps for variable-dependent sub-stepping
    ! NOTE: if we get to here then we are accepting then dtSum should dt

    ! save the fluxes as averages
    do iVar=1,size(flux_meta)
      if(count(fluxMask%var(iVar)%dat)>0) flux_data%var(iVar)%dat(:) = flux_mean%var(iVar)%dat(:)
    enddo

    ! save the energy fluxes as averages
    flux_data%var(iLookFLUX%scalarCanopyEvaporation)%dat(1) = sumCanopyEvaporation /dt      ! canopy evaporation/condensation (kg m-2 s-1)
    flux_data%var(iLookFLUX%scalarLatHeatCanopyEvap)%dat(1) = sumLatHeatCanopyEvap /dt      ! latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
    flux_data%var(iLookFLUX%scalarSenHeatCanopy)%dat(1)     = sumSenHeatCanopy     /dt      ! sensible heat flux from the canopy to the canopy air space (W m-2)

    ! associate block for indx_data components
    associate(&
     nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in):    [i4b]    number of soil layers
     nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)               ,& ! intent(in):    [i4b]    total number of layers
     ! get indices for balances
     ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy air space energy state variable
     ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy energy state variable
     ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy hydrology state variable (mass)
     ixAqWat                 => indx_data%var(iLookINDEX%ixAqWat)%dat(1)               ,& ! intent(in):    [i4b]    index of water storage in the aquifer
     ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat            ,& ! intent(in):    [i4b(:)] index in the state subset for hydrology state variables in the soil domain
     ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in):    [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
     ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in):    [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
     nSnowSoilNrg            => indx_data%var(iLookINDEX%nSnowSoilNrg)%dat(1)          ,& ! intent(in):    [i4b]    number of energy state variables in the snow+soil domain
     nSnowSoilHyd            => indx_data%var(iLookINDEX%nSnowSoilHyd)%dat(1)           & ! intent(in):    [i4b]    number of hydrology state variables in the snow+soil domain
    &)

      ! save the soil compression diagnostics as averages
      diag_data%var(iLookDIAG%scalarSoilCompress)%dat(1) = sumSoilCompress/dt
      do iSoil=1,nSoil
        if(ixSoilOnlyHyd(iSoil)/=integerMissing)&
        diag_data%var(iLookDIAG%mLayerCompress)%dat(iSoil) = sumLayerCompress(iSoil)/dt
      end do
      deallocate(sumLayerCompress)

      ! save the balance diagnostics as averages
      if(ixCasNrg/=integerMissing) diag_data%var(iLookDIAG%balanceCasNrg)%dat(1) = sumBalance(ixCasNrg)/dt
      if(ixVegNrg/=integerMissing) diag_data%var(iLookDIAG%balanceVegNrg)%dat(1) = sumBalance(ixVegNrg)/dt
      if(nSnowSoilNrg>0) then
        do concurrent (ixLayer=1:nLayers,ixSnowSoilNrg(ixLayer)/=integerMissing)
          diag_data%var(iLookDIAG%balanceLayerNrg)%dat(ixLayer) = sumBalance(ixSnowSoilNrg(ixLayer))/dt
        end do
      endif
      if(ixVegHyd/=integerMissing) diag_data%var(iLookDIAG%balanceVegMass)%dat(1) = sumBalance(ixVegHyd)/dt
      if(nSnowSoilHyd>0) then
        do concurrent (ixLayer=1:nLayers,ixSnowSoilHyd(ixLayer)/=integerMissing)
          diag_data%var(iLookDIAG%balanceLayerMass)%dat(ixLayer) = sumBalance(ixSnowSoilHyd(ixLayer))/dt
        end do
      endif 
      if(ixAqWat/=integerMissing) diag_data%var(iLookDIAG%balanceAqMass)%dat(1) = sumBalance(ixAqWat)/dt

    end associate

    ! update error codes
    if (failedMinimumStep) then
      err=-20 ! negative = recoverable error
      message=trim(message)//'failed minimum step'
    end if
  ! end associate statements
  end associate globalVars
end subroutine varSubstep

end module varSubstep_module
