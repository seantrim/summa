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

module stateFilter_module

! data types
USE nrtype

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing   ! missing integer

! state variable type
USE globalData,only:iname_nrgCanair  ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy  ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy  ! named variable defining the mass of total water on the vegetation canopy
USE globalData,only:iname_liqCanopy  ! named variable defining the mass of liquid water on the vegetation canopy
USE globalData,only:iname_nrgLayer   ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer   ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer   ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer   ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer   ! named variable defining the liquid matric potential state variable for soil layers
USE globalData,only:iname_watAquifer ! named variable defining the water storage in the aquifer
  
! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookINDEX       ! named variables for structure elements

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_flagVec,                                               & ! data vector with variable length dimension (lgt)
                    var_ilength,                                               & ! data vector with variable length dimension (i4b)
                    out_type_stateFilter                                         ! classes for stateFilter objects


! safety: set private unless specified otherwise
implicit none
private

! named variables for the coupling method
integer(i4b),parameter,public  :: fullyCoupled=1   ! 1st try: fully coupled solution
integer(i4b),parameter,public  :: stateTypeSplit=2 ! 2nd try: separate solutions for each state type

! named variables for the state variable split
integer(i4b),parameter,public  :: nrgSplit=1       ! order in sequence for the energy operation
integer(i4b),parameter,public  :: massSplit=2      ! order in sequence for the mass operation

! named variables for the domain type split
integer(i4b),parameter,public  :: vegSplit=1       ! order in sequence for the vegetation split
integer(i4b),parameter,public  :: snowSplit=2      ! order in sequence for the snow split
integer(i4b),parameter,public  :: soilSplit=3      ! order in sequence for the soil split
integer(i4b),parameter,public  :: aquiferSplit=4   ! order in sequence for the aquifer split

! named variables for the solution method
integer(i4b),parameter,public  :: vector=1         ! vector solution method
integer(i4b),parameter,public  :: scalar=2         ! scalar solution method
integer(i4b),parameter,public  :: nSolutions=2     ! number of solution methods

! named variables for the switch between states and domains
integer(i4b),parameter,public  :: fullDomain=1     ! full domain (veg+snow+soil)
integer(i4b),parameter,public  :: subDomain=2      ! sub domain (veg, snow, soil, and aquifer separately)

! maximum number of possible splits
integer(i4b),parameter,public  :: nStateTypes=2    ! number of state types (energy, water)
integer(i4b),parameter,public  :: nDomains=4       ! number of domains (vegetation, snow, soil, and aquifer)

! class definitions

type, public :: split_select_type  ! class for selecting operator splitting methods
  ! opSplittin indices (in order)
  integer(i4b)             :: iSplit                      ! iteration counter for split_select_loop
  integer(i4b)             :: ixCoupling
  integer(i4b)             :: iStateTypeSplit
  integer(i4b)             :: ixStateThenDomain           ! 1=state type split; 2=domain split within a given state type 
  integer(i4b)             :: iDomainSplit
  integer(i4b)             :: ixSolution
  integer(i4b)             :: iStateSplit
  ! variables for specifying the split
  integer(i4b)             :: nState                      ! # of state variables
  integer(i4b)             :: nSubset                     ! number of selected state variables for a given split
  type(var_flagVec)        :: fluxMask                    ! integer mask defining model fluxes
  logical(lgt),allocatable :: stateMask(:)                ! mask defining desired state variables
  ! flags for splitting method control
  logical(lgt)             :: stateTypeSplitting,stateThenDomain,domainSplit,solution,stateSplit
 contains
  procedure :: initialize_flags             => split_select_initialize_flags             ! initialize flags that control operations
  procedure :: initialize_ixCoupling        => split_select_initialize_ixCoupling        ! initialize operator splitting indices
  procedure :: initialize_iStateTypeSplit   => split_select_initialize_iStateTypeSplit   ! initialize operator splitting indices
  procedure :: initialize_ixStateThenDomain => split_select_initialize_ixStateThenDomain ! initialize operator splitting indices
  procedure :: initialize_iDomainSplit      => split_select_initialize_iDomainSplit      ! initialize operator splitting indices
  procedure :: initialize_ixSolution        => split_select_initialize_ixSolution        ! initialize operator splitting indices
  procedure :: initialize_iStateSplit       => split_select_initialize_iStateSplit       ! initialize operator splitting indices

  procedure :: get_stateMask                => split_select_compute_stateMask            ! compute stateMask and nSubset and load into class object

  procedure :: advance_iSplit               => split_select_advance_iSplit               ! advance coupling iterator
  procedure :: advance_ixCoupling           => split_select_advance_ixCoupling           ! advance coupling iterator
  procedure :: advance_iStateTypeSplit      => split_select_advance_iStateTypeSplit      ! advance stateTypeSplitting iterator
  procedure :: advance_ixStateThenDomain    => split_select_advance_ixStateThenDomain    ! advance stateThenDomain iterator
  procedure :: advance_iDomainSplit         => split_select_advance_iDomainSplit         ! advance domainSplit iterator
  procedure :: advance_ixSolution           => split_select_advance_ixSolution           ! advance solution iterator
  procedure :: advance_iStateSplit          => split_select_advance_iStateSplit          ! advance stateSplit iterator
  
  procedure :: logic_exit_stateTypeSplitting => split_select_logic_exit_stateTypeSplitting ! get logical for branch
  procedure :: logic_exit_stateThenDomain    => split_select_logic_exit_stateThenDomain    ! get logical for branch
  procedure :: logic_exit_domainSplit        => split_select_logic_exit_domainSplit        ! get logical for branch
  procedure :: logic_exit_solution           => split_select_logic_exit_solution           ! get logical for branch
  procedure :: logic_exit_stateSplit         => split_select_logic_exit_stateSplit         ! get logical for branch

  procedure :: logic_initialize_stateTypeSplitting => split_select_logic_initialize_stateTypeSplitting ! get logical for branch
  procedure :: logic_initialize_stateThenDomain    => split_select_logic_initialize_stateThenDomain    ! get logical for branch
  procedure :: logic_initialize_domainSplit        => split_select_logic_initialize_domainSplit        ! get logical for branch
  procedure :: logic_initialize_solution           => split_select_logic_initialize_solution           ! get logical for branch
  procedure :: logic_initialize_stateSplit         => split_select_logic_initialize_stateSplit         ! get logical for branch

  procedure :: logic_finalize_stateTypeSplitting => split_select_logic_finalize_stateTypeSplitting     ! get logical for branch
  procedure :: logic_finalize_stateThenDomain    => split_select_logic_finalize_stateThenDomain        ! get logical for branch
  procedure :: logic_finalize_domainSplit        => split_select_logic_finalize_domainSplit            ! get logical for branch
  procedure :: logic_finalize_solution           => split_select_logic_finalize_solution               ! get logical for branch
  procedure :: logic_finalize_stateSplit         => split_select_logic_finalize_stateSplit             ! get logical for branch
end type split_select_type

contains

! ****** Class procedures for split_select_type class ******

subroutine split_select_initialize_flags(split_select)
 ! *** Initialize flags for opSplittin split methods ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % stateTypeSplitting=.false. 
 split_select % stateThenDomain=.false.
 split_select % domainSplit=.false.
 split_select % solution=.false.
 split_select % stateSplit=.false.
end subroutine split_select_initialize_flags

subroutine split_select_advance_iSplit(split_select)
 ! *** Advance index for coupling split method ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % iSplit = split_select % iSplit + 1
end subroutine split_select_advance_iSplit

subroutine split_select_advance_ixCoupling(split_select)
 ! *** Advance index for coupling split method ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % ixCoupling = split_select % ixCoupling + 1
end subroutine split_select_advance_ixCoupling

subroutine split_select_advance_iStateTypeSplit(split_select)
 ! *** Advance index for stateTypeSplit split method ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % iStateTypeSplit = split_select % iStateTypeSplit + 1
end subroutine split_select_advance_iStateTypeSplit

subroutine split_select_advance_ixStateThenDomain(split_select)
 ! *** Advance index for stateThenDomain split method ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % ixStateThenDomain = split_select % ixStateThenDomain + 1
end subroutine split_select_advance_ixStateThenDomain

subroutine split_select_advance_iDomainSplit(split_select)
 ! *** Advance index for domainSplit split method ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % iDomainSplit = split_select % iDomainSplit + 1
end subroutine split_select_advance_iDomainSplit

subroutine split_select_advance_ixSolution(split_select)
 ! *** Advance index for solution split method ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % ixSolution = split_select % ixSolution + 1
end subroutine split_select_advance_ixSolution

subroutine split_select_advance_iStateSplit(split_select)
 ! *** Advance index for stateSplit split method ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % iStateSplit = split_select % iStateSplit + 1
end subroutine split_select_advance_iStateSplit

subroutine split_select_initialize_ixCoupling(split_select)
 ! *** initialize operator splitting indices for split_select_type class ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % ixCoupling        = 1       
end subroutine split_select_initialize_ixCoupling

subroutine split_select_initialize_iStateTypeSplit(split_select)
 ! *** initialize operator splitting indices for split_select_type class ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % iStateTypeSplit        = 1       
end subroutine split_select_initialize_iStateTypeSplit

subroutine split_select_initialize_ixStateThenDomain(split_select)
 ! *** initialize operator splitting indices for split_select_type class ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % ixStateThenDomain        = 1       
end subroutine split_select_initialize_ixStateThenDomain

subroutine split_select_initialize_iDomainSplit(split_select)
 ! *** initialize operator splitting indices for split_select_type class ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % iDomainSplit        = 1       
end subroutine split_select_initialize_iDomainSplit

subroutine split_select_initialize_ixSolution(split_select)
 ! *** initialize operator splitting indices for split_select_type class ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % ixSolution        = 1       
end subroutine split_select_initialize_ixSolution

subroutine split_select_initialize_iStateSplit(split_select)
 ! *** initialize operator splitting indices for split_select_type class ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % iStateSplit        = 1       
end subroutine split_select_initialize_iStateSplit

logical(lgt) function split_select_logic_initialize_stateTypeSplitting(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_initialize_stateTypeSplitting=&
 &(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..false.).and.(split_select % domainSplit.eqv..false.).and.(split_select % stateThenDomain.eqv..false.).and.(split_select % stateTypeSplitting.eqv..false.)
end function split_select_logic_initialize_stateTypeSplitting

logical(lgt) function split_select_logic_exit_stateTypeSplitting(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_exit_stateTypeSplitting=&
 &(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..false.).and.(split_select % domainSplit.eqv..false.).and.(split_select % stateThenDomain.eqv..false.).and.(split_select % stateTypeSplitting.eqv..true.)
end function split_select_logic_exit_stateTypeSplitting

logical(lgt) function split_select_logic_initialize_stateThenDomain(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_initialize_stateThenDomain=&
 &(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..false.).and.(split_select % domainSplit.eqv..false.).and.(split_select % stateThenDomain.eqv..false.)
end function split_select_logic_initialize_stateThenDomain

logical(lgt) function split_select_logic_exit_stateThenDomain(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_exit_stateThenDomain=&
 &(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..false.).and.(split_select % domainSplit.eqv..false.).and.(split_select % stateThenDomain.eqv..true.)
end function split_select_logic_exit_stateThenDomain

logical(lgt) function split_select_logic_initialize_domainSplit(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_initialize_domainSplit=&
 &(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..false.).and.(split_select % domainSplit.eqv..false.)
end function split_select_logic_initialize_domainSplit

logical(lgt) function split_select_logic_exit_domainSplit(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_exit_domainSplit=&
 &(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..false.).and.(split_select % domainSplit.eqv..true.)
end function split_select_logic_exit_domainSplit

logical(lgt) function split_select_logic_initialize_solution(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_initialize_solution=(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..false.)
end function split_select_logic_initialize_solution

logical(lgt) function split_select_logic_exit_solution(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_exit_solution=(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..true.)
end function split_select_logic_exit_solution

logical(lgt) function split_select_logic_initialize_stateSplit(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_initialize_stateSplit=(split_select % stateSplit.eqv..false.)
end function split_select_logic_initialize_stateSplit

logical(lgt) function split_select_logic_exit_stateSplit(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_exit_stateSplit=(split_select % stateSplit.eqv..true.)
end function split_select_logic_exit_stateSplit

logical(lgt) function split_select_logic_finalize_stateSplit(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_finalize_stateSplit=&
 &(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..true.).and.(split_select % stateThenDomain.eqv..true.)
end function split_select_logic_finalize_stateSplit

logical(lgt) function split_select_logic_finalize_solution(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_finalize_solution=&
 &(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..false.).and.(split_select % stateThenDomain.eqv..true.)
end function split_select_logic_finalize_solution

logical(lgt) function split_select_logic_finalize_domainSplit(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_finalize_domainSplit=&
 &(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..false.).and.(split_select % domainSplit.eqv..false.).and.(split_select % stateThenDomain.eqv..true.)
end function split_select_logic_finalize_domainSplit

logical(lgt) function split_select_logic_finalize_stateThenDomain(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_finalize_stateThenDomain=&
 &(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..false.).and.(split_select % domainSplit.eqv..false.).and.(split_select % stateThenDomain.eqv..false.)
end function split_select_logic_finalize_stateThenDomain

logical(lgt) function split_select_logic_finalize_stateTypeSplitting(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_finalize_stateTypeSplitting=&
 &(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..false.).and.(split_select % domainSplit.eqv..false.).and.(split_select % stateThenDomain.eqv..false.).and.(split_select % stateTypeSplitting.eqv..false.)
end function split_select_logic_finalize_stateTypeSplitting

subroutine split_select_compute_stateMask(split_select,indx_data,err,cmessage,message,return_flag)
 ! *** Get the mask for the state subset ***
 class(split_select_type),intent(inout) :: split_select              ! class object for operator splitting selector
 type(var_ilength),intent(in)           :: indx_data                 ! indices for a local HRU
 integer(i4b),intent(out)               :: err                       ! intent(out): error code
 character(*),intent(out)               :: cmessage                  ! intent(out): error message
 character(*),intent(out)               :: message                   ! error message
 logical(lgt),intent(out)               :: return_flag               ! return flag
 ! local variables
 type(out_type_stateFilter)             :: out_stateFilter           ! number of selected state variables for a given split and error control

 err=0               ! initialize error code 
 return_flag=.false. ! initialize flag
 call stateFilter(indx_data,split_select,out_stateFilter)
 call out_stateFilter % finalize(err,cmessage)
 if (err/=0) then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if  ! error control
end subroutine split_select_compute_stateMask


! **********************************************************************************************************
! private subroutine stateFilter: get a mask for the desired state variables
! **********************************************************************************************************
 subroutine stateFilter(indx_data,split_select,out_stateFilter)
 USE indexState_module,only:indxSubset                               ! get state indices
 implicit none
 ! input
 type(var_ilength),intent(in)           :: indx_data                 ! indices for a local HRU
 ! input-output
 type(split_select_type),intent(inout)  :: split_select              ! class object for operator splitting selector
 ! output
 type(out_type_stateFilter),intent(out) :: out_stateFilter           ! number of selected state variables for a given split and error control
 ! local
 integer(i4b),allocatable               :: ixSubset(:)               ! list of indices in the state subset
 character(len=256)                     :: cmessage                  ! error message
 logical(lgt)                           :: return_flag               ! flag to indicate a return 
 ! ----------------------------------------------------------------------------------------------------------------------------------------------------
 ! data structures
 associate(ixCoupling => split_select % ixCoupling  ,&  ! intent(in): [i4b] index of coupling method (1,2)
           err        => out_stateFilter % err      ,&  ! intent(out): error code
           message    => out_stateFilter % cmessage  )  ! intent(out): error message
   
  err=0; message='stateFilter/'; return_flag=.false. ! initialize error control

  ! identify splitting option
  select case(ixCoupling)
   ! *** fully coupled ***
   case(fullyCoupled); call fullyCoupled_stateMask ! get stateMask for fully coupled method 
   ! *** splitting by state type ***
   case(stateTypeSplit) ! initial split by state type
    call stateTypeSplit_stateMask; if (return_flag) return ! get stateMask for state split method -- return if error
    ! check
   case default; err=20; message=trim(message)//'unable to identify coupling method'; return_flag=.true.; return
  end select  ! selecting solution method 

  ! initialize ixSubset
  allocate(ixSubset(1_i4b),STAT=err)
  if (err/=0) then; message=trim(message)//'allocation error in stateFilter for ixSubset'; return_flag=.true.; return; end if
  ixSubset = 0._rkind
 end associate

 call identify_scalar_solutions; if (return_flag) return ! identify scalar solutions -- return if error occurs

 ! get the number of selected state variables
 split_select % nSubset = count(split_select % stateMask)

contains

 subroutine fullyCoupled_stateMask
  ! *** Get fully coupled stateMask ***
  split_select % stateMask(:) = .true. ! use all state variables
 end subroutine fullyCoupled_stateMask

 subroutine stateTypeSplit_stateMask
  ! *** Get state type split stateMask ***
  return_flag=.false. ! initialize flag
  ! switch between full domain and sub domains
  associate(&
   ixStateThenDomain => split_select % ixStateThenDomain,& ! intent(in): [i4b] switch between full domain and sub domains
   err               => out_stateFilter % err           ,& ! intent(out): error code
   message           => out_stateFilter % cmessage       ) ! intent(out): error message
   select case(ixStateThenDomain)
     ! split into energy and mass
     case(fullDomain); call stateTypeSplit_fullDomain_stateMask; if (return_flag) return
     ! split into vegetation, snow, and soil
     case(subDomain); call stateTypeSplit_subDomain_stateMask; if (return_flag) return
     ! check
     case default
       err=20; message=trim(message)//'unable to identify the switch between full domains and sub domains'; return_flag=.true.; return
   end select 
  end associate
 end subroutine stateTypeSplit_stateMask

 subroutine stateTypeSplit_fullDomain_stateMask
  ! *** Get full domain stateMask ***
  return_flag=.false. ! initialize flag
  associate(iStateTypeSplit => split_select % iStateTypeSplit,& ! intent(in): [i4b] index of the state type split
            err             => out_stateFilter % err         ,& ! intent(out): error code
            message         => out_stateFilter % cmessage     ) ! intent(out): error message
   select case(iStateTypeSplit)
    case(nrgSplit);  call stateTypeSplit_fullDomain_nrgSplit_stateMask
    case(massSplit); call stateTypeSplit_fullDomain_massSplit_stateMask
    case default; err=20; message=trim(message)//'unable to identify split based on state type'; return_flag=.true.; return
   end select
  end associate
 end subroutine stateTypeSplit_fullDomain_stateMask

 subroutine stateTypeSplit_fullDomain_nrgSplit_stateMask
  ! *** Get state type full domain energy split stateMask ***
  associate(ixStateType     => indx_data%var(iLookINDEX%ixStateType)%dat) ! intent(in): [i4b(:)] indices defining the type of the state (ixNrgState...)
   split_select % stateMask(:) = (ixStateType==iname_nrgCanair .or. ixStateType==iname_nrgCanopy .or. ixStateType==iname_nrgLayer)
  end associate
 end subroutine stateTypeSplit_fullDomain_nrgSplit_stateMask

 subroutine stateTypeSplit_fullDomain_massSplit_stateMask
  ! *** Get state type full domain mass split stateMask ***
  associate(ixStateType     => indx_data%var(iLookINDEX%ixStateType)%dat) ! intent(in): [i4b(:)] indices defining the type of the state (ixNrgState...)
   split_select % stateMask(:) = (ixStateType==iname_liqCanopy .or. ixStateType==iname_liqLayer  .or. &
                            &  ixStateType==iname_lmpLayer  .or. ixStateType==iname_watAquifer)
  end associate
 end subroutine stateTypeSplit_fullDomain_massSplit_stateMask

 subroutine stateTypeSplit_subDomain_stateMask
  ! *** Get subdomain stateMask ***
  return_flag=.false. ! initialize flag
  ! define state mask
  associate(&
            iStateTypeSplit => split_select % iStateTypeSplit,& ! intent(in): [i4b] index of the state type split
            err             => out_stateFilter % err         ,& ! intent(out): error code
            message         => out_stateFilter % cmessage     ) ! intent(out): error message
   split_select % stateMask(:)=.false. ! initialize state mask
   select case(iStateTypeSplit)
    ! define mask for energy
    case(nrgSplit); call stateTypeSplit_subDomain_nrgSplit_stateMask; if (return_flag) return
    ! define mask for water
    case(massSplit); call stateTypeSplit_subDomain_massSplit_stateMask; if (return_flag) return
    ! check
    case default; err=20; message=trim(message)//'unable to identify the state type'; return_flag=.true.; return
   end select  ! (split based on state type)
  end associate
 end subroutine stateTypeSplit_subDomain_stateMask

 subroutine stateTypeSplit_subDomain_nrgSplit_stateMask
  ! *** Get subdomain energy split stateMask ***
  return_flag=.false. ! initialize flag
  associate(&
   iDomainSplit    => split_select % iDomainSplit,& ! intent(in): [i4b] index of the domain split
   err             => out_stateFilter % err      ,& ! intent(out): error code
   message         => out_stateFilter % cmessage  ) ! intent(out): error message
   select case(iDomainSplit)
    case(vegSplit);  call stateTypeSplit_subDomain_nrgSplit_vegSplit_stateMask       ! vegetation subdomain
    case(snowSplit); call stateTypeSplit_subDomain_nrgSplit_snowSplit_stateMask      ! snow subdomain
    case(soilSplit); call stateTypeSplit_subDomain_nrgSplit_soilSplit_stateMask      ! soil subdomain
    case(aquiferSplit) ! do nothing: no energy state variable for the aquifer domain ! aquifer subdomain 
    case default; err=20; message=trim(message)//'unable to identify model sub-domain'; return_flag=.true.; return
   end select
  end associate
 end subroutine stateTypeSplit_subDomain_nrgSplit_stateMask

 subroutine stateTypeSplit_subDomain_nrgSplit_vegSplit_stateMask
  ! *** Get state type subdomain energy vegetation split ***
  associate(&
   ixNrgCanair => indx_data%var(iLookINDEX%ixNrgCanair)%dat,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for energy states in canopy air space domain
   ixNrgCanopy => indx_data%var(iLookINDEX%ixNrgCanopy)%dat,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for energy states in the canopy domain
   ixNrgLayer  => indx_data%var(iLookINDEX%ixNrgLayer)%dat  ) ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for energy states in the snow+soil domain
     if (ixNrgCanair(1)/=integerMissing) split_select % stateMask(ixNrgCanair) = .true.  ! energy of the canopy air space
     if (ixNrgCanopy(1)/=integerMissing) split_select % stateMask(ixNrgCanopy) = .true.  ! energy of the vegetation canopy
     split_select % stateMask(ixNrgLayer(1)) = .true.  ! energy of the upper-most layer in the snow+soil domain
  end associate
 end subroutine stateTypeSplit_subDomain_nrgSplit_vegSplit_stateMask

 subroutine stateTypeSplit_subDomain_nrgSplit_snowSplit_stateMask
  associate(&
   nSnow      => indx_data%var(iLookINDEX%nSnow)%dat(1)  ,& ! intent(in): [i4b] number of snow layers
   ixNrgLayer => indx_data%var(iLookINDEX%ixNrgLayer)%dat ) ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for energy states in the snow+soil domain
  ! *** Get state type subdomain energy snow split ***
   if (nSnow>1) split_select % stateMask(ixNrgLayer(2:nSnow)) = .true. ! NOTE: (2:) because the top layer in the snow+soil domain included in vegSplit
  end associate
 end subroutine stateTypeSplit_subDomain_nrgSplit_snowSplit_stateMask

 subroutine stateTypeSplit_subDomain_nrgSplit_soilSplit_stateMask
  ! *** Get state type subdomain energy soil split ***
  associate(&
   nSnow      => indx_data%var(iLookINDEX%nSnow)%dat(1)  ,& ! intent(in): [i4b] number of snow layers
   nLayers    => indx_data%var(iLookINDEX%nLayers)%dat(1),& ! intent(in): [i4b] total number of layers
   ixNrgLayer => indx_data%var(iLookINDEX%ixNrgLayer)%dat ) ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for energy states in the snow+soil domain
   split_select % stateMask(ixNrgLayer(max(2,nSnow+1):nLayers)) = .true. ! NOTE: max(2,nSnow+1) gives second layer unless more than 2 snow layers
  end associate
 end subroutine stateTypeSplit_subDomain_nrgSplit_soilSplit_stateMask

 subroutine stateTypeSplit_subDomain_massSplit_stateMask
  ! *** Get subdomain mass split stateMask ***
  return_flag=.false. ! initialize flag
  associate(&
   iDomainSplit => split_select % iDomainSplit,& ! intent(in): [i4b] index of the domain split
   err          => out_stateFilter % err      ,& ! intent(out): error code
   message      => out_stateFilter % cmessage  ) ! intent(out): error message
   select case(iDomainSplit)
    case(vegSplit);     call stateTypeSplit_subDomain_massSplit_vegSplit_stateMask     ! vegetation subdomain
    case(snowSplit);    call stateTypeSplit_subDomain_massSplit_snowSplit_stateMask    ! snow subdomain
    case(soilSplit);    call stateTypeSplit_subDomain_massSplit_soilSplit_stateMask    ! soil subdomain
    case(aquiferSplit); call stateTypeSplit_subDomain_massSplit_aquiferSplit_stateMask ! aquifer subdomain 
    case default; err=20; message=trim(message)//'unable to identify model sub-domain'; return_flag=.true.; return
   end select
  end associate
 end subroutine stateTypeSplit_subDomain_massSplit_stateMask

 subroutine stateTypeSplit_subDomain_massSplit_vegSplit_stateMask
  ! *** Get mass state vegetation subdomain split stateMask  ***
  associate(ixHydCanopy => indx_data%var(iLookINDEX%ixHydCanopy)%dat) ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the canopy domain
   if (ixHydCanopy(1)/=integerMissing) split_select % stateMask(ixHydCanopy) = .true. ! hydrology of the vegetation canopy
  end associate
 end subroutine stateTypeSplit_subDomain_massSplit_vegSplit_stateMask

 subroutine stateTypeSplit_subDomain_massSplit_snowSplit_stateMask
  ! *** Get mass state snow subdomain split stateMask  ***
  associate(&
   nSnow      => indx_data%var(iLookINDEX%nSnow)%dat(1)  ,& ! intent(in): [i4b] number of snow layers
   ixHydLayer => indx_data%var(iLookINDEX%ixHydLayer)%dat ) ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the snow+soil domain
   split_select % stateMask(ixHydLayer(1:nSnow)) = .true.  ! snow hydrology
  end associate
 end subroutine stateTypeSplit_subDomain_massSplit_snowSplit_stateMask

 subroutine stateTypeSplit_subDomain_massSplit_soilSplit_stateMask
  ! *** Get mass state soil subdomain split stateMask  ***
  associate(&
   nSnow      => indx_data%var(iLookINDEX%nSnow)%dat(1)  ,& ! intent(in): [i4b] number of snow layers
   nLayers    => indx_data%var(iLookINDEX%nLayers)%dat(1),& ! intent(in): [i4b] total number of layers
   ixHydLayer => indx_data%var(iLookINDEX%ixHydLayer)%dat ) ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the snow+soil domain
   split_select % stateMask(ixHydLayer(nSnow+1:nLayers)) = .true.  ! soil hydrology
  end associate
 end subroutine stateTypeSplit_subDomain_massSplit_soilSplit_stateMask

 subroutine stateTypeSplit_subDomain_massSplit_aquiferSplit_stateMask
  ! *** Get mass state aquifer subdomain split stateMask  ***
  associate(ixWatAquifer => indx_data%var(iLookINDEX%ixWatAquifer)%dat) ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for water storage in the aquifer
   if (ixWatAquifer(1)/=integerMissing) split_select % stateMask(ixWatAquifer) = .true. ! aquifer storage
  end associate
 end subroutine stateTypeSplit_subDomain_massSplit_aquiferSplit_stateMask

 subroutine identify_scalar_solutions
  ! *** Identify scalar solutions ***
  return_flag=.false. ! initialize flag
  associate(ixAllState  => indx_data%var(iLookINDEX%ixAllState)%dat,& ! intent(in): [i4b(:)] list of indices for all model state variables (1,2,3,...nState)
            ixSolution  => split_select % ixSolution               ,& ! intent(in): [i4b] index of solution method (1,2)
            iStateSplit => split_select % iStateSplit              ,& ! intent(in): [i4b] index of the layer split
            err         => out_stateFilter % err                   ,& ! intent(out): error code
            message     => out_stateFilter % cmessage               ) ! intent(out): error message
   if (ixSolution==scalar) then
    ! get the subset of indices
    call indxSubset(ixSubset, ixAllState, split_select % stateMask, err, cmessage)
    if (err/=0) then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if
    ! get the mask
    split_select % stateMask(:) = .false.
    split_select % stateMask( ixSubset(iStateSplit) ) = .true.
    ! check
    if (count(split_select % stateMask)/=1) then
     message=trim(message)//'expect size=1 (scalar)'
     err=20; return_flag=.true.; return
    end if
   end if
  end associate

 end subroutine identify_scalar_solutions
 
end subroutine stateFilter

end module stateFilter_module
