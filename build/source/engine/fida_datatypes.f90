

module fida_datatypes

! data types
USE nrtype
USE, intrinsic :: iso_c_binding

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (dp)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (dp)
                    zLookup         ! data vector with variable length dimension (dp)

implicit none

 type eqnsData
   type(c_ptr)                       	:: ida_mem                			! IDA memory
   real(dp)             				:: dt                     			! time step
   integer(i4b)         				:: nSnow                  			! number of snow layers
   integer(i4b)         				:: nSoil                  			! number of soil layers
   integer(i4b)         				:: nLayers                			! total number of layers
   integer              				:: nState                 			! total number of state variables
   integer(i4b)         				:: ixMatrix               			! form of matrix (dense or banded)
   logical(lgt)         				:: firstSubStep           			! flag to indicate if we are processing the first sub-step
   logical(lgt)         				:: computeVegFlux         			! flag to indicate if computing fluxes over vegetation
   logical(lgt)         				:: scalarSolution         			! flag to denote if implementing the scalar solution
   type(zLookup)         				:: lookup_data            			! lookup tables
   type(var_i)           				:: type_data              			! type of vegetation and soil
   type(var_d)           				:: attr_data              			! spatial attributes
   type(var_dlength)     				:: mpar_data              			! model parameters
   type(var_d)           				:: forc_data              			! model forcing data
   type(var_dlength)     				:: bvar_data              			! model variables for the local basin
   type(var_dlength)     				:: prog_data              			! prognostic variables for a local HRU
   type(var_ilength)     				:: indx_data              			! indices defining model states and layers
   type(var_dlength)     				:: diag_data              			! diagnostic variables for a local HRU
   type(var_dlength)     				:: flux_data              			! model fluxes for a local HRU
   type(var_dlength)     				:: deriv_data             			! derivatives in model fluxes w.r.t. relevant state variables
   real(dp)								:: scalarCanopyTempTrial  			! trial value of canopy temperature (K)
   real(dp)								:: scalarCanopyTempPrev   			! previous value of canopy temperature (K)
   real(dp)								:: scalarCanopyIceTrial				! trial value of canopy ice content (kg m-2)
   real(dp)								:: scalarCanopyIcePrev				! value of canopy ice content (kg m-2) at previous step
   real(dp)								:: scalarCanopyEnthalpyTrial  		! trial enthalpy of the vegetation canopy (J m-3)
   real(dp)								:: scalarCanopyEnthalpyPrev 		! previous enthalpy of the vegetation canopy (J m-3)
   real(qp), allocatable             	:: sMul(:)      		   			! state vector multiplier (used in the residual calculations)
   real(dp), allocatable             	:: dMat(:) 							! diagonal of the Jacobian matrix
   real(dp), allocatable              	:: fluxVec(:)             			! flux vector
   real(qp), allocatable              	:: resSink(:) 						! additional (sink) terms on the RHS of the state equation
   real(dp), allocatable   				:: atol(:)							! vector of absolute tolerances
   real(dp), allocatable   				:: rtol(:) 							! vector of relative tolerances
   real(dp), allocatable              	:: mLayerMatricHeadLiqTrial(:)
   real(dp), allocatable              	:: mLayerMatricHeadTrial(:)
   real(dp), allocatable              	:: mLayerMatricHeadPrev(:)
   real(dp), allocatable              	:: mLayerVolFracWatTrial(:)     	! trial value for volumetric fraction of total water (-)
   real(dp), allocatable              	:: mLayerVolFracWatPrev(:)
   real(dp), allocatable              	:: mLayerVolFracIceTrial(:)     	! trial value for volumetric fraction of total water (-)
   real(dp), allocatable              	:: mLayerVolFracIcePrev(:)
   real(dp), allocatable              	:: mLayerEnthalpyTrial(:)
   real(dp), allocatable              	:: mLayerEnthalpyPrev(:)
   real(dp), allocatable              	:: mLayerTempTrial(:)
   real(dp), allocatable              	:: mLayerTempPrev(:)
   real(dp), allocatable              	:: dBaseflow_dMatric(:,:) 			! derivative in baseflow w.r.t. matric head (s-1)
   integer(i4b)          				:: err                    			! error code
   character(len = 50)          		:: message                			! error message
 end type eqnsData
 
 
end module fida_datatypes





