module md_plant
  !////////////////////////////////////////////////////////////////
  !  Module contains (constrainable) model parameters.
  !  Model parameters adopted here are from LPX C3 grass PFT
  !  Litter and soil turnover parameters are divided by 365 to 
  !  convert from [1/yr] to [1/d].
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_classdefs
  use md_params_core, only: ndayyear, npft, nlu, lunat
  use md_interface_pmodel, only: myinterface

  implicit none

  private
  public plant_type, plant_fluxes_type, getpar_modl_plant, &
    initglobal_plant, params_plant, params_pft_plant, ftemp, &
    fmoist, add_seed, get_leaftraits, get_lai, get_fapar, init_annual_plant

  !----------------------------------------------------------------
  ! NON PFT-DEPENDENT PARAMETERS
  !----------------------------------------------------------------
  type params_plant_type
    real :: kbeer             ! canopy light extinction coefficient
    real :: r_root            ! Fine root-specific respiration rate (gC gC-1 d-1)
    real :: r_sapw            ! Sapwood-specific respiration rate (gC gC-1 d-1)
    real :: exurate           ! Fine root-specific C export rate (gC gC-1 d-1)
    real :: f_nretain         ! fraction of N retained at leaf abscission 
    real :: fpc_tree_max      ! maximum fractional plant coverage of trees
    real :: growtheff         ! growth respiration coefficient = yield factor [unitless]
    ! real :: cton_soil         ! C:N ratio of soil organic matter (consider this to be equal to that of microbial biomass)
    real :: frac_leaf         ! fraction of allocatable C to leaf 
  end type params_plant_type

  type( params_plant_type ) :: params_plant

  !----------------------------------------------------------------
  ! PFT-DEPENDENT PARAMETERS
  !----------------------------------------------------------------
  type params_pft_plant_type

    character(len=4) :: pftname    ! standard PFT name with 4 characters length
    integer :: lu_category         ! land use category associated with PFT
    logical, dimension(nlu) :: islu! islu(ipft,ilu) is true if ipft belongs to ilu
    logical :: grass               ! boolean for growth form 'grass'
    logical :: tree                ! boolean for growth form 'tree'
    logical :: nfixer              ! whether plant is capable of symbiotically fixing N
    logical :: c3                  ! whether plant follows C3 photosynthesis
    logical :: c4                  ! whether plant follows C4 photosynthesis
    real    :: sla                 ! specific leaf area (m2 gC-1)
    real    :: lma                 ! leaf mass per area (gC m-2)
    real    :: r_ntolma            ! constant ratio of structural N to C (LMA) (gN/gC)

    ! new for cnmodel
    real    :: k_decay_leaf_base   ! base leaf decay constant [year-1]
    real    :: k_decay_leaf_width  ! shape parameter for turnover function if LAI
    real    :: k_decay_sapw        ! sapwood decay constant [year-1]
    real    :: k_decay_root        ! root decay constant [year-1]
    real    :: k_decay_labl        ! labile pool decay constant [year-1]
    real    :: r_cton_root         ! C:N ratio in roots (gC/gN)
    real    :: r_ntoc_root         ! N:C ratio in roots (inverse of 'r_cton_root', gN/gC)
    real    :: ncw_min             ! y-axis intersection in the relationship of non-metabolic versus metabolic N per leaf area    
    real    :: r_n_cw_v            ! slope in the relationship of non-metabolic versus metabolic N per leaf area              
    real    :: r_ctostructn_leaf   ! constant ratio of C to structural N (mol C / mol N)

  end type params_pft_plant_type

  type(params_pft_plant_type), dimension(npft) :: params_pft_plant

  !----------------------------------------------------------------
  ! Pools and other variables with year-to-year memory
  !----------------------------------------------------------------
  type plant_type

    ! PFT index that goes along with this instance of 'plant'
    integer :: pftno

    ! canopy at individual-level
    integer :: nind             ! number of individuals (m-2)
    real :: fpc_grid            ! fractional projective cover
    real :: lai_ind             ! fraction of absorbed photosynthetically active radiation
    real :: fapar_ind           ! fraction of absorbed photosynthetically active radiation
    real :: acrown              ! crown area

    ! leaf traits
    real :: narea               ! total leaf N per unit leaf area (gN m-2)
    real :: narea_metabolic     ! metabolic leaf N per unit leaf area (gN m-2)
    real :: narea_structural    ! structural leaf N per unit leaf area (gN m-2)
    real :: lma                 ! leaf mass per area (gC m-2)
    real :: sla                 ! specific leaf area (m2 gC-1)
    real :: nmass               ! leaf N per unit leaf mass, g N / g-dry mass
    real :: r_cton_leaf         ! leaf C:N ratio [gC/gN] 
    real :: r_ntoc_leaf         ! leaf N:C ratio [gN/gC]

    ! pools
    type(orgpool) :: pleaf     ! leaf biomass [gC/ind.] (=lm_ind)
    type(orgpool) :: proot     ! root biomass [gC/ind.] (=rm_ind)
    type(orgpool) :: psapw     ! sapwood biomass [gC/ind.] (=sm_ind)
    type(orgpool) :: pwood     ! heartwood (non-living) biomass [gC/ind.] (=hm_ind)
    type(orgpool) :: plabl     ! labile pool, temporary storage of N and C [gC/ind.] (=bm_inc but contains also N) 

  end type plant_type


  !----------------------------------------------------------------
  ! Fluxes and other variables with no memory
  !----------------------------------------------------------------
  type plant_fluxes_type

    ! daily updated variables
    real :: dgpp              ! daily gross primary production [gC/m2/d]           
    real :: drd               ! daily dark respiration [gC/m2/d]
    real :: assim             ! daily assimilation (mol CO2 m-2 s-1)
    real :: dtransp           ! daily transpiration [mm]
    real :: dlatenth          ! daily latent heat flux [J m-2 d-1]

    real :: drleaf   ! daily total leaf respiration, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
    real :: drroot   ! root maintenance respiration, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
    real :: drsapw   ! sapwood maintenance respiration, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
    real :: drgrow   ! growth respiration (growth+maintenance resp. of all compartments), no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
    real :: dcex     ! labile C exudation for N uptake, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
    
    type(carbon)   :: dnpp     ! daily net primary production (gpp-ra, npp=bp+cex) [gC/m2/d]
    type(nitrogen) :: dnup     ! daily N uptake [gN/m2/d]

    real :: dnup_pas          ! daily N uptake by passsive uptake (transpiration) [gN/m2/d]
    real :: dnup_act          ! daily N uptake by active uptake [gN/m2/d]
    real :: dnup_fix          ! daily N uptake by plant symbiotic N fixation [gN/m2/d]
    real :: dnup_ret          ! daily N "uptake" by plant symbiotic N fixation [gN/m2/d]

    real :: vcmax25           ! acclimated Vcmax, normalised to 25 deg C (mol CO2 m-2 s-1)
    real :: jmax25            ! acclimated Jmax, normalised to 25 deg C (mol CO2 m-2 s-1)
    real :: vcmax             ! daily varying Vcmax (mol CO2 m-2 s-1)
    real :: jmax              ! daily varying Jmax (mol CO2 m-2 s-1)
    real :: gs_accl           ! acclimated stomatal conductance (xxx)
    real :: chi               ! ci:ca ratio (unitless)
    real :: iwue              ! intrinsic water use efficiency (A/gs = ca*(1-chi))
    real :: actnv_unitiabs    ! metabolic leaf N per unit absorbed light (g N m-2 mol-1)

    type(orgpool) :: dharv    ! daily total biomass harvest (g m-2 d-1)

  end type plant_fluxes_type

  !-----------------------------------------------------------------------
  ! Fixed parameters
  !-----------------------------------------------------------------------
  ! type( orgpool ), parameter :: seed = orgpool( carbon(5.0), nitrogen(0.0) )
  type( orgpool ), parameter :: seed = orgpool( carbon(5.0), nitrogen(0.12) )
  ! type( orgpool ), parameter :: seed = orgpool( carbon(100.0), nitrogen(1 .0) )


contains

  function get_fapar( lai ) result( fapar )
    !////////////////////////////////////////////////////////////////
    ! FOLIAGE PROJECTIVE COVER 
    ! = Fraction of Absorbed Photosynthetically Active Radiation
    ! Function returns fractional plant cover an individual
    ! Eq. 7 in Sitch et al., 2003
    !----------------------------------------------------------------
    ! arguments
    real, intent(in) :: lai

    ! function return variable
    real :: fapar

    fapar = ( 1.0 - exp( -1.0 * params_plant%kbeer * lai) )

  end function get_fapar


  function get_lai( pft, cleaf, meanmppfd, nv ) result( lai )
    !////////////////////////////////////////////////////////////////
    ! Calculates LAI as a function of leaf-C. This is not so straight
    ! forward due to the dependence of canopy-metabolic leaf N on LAI,
    ! and the dependence of canopy-structural leaf N and C on canopy-
    ! metabolic leaf N.
    !----------------------------------------------------------------
    use md_params_core, only: nmonth, c_molmass
    use md_lambertw, only: calc_wapr

    ! arguments
    integer, intent(in) :: pft
    real, intent(in) :: cleaf
    real, intent(in) :: meanmppfd
    real, intent(in) :: nv 

    ! function return variable
    real :: lai

    ! local variables
    real    :: alpha, beta, gamma ! variable substitutes
    real    :: maxnv
    real    :: arg_to_lambertw
    integer :: nerror


    if (cleaf > 0.0) then

      ! Monthly variations in metabolic N, determined by variations in meanmppfd and nv should not result in variations in leaf traits. 
      ! In order to prevent this, assume annual maximum metabolic N, part of which is deactivated during months with lower insolation (and Rd reduced.)
      maxnv = meanmppfd * nv

      alpha = maxnv * params_pft_plant(pft)%r_n_cw_v
      beta  = params_pft_plant(pft)%ncw_min
      gamma = cleaf / ( c_molmass * params_pft_plant(pft)%r_ctostructn_leaf ) 
      arg_to_lambertw = alpha * params_plant%kbeer / beta * exp( (alpha - gamma) * params_plant%kbeer / beta )
      lai = 1.0 / (beta * params_plant%kbeer ) * &
        ( -alpha * params_plant%kbeer + &
          gamma * params_plant%kbeer + &
          beta * calc_wapr( arg_to_lambertw, 0, nerror, 9999 ) &
        )
    else

      lai = 0.0

    end if
    
  end function get_lai


  function get_leaf_n_metabolic_canopy( mylai, meanmppfd, nv, myfapar ) result( mynleaf_metabolic )
    !////////////////////////////////////////////////////////////////
    ! Calculates metabolic leaf N at canopy-level, determined by 
    ! light conditions (meanmppfd) and the Rubisco-N per unit absorbed
    ! light.
    !----------------------------------------------------------------
    use md_params_core, only: nmonth

    ! arguments
    real, intent(in)                    :: mylai
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv
    real, intent(in), optional          :: myfapar

    ! function return variable
    real :: mynleaf_metabolic  ! mol N m-2-ground

    ! local variables
    real :: maxnv

    ! Metabolic N is predicted and is optimised at a monthly time scale. 
    ! Leaf traits are calculated based on metabolic N => cellwall N => cellwall C / LMA
    ! Leaves get thinner at the bottom of the canopy => increasing LAI through the season comes at a declining C and N cost
    ! Monthly variations in metabolic N, determined by variations in meanmppfd and nv should not result in variations in leaf traits. 
    ! In order to prevent this, assume annual maximum metabolic N, part of which is deactivated during months with lower insolation (and Rd reduced.)
    maxnv = maxval( meanmppfd(:) * nv(:) )

    if (present(myfapar)) then
      mynleaf_metabolic = maxnv * myfapar
    else
      mynleaf_metabolic = maxnv * get_fapar( mylai )
    end if

  end function get_leaf_n_metabolic_canopy


  subroutine get_leaftraits( plant, meanmppfd, nv )
    !////////////////////////////////////////////////////////////////
    ! Calculates leaf traits based on (predicted) metabolic Narea and
    ! (prescribed) parameters that relate structural to metabolic
    ! Narea and Carea to structural Narea:
    ! Narea_metabolic  = predicted
    ! Narea_structural = rN:C_struct * LMA
    !----------------------------------------------------------------
    use md_params_core, only: c_content_of_biomass, nmonth, n_molmass, c_molmass

    ! arguments
    type( plant_type ), intent(inout) :: plant
    real, intent(in) :: meanmppfd
    real, intent(in) :: nv

    ! local variables
    real :: narea_metabolic_canopy   ! g N m-2-ground

    ! canopy-level, in units of gN / m2-ground 
    narea_metabolic_canopy  = n_molmass * plant%fapar_ind * meanmppfd * nv

    ! leaf-level, in units of gN / m2-leaf 
    ! assume narea_metabolic is representative of the outer canopy, therefore divide by 1.0 (or just leave)
    plant%narea_metabolic  = narea_metabolic_canopy / 1.0
    plant%narea_structural = params_pft_plant(plant%pftno)%r_ntolma * params_pft_plant(plant%pftno)%lma
    plant%narea            = plant%narea_metabolic + plant%narea_structural
    plant%lma              = params_pft_plant(plant%pftno)%lma

    ! additional traits
    plant%nmass            = plant%narea / ( plant%lma / c_content_of_biomass )
    plant%r_cton_leaf      = params_pft_plant(plant%pftno)%lma / plant%narea
    plant%r_ntoc_leaf      = 1.0 / plant%r_cton_leaf

  end subroutine get_leaftraits


  subroutine getpar_modl_plant()
    !////////////////////////////////////////////////////////////////
    !  Subroutine reads model parameters from input file.
    !  It was necessary to separate this SR from module md_plant
    !  because this SR uses module md_waterbal, which also uses
    !  _plant.
    ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
    ! contact: b.stocker@imperial.ac.uk
    !----------------------------------------------------------------    
    ! local variables
    integer :: pft
    integer :: npft_site

    !----------------------------------------------------------------
    ! NON-PFT DEPENDENT PARAMETERS
    !----------------------------------------------------------------
    ! canopy light extinction coefficient for Beer's Law
    params_plant%kbeer = myinterface%params_calib%kbeer

    ! fraction of N retained at leaf abscission 
    params_plant%f_nretain = myinterface%params_calib%f_nretain
    
    ! maximum fractional plant coverage of trees (sum of all tree PFTs)
    params_plant%fpc_tree_max = myinterface%params_calib%fpc_tree_max

    ! growth efficiency = yield factor, central value: 0.6, range: 0.5-0.7; Zhang et al. (2009), see Li et al., 2014
    params_plant%growtheff = myinterface%params_calib%growtheff

    ! Fine-root mass specific respiration rate (gC gC-1 year-1)
    ! Central value: 0.913 year-1 (Yan and Zhao (2007); see Li et al., 2014)
    params_plant%r_root = myinterface%params_calib%r_root

    ! Sapwood specific respiration rate (gC gC-1 year-1)
    ! Central value: 0.044 year-1 (Yan and Zhao (2007); see Li et al., 2014)
    ! (= 0.044 nmol mol-1 s-1; range: 0.5–10, 20 nmol mol-1 s-1 (Landsberg and Sands (2010))
    params_plant%r_sapw = myinterface%params_calib%r_sapw

    ! C export rate per unit root mass
    params_plant%exurate = myinterface%params_calib%exurate

    ! ! C:N ratio of soil organic matter [1]
    ! params_plant%cton_soil = myinterface%params_calib%cton_soil


    !----------------------------------------------------------------
    ! PFT DEPENDENT PARAMETERS
    ! read parameter input file and store values in single array
    ! important: Keep this order of reading PFT parameters fixed.
    !----------------------------------------------------------------
    pft = 0
    if ( myinterface%params_siml%lTrE ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'tre' )
    end if

    if ( myinterface%params_siml%lTNE ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'tne' )
    end if

    if ( myinterface%params_siml%lTrD ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'trd' )
    end if

    if ( myinterface%params_siml%lTND ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'tnd' )
    end if

    if ( myinterface%params_siml%lGr3 ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'gr3' )
    end if

    if ( myinterface%params_siml%lGN3 ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'gn3' )
    end if

    if ( myinterface%params_siml%lGr4 ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'gr4' )
    end if

    npft_site = pft
    ! if (npft_site==0) stop 'PLANT:GETPAR_MODL_PLANT: PFT name not valid. See run/<simulationname>.sofun.parameter'

  end subroutine getpar_modl_plant


  function getpftparams( pftname ) result( out_getpftparams )
    !----------------------------------------------------------------
    ! Read PFT parameters from respective file, given the PFT name
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: pftname

    ! local variables
    real :: lu_category_prov = 0   ! land use category associated with PFT (provisional)

    ! function return variable
    type( params_pft_plant_type ) :: out_getpftparams

    ! standard PFT name
    out_getpftparams%pftname = pftname

    ! PFT names
    ! Gr3 : C3 grass                          
    ! Gr4 : C4 grass     
    if (trim(pftname)=='gr3') then
      out_getpftparams%grass   = .true.
      out_getpftparams%tree    = .false.
      out_getpftparams%c3      = .true.
      out_getpftparams%c4      = .false.
      out_getpftparams%nfixer  = .false.
    else if (trim(pftname)=='gn3') then
      out_getpftparams%grass   = .true.
      out_getpftparams%tree    = .false.
      out_getpftparams%c3      = .true.
      out_getpftparams%c4      = .false.
      out_getpftparams%nfixer  = .true.
    else if (trim(pftname)=='gr4') then
      out_getpftparams%grass   = .true.
      out_getpftparams%tree    = .false.
      out_getpftparams%c3      = .false.
      out_getpftparams%c4      = .true.
      out_getpftparams%nfixer  = .false.
    else if (trim(pftname)=='tre') then
      out_getpftparams%grass   = .false.
      out_getpftparams%tree    = .true.
      out_getpftparams%c3      = .true.
      out_getpftparams%c4      = .false.
      out_getpftparams%nfixer  = .false.
    else if (trim(pftname)=='tne') then
      out_getpftparams%grass   = .false.
      out_getpftparams%tree    = .true.
      out_getpftparams%c3      = .true.
      out_getpftparams%c4      = .false.
      out_getpftparams%nfixer  = .true.
    else if (trim(pftname)=='tnd') then
      out_getpftparams%grass   = .false.
      out_getpftparams%tree    = .true.
      out_getpftparams%c3      = .true.
      out_getpftparams%c4      = .false.
      out_getpftparams%nfixer  = .true.
    end if      

    ! land use category associated with PFT (provisional) 
    if (lu_category_prov==1) then
      out_getpftparams%lu_category = lunat
      out_getpftparams%islu(lunat) = .true.
    else
      out_getpftparams%islu(lunat) = .false.
    end if

    ! ! leaf mass per area (gC m-2)
    ! out_getpftparams%lma = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'lma' )
    ! out_getpftparams%sla = 1.0 / out_getpftparams%lma

    ! ! constant ratio of leaf structural N to LMA
    ! out_getpftparams%r_ntolma = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'r_ntolma' )

    ! leaf decay constant, read in as [years-1], central value: 0.0 yr-1 for deciduous plants
    out_getpftparams%k_decay_leaf_base = myinterface%params_calib%k_decay_tissue / ndayyear 

    ! shape parameter for turnover function if LAI
    out_getpftparams%k_decay_leaf_width = myinterface%params_calib%k_decay_leaf_width

    ! sapwood decay constant [days], read in as [years-1], central value: xxx
    out_getpftparams%k_decay_sapw =  myinterface%params_calib%k_decay_sapw / ndayyear 

    ! root decay constant [days], read in as [years-1], central value: 1.04 (Shan et al., 1993; see Li et al., 2014)  
    out_getpftparams%k_decay_root = myinterface%params_calib%k_decay_tissue / ndayyear 

    ! root C:N and N:C ratio (gC/gN and gN/gC)
    out_getpftparams%r_cton_root = myinterface%params_calib%r_cton_root
    out_getpftparams%r_ntoc_root = 1.0 / out_getpftparams%r_cton_root

    ! y-axis intersection in the relationship of non-metabolic versus metabolic N per leaf area
    out_getpftparams%ncw_min = myinterface%params_calib%ncw_min

    ! slope in the relationship of non-metabolic versus metabolic N per leaf area
    out_getpftparams%r_n_cw_v = myinterface%params_calib%r_n_cw_v

    ! constant ratio of C to structural N
    out_getpftparams%r_ctostructn_leaf = myinterface%params_calib%r_ctostructn_leaf

  end function getpftparams


  subroutine add_seed( plant )
    !//////////////////////////////////////////////////////////////////
    ! To initialise plant pools, add "sapling" mass
    !------------------------------------------------------------------
    use md_classdefs

    ! arguments
    type(plant_type), intent(inout) :: plant

    plant%plabl = orgplus( plant%plabl, seed )

  end subroutine add_seed


  subroutine initglobal_plant( plant )
    !////////////////////////////////////////////////////////////////
    !  Initialisation of all _pools on all gridcells at the beginning
    !  of the simulation.
    !  June 2014
    !  b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    ! argument
    type( plant_type ), dimension(npft), intent(inout) :: plant

    ! local variables
    integer :: pft

    !-----------------------------------------------------------------------------
    ! derive which PFTs are present from fpc_grid (which is prescribed)
    !-----------------------------------------------------------------------------
    do pft=1,npft
      call initpft( plant(pft) )
      plant(pft)%pftno = pft
    end do

  end subroutine initglobal_plant


  subroutine init_annual_plant( plant_fluxes )
    !////////////////////////////////////////////////////////////////
    ! Set (iterative) annual sums to zero
    !----------------------------------------------------------------
    ! arguments
    type(plant_fluxes_type), dimension(npft), intent(inout) :: plant_fluxes

    ! plant_fluxes(:)%dharv = 0.0

  end subroutine init_annual_plant


  subroutine initpft( plant )
    !////////////////////////////////////////////////////////////////
    !  Initialisation of specified PFT on specified gridcell
    !  June 2014
    !  b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    ! argument
    type( plant_type ), intent(inout) :: plant

    plant%fpc_grid  = 0.0
    plant%lai_ind   = 0.0
    plant%fapar_ind = 0.0
    plant%acrown    = 0.0

    ! canpopy state variables
    plant%narea            = 0.0
    plant%narea_metabolic  = 0.0
    plant%narea_structural = 0.0
    plant%lma              = 0.0
    plant%sla              = 0.0
    plant%nmass            = 0.0
    plant%r_cton_leaf      = 0.0
    plant%r_ntoc_leaf      = 0.0

    plant%pleaf = orgpool(carbon(0.0), nitrogen(0.0))
    plant%proot = orgpool(carbon(0.0), nitrogen(0.0))
    plant%psapw = orgpool(carbon(0.0), nitrogen(0.0))
    plant%pwood = orgpool(carbon(0.0), nitrogen(0.0))
    plant%plabl = orgpool(carbon(0.0), nitrogen(0.0))
    
  end subroutine initpft


  function ftemp( temp, method, ref_temp )
    !////////////////////////////////////////////////////////////////
    ! Generic temperature response function
    !----------------------------------------------------------------
    ! arguments
    real, intent(in)             :: temp ! temperature [in Degrees Celsius]
    character(len=*), intent(in) :: method
    real, intent(in), optional   :: ref_temp

    ! local variables
    real                         :: ref_temp_local  ! local copy of ref_temp

    ! for lloyd and taylor method
    real, parameter :: E0 = 308.56      ! Activation Energy
    real, parameter :: T0 = 227.13      ! calibration temperature [K]
    real, parameter :: Tzero = 273.15   ! 0°C = 273.15 K 

    ! function return variable
    real :: ftemp

    ! set default reference temperature to 10 deg C
    if (present(ref_temp)) then
     ref_temp_local = ref_temp
    else
     ref_temp_local = 10.0
    endif

    select case (method)

      case ("lloyd_and_taylor")
        !----------------------------------------------------------------
        ! LLOYD AND TAYLOR
        ! Temperature response function is a modified Q10 relationship
        ! (Lloyd & Taylor 1994)
        !----------------------------------------------------------------
        if (temp.ge.-40.0) then 
          ! avoid numerical errors
          ftemp = exp(E0*((1.0/(ref_temp_local+Tzero-T0))-(1.0/(temp+Tzero-T0))))
        else
          ! set temperature response to a constant at value of -40°C
          ftemp = exp(E0*((1.0/(ref_temp_local+Tzero-T0))-(1.0/(-40.0+Tzero-T0))))
        end if

      case default

        stop 'FTEMP: select valid method'

    end select

    return

  end function ftemp


  function fmoist( moist, method )
    !////////////////////////////////////////////////////////////////
    ! Generic moisture response function
    !----------------------------------------------------------------
    ! arguments
    real, intent(in)             :: moist ! temperature [in Degrees Celsius]
    character(len=*), intent(in) :: method

    ! function return variable
    real :: fmoist

    select case (method)

      case ("foley")
        !----------------------------------------------------------------
        ! FOLEY
        ! Calculates decomposition rate modifier for a given water fraction
        ! according to Foley 1995
        !----------------------------------------------------------------
        fmoist = 0.25 + 0.75 * moist

      case default

        stop 'FMOIST: select valid method'

    end select

    return

  end function fmoist

end module md_plant
