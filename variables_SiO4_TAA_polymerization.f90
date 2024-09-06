module variables
  implicit none
  
  !dimensions of the lattice
  integer :: lx,ly,lz

  !tag_site: array for site number <--> site label transformation
  integer,dimension(:),allocatable :: tag_site

  !arrays storing the neighbor list:
  !n1list: first neighbor list
  !n2list: second neighbor list
  !n3list: third neighbor listweq
  integer,dimension(:,:),allocatable :: n1list,n2list,n3list,n4list,n5list

  !nsites:total sites
  !natom: total atoms 
  !allmonomers: all the monomers on the lattice
  integer :: nsites,natom,allmolecules,unitcells

  !initial configuration flag,restart code flag, MC step, initial MC step
  integer :: initial_config_flag,step,initial_step
  !forerunnung string for configuration files
  character*10 :: snapshot_
  
  !coordinates of the particles 
  integer,dimension(:),allocatable :: rx,ry,rz
  
  !arrays storing the occupancy of sites
  integer,dimension(:),allocatable :: occupancy

  !occupancies parameters
  integer,parameter :: occSI = 7    !occupancy of Si in SI
  integer,parameter :: occSN = 2    !occupancy of Si in SN
  integer,parameter :: occW = 0     !occupancy of water(W)
  integer,parameter :: occTAA = 13   !occupancy of TAA site
  integer,parameter :: occISIO = 113 !occupancy of oxygen in SI
  integer,parameter :: occSNO = 17   !occupancy of oxygen in SN
  integer,parameter :: occSIO = 51   !occupancy of oxygen in SI
  !bridging O- between SI and SN
  integer,parameter :: occSIOSNO = occSIO + occSNO
  !bridging O- between SI and SI
  integer,parameter :: occSIOSIO = occSIO + occSIO
  !bridging O- between SN and SN
  integer,parameter :: occSNOSNO = occSNO + occSNO
  
  !array storing the identity of molecules on sites
  !1:SI,2:SN,3:TAA,4:Water(W)
  integer,dimension(:),allocatable :: spin
  integer,parameter :: spinSN = 2,spinSI = 1,spinW = 4,spinTAA = 3

  !noccupy:array storing the site of Si molecules
  !head:target of all the sites
  !resite: transforms isite to monomer no.
  integer,dimension(:),allocatable :: noccupy, head, resite

  !clabel: stores the labels of the clusters
  !csize: stores the size of the corresponding cluster
  integer,dimension(:),allocatable :: clabel,csize
  integer :: max_size,cluster_num,cluster_threshold
  real*8 :: avg_size
  
  !array storing the orientation detail of the monomer unit
  integer,dimension(4,8) :: Si_O,Si_notO
  
  !variables for the no of sweeps,equilibrium sweeps
  !and snapshot taking variables
  integer :: nsweeps,neqsweeps,nprint,nsnapshot,point,snapshot
  
  !temperature,total energy variable
  real*8 :: tstar
  
  !variable to store the charge per SI
  real :: chargeperSI,sublat_ordering,frac_SITAA,frac_SISNTAA
  
  !Qn distribution variables
  integer,dimension(0:4) :: Qn,Qni,Qnn
  
  !rings size distribution variables
  !integer :: rings3,rings4
  
  integer :: seed!seed for the random number generator
  integer,parameter :: nc=8!coordination number
  real*8 :: c,cn!degree of condensation
  
  !arrays for storing concentrations and the number of molecules
  !indices same as the spin number
  real*8,dimension(4) :: xi
  integer,dimension(4) :: ni
  real :: nTEOS,nTAAOH,nH2O,pH
  
  !penalties on 3 and 4 membered rings
  real*8 :: pen3, pen4
  
  !interaction energy strengths
  !eISITAA: interaction energy b/w ISI and TAA
  !eSNTAA: interaction energy b/w SN and TAA
  !eSNSN:     "         "    "   SN and SN
  !eSISN      "         "    "   SI and SN
  !eSISI      "         "    "   SI and S
  !eTAATAA        "         "    "   TAA and TAA
  real*8 :: eISITAA,eSNTAA,eSNSN,eSISN,eSISI,eTAATAA

  !variable for the temperature profile change
  real*8 :: tstar_initial,tstar_final
  real*8 :: up_rate,down_rate
  real :: m0,m1,m2,m3
  
  !final averages of various variables
  real :: final_avg_energy=0.00
  real :: final_avg_max_size=0.00,final_avg_cluster_num=0.00,final_avg_avg_size=0.00
  
  !time variables
  real :: start_time,end_time,runtime,time_limit
  
  !book keeping method variables
  integer,dimension(:),allocatable :: SIlink,TAAlink
  
  integer :: rotation_attempted=0, rotation_accepted=0
 
end module variables
  
  
