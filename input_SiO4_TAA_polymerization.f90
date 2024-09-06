module input_SiO4_TAA  
  use variables     
  implicit none
contains
  subroutine input
    implicit none
    integer :: ioerr,status
    integer :: i,j,manual_concen_flag
    real :: a,b,c,alpha

    point = 11
    cluster_threshold = 15
    
    open(101,file='input.dat',status='unknown',iostat=ioerr)
    if(ioerr.ne.0) then
       print*,"Error: Reading input.dat file"
       stop
    end if
    
    !read lattice dimensions
    read(101,*,iostat=status) lx, ly, lz
    if(status.ne.0)then
       print*,'error reading lattice dimensions'
       stop
    end if

    !total number of sites and unit cells
    nsites = 2*lx*ly*lz
    unitcells = lx*ly*lz

    !read the temperature 
    read(101,*,iostat=status) tstar
    if(tstar.lt.0.0d0.or.status.ne.0)then
       print*,'error reading tempersture'
       stop
    end if

    !read the no. of various sweeps
    read(101,*,iostat=status) nsweeps,nprint,nsnapshot,neqsweeps
    if(status.ne.0)then
       print*,'error reading sweeps'
       stop
    end if
    !neqsweeps = int(0.7 * nsweeps)

    !read the concentrations
    read(101,*,iostat=status)nTEOS,nTAAOH,nH2O,manual_concen_flag
    if(status.ne.0)then
       print*,'error reading composition values'
       stop
    end if
    !if manual_concen_flag is zero then
    !read the numbers in the next line as nTEOS,nTAAOH,nH2O
    !if the flag is non zero then
    !nTEOS => nSN
    !nTAAOH => nSI
    !nH2O => nTAA
    !nSN,nSI,nTAA are the molecule number of the respective species
            
    if(manual_concen_flag.eq.0)then!if manual concentration flag is zero

       !do the usual calculation of pH, ni, xi
       !call molar_pH(nTEOS,nTAAOH,nH2O,a,b,c,alpha,pH)
       call molar_to_pH(alpha)
       call molar_fractions(alpha,xi(1),xi(2),xi(3),xi(4))
       call fractions_number(xi(1),xi(2),xi(3),xi(4),ni(1),ni(2),ni(3),ni(4))
       
    elseif(manual_concen_flag.eq.1)then!if the manual concentration flag is non zero

       !negelct the nH2O read from the file. only the first two values of nTEOS and nTAAOH
       !are important. Refer to job_queue_molecules.sh for more information.
       ni(1) = int(0)
       ni(2) = nTEOS
       ni(3) = nTAAOH
       ni(4) = nsites - (ni(1)+ni(2)+ni(3))
       xi = 1.0d0*ni/real(unitcells)
       xi(4) = 1.00000d0 - (xi(1)+xi(2)+xi(3))
       
    else
       print*,'ERROR: incorrect manual_config_flag :: subroutine input'
       stop
    end if!endif manual concentration flag is zero
        
    !particle numbers for manual_check
    !make sure that the system size is 60,60,60
    ni(1) = 0
    ni(2) = 2
    ni(3) = 0
    ni(4) = nsites - ni(2) - ni(1) - ni(3)
        
    read(101,*,iostat=status) pen3, pen4
    if(status.ne.0) then
       print*, 'error: reading ring penalties'
       stop
    end if
    
    read(101,*,iostat=status) eISITAA,eSNSN,eSISN
    if(status.ne.0) then
       print*, 'error: reading ring penalties'
       stop
    end if
    
    read(101,*,iostat=status) tstar_final,m0,m1,m2
    if(status.ne.0)then
       print*,'error: reading temperature profile data'
       stop
    end if
    tstar_initial=tstar
    up_rate = real((tstar_final-tstar_initial)/(m1-m0))
    down_rate = 8.000

    read(101,*,iostat=status) initial_config_flag,time_limit
    if(status.ne.0)then
       print*,'error: reading initial configuration flag'
       stop
    end if
    
    close(101)

    !total atoms on the lattice
    natom = ni(1)*5 + ni(2)*5 + ni(3)*5
    !all the monomers/molecules on the lattice
    allmolecules = ni(1) + ni(2) + ni(3)
    !allmonomers = allmolecules - ni(3)
        
    allocate(occupancy(nsites))
    allocate(tag_site(nsites * 4))
    
    allocate(noccupy(allmolecules))
    allocate(resite(nsites))
    allocate(head(nsites))
    allocate(spin(nsites))
    allocate(clabel(nsites))
    allocate(csize(nsites))
    
    allocate(n1list(8,nsites))
    allocate(n2list(6,nsites))
    allocate(n3list(12,nsites))
    allocate(n4list(24,nsites))
    allocate(n5list(8,nsites))    

    allocate(rx(nsites))
    allocate(ry(nsites))
    allocate(rz(nsites))
    
    allocate(SIlink(nsites))
    allocate(TAAlink(nsites))
    
  end subroutine input
  !------------------------------------------------------------
  !--- molar_pH(x,y,z,pH): converts molar ratios to pH values
  !------------------------------------------------------------
  subroutine molar_pH(x,y,z,a,b,c,alpha,pH)
    implicit none
    
    real,intent(in) :: x,y,z
    real,intent(out) :: a,b,c,alpha,pH
    real :: alpha1,alpha2,temp
    real :: tempa,tempb,tempc
    real,parameter :: k=1.75e+06,nw=55.555,pKw=13.8

    a = x/(2*x + 5*y + z)
    b = y/(2*x + 5*y + z)
    c = z/(2*x + 5*y + z)
    
    tempa = x/(2*x + 5*y + z)
    tempb = y/(2*x + 5*y + z)
    tempc = z/(2*x + 5*y + z)
    
    alpha1 = (k*tempa+k*tempb+tempc + sqrt((k*tempa+k*tempb+tempc)**2 - 4*k*tempa*tempb*(k-1)))/(2*(k-1))
    alpha2 = (k*tempa+k*tempb+tempc - sqrt((k*tempa+k*tempb+tempc)**2 - 4*k*tempa*tempb*(k-1)))/(2*(k-1))
    
    temp = min(tempa,tempb)
    
    if((alpha1.lt.temp).and.(alpha2.gt.temp))then
       alpha = alpha1
    elseif((alpha1.gt.temp).and.(alpha2.lt.temp))then
       alpha = alpha2
    elseif(y.eq.real(0))then!@ iso-electric point
       !do nothing
    else
       print*,'Error in calculations'
       print*,'tempa',tempa
       print*,'tempb',tempb
       print*,'tempc',tempc
       print*,'alpha1',alpha1
       print*,'alpha2',alpha2
       print*,'alpha',alpha
       print*,'tempa-alpha',tempa-alpha
       print*,'tempb-alpha',tempb-alpha
       stop
    end if
    
    if(y.gt.real(0))then
       temp = temp - alpha
       pH = pKw + log10(nw * temp)
    else!iso-electric point
       pH = 2.5
    end if
    
    return
  end subroutine molar_pH
  !------------------------------------------------------------                       
  !--- molar_pH(x,y,z,pH): converts molar ratios to pH values                         
  !------------------------------------------------------------                       
  subroutine molar_to_pH(alpha)
    use variables
    implicit none

    real,intent(out) :: alpha
    real :: y, x, z!mole proportions                                                  
    real :: a1,b1,c1!quadratic equation coefficients                                  
    real*8,parameter :: c0=real(55.55),K=real(1.75*10**6),pKw = real(14)
    real :: alpha1,alpha2,beta,pOH

    y = nTEOS;x=nTAAOH;z=nH2O

    if(x.eq.real(0))then
       pH = real(2.5)
       alpha = real(0)
       return
    end if


    !concentrated system                                                              
    !quadratic equation is: (K-1)alpha^2 - (x*K + y*K + 4*y + z)alpha + x*y*k = 0     
    a1 = K-real(1)
    b1 = real(-1)*(x*K + y*K + real(4)*y + z)
    c1 = x*y*K

    !dilute system                                                                    
    !quadratic equation is: Kalpha^2 - (x*K+y*K+1) + K*x*y = 0                        
    !a1 = K                                                                           
    !b1 = real(-1)*(x*K + y*K + 1)                                                    
    !c1 = x*y*K                                                                       

    alpha1 = (-b1 + sqrt(b1**2 - real(4)*a1*c1))/real(2*a1)
    alpha2 = (-b1 - sqrt(b1**2 - real(4)*a1*c1))/real(2*a1)

    alpha = min(alpha1,alpha2)

    pH = pKw + log10((x - alpha) * c0)

    if((alpha.lt.0).or.(alpha.gt.x))then
       print*,'alpha1',alpha1
       print*,'alpha2',alpha2
       print*,'x',x
       print*,'y',y
       print*,'z',z
       stop
    end if

    return
  end subroutine molar_to_pH
  !------------------------------------------------------------
  !--- molar_fractions(y,x,z,alpha,xSI,xSN,xSDA,xH2O): converts molar ratios to fractions
  !------------------------------------------------------------
  subroutine molar_fractions(alpha,xSI,xSN,xSDA,xH2O)
    implicit none
    
    real,intent(in) :: alpha
    real :: y,x,z
    real*8,intent(out) :: xSI,xSN,xSDA,xH2O
    
    y = nTEOS;x=nTAAOH;z=nH2O
    
    if(x.eq.real(0))then
       
       xSN = (y)/(real(5)*y + z)
       xSI = real(0)
       xSDA = real(0)
       xH2O = real(1) - xSN - xSI - xSDA
       
    else
       
       xSN = (y-alpha)/(x + real(5)*y + z)
       xSI = (alpha)/(x + real(5)*y + z)
       xSDA = (x)/(x + real(5)*y + z)
       xH2O = real(1) - xSN - xSI - xSDA
       
    end if
    
    return
  end subroutine molar_fractions
  !------------------------------------------------------------
  !--- fractions_number(xSI,xSN,xSDA,xH2O,nSI,nSN,nSDA,nH2O): converts fractions to numbers
  !------------------------------------------------------------
  subroutine fractions_number(xSI,xSN,xSDA,xH2O,nSI,nSN,nSDA,nH2O)
    implicit none 

    integer,intent(out) :: nSI,nSN,nSDA,nH2O
    real*8,intent(in) :: xSI,xSN,xSDA,xH2O
    integer :: nsites
    real*8 :: temp_SI,temp_SN,temp_SDA
    
    nsites = real(lx) * real(ly) * real(lz)
    temp_SI = xSI * nsites
    nSI = nint(temp_SI)
    temp_SN = xSN * nsites
    nSN = nint(temp_SN)
    temp_SDA = xSDA * nsites
    nSDA = nint(temp_SDA)
    nH2O = nsites - nSI - nSN - nSDA
        
    return
  end subroutine fractions_number
  !*************************************************************
end module input_SiO4_TAA
