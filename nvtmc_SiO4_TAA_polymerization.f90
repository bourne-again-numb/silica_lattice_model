program NVT_SiO4_TAA
  use variables
  use input_SiO4_TAA
  implicit none
  
  integer :: opstatus,count
  logical :: restart_file_exists
  logical,parameter :: restart_flag = .true.
  
  !if restart_flag = .true. then make the code restartable
  !else don't
  
  !initialize random seed and start the system clock
  call system_clock(count)
  seed = - mod(1.0d0*count,1.0d4)
  !seed = -2
  call cpu_time(start_time)

  inquire(file="system_state.cfg",exist=restart_file_exists)
  
  !***********************RESTARTING THE PROGRAM FLAG********************************
  if(restart_file_exists.and.restart_flag)then!if the file restarting file exists
     
     !create output files
     call create_files_append
     !read all the vaules for the variables
     call input
     !setup the BCC lattice and its neighbor list
     call lattice
     
     !read the state of the system from system_state.cfg
     call read_state
     
     !start the NVT MC algorithm
     call nvt_mc
     call cpu_time(end_time)
     runtime = (end_time-start_time)/3600.0d0
     
     !Print simulation summary
     write(*,8000),'Program runtime(hrs):',runtime
     write(1008,8000)'Program runtime(hrs):',runtime
     write(*,*);write(*,*)
     
     !output the final state of the syste
     call output_state
     
     !deallocate all the allocated memory
     call freeandclose
     !************************************************************
  else!if the file restarting file doesn't exist
     
     initial_step = 1
     snapshot = 0
     
     !create output files
     call create_files_onetime
     !read all the vaules for the variables
     call input
     !setup the BCC lattice and its neighbor list
     call lattice
          
     !INITIAL CONFIGURATION OF THE LATTICE
     if(initial_config_flag.eq.0)then
        call random_configuration!random configuration
        !call manual_check
     elseif(initial_config_flag.eq.1)then
        if((xi(2)+xi(4)).eq.1.00)then
           call initial_cristo_config!beta-crystabolite configuration
        else
           print*,'BUG:: incorrect concen for crystobalite configurations'
           call freeandclose
           stop
        end if
        !elseif(initial_config_flag.eq.2)then
        !   call read_configuration
     else
        print*,'Error: incorrect value for initial_configuration_flag :: main program'
        call freeandclose
        stop
     end if
     
     !Print simulation input conditions
     write(*,1000),'-----MC simulation starts-----'
     write(*,2000),'System size:',lx,ly,lz
     write(*,3000),'Concentration(SI/SN/TAA/W)//pH:',xi(1),xi(2),xi(3),xi(4),pH
     write(*,4000),'Molecules(SI/SN/TAA/W):',ni(1),ni(2),ni(3),ni(4)
     write(*,6000),'Temperature:',tstar
     write(*,5000),'Ring penalties(3-rings/4rings):',pen3,pen4
     write(*,1000),'Interaction strengths:'
     write(*,6000),'SI-SN:',eSISN
     write(*,6000),'SN-SN:',eSNSN
     write(*,6000),'ISI-TAA:',eISITAA
     write(*,6000),'SI-SI:',eSISI
     write(*,6000),'TAA-TAA:',eTAATAA
     write(*,*)
     write(1008,1000)'MC simulation starts';flush(1008)
     write(1008,2000)'System size:',lx,ly,lz;flush(1008)
     write(1008,3000),'Concentration(SI/SN/TAA/W)//pH:',xi(1),xi(2),xi(3),xi(4),pH;flush(1008)
     write(1008,4000)'Molecules(SI/SN/TAA/W):',ni(1),ni(2),ni(3),ni(4);flush(1008)
     write(1008,6000)'Temperature:',tstar;flush(1008)
     write(1008,5000)'Ring penalties(3-rings/4rings):',pen3,pen4;flush(1008)
     write(1008,1000)'Interaction strengths:';flush(1008)
     write(1008,6000)'SI-SN:',eSISN;flush(1008)
     write(1008,6000)'SN-SN:',eSNSN;flush(1008)
     write(1008,6000)'ISI-TAA:',eISITAA;flush(1008)
     write(1008,6000),'SI-SI:',eSISI;flush(1008)
     write(1008,6000),'TAA-TAA:',eTAATAA;flush(1008)
     write(1008,*);flush(1008)
     
     !start the NVT MC algorithm
     call nvt_mc
     call cpu_time(end_time)
     runtime = (end_time-start_time)/3600.00

     !Print simulation summary
     write(*,8000),'Program runtime(hrs):',runtime
     write(1008,8000)'Program runtime(hrs):',runtime
     write(*,*);write(*,*)
     
     !output the final state of the syste
     call output_state
               
     !deallocate all the allocated memory
     call freeandclose
     
  end if!if the file restarting file exists
  !***********************************************************************
  
1000 format(x,A)
2000 format(x,A,i8,i8,i8)
3000 format(x,A,2x,f6.4,2x,f6.4,2x,f6.4,2x,f6.4,2x,f6.4)
4000 format(x,A,x,i8,x,i8,x,i8,x,i8)
5000 format(x,A,2x,f6.4,2x,f6.4)
6000 format(x,A,x,f7.4)
7000 format(x,A,i8,A,i8)
8000 format(A,es10.2)
9000 format(x,i10,2x,f6.4,2x,f6.4,2x,f6.4,2x,f6.4,2x,f7.4,2x,f7.4,2x,&
          f7.4,2x,f7.4,2x,f10.4,2x,f7.4)
contains
  !---------------------------------------------------------
  !--- CREATE_FILES: opens up the output files
  !---------------------------------------------------------
  subroutine create_files_append
    implicit none
    
    open(1001,file='q.c_mcsteps.csv',status='unknown',access='append',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'ERROR: opening Qn distribution file'
       call freeandclose
       stop
    end if
    
    open(1002,file='rings_dist.csv',status='unknown',access='append',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'ERROR: opening Qn distribution file'
       call freeandclose
       stop
    end if
    
    open(1003,file='lattice_energy.csv',status='unknown',access='append',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'ERROR: opening lattice_energy.csv'
       call freeandclose
       stop
    end if
    
    open(1006,file='vmd.xyz',status='unknown',access='append',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'Error: opening vmd.xyz file'
       call freeandclose
       stop
    end if
    
    open(1007,file='cluster_stat.csv',status='unknown',access='append',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'Error: opening cluster_stat.csv file'
       call freeandclose
       stop
    end if
    
    open(1008,file='screen_output.log',status='unknown',access='append',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'Error: opening screen_output.log file'
       call freeandclose
       stop
    end if
    
    open(1009,file='temperature_profile.csv',status='unknown',access='append',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'Error: opening temperature_profile.csv file'
       call freeandclose
       stop
    end if
    
    open(1012,file='cluster_order.csv',status='unknown',access='append',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'Error: opening cluster_order file'
       call freeandclose
       stop
    end if
    
    open(1014,file='summary.csv',status='unknown',access='append',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'Error: opening summary file'
       call freeandclose
       stop
    end if
    
    open(1015,file='qn.c_mcsteps.csv',status='unknown',access='append',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'Error: opening qn.c_mcsteps file'
       call freeandclose
       stop
    end if
    
    open(1018,file='size_range.csv',status='unknown',access='append',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'Error: opening size_range.csv file'
       call freeandclose
       stop
    end if    

  end subroutine create_files_append
  !---------------------------------------------------------
  !--- CREATE_FILES: opens up the output files
  !---------------------------------------------------------
  subroutine create_files_onetime
    implicit none
    
    open(1001,file='q.c_mcsteps.csv',status='unknown',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'ERROR: opening Qn distribution file'
       call freeandclose
       stop
    end if
    
    open(1002,file='rings_dist.csv',status='unknown',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'ERROR: opening Qn distribution file'
       call freeandclose
       stop
    end if
    
    open(1003,file='lattice_energy.csv',status='unknown',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'ERROR: opening lattice_energy.csv'
       call freeandclose
       stop
    end if
    
    open(1006,file='vmd.xyz',status='unknown',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'Error: opening vmd.xyz file'
       call freeandclose
       stop
    end if
    
    open(1007,file='cluster_stat.csv',status='unknown',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'Error: opening cluster_stat.csv file'
       call freeandclose
       stop
    end if
    
    open(1008,file='screen_output.log',status='unknown',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'Error: opening screen_output.log file'
       call freeandclose
       stop
    end if
    
    open(1009,file='temperature_profile.csv',status='unknown',access='append',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'Error: opening temperature_profile.csv file'
       call freeandclose
       stop
    end if
    
    open(1012,file='cluster_order.csv',status='unknown',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'Error: opening cluster_order file'
       call freeandclose
       stop
    end if
    
    open(1014,file='summary.csv',status='unknown',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'Error: opening summary file'
       call freeandclose
       stop
    end if
    
    open(1015,file='qn.c_mcsteps.csv',status='unknown',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'Error: opening qn.c_mcsteps file'
       call freeandclose
       stop
    end if
    
    open(1018,file='size_range.csv',status='unknown',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'Error: opening size_range.csv file'
       call freeandclose
       stop
    end if    

  end subroutine create_files_onetime
  !---------------------------------------------------------
  !--- LATTICE: CONSTRUSTS THE NEIGHBOR LIST
  !---------------------------------------------------------
  subroutine lattice
    use variables
    implicit none
    
    integer :: ix, iy, iz, itag, icheck, isite
    integer :: ix1,ix2,iy1,iy2,iz1,iz2
    integer :: ixs1,ixs2,iys1,iys2,izs1,izs2
    integer :: ixt1,ixt2,iyt1,iyt2,izt1,izt2
    
    !initialize the coordinates array and the tag to site transformation array
    isite=0
    do ix=1,2*lx
       do iy=1,2*ly
          do iz=1,2*lz
             icheck=mod(ix,2)+mod(iy,2)+mod(iz,2)
             if(icheck.eq.0.or.icheck.eq.3) then
                isite = isite + 1
                itag = (ix-1)*2*ly*2*lz + (iy-1)*2*lz + iz
                rx(isite) = ix
                ry(isite) = iy
                rz(isite) = iz
                tag_site(itag) = isite
             end if
          end do
       end do
    end do
            
    !Potential Neighbor List Errors
    if(isite.ne.nsites) then
       print*,"Error:nlist(site.neq.sites)::subroutine lattice",isite,nsites
       call freeandclose
       stop
    end if
    
    !assign the first,second and third nearest neighbors
    isite=0
    do ix = 1,2*lx
       do iy = 1,2*ly
          do iz = 1,2*lz
             icheck = mod(ix,2) + mod(iy,2) + mod(iz,2)
             if(icheck.eq.0.or.icheck.eq.3) then
                isite = isite + 1
                itag = (ix-1)*2*ly*2*lz + (iy-1)*2*lz + iz

                ! first neighbor periodic boundary conditions
                ix1 = ix - 1
                if(ix1.lt.1) ix1 = 2*lx
                ix2 = ix + 1
                if(ix2.gt.2*lx) ix2 = 1
                iy1 = iy - 1
                if(iy1.lt.1) iy1 = 2*ly
                iy2 = iy + 1
                if(iy2.gt.2*ly) iy2 = 1
                iz1 = iz - 1
                if(iz1.lt.1) iz1 = 2*lz
                iz2 = iz + 1
                if(iz2.gt.2*lz) iz2 = 1

                ! first neighbors: within a dist. of 'a/(2)^0.5' frm a site
                n1list(1,isite) = tag_site((ix1-1)*2*ly*2*lz + (iy1-1)*2*lz + iz2)
                n1list(2,isite) = tag_site((ix2-1)*2*ly*2*lz + (iy1-1)*2*lz + iz2)
                n1list(3,isite) = tag_site((ix2-1)*2*ly*2*lz + (iy2-1)*2*lz + iz2)
                n1list(4,isite) = tag_site((ix1-1)*2*ly*2*lz + (iy2-1)*2*lz + iz2)
                n1list(5,isite) = tag_site((ix1-1)*2*ly*2*lz + (iy1-1)*2*lz + iz1)
                n1list(6,isite) = tag_site((ix2-1)*2*ly*2*lz + (iy1-1)*2*lz + iz1)
                n1list(7,isite) = tag_site((ix2-1)*2*ly*2*lz + (iy2-1)*2*lz + iz1)
                n1list(8,isite) = tag_site((ix1-1)*2*ly*2*lz + (iy2-1)*2*lz + iz1)
                
                ixs1 = ix-2
                ixs2 = ix+2
                iys1 = iy-2
                iys2 = iy+2
                izs1 = iz-2
                izs2 = iz+2
                !second/third/fifth neighbor periodic boundary conditions
                if(ixs1.lt.1) ixs1 = ixs1 + 2*lx
                if(ixs2.gt.2*lx) ixs2 = ixs2 - 2*lx
                if(iys1.lt.1) iys1 = iys1 + 2*ly
                if(iys2.gt.2*ly) iys2 = iys2 - 2*ly
                if(izs1.lt.1) izs1 = izs1 + 2*lz
                if(izs2.gt.2*lz) izs2 = izs2 - 2*lz
                !second neighbor: within a dist. of 'a' from a site
                n2list(1,isite) = tag_site((ixs1-1)*2*ly*2*lz + (iy-1)*2*lz + iz)
                n2list(2,isite) = tag_site((ixs2-1)*2*ly*2*lz + (iy-1)*2*lz + iz)
                n2list(3,isite) = tag_site((ix-1)*2*ly*2*lz + (iys1-1)*2*lz + iz)
                n2list(4,isite) = tag_site((ix-1)*2*ly*2*lz + (iys2-1)*2*lz + iz)
                n2list(5,isite) = tag_site((ix-1)*2*ly*2*lz + (iy-1)*2*lz + izs1)
                n2list(6,isite) = tag_site((ix-1)*2*ly*2*lz + (iy-1)*2*lz + izs2)
                
                !third neigbors within a distance of '2^0.5*a' from a site
                n3list(1,isite) = tag_site((ix-1)*2*ly*2*lz + (iys1-1)*2*lz + izs2)
                n3list(2,isite) = tag_site((ixs2-1)*2*ly*2*lz + (iys1-1)*2*lz + iz)
                n3list(3,isite) = tag_site((ix-1)*2*ly*2*lz + (iys1-1)*2*lz + izs1)
                n3list(4,isite) = tag_site((ixs1-1)*2*ly*2*lz + (iys1-1)*2*lz + iz)
                n3list(5,isite) = tag_site((ixs1-1)*2*ly*2*lz + (iy-1)*2*lz + izs2)
                n3list(6,isite) = tag_site((ixs2-1)*2*ly*2*lz + (iy-1)*2*lz + izs2)
                n3list(7,isite) = tag_site((ixs2-1)*2*ly*2*lz + (iy-1)*2*lz + izs1)
                n3list(8,isite) = tag_site((ixs1-1)*2*ly*2*lz + (iy-1)*2*lz + izs1)
                n3list(9,isite) = tag_site((ixs1-1)*2*ly*2*lz + (iys2-1)*2*lz + iz)
                n3list(10,isite) = tag_site((ix-1)*2*ly*2*lz + (iys2-1)*2*lz + izs2)
                n3list(11,isite) = tag_site((ixs2-1)*2*ly*2*lz + (iys2-1)*2*lz + iz)
                n3list(12,isite) = tag_site((ix-1)*2*ly*2*lz + (iys2-1)*2*lz + izs1)
                
                ixt1 = ix - 3
                ixt2 = ix + 3
                iyt1 = iy - 3
                iyt2 = iy + 3
                izt1 = iz - 3
                izt2 = iz + 3
                !fourth neighbor list boundary condition
                if(ixt1.lt.1) ixt1 = ixt1 + 2*lx
                if(ixt2.gt.2*lx) ixt2 = ixt2 - 2*lx
                if(iyt1.lt.1) iyt1 = iyt1 + 2*ly
                if(iyt2.gt.2*ly) iyt2 = iyt2 - 2*ly
                if(izt1.lt.1) izt1 = izt1 + 2*lz
                if(izt2.gt.2*lz) izt2 = izt2 - 2*lz
                
                !fourth neighbor list
                n4list(1,isite) = tag_site((ixt1-1)*2*ly*2*lz + (iy1-1)*2*lz + iz1)
                n4list(2,isite) = tag_site((ixt1-1)*2*ly*2*lz + (iy2-1)*2*lz + iz1)
                n4list(3,isite) = tag_site((ixt1-1)*2*ly*2*lz + (iy2-1)*2*lz + iz2)
                n4list(4,isite) = tag_site((ixt1-1)*2*ly*2*lz + (iy1-1)*2*lz + iz2)
                n4list(5,isite) = tag_site((ixt2-1)*2*ly*2*lz + (iy1-1)*2*lz + iz1)
                n4list(6,isite) = tag_site((ixt2-1)*2*ly*2*lz + (iy2-1)*2*lz + iz1)
                n4list(7,isite) = tag_site((ixt2-1)*2*ly*2*lz + (iy2-1)*2*lz + iz2)
                n4list(8,isite) = tag_site((ixt2-1)*2*ly*2*lz + (iy1-1)*2*lz + iz2)
                n4list(9,isite) = tag_site((ix1-1)*2*ly*2*lz + (iyt1-1)*2*lz + iz1)
                n4list(10,isite) = tag_site((ix1-1)*2*ly*2*lz + (iyt1-1)*2*lz + iz2)
                n4list(11,isite) = tag_site((ix2-1)*2*ly*2*lz + (iyt1-1)*2*lz + iz2)
                n4list(12,isite) = tag_site((ix2-1)*2*ly*2*lz + (iyt1-1)*2*lz + iz1)
                n4list(13,isite) = tag_site((ix1-1)*2*ly*2*lz + (iyt2-1)*2*lz + iz1)
                n4list(14,isite) = tag_site((ix1-1)*2*ly*2*lz + (iyt2-1)*2*lz + iz2)
                n4list(15,isite) = tag_site((ix2-1)*2*ly*2*lz + (iyt2-1)*2*lz + iz2)
                n4list(16,isite) = tag_site((ix2-1)*2*ly*2*lz + (iyt2-1)*2*lz + iz1)
                n4list(17,isite) = tag_site((ix1-1)*2*ly*2*lz + (iy1-1)*2*lz + izt1)
                n4list(18,isite) = tag_site((ix2-1)*2*ly*2*lz + (iy1-1)*2*lz + izt1)
                n4list(19,isite) = tag_site((ix2-1)*2*ly*2*lz + (iy2-1)*2*lz + izt1)
                n4list(20,isite) = tag_site((ix1-1)*2*ly*2*lz + (iy2-1)*2*lz + izt1)
                n4list(21,isite) = tag_site((ix1-1)*2*ly*2*lz + (iy1-1)*2*lz + izt2)
                n4list(22,isite) = tag_site((ix2-1)*2*ly*2*lz + (iy1-1)*2*lz + izt2)
                n4list(23,isite) = tag_site((ix2-1)*2*ly*2*lz + (iy2-1)*2*lz + izt2)
                n4list(24,isite) = tag_site((ix1-1)*2*ly*2*lz + (iy2-1)*2*lz + izt2)
                
                !fifth neighbors: within a distance of '3^0.5*a' from a site
                n5list(1,isite) = tag_site((ixs1-1)*2*ly*2*lz + (iys1-1)*2*lz + izs2)
                n5list(2,isite) = tag_site((ixs2-1)*2*ly*2*lz + (iys1-1)*2*lz + izs2)
                n5list(3,isite) = tag_site((ixs2-1)*2*ly*2*lz + (iys2-1)*2*lz + izs2)
                n5list(4,isite) = tag_site((ixs1-1)*2*ly*2*lz + (iys2-1)*2*lz + izs2)
                n5list(5,isite) = tag_site((ixs1-1)*2*ly*2*lz + (iys1-1)*2*lz + izs1)
                n5list(6,isite) = tag_site((ixs2-1)*2*ly*2*lz + (iys1-1)*2*lz + izs1)
                n5list(7,isite) = tag_site((ixs2-1)*2*ly*2*lz + (iys2-1)*2*lz + izs1)
                n5list(8,isite) = tag_site((ixs1-1)*2*ly*2*lz + (iys2-1)*2*lz + izs1)
            end if
          end do
       end do
    end do
    
    !Potential Neighbor List Errors
    if(isite.ne.nsites) then
       print*,"Error:nlist(site.neq.sites)::subroutine lattice",isite,nsites
       call freeandclose
       stop
    end if

    Si_O(1,1) = 1
    Si_O(2,1) = 3
    Si_O(3,1) = 6
    Si_O(4,1) = 8
    
    Si_O(1,2) = 2
    Si_O(2,2) = 4
    Si_O(3,2) = 5 
    Si_O(4,2) = 7
    
    Si_O(1,3) = 3
    Si_O(2,3) = 8
    Si_O(3,3) = 1
    Si_O(4,3) = 6

    Si_O(1,4) = 4
    Si_O(2,4) = 2
    Si_O(3,4) = 5
    Si_O(4,4) = 7
    
    Si_O(1,5) = 5
    Si_O(2,5) = 2
    Si_O(3,5) = 4
    Si_O(4,5) = 7
    
    Si_O(1,6) = 6
    Si_O(2,6) = 1
    Si_O(3,6) = 3
    Si_O(4,6) = 8
    
    Si_O(1,7) = 7
    Si_O(2,7) = 2
    Si_O(3,7) = 4
    Si_O(4,7) = 5
    
    Si_O(1,8) = 8
    Si_O(2,8) = 3
    Si_O(3,8) = 1
    Si_O(4,8) = 6
    
    Si_notO(1,1) = 2
    Si_notO(2,1) = 4
    Si_notO(3,1) = 5
    Si_notO(4,1) = 7
    
    Si_notO(1,2) = 3
    Si_notO(2,2) = 8
    Si_notO(3,2) = 1
    Si_notO(4,2) = 6
    
    Si_notO(1,3) = 4
    Si_notO(2,3) = 2
    Si_notO(3,3) = 5
    Si_notO(4,3) = 7
    
    Si_notO(1,4) = 1
    Si_notO(2,4) = 3
    Si_notO(3,4) = 6
    Si_notO(4,4) = 8
    
    Si_notO(1,5) = 6
    Si_notO(2,5) = 8
    Si_notO(3,5) = 1
    Si_notO(4,5) = 3
    
    Si_notO(1,6) = 7
    Si_notO(2,6) = 5
    Si_notO(3,6) = 2
    Si_notO(4,6) = 4
    
    Si_notO(1,7) = 8
    Si_notO(2,7) = 6
    Si_notO(3,7) = 1
    Si_notO(4,7) = 3
    
    Si_notO(1,8) = 2
    Si_notO(2,8) = 5
    Si_notO(3,8) = 4
    Si_notO(4,8) = 7

    return
  end subroutine lattice
  !---------------------------------------------------------
  !--- LATTICE_INITIALIZATION: SETUP THE INITIAL CONDITION
  !--------------------------------------------------------- 
  subroutine random_configuration
    use variables
    implicit none
    
    integer :: isite,ivertex,jsite,jocc,jlink,ihead
    logical :: can_insert            
    integer :: molecule,bond,occfn,i
    integer,dimension(4) :: fn,fnsite

    !set the occupancy and spin of sites to that od water: 
    !all sites are occupied by water
    !the sites of lattice initially have no pointer
    occupancy = occW
    spin = spinW
    head=0
    noccupy = 0
    resite = 0
    SIlink = -1;TAAlink = -1
    
    !-----FIRST: place all the SN on the sites
    molecule=1
    do while(molecule.le.ni(2))
       !select a random site among nsites
       isite = int(1 + ran3(seed)*real(nsites))
       !select a random vertex of the monomer
       ivertex = int(1 + ran3(seed)*real(2))
       
       !check for isite insertion and get the logical flag can_insert
       call can_insert_SN(isite,ivertex,can_insert)
       
       !if insertion of the monomer is possible
       if(can_insert) then

          !update the noccupy,ioxy,head and spin arrays
          noccupy(ni(1)+molecule) = isite
          head(isite) = ivertex
          spin(isite) = spinSN
          resite(isite) = ni(1) + molecule

          !change the occupancy of the site to that of SN
          occupancy(isite) = occSN
          !update the occupancy of the vertices of the monomer
          fn(1) = ivertex
          fnsite(1) = n1list(fn(1),isite)
          occupancy(fnsite(1)) = occupancy(fnsite(1)) + occSNO
          do bond = 2,4
             fn(bond) = Si_O(bond,fn(1))
             fnsite(bond) = n1list(fn(bond),isite)
             occupancy(fnsite(bond)) = occupancy(fnsite(bond)) + occSNO
          end do
          
          !proceed to inserting next molecule
          molecule = molecule+1          

       end if  !end checking can_insert if statement

    end do     !end loop over all the SN molecules

    if(molecule-1.ne.ni(2))then
       print*,'BUG: molecule-1.ne.ni(2) :: subroutine initial_configuration'
       call freeandclose
       stop
    end if
    !-----
    !-----SECOND:place all the SI molecules on the lattice
    molecule = 1
    do while(molecule.le.ni(1))
       ! select a random site
       isite = int(1 + ran3(seed)*real(nsites))
       !select a random orientation of the O-
       ivertex = int(1 + ran3(seed)*real(nc))

       !get the logical variable can_insert
       call can_insert_SI(isite,ivertex,can_insert)
       
       if(can_insert) then

          !change the occupancy of the site to that of SN
          occupancy(isite) = occSI
          !update the occupancy of the vertices of the monomer
          fn(1) = ivertex
          fnsite(1) = n1list(fn(1),isite)
          occupancy(fnsite(1)) = occISIO
          do bond = 2,4
             fn(bond) = Si_O(bond,fn(1))
             fnsite(bond) = n1list(fn(bond),isite)
             occupancy(fnsite(bond)) = occupancy(fnsite(bond)) + occSIO
          end do
          
          !update the noccupy,ioxy,head and spin  arrays
          noccupy(molecule) = isite
          head(isite) = ivertex
          spin(isite) = spinSI
          resite(isite) = molecule
          
          !link the molecule to any 4th neighbouting TAA
          ihead = n1list(ivertex,isite)
          do i = 1,6
             jsite = n2list(i,ihead);jlink = TAAlink(jsite);jocc = occupancy(jsite)
             if((jlink.lt.0).and.(jocc.eq.occTAA))then
                SIlink(ihead) = jsite
                TAAlink(jsite) = ihead
                exit
             end if
          end do

          !proceed to inserting next molecule
          molecule = molecule + 1 
          
       end if!can_insert=yes if statement ends
       
    end do!end loop over all molecules of SI
    
    if(molecule-1.ne.ni(1))then
       print*,'BUG: molecule-1.ne.ni(1) :: subroutine initial_configuration'
       call freeandclose
       stop
    end if
    !-----
    !-----THIRD: PLACE ALL THE TAA MOLECULES ON THE LATTICE
    molecule = 1
    do while(molecule.le.ni(3))
       !select a random site
       isite = int(1+ran3(seed)*real(nsites))
       
       !check if the TAA molecule can be inserted in isite
       call can_insert_TAA(isite,can_insert)
       
       if(can_insert)then!if the molecule can be inserted
          
          !change the occupancy of isite
          occupancy(isite) = occTAA
          !change the spin of isite
          spin(isite) = spinTAA
          !update the noccupy
          noccupy(ni(1)+ni(2)+molecule) = isite
          !update the resite arrays
          resite(isite) = molecule+ni(1)+ni(2)  
          
          !link the molecule to any 4th neighbouring ISIO atom
          do i = 1,6
             jsite = n2list(i,isite);jlink = SIlink(jsite);jocc = occupancy(jsite)
             if((jlink.lt.0).and.(jocc.eq.occISIO))then
                TAAlink(isite) = jsite
                SIlink(jsite) = isite
                exit
             end if
          end do
          
          molecule = molecule + 1
       end if!endif: the molecule can be inserted
    end do
    
    if(molecule-1.ne.ni(3))then
       print*,'BUG: molecule-1.ne.ni(3) :: subroutine initial_configuration'
       call freeandclose
       stop
    end if
    
    !output the psf file
    call psf_out
    
    !output initial configuration
    call config_out_atom
    
  end subroutine random_configuration
  !---------------------------------------------------------
  !--- INITIAL_CRITOBALITE CONFIGURATION
  !--- only apply this subroutine with SI/TAA = 0
  !--------------------------------------------------------- 
  subroutine initial_cristo_config
    use variables
    implicit none
    
    integer :: i,j,k,bond,opsttus,k1,k2
    integer,dimension(8) :: x,y,z,fn
    integer :: isite,m,ivertex,molecule
    integer,dimension(4) :: fn1,fnsite
    logical :: can_insert
    
    !setup the initial lattice state
    spin = spinW      !all water on site
    occupancy = occW  !all sites occupied by water
    head = 0          !no pointer for any site
    k1 = 1            !start with molecule 1
        
    !Only place neutral silica there
    do i = 1, lx / 4 
       !do i = 1, lx / 4 
       do j = 1, ly / 4
          do k = 1, lz / 4
         
             !unit cell
             !Point 1 (7, 7, 7)
             x(1) = i * 8 - 1
             y(1) = j * 8 - 1
             z(1) = k * 8 - 1
             fn(1) = 1
             
             !Point 2 (7, 3, 3)
             x(2) = i * 8 - 1
             y(2) = j * 8 - 5
             z(2) = k * 8 - 5
             fn(2) = 1
             
             !Point 3 (3, 3, 7)
             x(3) = i * 8 - 5
             y(3) = j * 8 - 5
             z(3) = k * 8 - 1
             fn(3) = 1
             
             !Point 4 (3, 7, 3)
             x(4) = i * 8 - 5
             y(4) = j * 8 - 1
             z(4) = k * 8 - 5
             fn(4) = 1
             
             !Point 5 (1, 1, 1)
             x(5) = i * 8 - 7
             y(5) = j * 8 - 7
             z(5) = k * 8 - 7
             fn(5) = 2
             
             !Point 6 (1, 5, 5)
             x(6) = i * 8 - 7
             y(6) = j * 8 - 3
             z(6) = k * 8 - 3
             fn(6) = 2
             
             !Point 7 (5, 1, 5)
             x(7) = i * 8 - 3
             y(7) = j * 8 - 7
             z(7) = k * 8 - 3
             fn(7) = 2
             
             !!	Point 8 (5, 5, 1)
             x(8) = i * 8 - 3
             y(8) = j * 8 - 3
             z(8) = k * 8 - 7
             fn(8) = 2
             
             do m = 1, 8
                
                isite = tag_site((x(m)-1)*(2*ly)*(2*lz) + (y(m)-1) * (2*lz)+z(m))
                noccupy(k1) = isite
                resite(isite) = k1
                
                call can_insert_SN(isite, fn(m), can_insert)
                if(can_insert)then
                   call insert_site(isite,fn(m),spinSN,k1)
                else
                   write(*,*) 'Initialize beta-cristobalite fail'
                   call freeandclose
                   stop 
                endif
                
                !Increase occupied site counters
            	
                k1 = k1 + 1
                
                if(k1.eq.(ni(2)+1)) then
                   !write the psf file
                   call psf_out
                   !output initial configuration
                   call config_out_atom
                   return
                end if

             enddo
             !Coordinates for silica, on a bigger BCC lattice	
             
          enddo	!!	2*lz
          
       enddo	!!	2*ly
       
    enddo	!	2*lx
    
    return
  end subroutine initial_cristo_config
  !---------------------------------------------------------
  !--- PSF_OUT: GENERATES THE PSF FILE FOR VMD: 
  !--------------------------------------------------------- 
  subroutine psf_out
    use variables
    implicit none
    
    integer :: bond,i,j,k1,k2,sumoccupy,opstatus
    
    open(1005,file='config.dat.psf',status='unknown',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'Error: opening .psf file :: subroutine initial_configuration'
       call freeandclose
       stop
    end if
    
    !bond = ni(1)*4 + ni(2)*4 + ni(3)*4
    !sumoccupy = ni(1)*5 + ni(2)*5 + ni(3)*5
    bond = ni(1)*4 + ni(2)*4
    sumoccupy = ni(1)*5 + ni(2)*5 + ni(3)
    
    write(1005,*) 'PSF   '
    write(1005,*) '      '
    write(1005,700) 1, '!NTITLE'
    write(1005,*) 'REMARKS original generated silica nanoparticle psf file'
    write(1005,*) '       '
    write(1005,700) sumoccupy,'!NATOM'
    
    !ionic silica
    j = 0
    do i=1,ni(1)
       j = j + 1
       write(1005,200) j
       j = j + 1
       write(1005,800) j
       do k1 = 1,3
          j = j + 1
          write(1005,400) j
       end do
    end do
    
    !neutral silica
    do i=1,ni(2)
       j = j + 1
       write(1005,100) j
       do k1 = 1,4
          j = j + 1
          write(1005,400) j
       end do
    end do
    
    !TAA
    do i=1,ni(3)
       j = j + 1
       write(1005,300) j
       !do k1=1,4
       !   j = j + 1
       !   write(1005,300)j
       !end do
    end do
    
    !if(j.ne.natom)then
    if(j.ne.(5*ni(1)+5*ni(2)+ni(3)))then
       print*,'BUG:j.ne.natom::subroutine psf_out'
       
       call freeandclose
       stop
    end if
    
    write(1005,*) '       '
    write(1005,700) bond, '!NBOND: bonds' 
    
    !for ionic silica
    k2 = 0
    do i=1,ni(1)
       write(1005,600) k2+1,k2+2,k2+1,k2+3,k2+1,k2+4,k2+1,k2+5
       k2 = k2 + 5
    end do
    
    !for neutral silica
    do i=1,ni(2)
       write(1005,600) k2+1,k2+2,k2+1,k2+3,k2+1,k2+4,k2+1,k2+5
       k2 = k2 + 5
    end do
    
    !for TAA
    do i=1,ni(3)
       write(1005,600) k2+1
       k2 = k2 + 1
       !write(1005,600) k2+1,k2+2,k2+1,k2+3,k2+1,k2+4,k2+1,k2+5
       !k2 = k2 + 5
    end do
    
    !if(k2.ne.natom)then
    if(k2.ne.(5*ni(1)+5*ni(2)+ni(3)))then
       print*,'BUG:k2.ne.natom::subroutine psf_out'
       print*,'k2,natom',k2,natom
       call freeandclose
       stop
    end if
    
    write(1005,*) '       '
    write(1005,*) '       '
    write(1005,700) 0, '!NTHETA: angles'
    write(1005,*) '       '
    write(1005,700) 0, '!NPHI: dihedrals'
    write(1005,*) '       '
    write(1005,700) 0, '!NIMPHI: impropers'
    write(1005,*) '       '
    write(1005,700) 0, '!NCRTERM: cross-terms'
    
    close(1005)
    
100 format(I8,x,'CHAINS 1',x,'CHS',x,'CHS',x,'CHS',x,'0.000000 1.000000     0')
200 format(I8,x,'CHAINO 1',x,'CHO',x,'CHO',x,'CHO',x,'0.000000 1.000000     0')    
300 format(I8,x,'CHAINN 1',x,'CHN',x,'CHN',x,'CHN',x,'0.000000 1.000000     0')
400 format(I8,x,'CHAINC 1',x,'CHC',x,'CHC',x,'CHC',x,'0.000000 1.000000     0')
500 format(I8,x,'CHAINP 1',x,'CHP',x,'CHP',x,'CHP',x,'0.000000 1.000000     0')
800 format(I8,x,'CHAINH 1',x,'CHH',x,'CHH',x,'CHH',x,'0.000000 1.000000     0')
600 format(I8, I8, I8, I8, I8, I8, I8, I8)
700 format(I8,x,A) 
    
    return
  end subroutine psf_out
  !---------------------------------------------------------
  !--- CAN_INSERT_SN: CHECKS WHETHER SN CAN BE INSETED AT A SITE
  !---------------------------------------------------------   
  subroutine can_insert_SN(isite,ivertex,can_insert)
    use variables
    implicit none
    
    integer,intent(in) :: isite, ivertex
    logical,intent(out) :: can_insert
    integer :: occfn,neighs,overlap
    integer :: jsite,sn_neigh
    integer :: bridge,i,j,bond,snbond         
    integer,dimension(4) :: fn,fnsite,sn,snsite
    
    !assume insertion is possible
    can_insert = .false.
    overlap = 0
    
    if(occupancy(isite).eq.occW)then!if isite has water on it
       
       !There's no need to uncomment the following lines,
       !since, these conditions would be taken care by the constraints
       !imposed on OH and their 1st,2nd and 3rd neighbors
       
       !check first,second and third neighbors for TAA occupancy
       do i=1,8!checking first neighbors
          occfn = occupancy(n1list(i,isite))
          if(occfn.eq.occTAA)then
             can_insert = .false.
             return
          end if
       end do!endo: checking first neighbors
       !do i = 1,6!checking first neighbors
       !   occfn = occupancy(n2list(i,isite))
       !   if(occfn.eq.occTAA)then
       !      can_insert = .false.
       !      return
       !   end if
       !end do!enddo: checking first neighbors
       !do i = 1,12!checking third neighbors
       !   occfn = occupancy(n3list(i,isite))
       !   if(occfn.eq.occTAA)then
       !      can_insert = .false.
       !      return
       !   end if
       !end do!enddo: checking third neighbors
       
       fn(1) = ivertex
       fnsite(1) = n1list(fn(1),isite)
       if((occupancy(fnsite(1)).eq.occSNO).or.((occupancy(fnsite(1)).eq.occSIO)))&
            overlap = overlap + 1
       do bond = 2,4
          fn(bond) = Si_O(bond,fn(1))
          fnsite(bond) = n1list(fn(bond),isite)
          if((occupancy(fnsite(bond)).eq.occSNO).or.((occupancy(fnsite(bond)).eq.occSIO)))&
               overlap = overlap + 1
       end do
       
       !check if there are any water,SN,SNO or SIO molecules on the fnsite() sites
       do bond = 1,4
          occfn = occupancy(fnsite(bond))
          if((occfn.eq.occW).or.(occfn.eq.occSNO).or.(occfn.eq.occSIO))then
             can_insert = .true.
          else
             can_insert = .false.
             return
          end if
       end do
       
       do bond = 1,4!check first,second and third neighbors of OH for TAA occupancy
          
          do i=1,8!checking first neighbors
             occfn = occupancy(n1list(i,fnsite(bond)))
             if(occfn.eq.occTAA)then
                can_insert = .false.
                return
             end if
          end do!endo: checking first neighbors
          !do i = 1,6!checking first neighbors
          !   occfn = occupancy(n2list(i,fnsite(bond)))
          !   if(occfn.eq.occTAA)then
          !      can_insert = .false.
          !      return
          !   end if
          !end do!enddo: checking first neighbors
          !do i = 1,12!checking third neighbors
          !   occfn = occupancy(n3list(i,fnsite(bond)))
          !   if(occfn.eq.occTAA)then
          !      can_insert = .false.
          !      return
          !   end if
          !end do!enddo: checking third neighbors
          
       end do!enddo: check first,second and third neighbors of OH for TAA occupancy
             
       !-----CHECK FOR 2 MEMBERED RING FORMATION
       if(overlap.ge.2)then!if overlap>=2 chance of 2 membered ring

          do sn_neigh = 1,6!loop over 2 neighbors of isite
             
             jsite = n2list(sn_neigh,isite)
             !if there SN on neighboring site
             if((occupancy(jsite).eq.occSN).or.(occupancy(jsite).eq.occSI))then
                
                !find the bonded OH of jsite
                sn(1) = head(jsite)
                snsite(1) = n1list(sn(1),jsite)
                do snbond = 2,4
                   sn(snbond) = Si_O(snbond,sn(1))
                   snsite(snbond) = n1list(sn(snbond),jsite)
                end do
                
                !check for any bridging oxygens
                bridge=0
                do i=1,4
                   do j=1,4
                      if((fnsite(i).eq.snsite(j)).and.((occupancy(fnsite(i)).eq.occSNO).or.&
                           (occupancy(fnsite(i)).eq.occSIO)))then 
                         bridge = bridge + 1
                      end if
                   end do
                end do

                !if there are more than 1 bridging oxygens
                if(bridge.gt.1) then
                   can_insert = .false.
                   return
                end if
                
             end if!if there SN on neighboring site

          end do!enddo:loop over 2 neighbors of isite

       end if!if overlap.ge.2 chance of 2 membered ring
       !-----

    end if!if isite has water on it

    return
  end subroutine can_insert_SN
  !----------------------------------------------------------
  !--- CAN_INSERT_SI: CHECKs WHETHER SI CAN BE INSERTED AT A SITE
  !---------------------------------------------------------- 
  subroutine can_insert_SI(isite,ivertex,can_insert)
    
    integer,intent(in) :: isite, ivertex
    logical,intent(out) :: can_insert
    integer :: jsite,jocc,neighs,overlap,ihead
    integer :: bridge,i,j,bond,snbond,sn_neigh,occfn,TAA_num
    integer,dimension(4) :: fn,fnsite,sn,snsite
    
    !assume insertion is possible
    can_insert = .false.
    overlap = 0
    
    !check if there's water on the site
    if(occupancy(isite).eq.occW)then!if isite has water on it

       !There's no need to uncomment the following lines,
       !since, these conditions would be taken care by the constraints
       !imposed on OH and their 1st,2nd and 3rd neighbors
       
       !check first,second and third neighbors for TAA occupancy
       do i=1,8!checking first neighbors
          occfn = occupancy(n1list(i,isite))
          if(occfn.eq.occTAA)then
             can_insert = .false.
             return
          end if
       end do!endo: checking first neighbors
       !do i = 1,6!checking first neighbors
       !   occfn = occupancy(n2list(i,isite))
       !   if(occfn.eq.occTAA)then
       !      can_insert = .false.
       !      return
       !   end if
       !end do!enddo: checking first neighbors
       !do i = 1,12!checking third neighbors
       !   occfn = occupancy(n3list(i,isite))
       !   if(occfn.eq.occTAA)then
       !      can_insert = .false.
       !      return
       !   end if
       !end do!enddo: checking third neighbors
       
       fn(1) = ivertex
       fnsite(1) = n1list(fn(1),isite)
       ihead = fnsite(1)
       if(occupancy(fnsite(1)).eq.occW)then
       else
          can_insert = .false.
          return
       end if
       do bond = 2,4
          fn(bond) = Si_O(bond,fn(1))
          fnsite(bond) = n1list(fn(bond),isite)
          if(occupancy(fnsite(bond)).eq.occSNO)&
               overlap = overlap + 1
       end do
       
       do bond = 1,4!check first,second and third neighbors of OH for TAA occupancy
          
          do i=1,8!checking first neighbors
             occfn = occupancy(n1list(i,fnsite(bond)))
             if(occfn.eq.occTAA)then
                can_insert = .false.
                return
             end if
          end do!endo: checking first neighbors
          !do i = 1,6!checking first neighbors
          !   occfn = occupancy(n2list(i,fnsite(bond)))
          !   if(occfn.eq.occTAA)then
          !      can_insert = .false.
          !      return
          !   end if
          !end do!enddo: checking first neighbors
          !do i = 1,12!checking third neighbors
          !   occfn = occupancy(n3list(i,fnsite(bond)))
          !   if(occfn.eq.occTAA)then
          !      can_insert = .false.
          !      return
          !   end if
          !end do!enddo: checking third neighbors
          
       end do!enddo: check first,second and third neighbors of OH for TAA occupancy
       
       !check if there are any water,SN,SI,SNO or SIO molecules on the 
       !fnsite() sites
       !check for the location of the pointer variable

       if(occupancy(fnsite(1)).eq.occW)then!if the pointer has water on its site
          do bond = 2,4
             jocc = occupancy(fnsite(bond))
             if((jocc.eq.occW).or.(jocc.eq.occSNO))then
                can_insert = .true.
             else
                can_insert = .false.
                return
             end if
          end do
          !-----CHECK FOR 2 MEMBERED RING FORMATION
          if(overlap.ge.2)then!if overlap.ge.2 chance of 2 membered ring
             
             do sn_neigh = 1,6!loop over 2 neighbors of isite
                
                jsite = n2list(sn_neigh,isite)
                !if there SN on neighboring site
                if(occupancy(jsite).eq.occSN.or.(occupancy(jsite).eq.occSI))then
                   
                   !find the bonded OH of jsite
                   sn(1) = head(jsite)
                   snsite(1) = n1list(sn(1),jsite)
                   do snbond = 2,4
                      sn(snbond) = Si_O(snbond,sn(1))
                      snsite(snbond) = n1list(sn(snbond),jsite)
                   end do
                   
                   !check for any bridging oxygens
                   bridge=0
                   do i=1,4
                      do j=1,4
                         if((fnsite(i).eq.snsite(j)).and.(occupancy(fnsite(i)).eq.occSNO))then 
                            bridge = bridge + 1
                         end if
                      end do
                   end do
                   
                   !if there are more than 1 bridging oxygens
                   if(bridge.gt.1) then
                      can_insert = .false.
                      return
                   end if
                   
                end if!if there SN on neighboring site
                
             end do!enddo:loop over 2 neighbors of isite
             
          end if!if overlap.ge.2 chance of 2 membered ring
          !-----
       else
          can_insert = .false.
          return
       end if!endif: the pointer has water on its site
       
    end if!endif: isite has water on it
    
    return
  end subroutine can_insert_SI
  !---------------------------------------------------------
  !--- CAN_INSERT_TAA(isite,can_insert): FIND THE Qn DISTRIBUTION
  !--------------------------------------------------------- 
  subroutine can_insert_TAA(isite,can_insert)
    use variables
    implicit none
    
    integer,intent(in) :: isite
    logical,intent(out) :: can_insert
    integer :: i,jocc,occfn
    
    can_insert = .false.
    
    !check the occupancy of isite
    if(occupancy(isite).eq.occW)then!if there's water on isite
       
       !check first neighbors any occupancy other than water
       do i=1,8
          jocc = occupancy(n1list(i,isite))
          if(jocc.eq.occW)then
             can_insert = .true.
          else
             can_insert = .false.
             return
          end if
       end do
       
       !check first,second and third neighbors for TAA occupancy
       !do i=1,8!checking first neighbors
       !   occfn = occupancy(n1list(i,isite))
       !   if(occfn.ne.occW)then
       !      can_insert = .false.
       !      return
       !   else
       !      can_insert = .true.
       !   end if
       !end do!endo: checking first neighbors
       !do i = 1,6!checking first neighbors
       !   occfn = occupancy(n2list(i,isite))
       !   if(occfn.ne.occW)then
       !      can_insert = .false.
       !      return
       !   end if
       !end do!enddo: checking first neighbors
       !do i = 1,12!checking third neighbors
       !   occfn = occupancy(n3list(i,isite))
       !   if(occfn.ne.occW)then
       !      can_insert = .false.
       !      return
       !   else
       !      can_insert = .true.
       !   end if
       !end do!enddo: checking third neighbors

    else

       can_insert = .false.
       return

    end if!endif there's water on isite       
    
    return
  end subroutine can_insert_TAA
  !---------------------------------------------------------
  !--- Qn_DIST(Q0,Q1,Q2,Q3,Q4): FIND THE Qn DISTRIBUTION
  !--------------------------------------------------------- 
  subroutine Qn_dist
    use variables
    implicit none
    
    integer :: isite,ivertex,ispin,bond,occfn,centerocc
    integer :: overlap,monomer
    integer,dimension(4) :: fn,fnsite
    
    !set the accumulators to 0
    Qn = 0;Qni = 0;Qnn = 0
    
    do monomer=1,allmolecules!loop over all monomers
       
       isite = noccupy(monomer)
       ivertex = head(isite)
       ispin = spin(isite)
       centerocc = occupancy(isite)
       
       if((ispin.eq.spinSN).or.(ispin.eq.spinSI))then!if the molecule is SN/SI
          fn(1) = ivertex
          fn(2) = Si_O(2,ivertex)
          fn(3) = Si_O(3,ivertex)
          fn(4) = Si_O(4,ivertex)
          
          overlap = 0
          do bond=1,4
             fnsite(bond) = n1list(fn(bond),isite)
             occfn = occupancy(fnsite(bond))
             if((occfn.eq.occSNOSNO).or.(occfn.eq.occSIOSNO))then
                overlap = overlap + 1
             end if
          end do !end loop over hydroxide molecules
          
          !update the Qn variables
          if(overlap.eq.0) then
             Qn(0) = Qn(0) + 1
          elseif(overlap.eq.1) then
             Qn(1) = Qn(1) + 1
          elseif(overlap.eq.2) then
             Qn(2) = Qn(2) + 1
          elseif(overlap.eq.3) then
             Qn(3) = Qn(3) + 1
          elseif(overlap.eq.4) then
             Qn(4) = Qn(4) + 1
          else
             print*,'BUG: More than 4 bridging O to Si : subroutine Qn_distribution (Qn)'
             call freeandclose
             stop
          end if

       end if!endif the molecule is SN/SI
       
       !Qn distribution only for neutral silica
       if(ispin.eq.spinSN)then
          fn(1) = ivertex
          fn(2) = Si_O(2,ivertex)
          fn(3) = Si_O(3,ivertex)
          fn(4) = Si_O(4,ivertex)
          
          overlap = 0
          do bond=1,4
             fnsite(bond) = n1list(fn(bond),isite)
             occfn = occupancy(fnsite(bond))
             if((occfn.eq.occSNOSNO).or.(occfn.eq.occSIOSNO))then
                overlap = overlap + 1
             end if
          end do !end loop over hydroxide molecules
          
          !update the Qn variables
          if(overlap.eq.0) then
             Qnn(0) = Qnn(0) + 1
          elseif(overlap.eq.1) then
             Qnn(1) = Qnn(1) + 1
          elseif(overlap.eq.2) then
             Qnn(2) = Qnn(2) + 1
          elseif(overlap.eq.3) then
             Qnn(3) = Qnn(3) + 1
          elseif(overlap.eq.4) then
             Qnn(4) = Qnn(4) + 1
          else
             print*,'BUG: More than 4 bridging O to Si : subroutine Qn_distribution (Qnn)'
             call freeandclose
             stop
          end if
       end if

    end do!end loop over all the monomers
    
    return
  end subroutine Qn_dist
  !---------------------------------------------------------
  !---returnneigh: returns neighboring Si to a site
  !---------------------------------------------------------
  subroutine returnneigh(isite,isites,arethereSi)
    use variables
    implicit none
    
    integer,intent(in) :: isite
    !integer,dimension(:),allocatable,intent(out) :: neighSi
    logical,intent(out) :: arethereSi
    integer :: ispin,ivertex,bond,fnbond,OH,Obond,fnn,i,occfn
    integer :: neighboringSi,index
    integer,dimension(4) :: fn, fnsite, sn,snsite
    integer,dimension(4) :: isites
    logical :: connect
    
    ispin = spin(isite)
    ivertex = head(isite)
    
    !if(allocated(neighSi)) deallocate(neighSi)
    
    !set the accumulator to zero
    neighboringSi = 0
    isites=-1
    arethereSi = .false.
    
    if((ispin.eq.spinSI).or.(ispin.eq.spinSN))then
       !find the neighbors of isite
       fn(1) = ivertex
       fnsite(1) = n1list(fn(1),isite)
       do bond = 2,4
          fn(bond) = Si_O(bond,ivertex)
          fnsite(bond) = n1list(fn(bond),isite)
       end do
    else
       print*,'BUG : molecule other than SI/SN :: subroutine returnneigh'
       print*,'identity:',ispin
       call freeandclose
       stop
    end if
    
    do bond = 1,4!find the db occ. O on nearest sites to isite

       OH = fnsite(bond)
       occfn = occupancy(OH)

       if((occfn.eq.occSNOSNO).or.(occfn.eq.occSIOSNO))then!if there are double occupies Os
          
          do Obond = 1,8!loop over the first neighbors of OH

             fnn = n1list(Obond,OH)!fnn:nearest neighboring Si to isite

             if(((occupancy(fnn).eq.occSN).or.(occupancy(fnn).eq.occSI)).and.&
                  (fnn.ne.isite))then!if there's Si neighboring OH

                !find the neighbors of fnn
                sn(1) = head(fnn)
                snsite(1) = n1list(sn(1),fnn)
                do fnbond = 2,4
                   sn(fnbond) = Si_O(fnbond,sn(1))
                   snsite(fnbond) = n1list(sn(fnbond),fnn)
                end do

                !check for a connection
                connect = .false.
                do fnbond = 1,4
                   if(snsite(fnbond).eq.OH)then
                      connect = .true.
                   end if
                end do

                if(connect)then
                   neighboringSi = neighboringSi + 1
                   isites(neighboringSi) = fnn
                end if

             end if!endif: there's SN/SI neighboring OH
             
          end do!loop over the first neighbors of OH

       end if!endif there are double occupies Os

    end do!enddo:find the db occ. O on nearest sites to isite
    
    index = 0
    if(neighboringSi.ge.1)then!isite is a part of a ring system
       arethereSi = .true.
       !allocate(neighSi(neighboringSi))
       !do i = 1,size(isites)
       !   if(isites(i).ne.0)then
       !      index = index + 1
       !      neighSi(index) = isites(i)
       !   end if
       !end do
    end if!endif:isite is a part of a ring system
    
    return
  end subroutine returnneigh
  !---------------------------------------------------------
  !---rings3and4: 3 and 4 membered ring calculation
  !--------------------------------------------------------- 
  subroutine rings3and4(rings3,rings4)
    use variables
    implicit none

    integer,intent(out) :: rings3,rings4    
    integer :: rings3_temp,rings4_temp,molecule
    integer :: isite,ivertex,iocc,ispin,bond1,bond2
    integer :: jsite,jvertex,jocc
    integer :: ksite,kvertex,kocc
    integer :: lsite,lvertex,locc
    integer :: OH11,OH12,occOH11,occOH12,OH21,occOH21
    integer :: bondOH11,bondOH12,bondj,bondk,bondl,p,q,bondOH21
    integer,dimension(4) :: ifn,jfn,kfn,lfn,ifnsite,jfnsite,kfnsite,lfnsite
    logical :: connectedij,connectedik,connectedjl,connectedkl
    
    !set accumulators to zero
    rings3_temp = 0
    rings4_temp = 0
    
    do molecule = 1,allmolecules!loop over all molecules
       
       isite = noccupy(molecule)
       ivertex = head(isite)
       ispin = spin(isite)
       
       if((ispin.eq.spinSN).or.(ispin.eq.spinSI))then!if the molecule is SN/SI
          
          do bond1 = 1,4
             ifn(bond1) = Si_O(bond1,ivertex)
             ifnsite(bond1) = n1list(ifn(bond1),isite)              
          end do
          
          do bond1 = 1,4!loop over the OH of isite
             OH11 = ifnsite(bond1);occOH11 = occupancy(OH11)
             do bond2 = 1,4!loop over the OH of isite
                OH12 = ifnsite(bond2);occOH12 = occupancy(OH12)
                
                !if there are more then 2 brigdes
                if(((occOH11.eq.occSNOSNO).or.(occOH11.eq.occSIOSNO)).and.&
                     ((occOH12.eq.occSNOSNO).or.(occOH12.eq.occSIOSNO)).and.(OH11.ne.OH12))then
                   
                   !find if there are any SN/SI connected to OH11
                   do bondOH11 = 1,8!loop over first neighbors of OH11
                      jsite = n1list(bondOH11,OH11);jocc = occupancy(jsite);jvertex = head(jsite)
                      
                      if((jvertex.gt.0).and.(jsite.ne.isite))then!if jsite has SN/SI
                         
                         !find the first neighbors of jsite
                         do bondj = 1,4
                            jfn(bondj) = Si_O(bondj,jvertex)
                            jfnsite(bondj) = n1list(jfn(bondj),jsite)
                         end do
                         
                         !check if there's a bridge between jsite and isite
                         connectedij = .false.
                         do bondj = 1,4
                            if(OH11.eq.jfnsite(bondj)) connectedij = .true.
                         end do
                         
                         if(connectedij)then!if there's a bridge between isite and jsite
                            
                            !find the neighboring SI/SN to OH12
                            do bondOH12 = 1,8!loop over first neighbors of OH12
                               ksite = n1list(bondOH12,OH12);kocc = occupancy(ksite);kvertex = head(ksite)
                               
                               if((kvertex.gt.0).and.(ksite.ne.isite))then!if ksite has SN/SI
                                  
                                  !find the first neighbors of ksite
                                  do bondk = 1,4
                                     kfn(bondk) = Si_O(bondk,kvertex)
                                     kfnsite(bondk) = n1list(kfn(bondk),ksite)
                                  end do
                                  
                                  !check if there's a brigde between ksite and isite
                                  connectedik = .false.
                                  do bondk = 1,4
                                     if(kfnsite(bondk).eq.OH12)connectedik = .true.
                                  end do
                                  
                                  if(connectedik)then!if there's a bridge between isite and ksite
                                     
                                     !---CALCULATE THE 3 MEMBERED RINGS---
                                     !older version
                                     !check if there are any common OH to ksite and jsite
                                     !do p=1,4
                                     !   do q=1,4
                                     !      if(kfnsite(q).eq.jfnsite(p))rings3_temp = rings3_temp + 1
                                     !   end do
                                     !end do
                                     !------------------------------------
                                     
                                     !---CALCULATE THE 4 MEMBERED RINGS
                                     do bondj = 1,4!loop over first neighbors of jsite
                                        OH21 = jfnsite(bondj);occOH21 = occupancy(OH21)
                                        
                                        if((OH21.ne.OH11).and.((occOH21.eq.occSNOSNO).or.&
                                             (occOH21.eq.occSIOSNO)))then!if there's another bridge of jsite
                                           
                                           do bondOH21 = 1,8!loop over first neighbors of OH21
                                              lsite = n1list(bondOH21,OH21);locc = occupancy(lsite);lvertex = head(lsite)
                                              
                                              if((lvertex.gt.0).and.(lsite.ne.jsite).and.(lsite.ne.ksite))then!if there's SN/SI on lsite
                                                 
                                                 !find the OH of lsite
                                                 do bondl = 1,4
                                                    lfn(bondl) = Si_O(bondl,lvertex)
                                                    lfnsite(bondl) = n1list(lfn(bondl),lsite)
                                                 end do
                                                 
                                                 !loop through the first neighbors of OH21 to find the bridge between lsite and jsite
                                                 connectedjl = .false.
                                                 do bondl = 1,4
                                                    if(lfnsite(bondl).eq.OH21)connectedjl = .true.
                                                 end do
                                                 
                                                 if(connectedjl)then
                                                    !check if there's any common OH between lsite and ksite
                                                    do p = 1,4
                                                       do q=1,4
                                                          if((kfnsite(p).eq.lfnsite(q)).and.(kfnsite(p).ne.OH12))&
                                                               rings4_temp = rings4_temp + 1
                                                       end do
                                                    end do
                                                 end if
                                                 
                                              elseif((lvertex.gt.0).and.(lsite.ne.jsite).and.(lsite.eq.ksite))then!if there's SN/SI on lsite
                                                 
                                                 !---CALCULATE 3 MEMBERED RINGS---
                                                 !newer version
                                                 !loop over OH of lsite and check if it is bonded to OH12
                                                 do bondk = 1,4
                                                    if(kfnsite(bondk).eq.OH12)rings3_temp = rings3_temp + 1
                                                 end do
                                                 !--------------------------------
                                                 
                                              end if!endif: there's SN/SI on lsite
                                              
                                           end do!enddo: loop over first neighbors of OH21
                                           
                                        end if!endif: there's another bridge of jsite site
                                        
                                     end do!enddo: loop over first neighbors of jsite
                                     !---------------------------------
                                     
                                  end if!endif: there's a bridge between isite and ksite
                                  
                               end if!endif: ksite has SN/SI
                               
                            end do!enddo: loop over first neighbors of OH12
                            
                         end if!endif: there's a bridge between isite and jsite
                         
                      end if!endif jsite has SN/SI
                      
                   end do!enddo: loop over first neighbors of OH11
                   
                end if!endif: there are more then 2 brigdes
                
             end do!enddo: loop over the OH of isite
          end do!enndo: loop over the OH of isite
          
       end if!endif: the molecule is SN/SI
       
    end do!enddo: loop over all molecules
    
    rings3 = rings3_temp/6
    rings4 = rings4_temp/8
    
    return
  end subroutine rings3and4
  !---------------------------------------------------------
  !---site_rings3and4
  !--------------------------------------------------------- 
  subroutine site_rings3and4(isite,ivertex,rings3,rings4)
    use variables
    implicit none
    
    integer,intent(in) :: isite,ivertex
    integer,intent(out) :: rings3,rings4
    integer :: rings3_temp,rings4_temp
    integer :: bond1,bond2,iocc
    integer :: jsite,jvertex,jocc
    integer :: ksite,kvertex,kocc
    integer :: lsite,lvertex,locc
    integer :: OH11,OH12,occOH11,occOH12,OH21,occOH21
    integer :: bondOH11,bondOH12,bondj,bondk,bondl,p,q,bondOH21
    integer,dimension(4) :: ifn,jfn,kfn,lfn,ifnsite,jfnsite,kfnsite,lfnsite
    logical :: connectedij,connectedik,connectedjl,connectedkl
    
    !set the accumulators to zero
    rings3_temp = 0
    rings4_temp = 0
    
    do bond1 = 1,4
       ifn(bond1) = Si_O(bond1,ivertex)
       ifnsite(bond1) = n1list(ifn(bond1),isite)              
    end do
    
    do bond1 = 1,4!loop over the OH of isite
       OH11 = ifnsite(bond1);occOH11 = occupancy(OH11)
       do bond2 = 1,4!loop over the OH of isite
          OH12 = ifnsite(bond2);occOH12 = occupancy(OH12)
          
          !if there are more then 2 brigdes
          if(((occOH11.eq.occSNOSNO).or.(occOH11.eq.occSIOSNO)).and.&
               ((occOH12.eq.occSNOSNO).or.(occOH12.eq.occSIOSNO)).and.(OH11.ne.OH12))then
             
             !find if there are any SN/SI connected to OH11
             do bondOH11 = 1,8!loop over first neighbors of OH11
                jsite = n1list(bondOH11,OH11);jocc = occupancy(jsite);jvertex = head(jsite)
                
                if((jvertex.gt.0).and.(jsite.ne.isite))then!if jsite has SN/SI
                   
                   !find the first neighbors of jsite
                   do bondj = 1,4
                      jfn(bondj) = Si_O(bondj,jvertex)
                      jfnsite(bondj) = n1list(jfn(bondj),jsite)
                   end do
                   
                   !check if there's a bridge between jsite and isite
                   connectedij = .false.
                   do bondj = 1,4
                      if(OH11.eq.jfnsite(bondj)) connectedij = .true.
                   end do
                   
                   if(connectedij)then!if there's a bridge between isite and jsite
                      
                      !find the neighboring SI/SN to OH12
                      do bondOH12 = 1,8!loop over first neighbors of OH12
                         ksite = n1list(bondOH12,OH12);kocc = occupancy(ksite);kvertex = head(ksite)
                         
                         if((kvertex.gt.0).and.(ksite.ne.isite))then!if ksite has SN/SI
                            
                            !find the first neighbors of ksite
                            do bondk = 1,4
                               kfn(bondk) = Si_O(bondk,kvertex)
                               kfnsite(bondk) = n1list(kfn(bondk),ksite)
                            end do
                            
                            !check if there's a brigde between ksite and isite
                            connectedik = .false.
                            do bondk = 1,4
                               if(kfnsite(bondk).eq.OH12)connectedik = .true.
                            end do
                            
                            if(connectedik)then!if there's a bridge between isite and ksite
                               
                               !---CALCULATE THE 3 MEMBERED RINGS---
                               !older version
                               !check if there are any common OH to ksite and jsite
                               !do p=1,4
                               !   do q=1,4
                               !      if(kfnsite(q).eq.jfnsite(p))rings3_temp = rings3_temp + 1
                               !   end do
                               !end do
                               !------------------------------------
                               
                               !---CALCULATE THE 4 MEMBERED RINGS
                               do bondj = 1,4!loop over first neighbors of jsite
                                  OH21 = jfnsite(bondj);occOH21 = occupancy(OH21)
                                  
                                  if((OH21.ne.OH11).and.((occOH21.eq.occSNOSNO).or.&
                                       (occOH21.eq.occSIOSNO)))then!if there's another bridge of jsite
                                     
                                     do bondOH21 = 1,8!loop over first neighbors of OH21
                                        lsite = n1list(bondOH21,OH21);locc = occupancy(lsite);lvertex = head(lsite)
                                        
                                        if((lvertex.gt.0).and.(lsite.ne.jsite).and.(lsite.ne.ksite))then!if there's SN/SI on lsite
                                           
                                           !find the OH of lsite
                                           do bondl = 1,4
                                              lfn(bondl) = Si_O(bondl,lvertex)
                                              lfnsite(bondl) = n1list(lfn(bondl),lsite)
                                           end do
                                           
                                           !loop through the first neighbors of OH21 to find the bridge between lsite and jsite
                                           connectedjl = .false.
                                           do bondl = 1,4
                                              if(lfnsite(bondl).eq.OH21)connectedjl = .true.
                                           end do
                                           
                                           if(connectedjl)then
                                              !check if there's any common OH between lsite and ksite
                                              do p = 1,4
                                                 do q=1,4
                                                    if((kfnsite(p).eq.lfnsite(q)).and.(kfnsite(p).ne.OH12))&
                                                         rings4_temp = rings4_temp + 1
                                                 end do
                                              end do
                                           end if
                                           
                                        elseif((lvertex.gt.0).and.(lsite.ne.jsite).and.(lsite.eq.ksite))then!if there's SN/SI on lsite
                                              
                                           !---CALCULATE 3 MEMBERED RINGS---
                                           !newer version
                                           !loop over OH of lsite and check if it is bonded to OH12
                                           do bondk = 1,4
                                              if(kfnsite(bondk).eq.OH12)rings3_temp = rings3_temp + 1
                                           end do
                                           !--------------------------------
                                           
                                        end if!endif: there's SN/SI on lsite
                                        
                                     end do!enddo: loop over first neighbors of OH21
                                                                          
                                  end if!endif: there's another bridge of jsite site
                                  
                               end do!enddo: loop over first neighbors of jsite
                               !---------------------------------
                               
                            end if!endif: there's a bridge between isite and ksite
                            
                         end if!endif: ksite has SN/SI
                                                  
                      end do!enddo: loop over first neighbors of OH12
                      
                   end if!endif: there's a bridge between isite and jsite
                   
                end if!endif jsite has SN/SI

             end do!enddo: loop over first neighbors of OH11
             
          end if!endif: there are more then 2 brigdes
          
       end do!enddo: loop over the OH of isite
    end do!enndo: loop over the OH of isite
    
    rings3 = rings3_temp/2
    rings4 = rings4_temp/2
    
    return
  end subroutine site_rings3and4
  !---------------------------------------------------------
  !---site_energy(isite,ivertex,overlep)
  !--------------------------------------------------------- 
  subroutine site_energy(isite,ivertex,ispin,rings3,rings4,en)
    use variables
    implicit none
    
    integer,intent(in) :: isite,ivertex,ispin
    real*8,intent(out) :: en
    integer,intent(out) :: rings3,rings4
    integer :: i,bond,jocc,ihead,jsite
    integer :: interSNSN,interSISN,interISITAA
    integer :: lrx,lry,lrz,jrx,jry,jrz,mrx,mry,mrz,mtag
    integer,dimension(4) :: fn, fnsite
    integer :: lsite, msite
    logical :: linked = .false.
    
    !set the accumulators to zero
    en = real(0)
    interSNSN = 0;interSISN = 0;interISITAA = 0
    rings3 = 0;rings4 = 0
    
    !if isite has SN or SI, then calculate its neighbors
    if((ispin.eq.spinSI).or.(ispin.eq.spinSN))then
       fn(1) = ivertex
       fnsite(1) = n1list(fn(1),isite)
       ihead = fnsite(1)
       do bond = 2,4
          fn(bond) = Si_O(bond,fn(1))
          fnsite(bond) = n1list(fn(bond),isite)
       end do
    end if
   
    if(ispin.eq.spinSN)then!if isite has SN
       
       !loop over OHs and search for a bridge with SN or SI
       do bond = 1,4
          jocc = occupancy(fnsite(bond))
          if(jocc.eq.occSNOSNO)then
             interSNSN = interSNSN + 1!SN SN condensation
          elseif(jocc.eq.occSIOSNO)then
             interSISN = interSISN + 1!SI SN condensation
          end if
       end do

    elseif(ispin.eq.spinSI)then!if isite has SI
       
       !loop over OHs and search for a bridge with SN
       do bond = 2,4
          jocc = occupancy(fnsite(bond))
          if(jocc.eq.occSIOSNO)interSISN = interSISN + 1
       end do

       !search for SI-TAA linkage
       lsite = SIlink(ihead)
       if(lsite.gt.0)then
          msite = TAAlink(lsite)
          !if(msite.eq.ihead)interISITAA = interISITAA + 1
          if(msite.eq.ihead)interISITAA = 1
       end if
                    
    elseif(ispin.eq.spinTAA)then!if isite has TAA
       
       !search for SI-TAA linkage
       msite = TAAlink(isite)
       if(msite.gt.0)then
          lsite = SIlink(msite)
          !if(lsite.eq.isite)interISITAA = interISITAA + 1
          if(lsite.eq.isite)interISITAA = 1
          !interISITAA = 1
       end if
       
    else!otherwise
       print*,'BUG :: isite has neither SI/SN/TAA :: subroutine site_energy'
       print*,'jsite:',isite
       print*,'jvertex:',ivertex
       print*,'jspin:',ispin
       call freeandclose
       stop
    end if!endif: isite has SN
    
    !calculate the 3 rings and 4 rings the site is associated with
    if((ispin.eq.spinSN).or.(ispin.eq.spinSI))then
       call site_rings3and4(isite,ivertex,rings3,rings4)
    else
       rings3=0;rings4=0
    end if

    en = eSNSN*real(interSNSN) + eSISN*real(interSISN) + eISITAA*real(interISITAA) + &
         pen3*real(rings3) + pen4*real(rings4)
    
    return
  end subroutine site_energy
  !---------------------------------------------------------
  !--- total_energy: CALCULATE TOTAL ENERGY OF THE LATTICE
  !---------------------------------------------------------
  subroutine total_energy(tot_en)
    use variables
    implicit none
        
    real*8,intent(out) :: tot_en
    integer :: isite,ivertex,ispin,ihead
    integer :: i,molecule,bond,jocc,jsite
    integer :: rings3,rings4
    integer :: interSNSN,interSISN,interISITAA
    integer :: lrx,lry,lrz,jrx,jry,jrz,mrx,mry,mrz,mtag
    integer,dimension(4) :: fn,fnsite
    integer :: lsite, msite    
    
    call rings3and4(rings3,rings4)
    
    !set the accumulators to zero
    tot_en = real(0)
    interSNSN = 0;interSISN = 0;interISITAA = 0
        
    do molecule = 1,allmolecules!loop over all molecules
       
       isite = noccupy(molecule)
       ispin = spin(isite)
       ivertex = head(isite)
       msite = -1;lsite = -1
       
       if((ispin.eq.spinSI).or.(ispin.eq.spinSN))then
          fn(1) = ivertex
          fnsite(1) = n1list(fn(1),isite)
          ihead = fnsite(1)
          do bond = 2,4
             fn(bond) = Si_O(bond,fn(1))
             fnsite(bond) = n1list(fn(bond),isite)
          end do
       end if
       
       if(ispin.eq.spinSN)then!if isite has SN
          
          !loop over OHs and search for a bridge with SN or SI
          do bond = 1,4
             jocc = occupancy(fnsite(bond))
             if(jocc.eq.occSNOSNO)then
                interSNSN = interSNSN + 1!SN SN condensation
             elseif(jocc.eq.occSIOSNO)then
                interSISN = interSISN + 1!SI SN condensation
             end if
          end do
          
       elseif(ispin.eq.spinSI)then!if isite has SI
          
          !interactions taken care by SN and TAA
          !loop over OHs and search for a bridge with SN
          !do bond = 2,4
          !   jocc = occupancy(fnsite(bond))
          !   if(jocc.eq.occSIOSNO)interSISN = interSISN + 1
          !end do
          
          !search for SI-TAA linkage
          !lsite = SIlink(ihead)
          !if(lsite.gt.0)then
          !   msite = TAAlink(lsite)
          !   !if(msite.eq.ihead)interISITAA = interISITAA + 1
          !   interISITAA = interISITAA + 1
          !end if
                    
       elseif(ispin.eq.spinTAA)then!if isite has TAA
          
          !search for SI-TAA linkage
          msite = TAAlink(isite)
          if(msite.gt.0)then
             lsite = SIlink(msite)
             !if(lsite.eq.isite)interISITAA = interISITAA + 1
             interISITAA = interISITAA + 1
          end if
                    
       else!otherwise
          print*,'BUG :: isite has neither SI/SN/TAA :: subroutine site_energy'
          print*,'jsite:',isite
          print*,'jvertex:',ivertex
          print*,'jspin:',ispin
          call freeandclose
          stop
       end if!endif: isite has SN
       
    end do!enddo:loop over all molecules
    
    !calculate the energy enforced by the 3/4 membered ring penalties
    !tot_en = real(0.5)*(eSNSN*real(interSNSN) + eSISN*real(interSISN) + &
    !     eISITAA*real(interISITAA)) + pen3*real(rings3) + pen4*real(rings4)
    
    tot_en = real(0.5)*eSNSN*real(interSNSN) + eSISN*real(interSISN) + &
         eISITAA*real(interISITAA) + pen3*real(rings3) + pen4*real(rings4)
        
    return
  end subroutine total_energy
  !---------------------------------------------------------
  !---insert_SN(isite,ivertex,ispin): INSERT SN at isite
  !---------------------------------------------------------   
  subroutine insert_site(isite,ivertex,ispin,molecule)
    use variables
    implicit none
    
    integer,intent(in) :: isite,ivertex,ispin,molecule
    integer :: jsite,jocc,i,jlink,ihead
    integer :: bond,occfn
    integer,dimension(4) :: fn,fnsite
    integer :: addressTAA,addressSI
    logical :: linked
    
    !update the occupancy of isite
    if(ispin.eq.spinSN)then

       occupancy(isite) = occSN
       fn(1) = ivertex
       fnsite(1) = n1list(fn(1),isite)
       occupancy(fnsite(1)) = occupancy(fnsite(1)) + occSNO
       do bond = 2,4
          fn(bond) = Si_O(bond,fn(1))
          fnsite(bond) = n1list(fn(bond),isite)
          occupancy(fnsite(bond)) = occupancy(fnsite(bond)) + occSNO
       end do

    elseif(ispin.eq.spinSI)then

       !insert the molecule
       occupancy(isite) = occSI
       fn(1) = ivertex
       fnsite(1) = n1list(fn(1),isite)
       occupancy(fnsite(1)) = occISIO
       do bond = 2,4
          fn(bond) = Si_O(bond,fn(1))
          fnsite(bond) = n1list(fn(bond),isite)
          occupancy(fnsite(bond)) = occupancy(fnsite(bond)) + occSIO
       end do
       
       !make the connection with one of the open neighbouring TAA molecule
       ihead = fnsite(1)
       addressTAA = -100
       !SIlink(ihead) = -1
       do i = 1,6
          jsite = n2list(i,ihead);jocc = occupancy(jsite);jlink = TAAlink(jsite)
          if((jlink.eq.-1).and.(jocc.eq.occTAA))then
             !SIlink(ihead) = jsite
             !TAAlink(jsite) = ihead
             addressTAA = jsite
             !exit
          end if
       end do
       
       !make the connection between TAA and O-
       if(addressTAA.eq.-100)then
          SIlink(ihead) = -1
       elseif(addressTAA.gt.0)then
          SIlink(ihead) = addressTAA
          TAAlink(addressTAA) = ihead
       else
          print*,'wrong adressTAA :: subroutine insert_site'
          print*,'addressTAA,',addressTAA
          call freeandclose
          stop
       end if
       
    elseif(ispin.eq.spinTAA)then
       
       occupancy(isite) = occupancy(isite) + occTAA
       
       !make the connection with one of the open neighboring ISIO atoms
       !TAAlink(isite) = -1
       addressSI = -100
       do i =1,6
          jsite = n2list(i,isite);jocc = occupancy(jsite);jlink = SIlink(jsite)
          if((jlink.eq.-1).and.(jocc.eq.occISIO))then
             addressSI = jsite
             !TAAlink(isite) = jsite
             !SIlink(jsite) = isite
             !exit
          end if
       end do
       
       !make the connection between TAA and O-
       if(addressSI.eq.-100)then
          TAAlink(isite) = -1
       elseif(addressSI.gt.0)then
          TAAlink(isite) = addressSI
          SIlink(addressSI) = isite
       else
          print*,'wrong addressSI :: subroutine insert_site'
          call freeandclose
          stop
       end if
       
    else
       print*,'BUG:molecule is not SN/SI/TAA :: subroutine insert_site'
       print*,'ispin:',ispin
       call freeandclose
       stop
    end if
    
    !update head,spin,noccupy,resite arrays
    if(ispin.eq.spinSI.or.ispin.eq.spinSN)then
       head(isite) = ivertex
    elseif(ispin.eq.spinTAA)then
       head(isite) = 0
    end if
    spin(isite) = ispin
    noccupy(molecule) = isite
    resite(isite) = molecule

    return
  end subroutine insert_site
  !---------------------------------------------------------
  !---remove_SN(isite,ivertex): REMOVE SN from isite
  !---------------------------------------------------------   
  subroutine remove_site(isite,ivertex,ispin,molecule)
    use variables
    implicit none
    
    integer,intent(in) :: isite, ivertex,ispin,molecule
    integer :: bond,occfn,ihead,jsite,jocc,jlink,i
    integer,dimension(4) :: fn,fnsite
    
    !find the 1st neighbors of isite and update their occupancy
    if(ispin.eq.spinSN)then

       occupancy(isite) = occupancy(isite) - occSN
       fn(1) = ivertex
       fnsite(1) = n1list(fn(1),isite)
       occupancy(fnsite(1)) = occupancy(fnsite(1)) - occSNO
       do bond = 2,4
          fn(bond) = Si_O(bond,fn(1))
          fnsite(bond) = n1list(fn(bond),isite)
          occupancy(fnsite(bond)) = occupancy(fnsite(bond)) - occSNO
       end do

    elseif(ispin.eq.spinSI)then
       
       !remove any established connection with TAA
       ihead = n1list(ivertex,isite)
       jsite = SIlink(ihead)
       if(jsite.gt.0)then
          jocc = occupancy(jsite)
          if((TAAlink(jsite).eq.ihead).and.(jocc.eq.occTAA))then
             SIlink(ihead) = -1
             TAAlink(jsite) = -1
          else
             print*,'Bug: incorrect connection between O- and TAA (SI removal): subroutine remove_site'
             print*,'ihead/jsite',ihead,jsite
             print*,'SIlink(ihead)/TAAlink(jsite)',Silink(ihead),TAAlink(jsite)
             print*,'jocc/occTAA',jocc,occTAA
             call freeandclose
             stop
          end if
       else
          SIlink(ihead) = -1
       end if
       
       !remove the molecule
       occupancy(isite) = occupancy(isite) - occSI
       fn(1) = ivertex
       fnsite(1) = n1list(fn(1),isite)
       occupancy(fnsite(1)) = occupancy(fnsite(1)) - occISIO
       do bond = 2,4
          fn(bond) = Si_O(bond,fn(1))
          fnsite(bond) = n1list(fn(bond),isite)
          occupancy(fnsite(bond)) = occupancy(fnsite(bond)) - occSIO
       end do

    elseif(ispin.eq.spinTAA)then
       
       !remove any established connection with ISIO
       jsite = TAAlink(isite)
       if(jsite.gt.0)then
          jocc = occupancy(jsite)
          if((SIlink(jsite).eq.isite).and.(jocc.eq.occISIO))then
             TAAlink(isite) = -1
             SIlink(jsite) = -1
          else
             print*,'Bug: incorrect connection between O- and TAA (TAA removal): subroutine remove_site'
             call freeandclose
             stop
          end if
       else
          TAAlink(isite) = -1
       end if
       
       !remove the molecule
       occupancy(isite) = occupancy(isite) - occTAA
       
    else
       print*,'BUG:molecule is not SN/SI/TAA :: subroutine remove_site'
       print*,'ispin:',ispin
       call freeandclose
       stop
    end if

    !update the spin and head arrays
    !note: resite and noccupy arrays are not updated
    spin(isite) = 0
    head(isite) = 0
    
    !note:noccupy/resite arrays are not updated
    !they are updated in the subroutine insert_site
    return
  end subroutine remove_site
  !---------------------------------------------------------
  !---LINKAGE: this subroutine restores the SI-TAA linkage if an MC move is rejected
  !---------------------------------------------------------   
  subroutine restore_linkage(isite,ispin,ivertex,oldSIlinkage,oldTAAlinkage,move_accepted)
    use variables
    implicit none
    
    !isite => has SI/TAA
    !jsite => has the corresponding TAA/SI
    
    integer,intent(in) :: isite,ispin,ivertex,oldSIlinkage,oldTAAlinkage
    logical,intent(in) :: move_accepted
    integer :: jsite, ihead
    
    if(move_accepted)then!if the move has been accepted
       !do nothing
    else!if the move has been rejected
       
       if(ispin.eq.spinSI)then!if isite has SI
          
          ihead = n1list(ivertex,isite)
          !jsite = SIlink(ihead)
          SIlink(ihead) = oldSIlinkage
          if(oldSIlinkage.gt.0)TAAlink(oldSIlinkage) = oldTAAlinkage
          
       elseif(ispin.eq.spinTAA)then!if isite has TAA
          
          !jsite = TAAlink(isite)
          TAAlink(isite) = oldTAAlinkage
          if(oldTAAlinkage.gt.0)SIlink(oldTAAlinkage) = oldSIlinkage
          
       else
          print*,'BUG: incorrect SI-TAA linkage :: subroutine restore_linkage'
          print*,'isite,ispin,ivertex,oldSIlinkage,oldTAAlinkage',isite,ispin,ivertex,oldSIlinkage,oldTAAlinkage
          call freeandclose
          stop
       end if!endif: if isite has SI   
       
    end if!endif: the move has been accepted
    
    return
  end subroutine restore_linkage
  !---------------------------------------------------------
  !---LINKAGE: this subroutine restores the SI-TAA linkage if an MC move is rejected
  !---------------------------------------------------------   
  subroutine restore_linkage1(isite,ispin,ivertex,oldSIlinkage,oldTAAlinkage,move_accepted)
    use variables
    implicit none
    
    !isite => has SI/TAA
    !jsite => has the corresponding TAA/SI
    
    integer,intent(in) :: isite,ispin,ivertex,oldSIlinkage,oldTAAlinkage
    logical,intent(in) :: move_accepted
    integer :: jsite, ihead
    
    if(move_accepted)then!if the move has been accepted
       !do nothing
    else!if the move has been rejected
       
       if(ispin.eq.spinSI)then!if isite has SI
          
          ihead = n1list(ivertex,isite)
          jsite = SIlink(ihead)
          if((oldSIlinkage.lt.0).and.(jsite.gt.0))then!if the SI was ealier not connected to any TAA but is now
             SIlink(ihead) = -1
             TAAlink(jsite) = -1
          end if!if the SI was ealier not connected to any TAA but is now 
          
       elseif(ispin.eq.spinTAA)then!if isite has TAA
          
          jsite = TAAlink(isite)
          if((oldTAAlinkage.lt.0).and.(jsite.gt.0))then
             TAAlink(isite) = -1
             SIlink(jsite) = -1
          end if
          
       else
          print*,'BUG: incorrect SI-TAA linkage :: subroutine restore_linkage'
          print*,'isite,ispin,ivertex,oldSIlinkage,oldTAAlinkage',isite,ispin,ivertex,oldSIlinkage,oldTAAlinkage
          call freeandclose
          stop
       end if!endif: if isite has SI   
       
    end if!endif: the move has been accepted
    
    return
  end subroutine restore_linkage1
  !---------------------------------------------------------
  !---nvt_mc: NVT MONTE CARLO ALGORITHM
  !--------------------------------------------------------- 
  subroutine nvt_mc
    use variables
    implicit none
    
    integer :: molecule,bond
    integer :: isite,ispin,ivertex,imolecule
    integer :: jsite,jspin,jvertex,jmolecule
    integer :: ksite,kspin,kvertex,kmolecule
    integer :: occfn,overlap_old,overlap_new
    integer :: rings3_old,rings3_new,rings3_temp,rings4_old,rings4_new,rings4_temp,rings3,rings4
    real*8 :: old_en,new_en,delta_en,tot_en,temp_en,norm_en
    real*8 :: AR1,AR2,AR3,rand
    real*8,dimension(0:4) :: norm_Qn
    
    !calculate total energy, Qn distribution and 3 and 4 ring 
    !distribution before simulation
    call Qn_dist
    call total_energy(tot_en)
    call rings3and4(rings3,rings4)
    call cluster_stat
    norm_Qn = real(Qn)/real(size(noccupy))
    norm_en = tot_en/real(nsites)
    c = (norm_Qn(1) + norm_Qn(2)*2.0d0 + norm_Qn(3)*3.0d0 + norm_Qn(4)*4.0d0)/4.0d0
    if(initial_step.eq.1)then
       write(1001,*) step,norm_Qn,c;flush(1001)
       write(1002,*) step,rings3,rings4;flush(1002)
       write(1003,*) step,tot_en,norm_en;flush(1003)
       write(1007,*) step,max_size,avg_size,cluster_num;flush(1007)
    end if

    do step = initial_step,nsweeps!do all MC steps
       
       do molecule = 1,allmolecules!loop over all molecules

          !*********************************START SWAP MOVES*********************************          
          imolecule = int(1 + ran3(seed)*real(allmolecules))
          isite = noccupy(imolecule)
          ivertex = head(isite)
          ispin = spin(isite)
          
          !select a random site
          jsite = int(1 + ran3(seed)*real(nsites))
          do while(jsite.eq.isite)
             jsite = int(1 + ran3(seed)*real(nsites))
          end do
          jvertex = head(jsite)
          jspin = spin(jsite)
          jmolecule = resite(jsite)
          
          if(occupancy(isite).eq.occSN)then
             
             if(occupancy(jsite).eq.occW)then
                call translate(isite,ivertex,ispin,imolecule,jsite,rings3,rings4,tot_en,AR1)
             elseif(occupancy(jsite).eq.occSI)then
                call swap(jsite,jvertex,jspin,jmolecule,isite,ivertex,ispin,imolecule,&
                     rings3,rings4,tot_en,AR2)
             elseif(occupancy(jsite).eq.occTAA)then
                call swap(isite,ivertex,ispin,imolecule,jsite,jvertex,jspin,jmolecule,&
                     rings3,rings4,tot_en,AR2)
             end if
             
          elseif(occupancy(isite).eq.occSI)then
             
             if(occupancy(jsite).eq.occW)then
                call translate(isite,ivertex,ispin,imolecule,jsite,rings3,rings4,tot_en,AR1)
             elseif(occupancy(jsite).eq.occTAA)then
                call swap(isite,ivertex,ispin,imolecule,jsite,jvertex,jspin,jmolecule,&
                     rings3,rings4,tot_en,AR2)
             elseif(occupancy(jsite).eq.occSN)then
                call swap(isite,ivertex,ispin,imolecule,jsite,jvertex,jspin,jmolecule,&
                    rings3,rings4,tot_en,AR2)
             end if

          elseif(occupancy(isite).eq.occTAA)then

             if(occupancy(jsite).eq.occW)then
                call translate(isite,ivertex,ispin,imolecule,jsite,rings3,rings4,tot_en,AR1)
             elseif(occupancy(jsite).eq.occSI)then
                call swap(isite,ivertex,ispin,imolecule,jsite,jvertex,jspin,jmolecule,&
                     rings3,rings4,tot_en,AR2)
             elseif(occupancy(jsite).eq.occSN)then
                call swap(jsite,jvertex,jspin,jmolecule,isite,ivertex,ispin,imolecule,&
                     rings3,rings4,tot_en,AR2)    
             end if
             
          else
             print*,'BUG :: molecule neither SI/SN/TAA :: subroutine nvt_mc'
             print*,'isite:',isite
             print*,'ivertex:',ivertex
             print*,'imolecule:',imolecule
             print*,'ispin:',ispin
             call freeandclose
             stop
          end if
          !*********************************END SWAP MOVES***********************************

          !*********************************START ROTATION MOVES*****************************
          !select a monomer amongst all
          kmolecule = int(1 + ran3(seed)*real(allmolecules))
          ksite = noccupy(kmolecule)
          kspin = spin(ksite)
          kvertex = head(ksite)
          if((occupancy(ksite).eq.occSI).or.(occupancy(ksite).eq.occSN))then!if ksite has SI/SN on it
             call rotate(ksite,kvertex,kspin,kmolecule,rings3,rings4,tot_en,AR3)
          elseif(kspin.eq.spinTAA)then
             !reject the move
          else
             print*,'BUG :: molecule for rotation neither SI/SN/TAA :: subroutine nvt_mc'
             print*,'isite:',ksite
             print*,'ivertex:',kvertex
             print*,'imolecule:',kmolecule
             print*,'ispin:',kspin     
             call freeandclose
             stop
          end if!endif: ksite has SI/SN on it
          !*********************************END ROTATION MOVES*******************************
          
       end do!endloop over all molecules
       
       !print the output of the program
       call print_output(AR1,AR2,AR3,rings3,rings4,tot_en)
       
       !change the temperature
       !call update_temperature(step)
       !call linkage
       
    end do!enddo all MC steps
    
    return
  end subroutine nvt_mc
  !---------------------------------------------------------
  !---translate: translates molecule from isite to jsite(water site)
  !--------------------------------------------------------- 
  subroutine translate(isite,ivertex,ispin,imolecule,jsite,rings3,rings4,tot_en,AR)
    use variables
    implicit none
    
    integer,intent(in) :: isite,ivertex,ispin,imolecule,jsite
    integer,intent(inout) :: rings3,rings4
    real*8,intent(inout) :: tot_en
    real*8,intent(out) :: AR
    integer :: rings3_old,rings3_new,rings4_old,rings4_new
    real*8 :: old_en,new_en,delta_en,rand
    logical :: can_insert
    
    !first check if the molecule can be inserted at jsite
    if(ispin.eq.spinSN)then
       call can_insert_SN(jsite,ivertex,can_insert)
    elseif(ispin.eq.spinSI)then
       call can_insert_SI(jsite,ivertex,can_insert)
    elseif(ispin.eq.spinTAA)then
       call can_insert_TAA(jsite,can_insert)
    else
       print*,'BUG :: molecule neither SI/SN/TAA :: subroutine translate'
       print*,'isite:',isite
       print*,'ivertex:',ivertex
       print*,'imolecule:',imolecule
       print*,'ispin:',ispin
    end if
    
    if(can_insert)then!if molecule at isite can be inserted at jsite
       
       !calculate the energy associated with the molecule at isite
       call site_energy(isite,ivertex,ispin,rings3_old,rings4_old,old_en)
       !remove molecule from isite and insert it at jsite
       call remove_site(isite,ivertex,ispin,imolecule)
       
       !insert the molecule at jsite and calculate the new energy associated with it
       call insert_site(jsite,ivertex,ispin,imolecule)
       call site_energy(jsite,ivertex,ispin,rings3_new,rings4_new,new_en)
       
       !calculate the change in energy and apply the metropolis criteria
       delta_en = new_en - old_en
       AR = exp(-1.000d0*delta_en/tstar)
       rand = ran3(seed)
       if(AR.gt.rand)then!if metropolis criteria is satisfied
          tot_en = tot_en + delta_en
          rings3 = rings3 + (rings3_new - rings3_old)
          rings4 = rings4 + (rings4_new - rings4_old)
       else!if metropolis criteria is not satisfied
          !remove the molecule from jsite and put it back on isite
          call remove_site(jsite,ivertex,ispin,imolecule)
          call insert_site(isite,ivertex,ispin,imolecule)
       end if!endif:metropolis criteria is satisfied

    end if!if molecule at isite can be inserted at jsite
    
    !safety check
    !call safetycheck_translate(tot_en,ispin,AR,rand,old_en,new_en,rings3,rings4)
    !call linkage
    
    return
  end subroutine translate
  !---------------------------------------------------------
  !---swap: swaps two molecules
  !---note: molecule at jsite would be removed first 
  !---------------------------------------------------------   
  subroutine swap(isite,ivertex,ispin,imolecule,jsite,jvertex,jspin,jmolecule,rings3,rings4,tot_en,AR)
    use variables
    implicit none
    
    integer,intent(in) :: isite,ivertex,ispin,imolecule
    integer,intent(in) :: jsite,jvertex,jspin,jmolecule
    integer,intent(inout) :: rings3,rings4
    real*8,intent(inout) :: tot_en
    real*8,intent(out) :: AR
    integer :: rings3_old,rings3_new,rings4_old,rings4_new,rings3_temp,rings4_temp
    real*8 :: old_en,new_en,delta_en,new_temp_en,old_temp_en,rand
    logical :: can_inserti, can_insertj
    !variables for checking double couting in energy
    integer :: ihead,jhead,hrx,hry,hrz,irx,iry,irz,jrx,jry,jrz
    real :: d_ij,d_ih,d_jh,d_ij_sq,d_ih_sq,d_jh_sq
    
        
    !calculate the site energy associated with the molecule located at isite
    call site_energy(jsite,jvertex,jspin,rings3_temp,rings4_temp,old_temp_en)
    !remove the molecule from isite
    call remove_site(jsite,jvertex,jspin,jmolecule)
    
    !check if the molecule at isite can be inserted at jsite
    if(ispin.eq.spinSN)then
       call can_insert_SN(jsite,ivertex,can_inserti)
    elseif(ispin.eq.spinSI)then
       call can_insert_SI(jsite,ivertex,can_inserti)
    elseif(ispin.eq.spinTAA)then
       call can_insert_TAA(jsite,can_inserti)
    else
       print*,'BUG :: imolecule neither SI/SN/TAA :: subroutine swap'
       print*,'isite:',isite
       print*,'ivertex:',ivertex
       print*,'imolecule:',imolecule
       print*,'ispin:',ispin
    end if
    
    if(can_inserti)then!if the molecule at isite can be inserted at jsite
       
       !calculate the energy associated with the molecule located at isite
       call site_energy(isite,ivertex,ispin,rings3_old,rings4_old,old_en)
       old_en = old_en + old_temp_en
       rings3_old = rings3_old + rings3_temp
       rings4_old = rings4_old + rings4_temp
       !now remove imolecule from isite
       call remove_site(isite,ivertex,ispin,imolecule)
       
       !now check if jmolecule can be inserted at isite
       if(jspin.eq.spinSN)then
          call can_insert_SN(isite,jvertex,can_insertj)
       elseif(jspin.eq.spinSI)then
          call can_insert_SI(isite,jvertex,can_insertj)
       elseif(jspin.eq.spinTAA)then
          call can_insert_TAA(isite,can_insertj)
       else
          print*,'BUG :: jmolecule neither SI/SN/TAA :: subroutine swap'
          print*,'jsite:',jsite
          print*,'jvertex:',jvertex
          print*,'jmolecule:',jmolecule
          print*,'jspin:',jspin
       end if

       if(can_insertj)then!if jmolecule can be inserted at isite

          !put imolecule on jsite and calculate its energy
          call insert_site(jsite,ivertex,ispin,imolecule)
          call site_energy(jsite,ivertex,ispin,rings3_temp,rings4_temp,new_temp_en)
          !put imolecule on isite and calculate its energy
          call insert_site(isite,jvertex,jspin,jmolecule)
          call site_energy(isite,jvertex,jspin,rings3_new,rings4_new,new_en)

          ! calculate the new energy
          new_en = new_en + new_temp_en
          rings3_new = rings3_new + rings3_temp
          rings4_new = rings4_new + rings4_temp
          
          !calculate the change in energy after the move and apply the metropolis criteria
          delta_en = new_en - old_en
          AR = exp(-1.000000d0*delta_en/tstar)
          rand = ran3(seed)
          if(AR.gt.rand)then!if the netropolis criteria is satisfied
             tot_en = tot_en + delta_en
             rings3 = rings3 + (rings3_new - rings3_old)
             rings4 = rings4 + (rings4_new - rings4_old)
          else!if the netropolis criteria is not satisfied
             !remove imolecule from jsite and jmolecule from isite
             call remove_site(jsite,ivertex,ispin,imolecule)
             call remove_site(isite,jvertex,jspin,jmolecule)
             !put the molecules back in their original places
             call insert_site(jsite,jvertex,jspin,jmolecule)
             call insert_site(isite,ivertex,ispin,imolecule)
          end if!endif: the netropolis criteria is satisfied
                    
       else!if jmolecule cannot be inserted at isite

          !put both the molecules on their original location
          call insert_site(jsite,jvertex,jspin,jmolecule)
          call insert_site(isite,ivertex,ispin,imolecule)

       end if!endif: jmolecule can be inserted at isite
       
    else!if the molecule at isite cannot be inserted at jsite
       !put the jmolecule back on jsite
       call insert_site(jsite,jvertex,jspin,jmolecule)
    end if!endif: the molecule at isite can be inserted at jsite
    
    !safety check
    !call safetycheck_swap(tot_en,isite,jsite,AR,rand,old_en,new_en,rings3,rings4)
    !call linkage
        
    return
  end subroutine swap
  !---------------------------------------------------------
  !---rotate: rotates the molecule in its position
  !---------------------------------------------------------  
  subroutine rotate(ksite,kvertex,kspin,kmolecule,rings3,rings4,tot_en,AR)
    use variables
    implicit none
    
    integer,intent(in) :: ksite,kvertex,kspin,kmolecule
    integer,intent(inout) :: rings3,rings4
    real*8,intent(inout) :: tot_en
    real*8,intent(out) :: AR
    integer :: kvertex_new, oldSIlinkage,oldTAAlinkage
    integer :: rings3_old,rings3_new,rings4_old,rings4_new
    real*8 :: old_en,new_en,delta_en,rand
    logical :: can_rotate,move_accepted
    
    real*8 :: tot_en_old,tot_en_new
    integer :: SIconnected, TAAconnected
        
    oldSIlinkage=0;oldTAAlinkage=0
    move_accepted = .false.
    rotation_attempted = rotation_attempted + 1
    
    if((step.eq.3).and.(ksite.eq.49052))then
       call safetycheck_linkage(SIconnected,TAAconnected)
    end if
    
    if(kspin.eq.spinSN)then
       kvertex_new = int(1 + ran3(seed) * real(2))
       if(kvertex_new.eq.kvertex)return
    elseif(kspin.eq.spinSI)then
       kvertex_new = int(1 + ran3(seed) * real(nc))
       oldSIlinkage = SIlink(n1list(kvertex,ksite))
       if(oldSIlinkage.gt.0)oldTAAlinkage = TAAlink(oldSIlinkage)
       if(kvertex_new.eq.kvertex)return
    else
       print*,'BUG :: molecule not SI/SN :: subroutine rotate'
       print*,'kspin:',kspin
       print*,'kvertex:',kvertex
       print*,'kmolecule:',kmolecule
       print*,'ksite:',ksite
       call freeandclose
       stop
    end if
    
    !calculate the initial energy associated with the molecule ksite and then remove it
    call site_energy(ksite,kvertex,kspin,rings3_old,rings4_old,old_en)
    call remove_site(ksite,kvertex,kspin,kmolecule)
    
    !check if the molecule can be inserted again at ksite after rotation
    !rotation is defined as the assignment of the new pointer variable
    if(kspin.eq.spinSN)then
       call can_insert_SN(ksite,kvertex_new,can_rotate)
    elseif(kspin.eq.spinSI)then
       call can_insert_SI(ksite,kvertex_new,can_rotate)
    end if
    
    if(can_rotate)then
       
       !rotate the molecule at ksite and calculate the new energy associated with the molecule
       call insert_site(ksite,kvertex_new,kspin,kmolecule)
       call site_energy(ksite,kvertex_new,kspin,rings3_new,rings4_new,new_en)

       !calculate the change in energy and apply the metropolis criteria
       delta_en = new_en - old_en
       AR = exp(-1.000d0*delta_en/tstar)
       rand = ran3(seed)
       if(AR.gt.rand)then!if the metropolis criteria is satisfied
          
          !update energy and rings
          tot_en = tot_en + delta_en
          rings3 = rings3 + (rings3_new - rings3_old)
          rings4 = rings4 + (rings4_new - rings4_old)
          move_accepted = .true.
          rotation_accepted = rotation_accepted + 1
          
       else!if the metropolis criteria is not satisfied
          
          !put the molecule back in its original orientation
          call remove_site(ksite,kvertex_new,kspin,kmolecule)
          call insert_site(ksite,kvertex,kspin,kmolecule)          
          
          return
       end if!endif: the metropolis criteria is satisfied
       
    else
       !put the molecule back in its original orientation
       call insert_site(ksite,kvertex,kspin,kmolecule)
       return
    end if
    !restore the linkage if the move is rejected
    if(kspin.eq.spinSI)call restore_linkage(ksite,kspin,kvertex,oldSIlinkage,oldTAAlinkage,move_accepted)
    
    !safety check
    !call total_energy(test)
    !call safetycheck_linkage(SIconnected,TAAconnected)
    call safetycheck_rotate(tot_en,ksite,kspin,kvertex,kvertex_new,AR,rand,old_en,new_en,delta_en,rings3,rings4)
    
    return
  end subroutine rotate
  !---------------------------------------------------------
  !---print_output: prints the output into a file
  !--------------------------------------------------------- 
  subroutine print_output(AR1,AR2,AR3,rings3,rings4,en)
    use variables
    implicit none
    
    integer :: istep,points
    integer :: rings3,rings4
    real :: deltapoint,ratio,avg_particle_size
    real :: dmin,dmax
    real*8 :: AR1,AR2,AR3,en,norm_en,temp1,temp2,temp3,temp4
    integer,dimension(0:4) :: old_Qn,old_Qnn
    real*8,dimension(0:4) :: norm_Qn,new_Qn,norm_Qnn,new_Qnn
    integer :: SIconnected,TAAconnected
    
    !safety checking
    !if(mod(step,2*nprint).eq.0)call safetycheck(rings3,rings4,en)
    call safetycheck(rings3,rings4,en)
    call safetycheck_linkage(SIconnected,TAAconnected)
    
    !write the output configuration
    if(mod(step,2*nprint).eq.0)call config_out_atom

    !output Qn distribution
    !points for log plot of Qn distribution b/w 100 and 1000
    points = 8;ratio = real(1)/real(points)
    if(step.le.10)then
       call Qn_dist
       norm_Qn = real(Qn)/real(size(noccupy))
       norm_Qnn = real(Qnn)/real(ni(2))
       c = (norm_Qn(1) + norm_Qn(2)*2.0d0 + norm_Qn(3)*3.0d0 + norm_Qn(4)*4.0d0)/4.0d0
       cn = (norm_Qnn(1) + norm_Qnn(2)*2.0d0 + norm_Qnn(3)*3.0d0 + norm_Qnn(4)*4.0d0)/4.0d0
       write(1001,*) step,norm_Qn,c;flush(1001)
       write(1015,*) step,norm_Qnn,cn;flush(1015)
    else
       deltapoint = 10**(ratio)
       istep = int(deltapoint**point)
       if(istep.eq.step)then
          old_Qn = Qn
          old_Qnn = Qnn
          call Qn_dist
          new_Qn = real(Qn + old_Qn)/real(2)!average the Qn with the old value to get a smooth curve
          new_Qnn = real(Qnn + old_Qnn)/real(2)!average the Qnn with the old value to get a smooth curve
          norm_Qn = real(new_Qn)/real(size(noccupy))
          norm_Qnn = real(new_Qnn)/real(ni(2))
          c = (norm_Qn(1) + norm_Qn(2)*2.0d0 + norm_Qn(3)*3.0d0 + norm_Qn(4)*4.0d0)/4.0d0
          cn = (norm_Qnn(1) + norm_Qnn(2)*2.0d0 + norm_Qnn(3)*3.0d0 + norm_Qnn(4)*4.0d0)/4.0d0
          write(1001,*) step,norm_Qn,c;flush(1001)
          write(1015,*) step,norm_Qnn,cn;flush(1015)
          point = point + 1
       end if
    end if

    if(mod(step,nprint).eq.0)then

       !output cluster statistics only when all SN are present
       call cluster_stat

       !find the negative charge per SI and sublattice ordering
       call negchargeperSI       
       !call sublattice_ordering
       call shell_order

       !calculate the average spatial dimensions of the particles
       !call particle_size(avg_particle_size)
       
       norm_en = en/real(nsites)
       write(*,8000) 'MC step:',step
       write(*,9000) 'Acceptance Ratios:',AR1,AR2,AR3
       !write(*,10000) 'Total Energy:',en
       !write(*,11000) 'Normalized Qn Distribution:',norm_Qn;flush(1008)
       write(*,10000) 'Energy per lattice site:',norm_en
       !move summary
       !write(*,1000) 'Move Summary'
       !write(*,2000) 'SN(translations,rotations):',SN_translate_success,SN_rotate_success
       !write(*,2000) 'SI(translations,rotations):',SI_translate_success,SI_rotate_success
       !write(*,8000) 'SI-SN swaps:',SISN_swap_success
       !write(*,8000) 'SI-TAA swaps:',SITAA_swap_success
       !write(*,8000) 'SN-TAA swaps:',SNTAA_swap_success    
       write(*,*)
    
       write(1008,8000) 'MC step:',step;flush(1008)
       write(1008,9000) 'Acceptance Ratios:',AR1,AR2,AR3;flush(1008)
       !write(1008,10000) 'Total Energy:',en;flush(1008)
       !write(1008,11000) 'Normalized Qn Distribution:',norm_Qn;flush(1008)
       write(1008,10000) 'Energy per lattice site:',norm_en;flush(1008)
       !move summary
       !write(1008,1000) 'Move Summary';flush(1008)
       !write(1008,2000) 'SN(translations,rotations):',SN_translate_success,SN_rotate_success;flush(1008)
       !write(1008,2000) 'SI(translations,rotations):',SI_translate_success,SI_rotate_success;flush(1008)
       !write(1008,2000) 'TAA(translations):',TAA_translate_success;flush(1008)
       !write(1008,8000) 'SI-SN swaps:',SISN_swap_success;flush(1008)
       !write(1008,8000) 'SI-TAA swaps:',SITAA_swap_success;flush(1008)
       !write(1008,8000) 'SN-TAA swaps:',SNTAA_swap_success;flush(1008)
       write(1008,*);flush(1008)
       
       write(1002,*) step,rings3,rings4;flush(1002)
       write(1003,*) step,en,norm_en;flush(1003)
       write(1007,*) step,max_size,avg_size,cluster_num;flush(1007)
       write(1012,*) step,frac_SITAA,frac_SISNTAA;flush(1012)
    end if
    
    !calculate the average values of various quantities
    if((step.gt.neqsweeps).and.(mod(step,nprint).eq.0))then
       final_avg_energy = final_avg_energy + norm_en
       final_avg_max_size = final_avg_max_size + max_size
       final_avg_cluster_num = final_avg_cluster_num + cluster_num
       final_avg_avg_size = final_avg_avg_size + avg_size
       
       temp1 = final_avg_energy/real((step-neqsweeps)/nprint)
       temp2 = final_avg_max_size/real((step-neqsweeps)/nprint)
       temp3 = final_avg_cluster_num/real((step-neqsweeps)/nprint)
       temp4 = final_avg_avg_size/real((step-neqsweeps)/nprint)
       
       write(1014,12000) step,temp1,temp2,temp4,temp3;flush(1014)

    end if
    
    !save the snapshot of the lattice as well as the spatial dimensions
    !the nanoparticles
    if(mod(step,nsnapshot).eq.0)then
       snapshot = snapshot + 1
       call take_snapshot
       
       !call size_range(dmin,dmax)
       !write(1018,*)step,dmin,dmax;flush(1018)

    end if
    
    call cpu_time(end_time)
    runtime = (end_time-start_time)/3600.0d0
    if(runtime.gt.time_limit)then

       !output the final state of the syste
       call output_state
       
       write(*,10000),'Program runtime(hrs):',runtime
       write(1008,10000)'Program runtime(hrs):',runtime
       write(*,*);write(*,*)
       
       call freeandclose
       stop

    end if
           
1000 format(x,A)
2000 format(x,A,x,i8,x,i8)
8000 format(x,A,x,i8)
9000 format(x,A,5x,es10.4,5x,es10.4,5x,es10.4)
!9000 format(x,A,5x,1e10.5,5x,1e10.5,5x,1e10.5)
10000 format(x,A,es12.4)
11000 format(x,A,5x,es10.4,5x,es10.4,5x,es10.4,5x,es10.4)
12000 format(x,i10,x,f10.4,x,f10.4,x,f10.4,x,f10.4)
    
    return
  end subroutine print_output
  !---------------------------------------------------------
  !---UPDATE_TEMPERATURE: updates a temperature according to a profile
  !--------------------------------------------------------- 
  subroutine update_temperature(step)
    implicit none
    
    integer,intent(in) :: step
    real :: frac

    frac = real(step)/real(nsweeps)
    
    if((frac.gt.m0).and.(frac.lt.m1))then
       tstar = tstar_initial + (frac - m0)*(up_rate)
    elseif((frac.gt.m1).and.(frac.lt.m2))then
       !do nothing
    elseif(frac.gt.m2)then
       tstar = tstar_initial + (tstar_final-tstar_initial)/exp(down_rate*real(frac-m2))
    end if
    
    write(1009,*)step,tstar
    
    return
  end subroutine update_temperature
  !---------------------------------------------------------
  !---cluster_stat(max_size):cluster counting subroutine
  !---------------------------------------------------------   
  subroutine cluster_stat
    use variables
    implicit none
    
    integer :: isite,ivertex,ispin
    integer :: p,q,monomer,neighs,scanned_neighs,cluster_size,same
    integer :: jsite,jlabel,pjlabel
    integer :: ksite,klabel,pklabel
    integer :: min_label,largest_label
    integer :: cluster_units
    logical :: arethereSi
    integer,dimension(4) :: mark
    !integer,dimension(:),allocatable :: neighSi
    integer,dimension(4) :: neighSi
    integer :: check_size,check_label
    
    csize = -1
    clabel = -1
    largest_label = 0
        
    do monomer = 1,allmolecules!loop over all monomers
       
       isite = noccupy(monomer)
       ivertex = head(isite)
       ispin = spin(isite)
       
       if((occupancy(isite).eq.occSN).or.(occupancy(isite).eq.occSI))then!if isite has SN/SI

          neighSi=0
          call returnneigh(isite,neighSi,arethereSi)
          
          if(arethereSi)then!if there are neighboring SI/SN to isite

             neighs = size(neighSi)
             !check for scanned neighbors
             mark=-1;scanned_neighs = 0
             do p=1,neighs
                if(neighSi(p).gt.0)then
                   jsite = neighSi(p);jlabel = clabel(jsite)
                   if(jlabel.gt.0)then
                      scanned_neighs = scanned_neighs + 1
                      mark(p) = jsite
                   end if
                end if
             end do
             
             if(scanned_neighs.eq.1)then!if there's one scanned neighbor
                
                !if there only one scanned neighbor
                !assign isite the same label as the neighbor
                !and increment the size of the proper label of the neighbor by 1
                do p=1,neighs
                   if(mark(p).gt.0) jsite = mark(p)
                end do
                jlabel = clabel(jsite)
                call classify(jlabel,pjlabel)
                clabel(isite) = jlabel
                csize(pjlabel) = csize(pjlabel) + 1             
                
             elseif(scanned_neighs.ge.2)then!if there's more than one scanned neighbors
                !----------------------------------------------------------------------
                !find out the minimum proper label of the scanned neighbors
                min_label = allmolecules
                do p=1,neighs!loop over neighbors
                   if(mark(p).gt.0)then!if the site is scanned

                      jsite = mark(p);jlabel = clabel(jsite)
                      call classify(jlabel,pjlabel)
                      if(min_label.gt.pjlabel)then
                         min_label = pjlabel
                      end if

                   end if!endif: the site is scanned
                end do!enddo: loop over neighbors
                
                !coalesce the cluster into a larger one
                cluster_size = 0
                do p =1,neighs !loop over all neighbors
                   if(mark(p).gt.0)then!if the site is scanned
                      
                      jsite = mark(p);jlabel = clabel(jsite)
                      call classify(jlabel,pjlabel)
                      !check if pjlabel is a repetitive label amongst other labels
                      same = 0!variable to check if the label has already been operated upon
                      do q = 1,p-1
                         if(mark(q).gt.0)then
                            ksite = mark(q);klabel = clabel(ksite)
                            call classify(klabel,pklabel)
                            if(pklabel.eq.pjlabel) same = 1
                         end if
                      end do
                      
                      if((same.eq.0).and.(pjlabel.ne.min_label))then
                         cluster_size = cluster_size + csize(pjlabel)
                      end if
                   end if!endif: the site is scanned
                end do!enddo: loop over all neighbors
                
                !update the size of min_label
                csize(min_label) = csize(min_label) + cluster_size + 1
                
                !loop over the neighbors and set their sizes to zero except the min_label
                do p=1,neighs
                   if(mark(p).gt.0)then
                      jsite = mark(p);jlabel = clabel(jsite)
                      call classify(jlabel,pjlabel)
                      if(pjlabel.ne.min_label) csize(pjlabel) = -min_label
                   end if
                end do
                clabel(isite) = min_label!set the label of isite to the min_label
                !----------------------------------------------------------------------                
             else!if there's no scanned neighbor
                
                largest_label = largest_label + 1
                clabel(isite) = largest_label
                csize(largest_label) = 1
                
             end if!endif there's one scanned neighbor
             
          else!if there are no neighboring SI/SN to isite
             
             largest_label = largest_label + 1
             clabel(isite) = largest_label
             csize(largest_label) = 1
             
          end if!endif: there are neighboring SI/SN to isite
                    
       end if!endif: isite has SN/SI
       !if(allocated(neighSi)) deallocate(neighSi)

    end do!enddo: loop over all monomers
    
    !safety check
    check_size = 0
    do p=1,size(csize)
       if(csize(p).gt.0)check_size = check_size + csize(p)
    end do
    if(check_size.ne.(allmolecules-ni(3)))then
       print*,'BUG: cluster calculation error :: subroutine cluster_stat'
       print*,'check_size:',check_size
       print*,'monomer:',allmolecules
       call freeandclose
       stop       
    end if

    !calculate the various cluster statistics
    avg_size = 0;cluster_num = 0;cluster_units=0
    do p=1,size(csize)
       if(csize(p).ge.cluster_threshold)then
          cluster_num = cluster_num + 1
          cluster_units = cluster_units + csize(p)
       end if
    end do
    
    if(cluster_num.gt.0)then
       avg_size = real(cluster_units)/real(cluster_num)
       max_size = maxval(csize)
    else
       avg_size = real(0)
       max_size = int(0)
    end if
    
    return
  end subroutine cluster_stat
  !---------------------------------------------------------
  !---average_cluster_dimension: calculates average particle sizes
  !--------------------------------------------------------- 
  subroutine average_cluster_dimension(d0)
    use variables
    implicit none
    
    real,intent(out) :: d0
    real :: di,tot_di,tot_clus_units
    integer :: ilabel

    !call the cluster statistics
    call cluster_stat
    
    d0 = real(0)
    tot_di = real(0)
    tot_clus_units = real(0)
    
    if(max_size.ge.cluster_threshold)then!if max_size is greater than zero
       
       do ilabel = 1,nsites!loop through all the cluster labels
          
          if(csize(ilabel).ge.cluster_threshold)then
             call particle_size(di,ilabel)
             !tot_di = tot_di + di*csize(ilabel)!normalized by cluster units
             tot_di = tot_di + di
             tot_clus_units = tot_clus_units + csize(ilabel)!normalized by equal weights
          end if
          
       end do!loop through all the cluster labels
       
    end if!endif: max_size is greater than zero
    
    if(max_size.ge.cluster_threshold)then
       d0 = tot_di/(real(cluster_num))!normalized by equal weights
       !d0 = tot_di/(real(tot_clus_units))!normalized by cluster units
    else
       d0 = real(0)
    end if
    
    return
  end subroutine average_cluster_dimension
  !---------------------------------------------------------
  !---particle_size(d0) :: calculates the spatial dimensions of the particles
  !--------------------------------------------------------- 
  subroutine particle_size(d0,label)
    use variables
    implicit none
    
    real,intent(out) :: d0
    integer,intent(in) :: label
    real, parameter :: d=real(0.16),a0 = d/sqrt(real(3))
    integer :: molecule,jsite,jlabel,label_size,pjlabel,jocc
    integer :: xcm,ycm,zcm,cmsite,cmtag
    real :: tempd0,jcmdistance
    
    call cluster_stat
        
    label_size = csize(label)
    xcm = 0;ycm = 0;zcm = 0
    d0 = real(0);tempd0 = real(0)
    !first find the center of mass of the cluster with label j
    do molecule = 1,allmolecules
       jsite = noccupy(molecule);jocc = occupancy(jsite);jlabel = clabel(jsite)
       if((jlabel.gt.0).and.((jocc.eq.occSN).or.(jocc.eq.occSI)))then
       
          call classify(jlabel,pjlabel)
          if((pjlabel.eq.label).and.(csize(pjlabel).eq.label_size))then
             xcm = xcm + rx(jsite)
             ycm = ycm + ry(jsite)
             zcm = zcm + rz(jsite)
          end if

       end if
    end do

    xcm = int(real(xcm)/real(label_size))
    ycm = int(real(ycm)/real(label_size))
    zcm = int(real(zcm)/real(label_size))
    
    !calculate the closest site to CM on the lattice
    !cmcheck=mod(xcm,2)+mod(ycm,2)+mod(zcm,2)
    if(mod(xcm,2).eq.0)then
       if(mod(ycm,2).eq.0)then
          !do nothing
       elseif(mod(ycm,2).eq.1)then
          ycm = ycm + 1
       end if
       if(mod(zcm,2).eq.0)then
          !do nothing
       elseif(mod(zcm,2).eq.1)then
          zcm = zcm + 1
       end if
    elseif(mod(xcm,2).eq.1)then
       if(mod(ycm,2).eq.0)then
          ycm = ycm + 1
       elseif(mod(ycm,2).eq.1)then
          !do nothing
       end if
       if(mod(zcm,2).eq.0)then
          zcm = zcm + 1
       elseif(mod(zcm,2).eq.1)then
          !do nothing
       end if
    else
       print*,'BUG: lattice defect :: subroutine particle_size'
       call freeandclose
       stop
    end if
    
    cmtag = (xcm-1)*2*ly*2*lz + (ycm - 1)*2*lz + zcm
    cmsite = tag_site(cmtag)
    
    !now calculate the distance of every site with cm and calculate the
    !radius of gyration
    do molecule = 1,allmolecules
       jsite = noccupy(molecule);jocc = occupancy(jsite);jlabel = clabel(jsite)
       if((jlabel.gt.0).and.((jocc.eq.occSN).or.(jocc.eq.occSI)))then
          
          call classify(jlabel,pjlabel)
          if((pjlabel.eq.label).and.(csize(pjlabel).eq.label_size))then
             call distance(jsite,cmsite,jcmdistance)
             tempd0 = tempd0 + jcmdistance**2
          end if
          
       end if
    end do
    
    d0 = sqrt(tempd0/real(label_size))
    d0 = real(2)*d0*a0
    
    return
  end subroutine particle_size
  !---------------------------------------------------------
  !---CLASSIFY(JLABEL): CALCULATES THE PROPER LABEL OF A SITE
  !--------------------------------------------------------- 
  subroutine classify(jlabel,pjlabel)
    use variables
    implicit none
    
    integer,intent(in) :: jlabel
    integer,intent(out) :: pjlabel
    integer :: t,r
    
    r = jlabel
    t = r
    t = -csize(t)
    
    if(t.lt.0)then
       pjlabel = r
       return
    else
       r=t
       t = -csize(t)
       if(t.lt.0)then
          pjlabel = r
          return
       else
          
          do while(t.gt.0)
             r = t
             t = -csize(t)
          end do

          pjlabel = r
          csize(jlabel) = -r
          return

       end if
    end if
    
    return
  end subroutine classify
  !---------------------------------------------------------
  !---distance: calculates the distance between two sites
  !---------------------------------------------------------   
  subroutine distance(jsite,ksite,jkdistance)
    use variables
    implicit none
    
    integer,intent(in) :: jsite, ksite
    real,intent(out) :: jkdistance
    
    jkdistance = sqrt(real( &
         (rx(jsite) - rx(ksite))**2 + &
         (ry(jsite) - ry(ksite))**2 + &
         (rz(jsite) - rz(ksite))**2 &         
         ))
    
    return
  end subroutine distance
  !---------------------------------------------------------
  !---negchargeperSI: calculates the charges SI in solution
  !--------------------------------------------------------- 
  subroutine negchargeperSI
    use variables
    implicit none
    
    integer :: isite,ivertex,ispin,molecule
    integer :: i,occfn,bond
    integer,dimension(4) :: fn,fnsite
    logical :: dissolved_SI
    
    chargeperSI = 0
    
    do molecule = 1,allmolecules!loop over all molecules
       
       isite = noccupy(molecule)
       ivertex = head(isite)
       ispin = spin(isite)
       
       if(ispin.eq.spinSI)then!if the molecule is TAA
          
          fn(1) = ivertex
          fnsite(1) = n1list(fn(1),isite)
          do bond = 2,4
             fn(bond) = Si_O(bond,fn(1))
             fnsite(bond) = n1list(fn(bond),isite)             
          end do
          
          !check if SI is free or bonded by SN or TAA in the solid phase
          dissolved_SI = .true.
          do bond = 1,4
             occfn = occupancy(fnsite(bond))
             if((occfn.eq.occISIO).or.(occfn.eq.occSIO))then
                !do nothing
             elseif((occfn.eq.occSIOSNO).or.(occfn.eq.occSIOSIO))then
                dissolved_SI = .false.
             end if
          end do
          
          if(dissolved_SI)then
             !do nothing
          else
             chargeperSI = chargeperSI + 1.00
          end if
          
       end if!endif: the molecule is TAA
              
    end do!enddo:loop over all molecules
    
    chargeperSI = chargeperSI/real(ni(1))
    
    return
  end subroutine negchargeperSI
  !---------------------------------------------------------
  !---sublattice_ordering: calculates the sublattice ordering
  !                        in the system
  !--------------------------------------------------------- 
  subroutine sublattice_ordering
    use variables
    implicit none
    
    integer :: molecule,isite,ivertex,ispin
    integer :: i,occfn,bond
    integer,dimension(4) :: fn,fnsite
    integer :: different_connection

    sublat_ordering = 0.00
    
    do molecule = 1,allmolecules!loop over all molecules
       
       isite = noccupy(molecule)
       ivertex = head(isite)
       ispin = spin(isite)
       
       if(ispin.eq.spinSI)then!if the molecule is SI
          
          fn(1) = ivertex
          fnsite(1) = n1list(fn(1),isite)
          do bond = 2,4
             fn(bond) = Si_O(bond,fn(1))
             fnsite(bond) = n1list(fn(bond),isite)             
          end do
          
          !check the no. of SI connected to SN
          different_connection = 0
          do bond = 1,4
             occfn = occupancy(fnsite(bond))
             if(occfn.eq.occSIOSNO) different_connection = different_connection + 1
          end do
          
          if(different_connection.eq.3)then
             sublat_ordering = sublat_ordering + 1.00
          elseif(different_connection.eq.4)then
             print*,'BUG: occupancy of OH and O- :: subroutine sublattice_ordering'
             call freeandclose
             stop
          end if
          
       end if!endif: the molecule is SI
       
    end do!enddo: loop over all molecules
    
    sublat_ordering = sublat_ordering/real(ni(1))
    
    return
  end subroutine sublattice_ordering
  !---------------------------------------------------------
  !---shell_ordering: checks the fraction of SI connected to SI&TAA and only to TAA
  !--------------------------------------------------------- 
  subroutine shell_order
    use variables
    implicit none
    
    integer :: bond,molecule,isite,iocc,ivertex
    integer :: occfn
    integer,dimension(4) :: fn,fnsite
    logical :: connectedtoSN,connectedtoTAA
    
    frac_SITAA = real(0)
    frac_SISNTAA = real(0)
    
    do molecule = 1,allmolecules!loop over all the molecules
       isite = noccupy(molecule)
       iocc = occupancy(isite)
       ivertex = head(isite)
       
       if(iocc.eq.occSI)then!if isite has SI
          
          !find the neighbors
          do bond = 1,4
             fn(bond) = Si_O(bond,ivertex)
             fnsite(bond) = n1list(fn(bond),isite)
          end do
          
          !check for any connection with SN          
          connectedtoSN = .false.
          do bond = 2,4
             occfn = occupancy(fnsite(bond))
             if(occfn.eq.occSIOSNO) connectedtoSN = .true.
          end do

          !check for any connection with TAA          
          connectedtoTAA = .false.
          do bond = 1,6
             occfn = occupancy(n2list(bond,fnsite(1)))
             if(occfn.eq.occTAA) connectedtoTAA = .true.
          end do
          
          if(connectedtoSN.and.connectedtoTAA)then
             frac_SISNTAA = frac_SISNTAA + real(1)
          elseif(connectedtoTAA.and.(connectedtoSN.eqv..false.))then
             frac_SITAA = frac_SITAA + real(1)
          end if
          
       end if!endif: isite has SI       
       
    end do!enddo: loop over all the molecules

    frac_SITAA = frac_SITAA/real(ni(1))
    frac_SISNTAA = frac_SISNTAA/real(ni(1))
    
    return
  end subroutine shell_order
  !---------------------------------------------------------
  !---configuration_output: outputs the config. in vmd.xyz
  !--------------------------------------------------------- 
  subroutine config_out_molecule
    use variables
    implicit none
    
    integer :: isite,ivertex,ispin
    integer :: i,sumoccupy
    
    sumoccupy = ni(1) + ni(2) + ni(3)
    
    write(1006,*)lx, ly, lz
    write(1006,*) sumoccupy
    
    do i=1,size(noccupy)
       isite = noccupy(i)
       ispin = spin(isite)
       ivertex = head(isite)
       write(1006,100) rx(isite),ry(isite),rz(isite),ispin,ivertex
    end do
    
100 format(I4, I4, I4, I4, I4)
 
    return
  end subroutine config_out_molecule
  !---------------------------------------------------------
  !---config_out_atom: output the configuration of atoms
  !--------------------------------------------------------- 
  subroutine config_out_atom
    use variables
    implicit none
    
    integer :: isite, ivertex, ispin, iocc
    integer :: monomer, bond
    integer,dimension(4) :: fn
    integer :: i, j
    integer :: x, y, z, nx, ny, nz
    
    !write(1006,*) natom;flush(1006)
    write(1006,*) 5*ni(1)+5*ni(2)+ni(3);flush(1006)
    write(1006,*) 'index, x, y, z';flush(1006)   
    
    j=0
    do monomer=1,allmolecules
       
       isite = noccupy(monomer)
       ivertex = head(isite)
       ispin = spin(isite)
       iocc = occupancy(isite)
       
       if(ispin.eq.spinSN)then
          !output SN configuration
          fn(1) = ivertex
          do bond = 2,4
             fn(bond) = Si_O(bond,fn(1))
          end do
          x = rx(isite)
          y = ry(isite)
          z = rz(isite)
          j = j + 1
          write(1006,900) j,x,y,z;flush(1006)
          do bond = 1,4
             j = j+1
             call output(iocc,fn(bond),x,y,z,nx,ny,nz)
             write(1006,900) j,nx,ny,nz;flush(1006)
          end do

       elseif(ispin.eq.spinSI)then
          !output SI configuration
          fn(1) = ivertex
          do bond = 2,4
             fn(bond) = Si_O(bond,fn(1))
          end do
          x = rx(isite)
          y = ry(isite)
          z = rz(isite)
          j = j + 1
          write(1006,900) j,x,y,z;flush(1006)
          do bond = 1,4
             j = j+1
             call output(iocc,fn(bond),x,y,z,nx,ny,nz)
             write(1006,900) j,nx,ny,nz;flush(1006)
          end do
       elseif(ispin.eq.spinTAA)then
          !output TAA configuration
          x = rx(isite)
          y = ry(isite)
          z = rz(isite)
          j = j + 1
          write(1006,900) j,x,y,z;flush(1006)
          !j = j+1
          !call output(iocc,4,x,y,z,nx,ny,nz)
          !write(1006,900) j,nx,ny,nz;flush(1006)
          !j = j+1
          !call output(iocc,2,x,y,z,nx,ny,nz)
          !write(1006,900) j,nx,ny,nz;flush(1006)
          !j = j+1
          !call output(iocc,11,x,y,z,nx,ny,nz)
          !write(1006,900) j,nx,ny,nz;flush(1006)
          !j = j+1
          !call output(iocc,9,x,y,z,nx,ny,nz)
          !write(1006,900) j,nx,ny,nz;flush(1006)
       else
          print*,'Error: molecule other than SN/SI/TAA :: subroutine config_out_atom'
          print*,'ispin',ispin
          call freeandclose
          stop          
       end if
       
    end do
    
    j = natom
    
    if(j.ne.natom)then
       print*,'Error: j.ne.allatoms :: subroutine config_out_atom'
       print*,'j:',j
       print*,'allatoms:',natom
       call freeandclose
       stop
    end if
    
900 format(I8, I8, I8, I8)     
    
    return
  end subroutine config_out_atom
  !---------------------------------------------------------
  !---take_snapshot(snapshot)
  !---------------------------------------------------------   
  subroutine take_snapshot
    use variables
    implicit none
    
    integer :: molecule, isite,ispin,ivertex,bond,ihead,monomer
    integer,dimension(4) :: fn,fnsite
    integer :: i, j
    integer :: x, y, z, nx, ny, nz
    character(len=16) :: filenum,prefix,extension,temp1,filename
    integer :: opstatus
    
    extension = '.xyz'
    prefix = 'snapshot_'

    write(filenum,'(I3.3)')snapshot    
    filename = trim('snap')//'_'//trim(filenum)//'.'//trim('cfg')

    open(1016,file=filename,status='unknown')
    
    write(1016,900)step,snapshot,point
    
    !output the configuration of the molecules
    do molecule = 1,allmolecules

       isite = noccupy(molecule)
       ivertex = head(isite)
       ispin = spin(isite)
       
       if((ispin.eq.spinSI).or.(ispin.eq.spinSN))then
          do bond = 1,4
             fn(bond) = Si_O(bond,ivertex)
             fnsite(bond) = n1list(fn(bond),isite)
          end do
          ihead = fnsite(1)
       end if

       write(1016,*)molecule,isite,ivertex,ispin
       if(ispin.eq.spinSN)then
          write(1016,*)
       elseif(ispin.eq.spinSI)then
          write(1016,*)SIlink(ihead)
       elseif(ispin.eq.spinTAA)then
          write(1016,*)TAAlink(isite)
       else
          print*,'BUG: molecule neither SN/SI/TAA: subroutine output_state'
       end if

    end do

200 format(i8,x,i8,x,i8)
300 format(f6.2)
400 format(f10.2,x,f10.2,x,f10.2)
500 format(i8,x,i8,x,i8,x,i8)
600 format(i8,x,i8,x,i8,x,i8)
700 format(f6.2,x,f6.2)
800 format(f6.4,x,f6.4,x,f6.4,x,f6.4)
900 format(i15)
  
    close(1016)
    return
  end subroutine take_snapshot
  !---------------------------------------------------------
  !---output_state: outputs the state of the system in the file system_state.cfg
  !--------------------------------------------------------- 
  subroutine output_state
    use variables
    implicit none
    
    integer :: opstatus
    integer :: molecule, isite, ispin ,ivertex, bond, ihead
    integer,dimension(4) :: fn,fnsite
    
    open(1017,file='system_state.cfg',status='unknown',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'Error: opening system_state.cfg file :: subroutine output_state'
       call freeandclose
       stop
    end if    

    write(1017,900)step,snapshot,point
    
    !output the configuration of the molecules
    do molecule = 1,allmolecules

       isite = noccupy(molecule)
       ivertex = head(isite)
       ispin = spin(isite)
       
       if((ispin.eq.spinSI).or.(ispin.eq.spinSN))then
          do bond = 1,4
             fn(bond) = Si_O(bond,ivertex)
             fnsite(bond) = n1list(fn(bond),isite)
          end do
          ihead = fnsite(1)
       end if

       write(1017,*)molecule,isite,ivertex,ispin
       if(ispin.eq.spinSN)then
          write(1017,*)
       elseif(ispin.eq.spinSI)then
          write(1017,*)SIlink(ihead)
       elseif(ispin.eq.spinTAA)then
          write(1017,*)TAAlink(isite)
       else
          print*,'BUG: molecule neither SN/SI/TAA: subroutine output_state'
       end if

    end do

200 format(i8,x,i8,x,i8)
300 format(f6.2)
400 format(f10.2,x,f10.2,x,f10.2)
500 format(i8,x,i8,x,i8,x,i8)
600 format(i8,x,i8,x,i8,x,i8)
700 format(f6.2,x,f6.2)
800 format(f6.4,x,f6.4,x,f6.4,x,f6.4)
900 format(i15)

    close(1017)
    return
  end subroutine output_state
  !---------------------------------------------------------
  !---read_state: outputs the state of the system in the file system_state.cfg
  !--------------------------------------------------------- 
  subroutine read_state
    use variables
    use input_SiO4_TAA
    implicit none
    
    integer :: opstatus, data, bond, ihead
    integer :: molecule, isite, ispin ,ivertex
    integer,dimension(4) :: fn,fnsite
    logical :: can_insert

    !setup the initial lattice
    occupancy = occW
    spin = spinW
    head=0
    noccupy = 0
    resite = 0
    SIlink = -1;TAAlink = -1
    
    open(1017,file='system_state.cfg',status='unknown',iostat=opstatus)
    if(opstatus.ne.0)then
       print*,'Error: opening system_state.cfg file :: subroutine output_state'
       call freeandclose
       stop
    end if
    
    read(1017,900)initial_step,snapshot,point
    
    !output the configuration of the molecules
    do data = 1,allmolecules

       read(1017,*)molecule,isite,ivertex,ispin       
       if((ispin.eq.spinSI).or.(ispin.eq.spinSN))then
          do bond = 1,4
             fn(bond) = Si_O(bond,ivertex)
             fnsite(bond) = n1list(fn(bond),isite)
          end do
          ihead = fnsite(1)
       end if

       if(ispin.eq.spinSN)then
          read(1017,*)
       elseif(ispin.eq.spinSI)then
          read(1017,*)SIlink(ihead)
       elseif(ispin.eq.spinTAA)then
          read(1017,*)TAAlink(isite)
       else
          print*,'BUG: molecule neither SN/SI/TAA: subroutine output_state'
       end if

       if(data.ne.molecule)then
          print*,'BUG: data.ne.molecule :: subroutine read_output'
          print*,'data/molecule',data/molecule
          call freeandclose
          stop
       end if
              
       if(ispin.eq.spinSN)then!if the molecule is SN
          call can_insert_SN(isite,ivertex,can_insert)
          if(can_insert)then
             call insert_site(isite,ivertex,ispin,molecule)
          else
             print*,'BUG: Reading SN from system_state.cfg :: subroutine read_state'
             call freeandclose
             stop
          end if
       elseif(ispin.eq.spinSI)then!if the molecule is SI
          call can_insert_SI(isite,ivertex,can_insert)
          if(can_insert)then
             call insert_site(isite,ivertex,ispin,molecule)
          else
             print*,'BUG: Reading SI from system_state.cfg :: subroutine read_state'
             call freeandclose
             stop
          end if
       elseif(ispin.eq.spinTAA)then!if the molecule is TAA
          call can_insert_TAA(isite,can_insert)
          if(can_insert)then
             call insert_site(isite,ivertex,ispin,molecule)
          else
             print*,'BUG: Reading TAA from system_state.cfg :: subroutine read_state'
             call freeandclose
             stop
          end if
       end if!endif: the molecule is SN
       
    end do
    
    !output the psf file
    call psf_out

    !output initial configuration
    call config_out_atom

200 format(i8,x,i8,x,i8)
300 format(f6.2)
400 format(f10.2,x,f10.2,x,f10.2)
500 format(i8,x,i8,x,i8,x,i8)
600 format(i8,x,i8,x,i8,x,i8)
700 format(f6.2,x,f6.2)
800 format(f6.4,x,f6.4,x,f6.4,x,f6.4)
900 format(i15)

    close(1017)    
    return
  end subroutine read_state
  !---------------------------------------------------------
  !---output(ivertex,x,y,z,nx,ny,nz)
  !--------------------------------------------------------- 
  subroutine output(iocc,ivertex, x, y, z, nx, ny, nz)
    use variables
    implicit none
    
    integer :: iocc,ivertex,x,y,z,nx,ny,nz
    
    if((iocc.eq.occSN).or.(iocc.eq.occSI))then!if isite has SN/SI
    
       if(ivertex.eq.1)then
          nx = x - 1
          ny = y - 1
          nz = z + 1
       elseif(ivertex.eq.2)then
          nx = x + 1
          ny = y - 1
          nz = z + 1
       elseif(ivertex.eq.3)then
          nx = x + 1
          ny = y + 1
          nz = z + 1
       elseif(ivertex.eq.4)then
          nx = x - 1
          ny = y + 1
          nz = z + 1
       elseif(ivertex.eq.5)then
          nx = x - 1
          ny = y - 1
          nz = z - 1
       elseif(ivertex.eq.6)then
          nx = x + 1
          ny = y - 1
          nz = z - 1
       elseif(ivertex.eq.7)then
          nx = x + 1
          ny = y + 1
          nz = z - 1
       elseif(ivertex.eq.8)then
          nx = x - 1
          ny = y + 1
          nz = z - 1
       end if

    elseif(iocc.eq.occTAA)then!if isite has TAA
       
       if(ivertex.eq.4)then
          nx = x - 1
          ny = y - 1
          nz = z 
       elseif(ivertex.eq.2)then
          nx = x + 1
          ny = y - 1
          nz = z 
       elseif(ivertex.eq.11)then
          nx = x + 1
          ny = y + 1
          nz = z
       elseif(ivertex.eq.9)then
          nx = x - 1
          ny = y + 1
          nz = z
       end if
       
    else
       print*,'Error:: ivertex must be between 1 and 8:: subroutine output'
       call freeandclose
       stop
       
    end if!if isite has SN/SI
    
    return
  end subroutine output
  !---------------------------------------------------------
  !---LCG RANDOM NUMBER GENERATOR
  !--------------------------------------------------------- 
  function ran3(idum)
    implicit none
    integer idum
    integer mbig,mseed,mz
    !real*8 mbig,mseed,mz
    real*8 ran3,fac
    parameter(mbig=1000000000,mseed=161803398,mz=0,fac=1./mbig)  
    !parameter (mbig=4000000.,mseed=1618033.,mz=0.0,fac=1./mbig)
    integer i,iff,ii,inext,inextp,k
    integer mj,mk,ma(55)
    !real*8 mj,mk,ma(55)
    save iff,inext,inextp,ma
    data iff  /0/
    if(idum.lt.0.or.iff.eq.0) then
       iff=1
       mj=abs(mseed-abs(idum))
       mj=mod(mj,mbig)
       ma(55)=mj
       mk=1
       do i = 1,54
          ii=mod(21*i,55)
          ma(ii) = mk
          mk = mj-mk
          if(mk.lt.mz) mk = mk + mbig
          mj = ma(ii)
       enddo
       do  k=1,4
          do  i = 1,55
             ma(i) = ma(i)-ma(1+mod(i+30,55))
             if(ma(i).lt.mz)ma(i) = ma(i) + mbig
          enddo
       enddo
       inext = 0
       inextp = 31
       idum = 1
    endif
    inext = inext + 1
    if(inext.eq.56) inext = 1
    inextp = inextp + 1
    if(inextp.eq.56) inextp = 1
    mj = ma(inext) - ma(inextp)
    if(mj.lt.mz)mj = mj + mbig
    ma(inext) = mj
    ran3 = mj*fac
    return
  end function ran3
  !---------------------------------------------------------
  !--- FREE_ALL: DEALLOCATE ALL MEMORY
  !---------------------------------------------------------
  subroutine freeandclose
    use variables
    implicit none

    !deallocate all arrays
    deallocate(occupancy)
    deallocate(tag_site)
    
    deallocate(noccupy)
    deallocate(head)
    deallocate(spin)
    deallocate(resite)
    deallocate(clabel)
    deallocate(csize)
    
    deallocate(n1list)
    deallocate(n2list)
    deallocate(n3list)
    deallocate(n4list)
    deallocate(n5list)

    deallocate(rx)
    deallocate(ry)
    deallocate(rz)
    
    deallocate(SIlink)
    deallocate(TAAlink)
    
    !close all opened files
    close(1001);close(1002);close(1003);close(1004)
    close(1006);close(1007);close(1008);close(1009)
    !close(1010);close(1011)
    close(1012);close(1013);close(1014);close(1015)
    close(1018)

    return
  end subroutine freeandclose
  !***************************************************************************
  !***SAEFTY CHECKING
  !***************************************************************************
  subroutine safetycheck(rings3,rings4,tot_en)
    use variables
    implicit none
    
    integer,intent(in) :: rings3,rings4
    real*8,intent(in) :: tot_en
    integer :: rings3_temp,rings4_temp,occen,occfn,neighs,overlap,bridge,snbond
    integer :: i,j,molecule,bond,isite,ispin,ivertex,ihead,iocc,jsite,sn_neigh,jocc
    integer :: connectedSI=0,connectedTAA=0
    real*8 :: temp_en,ratio,epsilon = real(0.001),old_en
    integer,dimension(4) :: fn,fnsite,sn,snsite
    
    do i=1,nsites
       occen = occupancy(i)
       if((occen.eq.occW).or.(occen.eq.occSN).or.(occen.eq.occSI).or.(occen.eq.occSNOSNO)&
            .or.(occen.eq.occISIO).or.(occen.eq.occSIOSNO).or.(occen.eq.occSNO).or.&
            (occen.eq.occSIO).or.(occen.eq.occTAA))then
       else
          print*,'BUG : site has invalid occupancy :: subroutine safetycheck'
          print*,i,rx(i),ry(i),rz(i),occen
          call freeandclose
          stop
       end if
    end do
    
    !check the proper connection nos. between SI and TAA
    do i = 1,nsites
       if(SIlink(i).gt.0)then
          if(occupancy(i).eq.occISIO)then
             connectedSI = connectedSI + 1
             isite = SIlink(i)
             if((TAAlink(isite).ne.i).or.(occupancy(TAAlink(isite)).ne.occISIO))then
                print*,'Bug: TAAlink(isite).ne.i: cubroutine safetycheck'
                call freeandclose
                stop
             end if
          else
             print*,'BUG: occupancy(i).eq.occISIO: subroutine safetycheck'
             call freeandclose
             stop
          end if
       end if
       
       if(TAAlink(i).gt.0)then
          if(occupancy(i).eq.occTAA)then
             connectedTAA = connectedTAA + 1
             isite = TAAlink(i)
             if((SIlink(isite).ne.i).or.(occupancy(SIlink(isite)).ne.occTAA))then
                print*,'Bug: SIlink(isite).ne.i: cubroutine safetycheck'
                call freeandclose
                stop
             end if
          else
             print*,'BUG: occupancy(i).eq.occTAA: subroutine safetycheck'
             call freeandclose
             stop
          end if
       end if
    end do

    if(connectedSI.eq.connectedTAA)then
       !do nothing
    else
       print*,'Bug: connectedSI and connectedTAA are not same: cubroutine safetycheck'
       print*,'connectedSI/connectedTAA',connectedSI,connectedTAA
       call freeandclose
       stop
    end if
    
    call total_energy(temp_en)
    call rings3and4(rings3_temp,rings4_temp)
    !ratio = abs(1-(temp_en/tot_en))
    !ratio = abs(tot_en-temp_en)/tot_en
    ratio = abs(real(tot_en/temp_en)-real(1))
    if((ratio.gt.epsilon).or.(rings3_temp.ne.rings3).or.(rings4_temp.ne.rings4))then
       print*,'BUG :: energy/rings3/rings4 calculation :: subroutine safetycheck'
       print*,'Energy from subroutine:',temp_en
       print*,'Energy from simulation:',tot_en
       print*,'ratio/epsilon:',ratio,epsilon
       print*,'rings3/rings4 from subroutine:',rings3_temp,rings4_temp
       print*,'rings3/rings4 from simulation:',rings3,rings4
       print*,
       call freeandclose
       stop
    end if
    
    do molecule=1,allmolecules!loop over all molecules

       isite = noccupy(molecule)
       ivertex = head(isite)
       ispin = spin(isite)
       iocc = occupancy(isite)
       
       if(iocc.eq.occSN)then
          
          overlap = 0
          if(occupancy(isite).eq.occSN)then!if: isite has TAA
             
             !check first,second and third neighbors for TAA occupancy
             do i=1,8!checking first neighbors
                occfn = occupancy(n1list(i,isite))
                if(occfn.eq.occTAA)then
                   print*,'BUG : incorrect occupancy(n1list(i,isite) 1st (SN): subroutine safetycheck'
                   call freeandclose
                   stop
                end if
             end do!endo: checking first neighbors
             !do i = 1,6!checking first neighbors
             !   occfn = occupancy(n2list(i,isite))
             !   if(occfn.eq.occTAA)then
             !      print*,'BUG : incorrect occupancy(n1list(i,isite) 2nd (SN): subroutine safetycheck'
             !      call freeandclose
             !      stop
             !   end if
             !end do!enddo: checking first neighbors
             !do i = 1,12!checking third neighbors
             !   occfn = occupancy(n3list(i,isite))
             !   if(occfn.eq.occTAA)then
             !      print*,'BUG : incorrect occupancy(n1list(i,isite) 3rd (SN): subroutine safetycheck'
             !      call freeandclose
             !      stop
             !   end if
             !end do!enddo: checking third neighbors
             
             fn(1) = ivertex
             fnsite(1) = n1list(fn(1),isite)
             if((occupancy(fnsite(1)).eq.occSNO).or.((occupancy(fnsite(1)).eq.occSIO)))&
                  overlap = overlap + 1
             do bond = 2,4
                fn(bond) = Si_O(bond,fn(1))
                fnsite(bond) = n1list(fn(bond),isite)
                if((occupancy(fnsite(bond)).eq.occSNO).or.((occupancy(fnsite(bond)).eq.occSIO)))&
                     overlap = overlap + 1
             end do
             
             !check if there are any water,SN,SNO or SIO molecules on the fnsite() sites
             do bond = 1,4
                occfn = occupancy(fnsite(bond))
                if((occfn.eq.occSNO).or.(occfn.eq.occSIOSNO).or.(occfn.eq.occSNOSNO))then
                   !do nothing
                else
                   print*,'BUG : incorrect OH occupancy (SN): subroutine safetycheck'
                   call freeandclose
                   stop
                end if
             end do
             
             do bond = 1,4!check first,second and third neighbors of OH for TAA occupancy
                
                do i=1,8!checking first neighbors
                   occfn = occupancy(n1list(i,fnsite(bond)))
                   if(occfn.eq.occTAA)then
                      print*,'BUG : TAA in 1st neighbor of fnsite(bond) (SN): subroutine safetycheck'
                      call freeandclose
                      stop
                   end if
                end do!endo: checking first neighbors
                !do i = 1,6!checking first neighbors
                !   occfn = occupancy(n2list(i,fnsite(bond)))
                !   if(occfn.eq.occTAA)then
                !      print*,'BUG : TAA in 2nd neighbor of fnsite(bond) (SN): subroutine safetycheck'
                !      call freeandclose
                !      stop
                !   end if
                !end do!enddo: checking first neighbors
                !do i = 1,12!checking third neighbors
                !   occfn = occupancy(n3list(i,fnsite(bond)))
                !   if(occfn.eq.occTAA)then
                !      print*,'BUG : TAA in 3rd neighbor of fnsite(bond) (SN): subroutine safetycheck'
                !      call freeandclose
                !      stop
                !   end if
                !end do!enddo: checking third neighbors
                
             end do!enddo: check first,second and third neighbors of OH for TAA occupancy
             
             !-----CHECK FOR 2 MEMBERED RING FORMATION
             if(overlap.ge.2)then!if overlap>=2 chance of 2 membered ring
                
                do sn_neigh = 1,6!loop over 2 neighbors of isite
                   
                   jsite = n2list(sn_neigh,isite)
                   !if there SN on neighboring site
                   if((occupancy(jsite).eq.occSN).or.(occupancy(jsite).eq.occSI))then
                      
                      !find the bonded OH of jsite
                      sn(1) = head(jsite)
                      snsite(1) = n1list(sn(1),jsite)
                      do snbond = 2,4
                         sn(snbond) = Si_O(snbond,sn(1))
                         snsite(snbond) = n1list(sn(snbond),jsite)
                      end do
                      
                      !check for any bridging oxygens
                      bridge=0
                      do i=1,4
                         do j=1,4
                            if((fnsite(i).eq.snsite(j)).and.((occupancy(fnsite(i)).eq.occSNO).or.&
                                 (occupancy(fnsite(i)).eq.occSIO)))then 
                               bridge = bridge + 1
                            end if
                         end do
                      end do
                      
                      !if there are more than 1 bridging oxygens
                      if(bridge.gt.1) then
                         print*,'BUG : Two ring formation (SN): subroutine safetycheck'
                         call freeandclose
                         stop
                      end if
                      
                   end if!if there SN on neighboring site
                   
                end do!enddo:loop over 2 neighbors of isite
                
             end if!if overlap.ge.2 chance of 2 membered ring
             !-----
             
          else
             
             print*,'BUG: incorrect occupancy on isite (SN): subroutine safetycheck'
             call freeandclose
             stop             
             
          end if!if isite has SN on it
                    
       elseif(iocc.eq.occSI)then!if: isite has SI
          
          overlap = 0
          !check if there's water on the site
          if(occupancy(isite).eq.occSI)then!if isite has SI on it
             
             !check first,second and third neighbors for TAA occupancy
             do i=1,8!checking first neighbors
                occfn = occupancy(n1list(i,isite))
                if(occfn.eq.occTAA)then
                   print*,'BUG : incorrect occupancy(n1list(i,isite) 1st (SI): subroutine safetycheck'
                   call freeandclose
                   stop
                end if
             end do!endo: checking first neighbors
             !do i = 1,6!checking first neighbors
             !   occfn = occupancy(n2list(i,isite))
             !   if(occfn.eq.occTAA)then
             !      print*,'BUG : incorrect occupancy(n1list(i,isite) 2nd (SI): subroutine safetycheck'
             !      call freeandclose
             !      stop
             !   end if
             !end do!enddo: checking first neighbors
             !do i = 1,12!checking third neighbors
             !   occfn = occupancy(n3list(i,isite))
             !   if(occfn.eq.occTAA)then
             !      print*,'BUG : incorrect occupancy(n1list(i,isite) 3rd (SI): subroutine safetycheck'
             !      call freeandclose
             !      stop
             !   end if
             !end do!enddo: checking third neighbors
             
             fn(1) = ivertex
             fnsite(1) = n1list(fn(1),isite)
             ihead = fnsite(1)
             if(occupancy(fnsite(1)).eq.occISIO)then
             else
                print*,'BUG : incorrect occupancy at O- (SI): subroutine safetycheck'
                call freeandclose
                stop
             end if
             do bond = 2,4
                fn(bond) = Si_O(bond,fn(1))
                fnsite(bond) = n1list(fn(bond),isite)
                if(occupancy(fnsite(bond)).eq.occSNO)&
                     overlap = overlap + 1
             end do
             
             do bond = 1,4!check first,second and third neighbors of OH for TAA occupancy
                
                do i=1,8!checking first neighbors
                   occfn = occupancy(n1list(i,fnsite(bond)))
                   if(occfn.eq.occTAA)then
                      print*,'BUG : TAA in 1st neighbor of fnsite(bond) (SI): subroutine safetycheck'
                      call freeandclose
                      stop
                   end if
                end do!endo: checking first neighbors
                !do i = 1,6!checking first neighbors
                !   occfn = occupancy(n2list(i,fnsite(bond)))
                !   if(occfn.eq.occTAA)then
                !      print*,'BUG : TAA in 2nd neighbor of fnsite(bond) (SI): subroutine safetycheck'
                !      call freeandclose
                !      stop
                !   end if
                !end do!enddo: checking first neighbors
                !do i = 1,12!checking third neighbors
                !   occfn = occupancy(n3list(i,fnsite(bond)))
                !   if(occfn.eq.occTAA)then
                !      print*,'BUG : TAA in 3rd neighbor of fnsite(bond) (SI): subroutine safetycheck'
                !      call freeandclose
                !      stop
                !   end if
                !end do!enddo: checking third neighbors
                
             end do!enddo: check first,second and third neighbors of OH for TAA occupancy
             
             !check if there are any water,SN,SI,SNO or SIO molecules on the 
             !fnsite() sites
             !check for the location of the pointer variable
             
             if(occupancy(fnsite(1)).eq.occISIO)then!if the pointer has water on its site
                do bond = 2,4
                   jocc = occupancy(fnsite(bond))
                   if((jocc.eq.occSIO).or.(jocc.eq.occSIOSNO))then
                      !do nothing
                   else
                      print*,'BUG : incorrect OH occupancy (SI): subroutine safetycheck'
                      call freeandclose
                      stop
                   end if
                end do
                !-----CHECK FOR 2 MEMBERED RING FORMATION
                if(overlap.ge.2)then!if overlap.ge.2 chance of 2 membered ring
                   
                   do sn_neigh = 1,6!loop over 2 neighbors of isite
                      
                      jsite = n2list(sn_neigh,isite)
                      !if there SN on neighboring site
                      if(occupancy(jsite).eq.occSN.or.(occupancy(jsite).eq.occSI))then
                         
                         !find the bonded OH of jsite
                         sn(1) = head(jsite)
                         snsite(1) = n1list(sn(1),jsite)
                         do snbond = 2,4
                            sn(snbond) = Si_O(snbond,sn(1))
                            snsite(snbond) = n1list(sn(snbond),jsite)
                         end do
                         
                         !check for any bridging oxygens
                         bridge=0
                         do i=1,4
                            do j=1,4
                               if((fnsite(i).eq.snsite(j)).and.(occupancy(fnsite(i)).eq.occSNO))then 
                                  bridge = bridge + 1
                               end if
                            end do
                         end do
                         
                         !if there are more than 1 bridging oxygens
                         if(bridge.gt.1) then
                            print*,'BUG : Two ring formation (SI): subroutine safetycheck'
                            call freeandclose
                            stop
                         end if
                         
                      end if!if there SN on neighboring site
                      
                   end do!enddo:loop over 2 neighbors of isite
                   
                end if!if overlap.ge.2 chance of 2 membered ring
                !-----
             else
                print*,'BUG : incorrect O- occupancy (SI): subroutine safetycheck'
                call freeandclose
                stop
             end if!endif: the pointer has water on its site
             
          else
             
             print*,'BUG: incorrect occupancy on isite (SI): subroutine safetycheck'
             call freeandclose
             stop
             
          end if!endif: isite has SI on it
          
       elseif(iocc.eq.occTAA)then!if: isite has TAA
          
          !check the occupancy of isite
          if(occupancy(isite).eq.occTAA)then!if there's water on isite
                          
             !check first,second and third neighbors for TAA occupancy
             do i=1,8!checking first neighbors
                occfn = occupancy(n1list(i,isite))
                if(occfn.ne.occW)then
                   print*,'BUG: 1st neighbor does not have water on it (TAA): subroutine safetycheck'
                   call freeandclose
                   stop
                end if
             end do!endo: checking first neighbors
             !do i = 1,6!checking first neighbors
             !   occfn = occupancy(n2list(i,isite))
             !   if(occfn.ne.occW)then
             !      print*,'BUG: 2nd neighbor does not have water on it (TAA): subroutine safetycheck'
             !      call freeandclose
             !      stop
             !   end if
             !end do!enddo: checking first neighbors
             !do i = 1,12!checking third neighbors
             !   occfn = occupancy(n3list(i,isite))
             !   if(occfn.ne.occW)then
             !      print*,'BUG: 3rd neighbor does not have water on it (TAA): subroutine safetycheck'
             !      call freeandclose
             !      stop
             !   end if
             !end do!enddo: checking third neighbors
             
          else
             
             print*,'BUG: isite does not have TAA (TAA): subroutine safetycheck'
             call freeandclose
             stop
             
          end if!endif there's water on isite  
          
       else
          
          print*,'BUG: molecule on the list is not SI/SN/TAA: subroutine safetycheck'
          call freeandclose
          stop
          
       end if!endif: isite has SN
       
    end do!enddo: loop over all molecules
    
  end subroutine safetycheck
  !---------------------------------------------------------------------------
  subroutine safetycheck_translate(tot_en,ispin,AR1,rand,old_en,new_en,rings3,rings4)
    use variables
    implicit none
    
    integer :: ispin,rings3,rings4,rings3_temp,rings4_temp
    real*8 :: tot_en,old_en,new_en,temp_en,delta_en,AR1,rand,ratio,epsilon = 0.00000001d0
    
    call total_energy(temp_en)
    call rings3and4(rings3_temp,rings4_temp)
    ratio = abs(1-(temp_en/tot_en))
    if((ratio.gt.epsilon).or.(rings3_temp.ne.rings3).or.(rings4_temp.ne.rings4))then
       print*,'BUG :: energy calculation :: subroutine safetycheck_translate'
       print*,'ispin:',ispin
       print*,'AR1:',AR1
       print*,'ran3(seed)',rand
       print*,'Energy from subroutine:',temp_en
       print*,'Energy from simulation:',tot_en
       print*,'old_en',old_en
       print*,'new_en:',new_en
       print*,'delta_en',new_en-old_en
       print*,'rings3/rings4 from subroutine:',rings3_temp,rings4_temp
       print*,'rings3/rings4 from simulation:',rings3,rings4
       print*,
       call freeandclose
       stop
    end if
    
    return
  end subroutine safetycheck_translate
  !---------------------------------------------------------------------------
  subroutine safetycheck_swap(tot_en,isite,jsite,AR1,rand,old_en,new_en,rings3,rings4)
    use variables
    implicit none
    
    integer :: isite,rings3,rings4,rings3_temp,rings4_temp,jsite
    real*8 :: tot_en,old_en,new_en,temp_en,delta_en,AR1,rand,ratio,epsilon = 0.00000001d0
    
    integer :: bond,ihead
    
    call total_energy(temp_en)
    call rings3and4(rings3_temp,rings4_temp)
    ratio = abs(1-(temp_en/tot_en))
    if((ratio.gt.epsilon).or.(rings3_temp.ne.rings3).or.(rings4_temp.ne.rings4))then
       print*,'BUG :: energy calculation :: subroutine safetycheck_swap'
       print*,'ispin/jspin:',spin(isite),spin(jsite)
       print*,'AR1:',AR1
       print*,'ran3(seed):',rand
       print*,'Energy from subroutine:',temp_en
       print*,'Energy from simulation:',tot_en
       print*,'old_en',old_en
       print*,'new_en',new_en
       print*,'delta_en',new_en-old_en
       print*,'rings3/rings4 from subroutine:',rings3_temp,rings4_temp
       print*,'rings3/rings4 from simulation:',rings3,rings4
       print*,
       
       !neighbor information
       print*,'---imolecule information---'
       print*,'ivertex/jocc/spin',head(isite),occupancy(isite),spin(isite)
       print*,'rx(isite),ry(isite),rz(isite)',rx(isite),ry(isite),rz(isite)
       print*,'first neighbors'
       do bond = 1,8
          print*,bond,occupancy(n1list(bond,isite)),rx(n1list(bond,isite)),&
               ry(n1list(bond,isite)),rz(n1list(bond,isite))
       end do
       print*,'second neighbors'
       do bond = 1,6
          print*,bond,occupancy(n2list(bond,isite)),rx(n2list(bond,isite)),&
               ry(n2list(bond,isite)),rz(n2list(bond,isite))
       end do
       
       !neighbor information
       print*,'---jmolecule information---'
       print*,'jvertex/jocc/jspin',head(jsite),occupancy(jsite),spin(jsite)
       print*,'rx(jsite),ry(jsite),rz(jsite)',rx(jsite),ry(jsite),rz(jsite)
       print*,'first neighbors'
       do bond = 1,8
          print*,bond,occupancy(n1list(bond,jsite)),rx(n1list(bond,jsite)),&
               ry(n1list(bond,jsite)),rz(n1list(bond,jsite))
       end do
       print*,'second neighbors'
       do bond = 1,6
          print*,bond,occupancy(n2list(bond,jsite)),rx(n2list(bond,jsite)),&
               ry(n2list(bond,jsite)),rz(n2list(bond,jsite))
       end do
       
       call freeandclose
       stop
    end if
    
    return
  end subroutine safetycheck_swap
  !---------------------------------------------------------------------------
  subroutine safetycheck_rotate(tot_en,ksite,kspin,kvertex,kvertex_new,AR2,rand,old_en,new_en,delta_en,rings3,rings4)
    use variables
    implicit none
    
    integer :: ksite,kspin,rings3,rings4,rings3_temp,rings4_temp,kvertex,kvertex_new
    real*8 :: tot_en,old_en,new_en,temp_en,delta_en,AR2,rand,ratio,epsilon = real(0.001),en_site_old,en_site_new
    integer :: bond,khead,isite,jsite
    integer :: rings3_site_old,rings4_site_old,rings3_site_new,rings4_site_new
    integer :: SIconnected,TAAconnected
    integer,dimension(4) :: fn,fnsite
    
    call total_energy(temp_en)
    call rings3and4(rings3_temp,rings4_temp)
    ratio = abs(1-(temp_en/tot_en))

    call safetycheck_linkage(SIconnected,TAAconnected)    
    if(SIconnected.eq.TAAconnected)then
       !do nothing
    else
       print*,'SIconnected,TAAconnected',SIconnected,TAAconnected
       call freeandclose
       stop
    end if
    
    if((ratio.gt.epsilon).or.(rings3_temp.ne.rings3).or.(rings4_temp.ne.rings4))then!if there's a bug

       call site_energy(ksite,kvertex,kspin,rings3_site_old,rings4_site_old,en_site_old)
       call site_energy(ksite,kvertex_new,kspin,rings3_site_new,rings4_site_new,en_site_new)
       
       print*,'BUG :: energy calculation :: subroutine safetycheck_rotate'
       print*,'ksite,kspin:',ksite,kspin
       print*,'kvertex_old,kvertex_new,head(ksite)',kvertex,kvertex_new,head(ksite)
       print*,'SIconnected,TAAconnected',SIconnected,TAAconnected
       print*,'Old:rings3_site,rings4_site,en_site',rings3_site_old,rings4_site_old,en_site_old
       print*,'New:rings3_site,rings4_site,en_site',rings3_site_new,rings4_site_new,en_site_new
       print*,'ratio,epsilon',ratio,epsilon
       print*,'AR2,ran3(seed)',AR2,ran3(seed)
       print*,'Energy from subroutine:',temp_en
       print*,'Energy from simulation:',tot_en
       print*,'old_en,new_en,delta_en',old_en,new_en,delta_en
       print*,'rings3/rings4 from subroutine:',rings3_temp,rings4_temp
       print*,'rings3/rings4 from simulation:',rings3,rings4
       print*,
       print*,'!----------------------------------------'
       print*,'!---INFORMATION BEFORE THE ROTATION MOVE'
       print*,'!----------------------------------------'
       khead = n1list(kvertex,ksite)
       print*,'khead,Silink(khead)',khead,SIlink(khead)
       do bond = 1,6
          isite = n2list(bond,khead)
          print*,'isite,occupancy,TAAlink,SIlink',isite,occupancy(isite),TAAlink(isite),SIlink(isite)
       end do

       print*
       print*,'!----------------------------------------'
       print*,'!---INFORMATION AFTER THE ROTATION MOVE'
       print*,'!----------------------------------------'
       khead = n1list(kvertex_new,ksite)
       print*,'khead,Silink(khead)',khead,SIlink(khead)
       do bond = 1,6
          isite = n2list(bond,khead)
          print*,'isite,occupancy,TAAlink,SIlink',isite,occupancy(isite),TAAlink(isite),SIlink(isite)
       end do
       

       call freeandclose
       stop
    end if!endif: there's a bug
    
    return
  end subroutine safetycheck_rotate
  !---------------------------------------------------------
  !---manual_check: CHECK FOR CORRECTNESS MANUALLY
  !--------------------------------------------------------- 
  subroutine safetycheck_linkage(SIconnected,TAAconnected)
    use variables
    implicit none

    integer,intent(out) :: SIconnected, TAAconnected    
    integer :: imolecule, isite, ivertex, ihead, ispin
    integer :: jsite, ksite, SISNBO,SNSNBO, bond, occfn
    integer,dimension(4) :: fn,fnsite
    real :: en
    
    SIconnected = 0; TAAconnected = 0
    SNSNBO = 0; SISNBO = 0
    
    !convention
    !SI(ihead) => TAA(jsite) and TAA(jsite) => SI(ksite)
    !for proper connection => ksite = ihead 
    !TAA(isite) => SI(jsite) and SI(jsite) => TAA(ksite)
    !for proper connection => ksite = isite
    
    do imolecule = 1,allmolecules!loop over all molecules
       isite = noccupy(imolecule)
       ispin = spin(isite)
       ivertex = head(isite)
       
       if(ispin.eq.spinSN)then!if molecule is SN
          
          fn(1) = ivertex
          fnsite(1) = n1list(fn(1),isite)
          do bond = 2,4
             fn(bond) = Si_O(bond,ivertex)
             fnsite(bond) = n1list(fn(bond),isite)
          end do
          
          do bond = 1,4
             occfn = occupancy(fnsite(bond))
             
             if(occfn.eq.occSNOSNO)then
                SNSNBO = SNSNBO + 1
             elseif(occfn.eq.occSIOSNO)then
                SISNBO = SISNBO + 1
             end if
          end do         
          
       elseif(ispin.eq.spinSI)then!if molecule is SI
          ihead = n1list(ivertex,isite)
          jsite = SIlink(ihead)
          
          if(jsite.gt.0)then
             !SIconnected = SIconnected + 1
             ksite = TAAlink(jsite)
             if(ksite.eq.ihead)then
                SIconnected = SIconnected + 1
             else
                print*,'BUG :: incorrect assignment of SI-TAA connection :: subroutine safetycheck_linkage(SI)'
                print*,'imolecule,ispin,ivertex,isite,ihead',imolecule,ispin,ivertex,isite,ihead
                call freeandclose
                stop
             end if
          end if
          
       elseif(ispin.eq.spinTAA)then!if molecule is TAA
          jsite = TAAlink(isite)
          if(jsite.gt.0)then
             !TAAconnected = TAAconnected + 1
             
             ksite = SIlink(jsite)
             if(ksite.eq.isite)then
                TAAconnected = TAAconnected + 1
             else
                print*,'BUG :: incorrect assignment of SI-TAA connection :: subroutine safetycheck_linkage(TAA)'
                print*,'imolecule,ispin,ivertex,isite,ihead',imolecule,ispin,ivertex,isite,ihead
                call freeandclose
                stop                
             end if
          end if
          
       end if!endif molecule is SN
    end do!enddo: loop over all molecules
    
    if(TAAconnected.eq.SIconnected)then
       !do nothing
    else
       print*,'BUG: incorrect SI-TAA connections :: subroutine safetycheck_linkage'
       print*,'SIconnected,TAAconnected',SIconnected,TAAconnected
       call freeandclose
       stop
    end if
    
    en = real(-1)*(eISITAA*SIconnected + real(0.5)*eSNSN*SNSNBO + eSISN*SISNBO)
    write(100000,*)'SIconnected,TAAconnected',SIconnected,TAAconnected
    
    return
  end subroutine safetycheck_linkage
  !---------------------------------------------------------
  !---manual_check: CHECK FOR CORRECTNESS MANUALLY
  !--------------------------------------------------------- 
  subroutine manual_check
    use variables
    implicit none
    integer,dimension(13) :: s
    integer,dimension(18) :: o
    integer,dimension(17) :: d
    integer :: ivertex,ispin,monomer,tag
    logical :: arethereSi,can_insert
    integer :: isite,jsite,rings3,rings4,i,j,overlap
    integer,dimension(32) :: neighSi1
    integer,dimension(:),allocatable :: neighSi
    real*8 :: en,en_temp
    !
    s(1) = tag_site(73326)
    s(2) = tag_site(73806)
    s(3) = tag_site(101886)
    s(4) = tag_site(102366)
    s(5) = tag_site(102846)
    s(6) = tag_site(130446)
    s(7) = tag_site(130926)
    s(8) = tag_site(131406)
    s(9) = tag_site(159006)
    s(10) = tag_site(159966)
    s(11) = tag_site(101648)
    s(12) = tag_site(158768)
    s(13) = tag_site(187566)
    !
    o(1) = tag_site(58807)
    o(2) = tag_site(59045)
    o(3) = tag_site(59287)
    o(4) = tag_site(59525)
    o(5) = tag_site(88565)
    o(6) = tag_site(117367)
    o(7) = tag_site(144967)
    o(8) = tag_site(145205)
    o(9) = tag_site(145927)
    o(10) = tag_site(173527)
    o(11) = tag_site(174487)
    o(12) = tag_site(174245)
    o(13) = tag_site(116169)
    o(14) = tag_site(87129)
    o(15) = tag_site(144249)
    o(16) = tag_site(173289)
    o(17) = tag_site(202087)
    o(18) = tag_site(201845)
    !
    d(1) = tag_site(88327)
    d(2) = tag_site(87847)
    d(3) = tag_site(87367)
    d(4) = tag_site(87605)
    d(5) = tag_site(88085)
    d(6) = tag_site(116887)
    d(7) = tag_site(116407)
    d(8) = tag_site(115927)
    d(9) = tag_site(116165)
    d(10) = tag_site(116645)
    d(11) = tag_site(117125)
    d(12) = tag_site(145447)
    d(13) = tag_site(144487)
    d(14) = tag_site(144725)
    d(15) = tag_site(145685)
    d(16) = tag_site(173047)
    d(17) = tag_site(173285)    
    !
    !set the occupancy and spin of sites to that of water: 
    !all sites are occupied by water
    !the sites of lattice initially have no pointer
    occupancy = occW
    spin = spinW
    head = 0
    !
    do i=1,size(s)
       occupancy(s(i)) = occSN
    end do
    do i=1,size(d)
       occupancy(d(i)) = occSNOSNO
    end do
    !
    do i=1,size(o)
       occupancy(o(i)) = occSNO
    end do
    !
    do i=1,size(s)
       noccupy(i) = s(i)
    end do
    !
    do i = 1,size(s)
       head(s(i)) = 1
    end do
    !
    do i=1,size(s)
       spin(s(i)) = spinSN
    end do
    
    call psf_out
    call config_out_atom
    
    do monomer=1,allmolecules
       isite = noccupy(monomer);ivertex=head(isite)
       call site_rings3and4(isite,head(isite),rings3,rings4)
       print*,monomer,rings3,rings4
    end do
    !    
    return    
  end subroutine manual_check
  !********************************************************************
end program NVT_SiO4_TAA
