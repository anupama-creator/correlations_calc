!RADIAL DISTRIBUTION FUNCTION CODE (WORK FOR BOTH NVT & NPT)
       program msd
!===============================================================================
! Gromacs trajectroy reading utility available on github
!===============================================================================
       use gmxfort_trajectory
       use gmxfort_utils        
       implicit none
       type(Trajectory) :: trj
!=================================================================================
! This section is for the trajectory file input variables
!=================================================================================
       character(len=50) :: filename  ! input file name
       integer :: ierror              ! status flag for I/O statements
       character(len=80) :: msg       ! error message
       integer :: inp,ninput          ! number of input files to work with
       integer,parameter :: maxinp = 150 ! maximum number of trajectory files as an input
       character(len=50) :: inputfilename(maxinp) ! input file names to read information for calculation
       character(len=100) outfile1,outfile2,outfile3  ! out file names to write the final calculation results
       integer :: nmol_diff  ! number of different molecules present in the system
       integer :: nmol_PGC,nmol_water,nmol_PFL,natom_PGC,natom_water,natom_PFL ! information of atoms and molecules number
       character(len=5) :: atname1 ! atom for which we are finding correlation
       integer :: no_atname1
       character(len=100) :: grofile ! grofile to read and store information for calculation
       integer :: total_atom,natms
       character(len=100) head       ! grofile header
       integer :: at_num             ! total number of atoms present in gro file
       character(len=5) :: at_name,mol_name ! atom names defined in system
       character(len=80) :: infile   ! input trajectory
       integer :: frame
       integer :: fr_start,fr_end  ! begining and end of the frame to work with
       real(kind=8) :: x1,x2,y1,y2,z1,z2 ! to generate the cylinder
       real(kind=8),allocatable :: rw(:,:,:),rPGC(:,:,:),rPFL(:,:,:)
       character(len=5),allocatable :: atnm_rw(:,:),atnm_rPGC(:,:),atnm_rPFL(:,:)
       real(kind=8),allocatable :: norm(:)
       integer :: nread ,nskip,nused,intv
       integer :: n1               ! frames found in a trajectory
       real(kind=8) :: cell(3,3),boxl(3),dr_oo(3) !   ! lattice vectors
       integer :: ntot,nsites      ! number of atoms found in frame
       real(kind=8),allocatable :: xx1(:),yy1(:),zz1(:) ! coordinate of atoms
       real(kind=8):: xxx,yyy,zzz  ! temporary coordinate store variable
       integer :: atom_cbw
       integer :: imol_PGC,iatom_PGC,imol_water,iatom_water,imol_PFL,iatom_PFL
       integer :: good_water
       integer :: rgw,n ! 
       real(kind=8),allocatable :: distant(:)
       real(kind=8) :: maxdist,maxdist_sq,val,lent,r_oo
       integer,allocatable :: gwi(:),rgwi(:)
       integer(kind=1),allocatable :: ht_value(:,:) ! hit value to identify if water is inside the channel or not
       real(kind=8),allocatable :: rw_msd(:,:,:) ! msd of those water that are present inside the channel
       real(kind=8),allocatable :: r2(:),r4(:),alpha(:),fs(:),inv_r2(:),nngp(:)
       real(kind=8) :: dr(3)
       integer :: j1,mk,ii,im,jm
       integer :: i,j,k
       integer :: nstep
       integer :: t
       real(kind=8) :: sum_sq_dis,sum_4_dis,dimens
       real(kind=8),allocatable :: msd_tau(:)
       real :: start,finish
       real(kind=8) :: rsq,Q,r,iota,it,d
       character(len=5),allocatable :: atname(:)
       call cpu_time(start)
!=======================================================================
!read the number of input files for the calculation
!-----------------------------------------------------------------------
      write(*,*) "How many input files you have"
       read(*,*) ninput
       if (ninput .lt. 2) then
               write(*,*) "Enter the input filename for the analysis"
               read(*,*) inputfilename(ninput)
       else
               write(*,*) "=====================NOTE==================="
               write(*,*) "More than 1 file for the analysis.  ",ninput
       write(*,*) "Enter the input file name."
       read(*,'(A)') filename
       open(unit=2, file=filename,status='old',action='read',iostat=ierror,iomsg=msg)
errorcheck: IF (ierror > 0) then
       write(*,1010) filename
1010   format ('Error: File ',A,'does not exist !')
       write(*,'(A)') trim(msg)
       else
       !no of trajectory file you want to analyze'
       read(2,*) ninput   ! (maxinp=150)
       if(ninput>maxinp)then
               write(*,*)'Too many inputs'
       else
               write(*,*) "Working with ", ninput, " files." 
       endif
       do i=1,ninput
       read(2,*) inputfilename(i) ! input filename to read the files to work with and other information
       enddo
       endif errorcheck
       endif
!========================================================================
!Now read the input for the calculations from individual files
!------------------------------------------------------------------------
input: do inp = 1,ninput
       open(unit=3,file=inputfilename(inp),action='read')
       read(3,*)outfile1 ! output file name1
       open(unit=4,form='formatted',file=outfile1,action='write')
!=======================================================================
! Read the input data for the calculation
!=======================================================================
      write(*,*)'Enter the number of different molecule types'
       read(3,*) nmol_diff
       write(*,*)'Enter the no of molecules of all types'
       read(3,*) nmol_PGC, nmol_water,nmol_PFL
       write(*,*)'Enter the no of atom in each of the diff molecules'
       read(3,*) natom_PGC, natom_water,natom_PFL
       write(*,*)'nmol_PGC=',nmol_PGC,'nmol_water=',nmol_water,'nmol_PFL=',nmol_PFL
       write(*,*)'natom_PGC=',natom_PGC,'natom_water=',natom_water,'natom_PFL=',natom_PFL
       write(*,*)'Enter the name of atom for the calculation'
       read(3,*) atname1
       write(*,*)'Enter the number of atom for correlation'
       read(3,*) no_atname1
       write(*,*) 'Enter the distance in nm under which water molecule is to be tagged for orientational.'
       read(3,*) maxdist
       maxdist_sq = maxdist*maxdist
       write(*,*)'enter the name of gro file to read atom information'
       read(3,*) grofile
       write(*,*)'Finished reading inputs'
!=======================================================================
       total_atom = nmol_PGC*natom_PGC + nmol_water*natom_water + nmol_PFL*natom_PFL
!=======================================================================
! open *.gro file to read system information 
!-----------------------------------------------------------------------
       open(unit=10,file=grofile,status='old',action='read',iostat=ierror)
openif :     if (ierror == 0) then
!open was ok. Read values.
                read(10,*)head
                read(10,*)at_num
                allocate(atname(at_num))
                do n=1,at_num
                read(10,200) at_name,at_num
                !if (adjustl(at_name) == atname1) then
               ! write(44,*)at_num
                atname(at_num)=at_name
                !write(44,*) atname(at_num)
                200 format (10x,a5,i5)
       ! endif
       enddo
       endif openif
!=======================================================================
! from this line of code below trajectory file was taken directly from the file
!=======================================================================
       write(*,*)'Working with trajectory file number.',inp
       read (3,*) infile
       read (3,*) fr_start,fr_end
       frame = (fr_end-fr_start)+1
       write(*,*) "---------working with------",frame,"frames --------------"
!=======================================================================
!Number of atoms was going to be read and passed to a variable
!=======================================================================
       allocate(rw(nmol_water,natom_water,3),rPGC(nmol_PGC,natom_PGC,3),rPFL(nmol_PFL,natom_PFL,3))
       allocate(atnm_rw(nmol_water,natom_water),atnm_rPGC(nmol_PGC,natom_PGC),atnm_rPFL(nmol_PFL,natom_PFL))
       allocate(gwi(nmol_water),rgwi(nmol_water))
       allocate(distant(nmol_PFL))
       nread = 0
!----------------------------------------------------------------------- 
        call trj%open(infile)
        n1 = trj%read_next(fr_end)                          
!=======================================================================
!here the code will loop over the provided number of frames you want to 
!work with as an input trajectory
!======================================================================= 
nframe:     do t = fr_start,fr_end!,nskip
             nread = nread + 1     ! tag to identify how many frames have been read
             cell= trj%box(t)
             nsites = trj%natoms()
       allocate (xx1(nsites),yy1(nsites),zz1(nsites)) ! allocate for successive frame
site:         do j1=1,nsites
              xxx=trj%frameArray(t)%xyz(1,j1)
              yyy=trj%frameArray(t)%xyz(2,j1)
              zzz=trj%frameArray(t)%xyz(3,j1)
              xx1(j1)=xxx
              yy1(j1)=yyy
              zz1(j1)=zzz
      enddo site
      atom_cbw=0 ! function to categorize different atoms and molecules
      ntot = 0
PGC:  do imol_PGC=1,nmol_PGC
      do iatom_PGC=1,natom_PGC
      atom_cbw=atom_cbw+1
      rPGC(imol_PGC,iatom_PGC,1)=xx1(atom_cbw)
      rPGC(imol_PGC,iatom_PGC,2)=yy1(atom_cbw)
      rPGC(imol_PGC,iatom_PGC,3)=zz1(atom_cbw)
      atnm_rPGC(imol_PGC,iatom_PGC)=atname(atom_cbw)
      enddo
      enddo PGC

wt:   do imol_water=1,nmol_water
      do iatom_water=1,natom_water
      atom_cbw=atom_cbw+1
      rw(imol_water,iatom_water,1)=xx1(atom_cbw)
      rw(imol_water,iatom_water,2)=yy1(atom_cbw)
      rw(imol_water,iatom_water,3)=zz1(atom_cbw)
      atnm_rw(imol_water,iatom_water)=atname(atom_cbw)
      enddo
      enddo wt

PFL:  do imol_PFL=1,nmol_PFL
      do iatom_PFL=1,natom_PFL
      atom_cbw=atom_cbw+1
      rPFL(imol_PFL,iatom_PFL,1)=xx1(atom_cbw)
      rPFL(imol_PFL,iatom_PFL,2)=yy1(atom_cbw)
      rPFL(imol_PFL,iatom_PFL,3)=zz1(atom_cbw)
      atnm_rPFL(imol_PFL,iatom_PFL)=atname(atom_cbw)
      enddo
      enddo PFL

       do i = 1,3
       boxl(i) = cell(i,i)
       enddo
!-------------------------------------------------------------------------------
!---------------------------------------------------------------------
!tag those molecule that satisfy the criteria
!---------------------------------------------------------------------
      good_water = 0
wat:  do im = 1,nmol_water
      if (adjustl(atnm_rw(im,1)) == atname1) then
PF:   do jm = 1,nmol_PFL
      !   ja = 1,natom_PFL                     ! to check if wat is present near any of pfl atoms
      if (adjustl(atnm_rPFL(jm,7))=="O1") then ! to check 7th atom is oxygen
              dr_oo(:) = rPFL(jm,7,:)-rw(im,1,:)
              dr_oo(:) = dr_oo(:)-boxl(:)*anint(dr_oo(:)/boxl(:)) ! fix for pbc
              r_oo = dot_product(dr_oo,dr_oo)
              distant(jm) = r_oo ! store the distance based on propofol identity
      endif
              enddo PF            ! for the reason we are looking water near pfl
              if (minval(distant) .lt. maxdist_sq) then ! 0.5**2 nm
              good_water = good_water+1
              gwi(good_water) = im  ! storing the identity of water molecule in an index
              endif
      endif
              enddo wat
  if (nread==1) then    ! necessary to identify that atoms is present in both the frames
          rgw = good_water ! total number of good water in first frame
          rgwi(1:good_water) = gwi(1:good_water) ! transfer the identity of water stored in first frame to rgwi
  else
! to check the increment of water if there is any from the first and consecutive frames after that
          do imol_water = 1,good_water ! for second and consecutive frame
          if (any(rgwi==gwi(imol_water))) then ! check of water id from previous store
          else
                  rgw = rgw + 1        ! increment the number of good water if there is any new wat comes near
                  rgwi(rgw) = gwi(imol_water) ! store the id of that wat mol
          endif
          enddo
  endif
                deallocate (xx1,yy1,zz1)
                enddo nframe
            write(*,*) "Water in 5 angstorm of Propofol molecule is ",rgw
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! indexing of water present inside the channel was done. Now storing the coordinate of those atoms
!-------------------------------------------------------------------------------
            allocate(ht_value(fr_start:fr_end,rgw))       ! allocate when required
            allocate(rw_msd(fr_start:fr_end,nmol_water,3))
            ht_value = 0
            atom_cbw = 0
            nread = 0
nframe1:     do t = fr_start,fr_end!,nskip
             nread = nread + 1
             allocate (xx1(nsites),yy1(nsites),zz1(nsites))
site1:        do j1=1,nsites
              xxx=trj%frameArray(t)%xyz(1,j1)
              yyy=trj%frameArray(t)%xyz(2,j1)
              zzz=trj%frameArray(t)%xyz(3,j1)
              xx1(j1)=xxx
              yy1(j1)=yyy
              zz1(j1)=zzz
      enddo site1
      atom_cbw=0 ! function to categorize different atoms and molecules
      ntot = 0
PGC1:  do imol_PGC=1,nmol_PGC
       do iatom_PGC=1,natom_PGC
      atom_cbw=atom_cbw+1
      rPGC(imol_PGC,iatom_PGC,1)=xx1(atom_cbw)
      rPGC(imol_PGC,iatom_PGC,2)=yy1(atom_cbw)
      rPGC(imol_PGC,iatom_PGC,3)=zz1(atom_cbw)
      atnm_rPGC(imol_PGC,iatom_PGC)=atname(atom_cbw)
      enddo
      enddo PGC1

wt1:               do imol_water=1,nmol_water
                  do iatom_water=1,natom_water
                        atom_cbw=atom_cbw+1
                        rw(imol_water,iatom_water,1)=xx1(atom_cbw)
                        rw(imol_water,iatom_water,2)=yy1(atom_cbw)
                        rw(imol_water,iatom_water,3)=zz1(atom_cbw)
                        atnm_rw(imol_water,iatom_water)=atname(atom_cbw)
   enddo
   enddo wt1
PFL1:               do imol_PFL=1,nmol_PFL
                  do iatom_PFL=1,natom_PFL
                        atom_cbw=atom_cbw+1
                        rPFL(imol_PFL,iatom_PFL,1)=xx1(atom_cbw)
                        rPFL(imol_PFL,iatom_PFL,2)=yy1(atom_cbw)
                        rPFL(imol_PFL,iatom_PFL,3)=zz1(atom_cbw)
                        atnm_rPFL(imol_PFL,iatom_PFL)=atname(atom_cbw)
   enddo
   enddo PFL1

   do i = 1,3
   boxl(i) = cell(i,i)
   enddo
!-------------------------------------------------------------------------------
! storing coordinates with time of those water molecules whose index is found in the defined region of space
!-------------------------------------------------------------------------------
      do imol_water = 1,rgw
         ii = rgwi(imol_water)
      do jm = 1,nmol_PFL
      !   ja = 1,natom_PFL                     ! to check if wat is present near any of pfl atoms
      if (adjustl(atnm_rPFL(jm,7))=="O1") then ! to check 7th atom is oxygen
              dr_oo(:) = rPFL(jm,7,:)-rw(ii,1,:)
              dr_oo(:) = dr_oo(:)-boxl(:)*anint(dr_oo(:)/boxl(:))
              r_oo = (sqrt(dot_product(dr_oo,dr_oo)))
              distant(jm) = r_oo              ! store the distance based on propofol identity
      endif
              enddo                            ! for the reason we are looking water near pfl
             if (minval(distant) .lt. maxdist) then
             ht_value(t,imol_water) = 1
              rw_msd(t,imol_water,:) = rw(ii,1,:)
     endif
     enddo
     deallocate(xx1,yy1,zz1)
     enddo nframe1
     deallocate(rw,atnm_rw,rPGC,atnm_rPGC,rPFL,atnm_rPFL,gwi,rgwi) ! deallocate all the variables not required for the calculation 
     deallocate(atname,distant)                                            ! to realese the memory
     allocate(msd_tau(0:fr_end-fr_start),norm(0:fr_end-fr_start))!,fs(0:fr_end-fr_start))           ! allocate the ones required for calculation
       norm(:) = 0.0
       msd_tau(:) = 0.0
!===============================================================================
!calculation of mean squared displacement
!===============================================================================
br2:    do imol_water = 1,nmol_water
        do i = fr_start,fr_end   ! origin
        k = 0                    ! origin+dt calculation on every origin interval
lp:     do j = 0,fr_end-fr_start ! dt
        k=i+j
        if (k.le.fr_end) then    ! for value of k larger than fr_end
!--------------------------------------------------------------------------------
       dr(:) = rw_msd(k,imol_water,:)-rw_msd(i,imol_water,:)
       rsq = dot_product(dr,dr)                                  ! squared displacement z contribution
       msd_tau(j) = msd_tau(j) + rsq                             ! msd(dt)
       norm(j) = norm(j) + 1.0d0
       endif
       enddo lp
       enddo 
       enddo br2
!===============================================================================
!
!===============================================================================
                do i = 0,fr_end-fr_start    ! dt
                if (norm(i) == 0.0d0) then
                msd_tau(i) = 0.0d0
                else
                        msd_tau(i) = msd_tau(i)/norm(i)
                endif
                write(4,1011) i*1.0d0, msd_tau(i)! in ps
1011  format(F20.10,F20.10)
                enddo
               call trj%close()
               deallocate(msd_tau,norm)
               enddo input
!-------------------------------------------------------------------------------
               call cpu_time(finish)
               print'("Time= ",f6.3,"seconds")',finish-start,finish,start
               end program msd
!program was closed
!-g -fcheck=all -Wall -fbacktrace command to check error
!---------------------------------------------------------------------
