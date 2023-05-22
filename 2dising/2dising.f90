! This program implements Metropolis MC for the 2D Ising model
! Reduced units, so J = kB = 1, H = 0

program twodising

  implicit none

  integer,parameter::double=selected_real_kind(15)

  integer :: nt,nti,i,j,iplus,iminus,jplus,jminus,magnetization,u_old,du,nran,single_seed
  real(kind=double) :: ran,accprob,m_per_spin,average_magnetization,expvalues_4,expvalues_8
  logical :: random_init

  integer,dimension(:),allocatable :: iran

  integer :: n            ! lattice size
  real(kind=double) :: t  ! temperature
  integer :: ntskip       ! thinning factor
  integer :: ntjob        ! number of MC passes

  integer,dimension(:,:),allocatable :: lattice

! read parameters from input file

  open(unit=10,file='options.dat')
  read(unit=10,fmt=*)
  read(unit=10,fmt=*) n
  read(unit=10,fmt=*) t
  read(unit=10,fmt=*) ntskip
  read(unit=10,fmt=*) ntjob
  read(unit=10,fmt=*) random_init
  read(unit=10,fmt=*) single_seed
  close(unit=10)

! seed the RNG from value in file

  call random_seed(size=nran)
  allocate(iran(nran))
  do i=1,nran
    iran(i) = single_seed**i
  end do
  call random_seed(put=iran)

! save values for acceptance test for efficiency

  expvalues_4 = exp(-4.0_8/t)
  expvalues_8 = exp(-8.0_8/t)

! allocate dynamic array for lattice

  allocate(lattice(0:n-1,0:n-1))

! start with a random lattice if asked for

  if (random_init) then
    do i=0,n-1
      do j=0,n-1
        call random_number(ran)
        if (ran < 0.5) then
          lattice(i,j) = 1
        else
          lattice(i,j) = -1
        end if
      end do
    end do
  else
    lattice = 1
  end if
  
! run the actual MC

  average_magnetization = 0.0
  
  open(unit=20,file='magnetization.dat')
  write(unit=20, fmt='(a)') '# time M m'
  
  do nt=1,ntjob

    do nti=1,n**2*ntskip
    
! select a random spin

      call random_number(ran)
      i = int(ran*n)
      call random_number(ran)
      j = int(ran*n)
      
! neighbor spins using PBC
      
      iplus = modulo(i+1,n)
      iminus = modulo(i-1,n)
      jplus = modulo(j+1,n)
      jminus = modulo(j-1,n)
      
      u_old = -1*lattice(i,j)*( lattice(i,jplus) + lattice(i,jminus) + lattice(iminus,j) + lattice(iplus,j) )
            
      du =  -2*u_old
      
      if (du <= 0) then
        lattice(i,j) = -1*lattice(i,j)
      else
        if (du==4) then
          accprob = expvalues_4
        elseif (du==8) then
          accprob = expvalues_8
        end if
        call random_number(ran)
        if (ran < accprob) then 
          lattice(i,j) = -1*lattice(i,j)
        end if
      end if
          
    end do
    
! magnetization = sum of spins

    magnetization = 0
    do i=0,n-1
      do j=0,n-1
        magnetization = magnetization + lattice(i,j)
      end do
    end do
    m_per_spin = real(magnetization) / n**2
    average_magnetization = average_magnetization + m_per_spin
    
    write(unit=*,fmt=*) nt, magnetization
    write(unit=20, fmt='(2i16,f12.8)') nt, magnetization, m_per_spin 

  end do
  
  close(unit=20)
  
  average_magnetization = average_magnetization / ntjob
  
  write(*,*) "# Average magnetization per spin: ", average_magnetization

  deallocate(lattice, iran)
  
end program twodising
