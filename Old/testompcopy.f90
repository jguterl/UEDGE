module dim
integer,parameter:: Nx=3,Ny=20
integer,dimension(Ny)::yindex
integer,dimension(Nx,Ny),save :: fngx
 integer::PARALLELCOMP=0
!$OMP threadprivate(fngx)
end module

program testomp
   use omp_lib
   use dim
   implicit none
   integer i,j,ithread
   integer, parameter :: Nthread = 4


   integer,dimension(Nx,Ny),save :: fngx_

  do j=1,Ny
  do i = 1, Nx
      fngx_(i,j)=i+j*j
      enddo
   enddo
yIndex(1:ny)=-1
!$OMP PARALLEL DO Default(shared) PRIVATE(i,ithread) schedule(static) if (PARALLELCOMP.gt.0)
   do j = 1, Ny
   write(*,*)'0-ithread ',omp_get_thread_num()
yIndex(j)=omp_get_thread_num()
      do i=1,Nx
      fngx(i,j)=i+j*j
      enddo
      write(*,*)'0-fngx ',omp_get_thread_num(),fngx(1:Nx,j)
   enddo
!$OMP END PARALLEL DO
write(*,*) 'yindex=',yindex
call copyfngx()
write(*,*) '------------'
do i = 1, Nx
      do j=1,Ny
      write(*,'(a,i2,a,i2,a,i3,a,i3)') 'fngx(',i,',',j,')=',fngx(i,j),'|',fngx_(i,j)
      enddo
   enddo
end program testomp

subroutine copyfngx()
use omp_lib
use dim
integer,allocatable::fngxcopy(:,:)
integer::idxymin=1,idxymax=Ny
integer::idxxmin=1,idxxmax=Nx
allocate(fngxcopy(Nx,Ny))
ithread=0
call CheckIndex(idxxmin)



!$OMP PARALLEL SHARED(fngxcopy)
!$OMP& firstprivate(,idxxmin,idxxmax,idxymin,idxymax) private(ithread) if (PARALLELCOMP.gt.0)
!$OMP CRITICAL

write(*,*) ithread,idxymin,idxymax,idxxmin,idxxmax
fngxcopy(idxxmin:idxxmax,idxymin:idxymax)=fngx(idxxmin:idxxmax,idxymin:idxymax)

!$OMP END CRITICAL
!$OMP END PARALLEL
write(*,*) fngxcopy
fngx=fngxcopy
deallocate(fngxcopy)
end subroutine copyfngx


subroutine GetRangeIndex(idxxmin,idxxmax,idxymin,idxymax,ithread)
use dim,only:nx,ny
use OMPPandf,only:OMPyindex,OMPxindex

integer,intent(out)::idxxmin,idxxmax,idxxmin,idxxmax
integer,intent(in)::ithread
integer:: j
logical:: FindMinx,Findminy
idxxmin=0
idxxmax=Nx+1
idxymin=0
idxymax=0
FindMiny=.false.
do j=1,Ny
if (yindex(j)==ithread) then
if (.not.FindMiny) then
idxymin=j
FindMin=.True.
endif
idxymax=j
endif
enddo

FindMinx=.false.
do i=1,Nx
if (xindex(j)==ithread) then
if (.not.FindMinx) then
idxxmin=j
FindMin=.True.
endif
idxxmax=j
endif
enddo

return
end subroutine

subroutine checkindexes()
end subroutine checkindexes()

