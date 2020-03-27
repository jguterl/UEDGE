! omp_test.f90
! gfortran -fopenmp omp_test.f90 -oomp_test
module M
   use ISO_C_BINDING
   implicit none
   private
   public GetCurrentThreadStackLimits
   interface
      subroutine GetCurrentThreadStackLimits(LowLimit,HighLimit) &
         bind(C,name='GetCurrentThreadStackLimits')
         import
         implicit none
!GCC$ ATTRIBUTES STDCALL:: GetCurrentThreadStackLimits
!DEC$ ATTRIBUTES STDCALL:: GetCurrentThreadStackLimits
         type(C_PTR) LowLimit
         type(C_PTR) HighLimit
      end subroutine GetCurrentThreadStackLimits
   end interface
   type, public :: T
      type(C_PTR) :: Address = C_NULL_PTR
      type(C_PTR) :: LowLimit = C_NULL_PTR
      type(C_PTR) :: HighLimit = C_NULL_PTR
   end type T
   type(T), allocatable, public :: Tarray(:)
   public S
   contains
      recursive subroutine S(i)
         use omp_lib
         integer i
         integer, target :: j
         j = OMP_GET_THREAD_NUM()
         Tarray(j)%Address = C_LOC(j)
         call GetCurrentThreadStackLimits(Tarray(j)%LowLimit,Tarray(j)%HighLimit)
      end subroutine S
end module M

program P
   use M
   use ISO_C_BINDING
   use omp_lib
   implicit none
   integer, parameter :: N = 100
   integer i
   integer, parameter :: Nthread = 10

   allocate(Tarray(0:Nthread-1))
   call OMP_SET_NUM_THREADS(Nthread)
   call OMP_SET_DYNAMIC(.FALSE.)
!$OMP PARALLEL DO     &
!$OMP DEFAULT(NONE)   &
!$OMP private(i)      &
!$OMP SHARED(Tarray)
   do i = 1, N
      call S(i)
   end do
!$OMP END PARALLEL DO
   do i = 0, Nthread-1
      write(*,*) 'Thread = ', i
      write(*,*) 'Address = ', transfer(Tarray(i)%address,0_C_INTPTR_T)
      write(*,*) 'LowLimit = ', transfer(Tarray(i)%LowLimit,0_C_INTPTR_T)
      write(*,*) 'HighLimit = ', transfer(Tarray(i)%HighLimit,0_C_INTPTR_T)
   end do
end program P
