#include "bbb.h"
#include "../com/com.h"
#include "../mppl.h"
#include "../sptodp.h"
!-------------------------------------------------------------------------------------------------
subroutine InitParallel
    Use Output
    Use ParallelJacOptions,only:OMPParallelJac,MPIParallelJac,ParallelJac
    Use HybridOptions,only:HybridOMPMPI
    implicit none
    if (OMPParallelJac==1 .and. MPIParallelJac==0) then
        call InitOMP
        ParallelJac=1
    elseif (OMPParallelJac==0 .and. MPIParallelJac==1) then
        call InitMPI
        ParallelJac=1
    elseif (OMPParallelJac==1 .and. MPIParallelJac==1) then
        HybridOMPMPI=1
        call InitMPI() ! MPI first to setup nnzmxperproc used by InitOMP
        call InitOMP()
        ParallelJac=1
    else
        write(iout,'(a)') '*OMPJac* Jacobian calculation: OMP not enabled'
        write(iout,'(a)') '*MPIJac* Jacobian calculation: MPI not enabled'
        ParallelJac=0
    endif
end subroutine InitParallel
!-------------------------------------------------------------------------------------------------
subroutine jac_calc_parallel(neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)

    Use Output
    Use ParallelJacOptions,only:OMPParallelJac,MPIParallelJac
    implicit none
    ! ... Input arguments:
    integer,intent(in):: neq      !      total number of equations (all grid points)
    real,intent(in)   :: t              ! physical time
    real,intent(in)   :: yl(*)          ! dependent variables
    real,intent(in)   :: yldot00(neq+2) ! right-hand sides evaluated at yl
    integer,intent(in):: ml, mu         ! lower and upper bandwidths
    integer,intent(in):: nnzmx          ! maximum number of nonzeros in Jacobian

    ! ... Output arguments:
    real,intent(out)   :: jac(nnzmx)     ! nonzero Jacobian elements
    integer,intent(out):: ja(nnzmx)   ! col indices of nonzero Jacobian elements
    integer,intent(out):: ia(neq+1)   ! pointers to beginning of each row in jac,ja

    ! ... Work-array argument:
    real wk(neq)     ! work space available to this subroutine

    if (OMPParallelJac==1 .and. MPIParallelJac==0) then
        call jac_calc_omp (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)
    elseif (OMPParallelJac==0 .and. MPIParallelJac==1) then
        call jac_calc_mpi (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)
    elseif (OMPParallelJac==1 .and. MPIParallelJac==1) then
        call jac_calc_hybrid (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)
    else
        call xerrab('Cannot call calc_jac_parallel OMP and MPI jac calc are not enabled')
    endif

end subroutine jac_calc_parallel
!-------------------------------------------------------------------------------------------------
#ifdef _OPENMP
subroutine InitOMP()
    Use Output
    Use OmpOptions,only:Nthreads,nnzmxperthread,omplenpfac,ompneq,OMPVerbose,OMPStamp
    Use MPIOptions,only:nnzmxperproc,MPIrank
    Use HybridOptions,only: HybridOMPMPI,Hybridstamp
    Use Jacobian,only:nnzmx
    Use Lsode, only:neq

    implicit none
    integer:: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
    character*8 :: MPIRankTag

    ! prepare MPI/Hybrid stamps for output
    if (HybridOMPMPI>0) then
        write(MPIRankTag,'(I4)') MPIrank
        write(OMPstamp,'(a,a,a)') '[',trim(adjustl(trim(MPIRankTag))),'] OMPJac* '
    else
        write(OMPstamp,'(a)') '*OMPJac* '
    endif

    !$omp parallel
    if (OMP_GET_THREAD_NUM().eq.0) then
        if (Nthreads.gt.OMP_GET_NUM_THREADS()) then
            if (OMPVerbose.gt.1) write(iout,*) OMPStamp,' Warning: # threads requested > # threads available'
            Nthreads=OMP_GET_NUM_THREADS()
            write(iout,*) OMPStamp,' Nthreads:', Nthreads
        endif
        if (Nthreads.le.0) then
            call xerrab('Nthread must be >0')
        endif
        if (OMPVerbose.gt.0) write(iout,'(a,a,i3)') OMPStamp,' Number of threads for omp calculations:',Nthreads
    endif
    !$omp END parallel

    if (Nthreads.gt.1) then
        if (HybridOMPMPI>0) then
            nnzmxperthread=ceiling(real(nnzmxperproc)/real(Nthreads-1)*omplenpfac)
        else
            nnzmxperthread=ceiling(real(nnzmx)/real(Nthreads-1)*omplenpfac)
        endif
    else
        if (HybridOMPMPI>0) then
            nnzmxperthread=nnzmxperproc
        else
            nnzmxperthread=nnzmx
        endif
    endif
    ompneq=neq
    call gchange('OmpJacobian',0)
end subroutine InitOMP
!-------------------------------------------------------------------------------------------------
subroutine OMPCollectJacobian(neq,nnzmx,rcsc,icsc,jcsc,nnzcumout)
    use Output
    use OmpJacobian,only:iJacCol,rJacElem,iJacRow,OMPivmin,OMPivmax,nnz,nnzcum
    use OmpOptions,only:Nthreads,OMPDebug,OMPVerbose,OMPStamp
    integer,intent(in):: neq
    integer,intent(in):: nnzmx          ! maximum number of nonzeros in Jacobian
    real,intent(out)   :: rcsc(nnzmx)     ! nonzero Jacobian elements
    integer,intent(out):: icsc(nnzmx)   ! col indices of nonzero Jacobian elements
    integer,intent(out):: jcsc(neq+1)   ! pointers to beginning of each row in jac,ja
    integer,intent(out):: nnzcumout
    integer ith
    nnzcum(1:Nthreads)=-1
    nnzcum(1)=nnz(1)-1

    do ith=2,Nthreads
        nnzcum(ith)=nnzcum(ith-1)+nnz(ith)-1
    enddo
    if (OMPDebug.gt.0) then
        write(iout,*) OMPStamp,' nnz:',nnz(1:Nthreads)
        write(iout,*) OMPStamp,' nnzcum:',nnzcum(1:Nthreads)
    endif
    if (OMPVerbose.gt.0) write(iout,'(a,i9)') '**** Number of non-zero Jacobian elems:',nnzcum(Nthreads)

    if (nnzcum(Nthreads).gt.nnzmx) then
        call xerrab(' Problem: nnzcum > nnzmx...')
    endif



    jcsc(OMPivmin(1):OMPivmax(1))= iJacRow(OMPivmin(1):OMPivmax(1),1)
    do ith=2,Nthreads
        jcsc(OMPivmin(ith):OMPivmax(ith))= iJacRow(OMPivmin(ith):OMPivmax(ith),ith)+nnzcum(ith-1)
    enddo

    rcsc(1:nnz(1)-1)= rJacElem(1:nnz(1)-1,1)
    icsc(1:nnz(1)-1)= iJacCol(1:nnz(1)-1,1)
    do ith=2,Nthreads
        rcsc(nnzcum(ith-1)+1:nnzcum(ith))=rJacElem(1:nnz(ith)-1,ith)
        icsc(nnzcum(ith-1)+1:nnzcum(ith))=iJacCol(1:nnz(ith)-1,ith)
    enddo
    nnzcumout=nnzcum(Nthreads)
end subroutine OMPCOllectJacobian
!-------------------------------------------------------------------------------------------------
subroutine jac_calc_omp (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)

    !   Calculate Jacobian matrix (derivatives with respect to each
    !   dependent variable of the right-hand side of each rate equation).
    !   Lower and upper bandwidths are used to select for computation
    !   only those Jacobian elements that may be nonzero.
    !   Estimates of Jacobian elements are computed by finite differences.
    !   The Jacobian is stored in compressed sparse row format.
    Use Output
    use Timing,only:istimingon,ttjstor,ttotjf,ttimpjf,ShowTime,&
    OMPTotTimeCollect,OMPTotTimeBuild,OMPTotJacCalc
    use Grid,only:ngrid,ig,ijac,ijactot
    use Jacobian_csc,only:rcsc,jcsc,icsc,yldot_pert,yldot_unpt
    use OmpOptions,only:OMPDebug,Nthreads,OMPVerbose,OMPDebug,nnzmxperthread,OMPCheckNaN,WriteJacobian,&
        OMPLoadBalance,OMPAutoBalance,OMPStamp,OMPBalanceStrength
    use OmpJacobian,only:iJacCol,rJacElem,iJacRow,OMPivmin,OMPivmax,nnz,nnzcum,OMPLoadWeight,OMPTimeLocalJac
    use UEpar, only: svrpkg
    use Math_problem_size,only:neqmx

    implicit none
    ! ... Input arguments:
    integer,intent(in):: neq      !      total number of equations (all grid points)
    real,intent(in)   :: t              ! physical time
    real,intent(in)   :: yl(*)          ! dependent variables
    real,intent(in)   :: yldot00(neq+2) ! right-hand sides evaluated at yl
    integer,intent(in):: ml, mu         ! lower and upper bandwidths
    integer,intent(in):: nnzmx          ! maximum number of nonzeros in Jacobian

    ! ... Output arguments:
    real,intent(out)   :: jac(nnzmx)     ! nonzero Jacobian elements
    integer,intent(out):: ja(nnzmx)   ! col indices of nonzero Jacobian elements
    integer,intent(out):: ia(neq+1)   ! pointers to beginning of each row in jac,ja

    ! ... Work-array argument:
    real wk(neq)     ! work space available to this subroutine
    integer,allocatable :: iJacConstructor(:,:)
    real,allocatable:: rJacConstructor(:,:)
    integer:: nnzcumout
    ! ... Functions
    logical tstguardc
    real(kind=4) gettime

    ! ... Local variables:
    real(kind=4) sec4, tsjstor, tsimpjf, dtimpjf,time0,time1,OMPTimeBuild,OMPTimeCollect,OMPTimeJacCalc
    integer:: i,thread,ith,iv,TID, OMP_GET_THREAD_NUM
    character(len = 80) ::  filename
    ! Calculate load distribution for threads
    if (OMPLoadBalance.ne.1 .and. OMPAutoBalance.ne.1) then
        OMPLoadWeight(1:Nthreads)=1.0
    endif
    if (OMPAutoBalance.eq.1) then
        !Check that time are not zero
        if (OMPBalanceStrength<=0) call xerrab('OMPBalanceStrength must be >0')
        if (minval(OMPTimeLocalJac).gt.0.0) then
            do i=1,Nthreads
                OMPLoadWeight(i)=OMPLoadWeight(i)*1/(OMPTimeLocalJac(i)/sum(OMPTimeLocalJac)*real(Nthreads))**OMPBalanceStrength
            enddo
        else
            OMPLoadWeight(1:Nthreads)=1.0
        endif
    endif
    !   Get the range of the iv index for each thread
    call OMPSplitIndex(1,neq,Nthreads,OMPivmin,OMPivmax,OMPLoadWeight)

    if (OMPVerbose.gt.0) then
        write(iout,*)' *OMPJac* neq=',neq,neqmx
        write(iout,*)' *OMPJac* Ivmin(ithread),Ivmax(ithread), OMPLoadWeight(ithread) ***'
        do ith=1,Nthreads
            write(iout,'(a,I3,a,I7,I7,f5.1)') '  *    ithread ', ith,':',OMPivmin(ith),OMPivmax(ith),OMPLoadWeight(ith)
        enddo
    endif


    !    iJacCol(1:nnzmxperthread,1:Nthreads)=0
    !    rJacElem(1:nnzmxperthread,1:Nthreads)=0.0
    !    iJacRow(1:neq,1:Nthreads)=0
    !    if (OMPDebug.gt.0) then
    !        write(iout,*) OMPStamp,' Jacobian arrays set to zero'
    !    endif
    OMPTimeJacCalc= gettime(sec4)

    !   Get initial value of system cpu timer.
    if (istimingon .eq. 1) tsjstor = gettime(sec4)

    !   Count Jacobian evaluations, both for total and for this case
    ijactot = ijactot + 1
    ijac(ig) = ijac(ig) + 1
    if (svrpkg.eq.'nksol') write(iout,*) ' Updating Jacobian, npe =  ',ijac(ig)

    !   Set up diagnostic arrays for debugging
    do iv = 1, neq
        yldot_unpt(iv) = yldot00(iv)  ! for diagnostic only
        yldot_pert(iv) = 0.
    enddo

    !   build jacobian ##############################################################
    OMPTimeBuild=gettime(sec4)
    nnz(1:Nthreads)=-1
    call OMPJacBuilder(neq, t, yl,yldot00, ml, mu,wk,iJacCol,rJacElem,iJacRow,nnz)

    OMPTimeBuild=gettime(sec4)-OMPTimeBuild
    if (istimingon .eq. 1) OMPTotTimebuild = OMPTimeBuild+OMPTotTimebuild
    if (OMPVerbose.gt.0) write(iout,*)OMPStamp,' Time to build jac:',OMPTimeBuild
    !   end build jacobian ##############################################################

    !   collect jacobian ##############################################################
    OMPTimeCollect=gettime(sec4)
    call OMPCollectJacobian(neq,nnzmx,rcsc,icsc,jcsc,nnzcumout)
    OMPTimeCollect=gettime(sec4)-OMPTimeCollect
    if (istimingon .eq. 1) OMPTotTimeCollect = OMPTimeCollect+OMPTotTimeCollect
    if (OMPVerbose.gt.0) write(iout,*)OMPStamp,' Time to collect jac:',OMPTimeCollect
    !   end collect jacobian ##############################################################

    jcsc(neq+1) = nnzcumout+1 ! This is set here out of OMPJAcCollect for compatibility with hybrid jac_calc

    !   for Debug purpose
    if (WriteJacobian.eq.1) then
        write(filename,'(a,3i3,a)') "jac_omp_",ijac(ig),".txt"
        call jac_write(filename,neq, rcsc, icsc, jcsc)
    endif

    if (OMPCheckNaN.gt.0) then
        do i=1,nnzmx
            if (isnan(rcsc(i))) then
                write(iout,*) 'rcsc is NaN at i=',i,rcsc(i)
            endif
            if (icsc(i).ne.icsc(i)) then
                write(iout,*) 'icsc is NaN at i=',i,icsc(i)
            endif
        enddo
        do i=1,neq+1
            if (jcsc(i).ne.jcsc(i)) then
                write(iout,*) 'jcsc is NaN at i=',i
            endif
        enddo
    endif

    !   Convert Jacobian from compressed sparse column to compressedsparse row format.
    time1=gettime(sec4)
    call csrcsc (neq, 1, 1, rcsc, icsc, jcsc, jac, ja, ia)

    !! ... If desired, calculate Jacobian elements for ion continuity
    !!     equation by using INEL routines, and combine elements with
    !!     those calculated above.
    !      if (iondenseqn .eq. "inel") then
    !         call iondens2 (neq, nnzmx, jac, ja, ia, rcsc, icsc, jcsc)
    !      endif

    !! ... If necessary, calculate Jacobian elements for rows corresponding
    !!     to impurity-density equations (for cells other than guard cells),
    !!     and combine elements with those calculated above.
    !      if (isimpon .eq. 3 .or. isimpon .eq. 4) then
    !         if (istimingon .eq. 1) tsimpjf = gettime(sec4)
    !         call impdens2 (neq, nnzmx, jac, ja, ia, rcsc, icsc, jcsc)
    !         if (istimingon .eq. 1) then
    !            dtimpjf = gettime(sec4) - tsimpjf
    !            ttimpjf = ttimpjf + dtimpjf
    !            ttotjf = ttotjf + dtimpjf
    !         endif
    !      endif

    if (istimingon .eq. 1) ttjstor = ttjstor + gettime(sec4) - tsjstor

    OMPTimeJacCalc=gettime(sec4)-OMPTimeJacCalc
    if (istimingon .eq. 1) OMPTotJacCalc = OMPTimeJacCalc+OMPTotJacCalc
    if (ShowTime.gt.0) write(iout,*)'**** Time in jac_calc:',gettime(sec4)-OMPTimeJacCalc

    return
end subroutine jac_calc_omp
!-------------------------------------------------------------------------------------------------
subroutine OMPJacBuilder(neq, t, yl,yldot00, ml,mu,wk,iJacCol,rJacElem,iJacRow,nnz)
    use OmpOptions,only:OMPDebug,OMPCopyArray,OMPCopyScalar,nthreads,nnzmxperthread,OMPStamp,OMPVerbose
    use OMPJacobian, only:OMPivmin,OMPivmax,OMPTimeLocalJac
    use OmpCopybbb
    use OmpCopycom
    use OmpCopyapi
    use omp_lib

    implicit none
    integer,intent(inout)::nnz(Nthreads)
    integer,intent(in):: neq      ! total number of equations (all grid points)
    integer,intent(in):: ml, mu   ! lower and upper bandwidths
    real,intent(in):: t           ! physical time
    real,intent(in) ::yl(*)       ! dependent variables
    real,intent(in) ::yldot00(neq+2) ! right-hand sides evaluated at yl
    real,intent(inout) :: wk(neq)
    integer,intent(out)::iJacCol(nnzmxperthread,nthreads)
    integer,intent(out):: iJacRow(neq,nthreads)
    real,intent(out):: rJacElem(nnzmxperthread,nthreads)
    real ::wkcopy(neq)
    real::ylcopy(neq+2)

    integer ::iJacColCopy(nnzmxperthread),iJacRowCopy(neq)
    integer ::ivmincopy(Nthreads),ivmaxcopy(Nthreads),nnzcopy(Nthreads)
    integer ::Nthreadscopy,nnzmxperthreadcopy
    real :: rJacElemCopy(nnzmxperthread)
    integer:: ith,tid,nnzlocal,ithcopy
    DOUBLE PRECISION :: TimeThread

    if (OMPDebug.gt.0)write(iout,*) OMPStamp,' Copying data....'

    if (OMPCopyArray.gt.0) then
        if (OMPDebug.gt.0)write(iout,*) OMPStamp,' Copying array....'
        call OmpCopyPointerbbb
        call OmpCopyPointercom
        call OmpCopyPointerapi
    endif
    if (OMPCopyScalar.gt.0) then
        if (OMPDebug.gt.0)write(iout,*) OMPStamp,' Copying scalar....'
        call OmpCopyScalarbbb
        call OmpCopyScalarcom
        call OmpCopyScalarapi
    endif
    !   We cannot use variables in the parallel construct declarations below when these variables are not in the scope of the subroutine
    Nthreadscopy=Nthreads
    nnzmxperthreadcopy=nnzmxperthread
    ivmincopy(1:Nthreads)=OMPivmin(1:Nthreads)
    ivmaxcopy(1:Nthreads)=OMPivmax(1:Nthreads)
    iJacColCopy(1:nnzmxperthread)=0
    rJacElemCopy(1:nnzmxperthread)=0.0
    iJacRowCopy(1:neq)=0
    ylcopy(1:neq+2)=yl(1:neq+2) ! a very barbarian use of yl(neq+1) is implemented as a switch in pandf... Error-prone!
    wkcopy(1:neq)=wk(1:neq) ! Could be set equal to zero as well. The worker wk is not an output...

    if (OMPDebug.gt.0) then
        write(iout,*) OMPStamp,' Starting parallel loop'
    endif
    tid=-1
    nnzlocal=-10000
    ! ivmincopy,ivmaxcopy,yldot00, neq an t  could be shared as well as well as
    !$omp parallel do default(shared)&
    !$omp& firstprivate(ithcopy,ivmincopy,ivmaxcopy,tid,nnzlocal,ylcopy,wkcopy,ml,mu,yldot00,t,neq)&
    !$omp& firstprivate(nnzmxperthreadcopy,nthreadscopy,iJacRowCopy,iJacColCopy,rJacElemCopy)&
    !$omp& private(TimeThread)

    loopthread: do ith=1,Nthreads !ith from 1 to Nthread, tid from 0 to Nthread-1
        Timethread = omp_get_wtime()
        tid=omp_get_thread_num()
        ithcopy=ith

        if (OMPDebug.gt.0) write(iout,*) OMPStamp,' Thread id:',tid,' <-> ith:',ithcopy
        ! we keep all these parameters as it is easier to debug LocalJacBuilder and deal with private/shared attributes

        call LocalJacBuilder(ivmincopy(ithcopy),ivmaxcopy(ithcopy),neq, t, ylcopy,yldot00,ml,mu,wkcopy,&
            iJacColcopy,rJacElemcopy,iJacRowcopy,ithcopy,nnzlocal,nnzmxperthreadcopy,nthreadscopy)
        if (OMPDebug.gt.0) write(iout,*) OMPStamp,',',tid,' nzlocal:',nnzlocal
        !$omp  critical
        iJacCol(1:nnzlocal,ithcopy)=iJacColCopy(1:nnzlocal)
        rJacElem(1:nnzlocal,ithcopy)=rJacElemCopy(1:nnzlocal)
        iJacRow(1:neq,ithcopy)=iJacRowCopy(1:neq)
        nnzcopy(ithcopy)=nnzlocal
        !$omp  end critical
        OMPTimeLocalJac(ithcopy)=omp_get_wtime() - Timethread
        if (OMPVerbose.gt.1) write(*,*) OMPStamp,' Time in thread #', tid,':',OMPTimeLocalJac(ithcopy)
    enddo loopthread
    !$omp  END PARALLEL DO


    nnz(1:Nthreads)=nnzcopy(1:Nthreads) !nnzcopy is not necssary as nnz would be shared anyway in the parallel construct

    if (OMPDebug.gt.0) then
        write(iout,*) OMPStamp,' End of parallel loop....'
    endif

end subroutine OMPJacBuilder
#endif
!-------------------------------------------------------------------------------------------------
subroutine LocalJacBuilder(ivmin,ivmax,neq, t, yl,yldot00, ml, mu, wk,iJacCol,rJacElem,iJacRow,ith,nnz,nnzmxperthread,nthreads)

    ! ... Calculate Jacobian matrix (derivatives with respect to each
    ! ... Calculate Jacobian matrix (derivatives with respect to each
    !     dependent variable of the right-hand side of each rate equation).
    !     Lower and upper bandwidths are used to select for computation
    !     only those Jacobian elements that may be nonzero.
    !     Estimates of Jacobian elements are computed by finite differences.
    !     The Jacobian is stored in compressed sparse row format.
    use Output
    use Dim, only:nx
    use Timing,only:istimingon,ttjstor,ttotjf,ttimpjf
    use Math_problem_size,only:neqmx,numvar
    use Indexes,only:igyl,iseqalg,idxphi
    use Variable_perturbation,only:delperturb,dylconst,isjacreset
    use Jacobian_clipping,only:jaccliplim,istopjac,irstop,icstop
    use Ynorm,only:suscal,sfscal
    use UEpar,only:isphion,isnewpot,svrpkg,isbcwdt
    use Model_choice,only:iondenseqn
    use Bcond,only:isextrnpf,isextrtpf,isextrngc,isextrnw,isextrtw
    use Time_dep_nwt,only:nufak,dtreal,ylodt,dtuse,dtphi
    use OmpOptions,only:iidebugprint,ivdebugprint
    use Jacobian_csc,only:yldot_pert

    implicit none
    integer,intent(in):: ith,nnzmxperthread,nthreads,ivmin,ivmax,neq
    integer,intent(inout)::nnz
    real,intent(in):: t           ! physical time
    real,intent(inout) ::yl(*)       ! dependent variables
    integer,intent(inout)::iJacCol(nnzmxperthread),iJacRow(neq)
    real,intent(inout):: rJacElem(nnzmxperthread)
    real,intent(in) ::yldot00(neq+2) ! right-hand sides evaluated at yl
    integer,intent(in):: ml, mu   ! lower and upper bandwidths

    ! ... Work-array argument:
    real,intent(inout) :: wk(neq)     ! work space available to this subroutine

    ! ... Functions:
    logical ::tstguardc
    real(kind=4) gettime
    !     real(kind=4) ranf

    ! ... Local variables:
    real ::yold, dyl, jacelem
    real(kind=4) sec4, tsjstor, tsimpjf, dtimpjf,time0
    integer:: ii, iv, jv,ii1, ii2, xc, yc, ix, iy,omp_get_thread_num,tid
    !       time0=gettime(sec4) #here

    ! ... Only perturb variables that are being solved for (for Daspk option)
    !      if (iseqon(iv) .eq. 0) goto 18

    ! ... Set beginning and ending indices of right-hand sides that might be
    !     perturbed.

    nnz=1
    loopiv: do iv=ivmin,ivmax
        ii1 = max(iv-mu, 1)
        ii2 = min(iv+ml, neq)
        ! ... Reset range if this is a potential perturbation with isnewpot=1
        !         if (isphion*isnewpot.eq.1 .and. mod(iv,numvar).eq.0) then
        if (isphion*isnewpot.eq.1) then
            ii1 = max(iv-4*numvar*nx, 1)      ! 3*nx may be excessive
            ii2 = min(iv+4*numvar*nx, neq)    ! 3*nx may be excessive
        endif
        ! ... Reset range if extrapolation boundary conditions are used
        if (isextrnpf+isextrtpf+isextrngc+isextrnw+isextrtw.gt.0) then
            ii1 = max(iv-2*numvar*(nx+3), 1)      ! guess to include extrap. bc
            ii2 = min(iv+2*numvar*(nx+3), neq)    ! guess to include extrap. bc
        endif

        ! ... Initialize all of those right-hand sides to their unperturbed
        !     values.

        do ii = ii1, ii2   ! a) below wk is reset, but only over limited range
            wk(ii) = yldot00(ii)
        enddo

        ! ... Set spatial-location indices for this dependent variable.
        xc = igyl(iv,1)
        yc = igyl(iv,2)

        ! ... Save value of dependent variable, then perturb it.
        !     The perturbation to the variable is proportional to parameter
        !     del and to a measure of the size of the variable.  That measure
        !     increases with the absolute value of the variable if it exceeds
        !     the typical size given by dylconst/suscal but can never be less
        !     than that typical size.

        yold = yl(iv)
        dyl = delperturb * (abs(yold) + dylconst / suscal(iv))
        yl(iv) = yold + dyl

        !Calculate right-hand sides near location of perturbed variable.
        !call convsr_vo (xc, yc, yl)  ! test new convsr placement
        !call convsr_aux (xc, yc, yl) ! test new convsr placement
        call pandf1 (xc, yc, iv, neq, t, yl, wk)

        !yl(iv) = yold - dyl    ! for 2nd order Jac
        !call pandf1 (xc, yc, iv, neq, t, yl, yldot0) ! for 2nd order Jac

        !Calculate possibly nonzero Jacobian elements for this variable,
        !and store nonzero elements in compressed sparse column format.

        iJacRow(iv) = nnz      ! sets index for first Jac. elem. of var. iv
                                       ! the index is set to start at 0

        loopii: do ii = ii1, ii2
            jacelem = (wk(ii) - yldot00(ii)) / dyl
            !jacelem = (wk/(ii) - yldot0(ii)) / (2*dyl)  ! for 2nd order Jac

            !Add diagonal 1/dt for nksol
            ifdiagonal:if (((svrpkg.eq."nksol") .or. (svrpkg.eq."petsc")) .and. iv.eq.ii) then
                if (iseqalg(iv)*(1-isbcwdt).eq.0) then
                    jacelem = jacelem - 1/dtuse(iv)
                endif
                ix = igyl(iv,1)
                iy = igyl(iv,2)
                if (idxphi(ix,iy)==iv .and. dtphi<1e10) then ! selects phi eqn
                    jacelem = jacelem - 1/dtphi
                endif
            endif ifdiagonal

            ! ...  Add a pseudo timestep to the diagonal #! if eqn is not algebraic
            if (svrpkg .ne. "cvode" .and. nufak .gt. 0) then
                if (iv.eq.ii .and. yl(neq+1).eq.1) jacelem = jacelem - nufak  !omit .and. iseqalg(iv).eq.0)
            !     .                   jacelem = jacelem - nufak*suscal(iv)/sfscal(iv)
            endif

            ! Debug
            debug: if (ii==iidebugprint.and.iv==ivdebugprint) then
                tid=omp_get_thread_num()
                write(iout,'(a,i3,a,i7,i7,i7, E20.12,E20.12,E20.12,E20.12,E20.12,E20.12)') '#', tid,&
                    ' : ',iv,nnz,ii,jacelem,wk(ii),yldot00(ii),sfscal(iv),jaccliplim,nufak
                !call DebugHelper('dumpomp.txt') ! See at the end of this file
            endif debug

            storage: if (abs(jacelem*sfscal(iv)) .gt. jaccliplim) then
                if (nnz .gt. nnzmxperthread) then
                    write(iout,*) 'nnz=',nnz,'nnzmxperthread=',nnzmxperthread
                    write(iout,*) 'ith',ith,'*** jac_calc -- More storage needed for Jacobian. Increase lenpfac or omplenpfac.'
                    call xerrab("")
                endif
                !              if (rdoff.ne.0.e0) jacelem=jacelem*(1.0e0+ranf()*rdoff)
                iJacCol(nnz)=ii
                rJacElem(nnz)=jacelem
                nnz = nnz + 1
            endif storage

            check: if (istopjac.gt.0 .and. ii.eq.irstop .and. iv.eq.icstop) then
                !               yldot_pert(ii) = wk(ii)      ! for diagnostic only
                if (istopjac == 2) then
                    yl(iv) = yold
                    call pandf1 (xc, yc, iv, neq, t, yl, wk)
                endif
                !               call convsr_vo (xc, yc, yl)  ! was one call to convsr
                !               call convsr_aux (xc, yc, yl)
                call remark("***** non-zero jac_elem at irstop,icstop")
                write(iout,*) 'irstop = ', irstop, ', icstop = ', icstop
                call xerrab("")
            endif check

        enddo loopii

        ! ... Restore dependent variable and plasma variables near its location.
        yl(iv) = yold
        !         call convsr_vo (xc, yc, yl)  ! was one call to convsr
        !         call convsr_aux (xc, yc, yl)

        call pandf1 (xc, yc, iv, neq, t, yl, wk)

        ! 18   continue
        !...  If this is the last variable before jumping to new cell, reset pandf

        reset: if( isjacreset.ge.1) then
            yldot_pert(1:neq)=wk(1:neq)
            ! 18   continue
            !...  If this is the last variable before jumping to new cell, reset pandf
            !JG this call to pandf1 can be safely ignored with ijacreset=0 (and save some time...)
            if (mod(iv,numvar).eq.0) then
                call pandf1 (xc, yc, iv, neq, t, yl, wk)
            endif

            do ii=1,neq
                if (yldot_pert(ii).ne.wk(ii)) then
                    write(iout,'(a,i5,e20.12,e20.12)') ' *** wk modified on second call to pandf1 at ii=', ii,yldot_pert(ii),wk(ii)
                    call xerrab('*** Stop ***')
                endif
                if (isnan(yldot_pert(ii))) then
                    write(iout,*) 'NaN at ii=',ii
                    call xerrab('*** Nan in wk array in jac_calc ***')
                endif

            enddo
        endif reset
    enddo loopiv
! ... End loop over dependent variables and finish Jacobian storage.
end subroutine LocalJacBuilder
!-------------------------------------------------------------------------------------------------
#ifdef MPIJAC
subroutine InitMPI
    Use MpiOptions,only:Nprocs,ComSize,MPIRank,mpilenpfac,mpineq,MPIVerbose,nnzmxperproc,MPIstamp
    Use HybridOptions,only:Hybridstamp,HybridOMPMPI
    use mpi
    use Output
    use LSode,only:neq
    Use Jacobian,only:nnzmx
    implicit none
    integer(kind=4):: ierr
    integer(kind=4)::irank,icomsize
    character*8 :: MPIRankTag
    ! define rank and size

    ! check the size of the common world
    call MPI_COMM_SIZE(MPI_COMM_WORLD, iComSize, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, iRank, ierr)
    MPIrank=irank       ! irank is integer(kind=4)
    ComSize=iComSize !

    ! prepare MPI/Hybrid stamps for output
    write(MPIRankTag,'(I4)') MPIrank
    write(MPIstamp,'(a,a,a)') '[',trim(adjustl(trim(MPIRankTag))),'] MPIJac*'
    if (HybridOMPMPI>0) then
        write(Hybridstamp,'(a,a,a)') '[',trim(adjustl(trim(MPIRankTag))),'] HybridJac*'
    endif

    if (MPIVerbose>1) write(iout,*) MPIStamp,' MPI enabled with ComSize=',ComSize
    !
    if (Nprocs.gt.ComSize) then
        if (MPIVerbose.gt.1) write(iout,*) MPIStamp,' # processors requested > # processors available.',&
            'Nprocmax:',ComSize
        Nprocs=ComSize
        if (MPIVerbose.gt.0) write(iout,*) MPIStamp,' Resetting Nproc to Nprocmax: Nproc=',Nprocs
    endif
    if (Nprocs.le.0) then
        call xerrab('Nprocs must be >0')
    endif

    if (Nprocs.gt.1) then
        nnzmxperproc=ceiling(real(nnzmx)/real(Nprocs-1)*mpilenpfac)
    else
        nnzmxperproc=nnzmx
    endif
    if (MPIVerbose.gt.0) write(iout,'(a,a,i10)') MPIStamp,' nnzmxperproc:',nnzmxperproc
    mpineq=neq
    call gchange('MpiJacobian',0)
    if (MPIVerbose.gt.0) write(iout,'(a,a)') MPIStamp,'  MPI variables allocated'

end subroutine InitMPI
!-------------------------------------------------------------------------------------------------
subroutine jac_calc_mpi (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)

    !   Calculate Jacobian matrix (derivatives with respect to each
    !   dependent variable of the right-hand side of each rate equation).
    !   Lower and upper bandwidths are used to select for computation
    !   only those Jacobian elements that may be nonzero.
    !   Estimates of Jacobian elements are computed by finite differences.
    !   The Jacobian is stored in compressed sparse row format.
    use Output
    use Timing,only:istimingon,ttjstor,ttotjf,ttimpjf,ShowTime,&
    MPITotJacCalc,MPITotTimeCollect,MPITotTimeBuild
    use Grid,only:ngrid,ig,ijac,ijactot
    use Jacobian_csc,only:rcsc,jcsc,icsc,yldot_pert,yldot_unpt
    use mpi
    use MPIOptions,only:Nprocs,MPIVerbose,MPIDebug,nnzmxperproc,MPICheckNaN,MPIWriteJacobian,MPIRank,&
        MPILoadBalance,MPIAutoBalance,MPIStamp,MPIBalanceStrength
    use MPIJacobian,only:MPIivmin,MPIivmax,MPIiJacCol,MPIrJacElem,MPIiJacRow,MPILoadWeight,MPITimeLocalJac
    use UEpar, only: svrpkg
    implicit none
    ! ... Input arguments:
    integer,intent(in):: neq      !      total number of equations (all grid points)
    real,intent(in)   :: t              ! physical time
    real,intent(in)   :: yl(*)          ! dependent variables
    real,intent(in)   :: yldot00(neq+2) ! right-hand sides evaluated at yl
    integer,intent(in):: ml, mu         ! lower and upper bandwidths
    integer,intent(in):: nnzmx          ! maximum number of nonzeros in Jacobian

    ! ... Output arguments:
    real,intent(out)   :: jac(nnzmx)     ! nonzero Jacobian elements
    integer,intent(out):: ja(nnzmx)   ! col indices of nonzero Jacobian elements
    integer,intent(out):: ia(neq+1)   ! pointers to beginning of each row in jac,ja

    ! ... Work-array argument:
    real wk(neq)     ! work space available to this subroutine

    ! ... Functions
    logical tstguardc
    real(kind=4) gettime

    ! ... Local variables:
    real factor
    real(kind=4) sec4, tsjstor, tsimpjf, dtimpjf,time0,time1,MPITimeBuild,MPITimeCollect,MPITimeJacCalc
    integer:: i,iv,iproc,nnz
    integer (kind=4) :: ierr
    character(len = 80) ::  filename
    !integer::iJacCol(1:nnzmxperproc)
    !real ::rJacElem(1:nnzmxperproc)
    if (MPILoadBalance.ne.1 .and. MPIAutoBalance.ne.1) then
        MPILoadWeight(0:Nprocs-1)=1.0
    endif
    if (MPIAutoBalance.eq.1) then
    if (MPIBalanceStrength<=0) call xerrab('MPIBalanceStrength must be >0')
        !Check that time are not zero
        if (minval(MPITimeLocalJac(0:Nprocs-1)).gt.0.0) then
            do i=0,Nprocs-1
                factor=1/(MPITimeLocalJac(i)/sum(MPITimeLocalJac(0:Nprocs-1))*real(Nprocs))**MPIBalanceStrength
                MPILoadWeight(i)=MPILoadWeight(i)*factor
            enddo
        else
            MPILoadWeight(0:Nprocs-1)=1.0
        endif
    endif
    if (MPIDebug.gt.0) write(*,*) MPIStamp,' Starting jac_calc_mpi'
    !   Get the range of the iv index for each thread
    call MPISplitIndex(neq,Nprocs,MPIivmin,MPIivmax,MPILoadWeight)

    if (MPIVerbose.gt.1) then
        write(iout,*)MPIStamp,'neq=',neq
        write(iout,*)MPIStamp,' MPIivmin | MPIivmax | MPILoadWeight | MPITimeLocalJac'
        do iproc=0,Nprocs-1
    write(iout,'(a7,I3,a,I7,I7,f5.1,f8.3)') 'rank:', iproc,':',MPIivmin(iproc),MPIivmax(iproc),&
    MPILoadWeight(iproc),MPITimeLocalJac(iproc)
        enddo
    endif
    MPITimeJacCalc= gettime(sec4)

    !   Get initial value of system cpu timer.
    if (istimingon .eq. 1) tsjstor = gettime(sec4)

    !   Count Jacobian evaluations, both for total and for this case
    ijactot = ijactot + 1
    ijac(ig) = ijac(ig) + 1
    if (svrpkg.eq.'nksol') write(iout,*) ' Updating Jacobian, npe =  ',ijac(ig)

    !   Set up diagnostic arrays for debugging
    do iv = 1, neq
        yldot_unpt(iv) = yldot00(iv)  ! for diagnostic only
        yldot_pert(iv) = 0.
    enddo

    !   build jacobian ##############################################################
    MPITimeBuild=gettime(sec4)

    call MPIJacBuilder(neq, t, yl,yldot00, ml, mu,wk,MPIiJacCol,MPIrJacElem,MPIiJacRow,nnz)

    MPITimeBuild=gettime(sec4)-MPITimeBuild
    if (istimingon .eq. 1) MPITotTimebuild = MPITimeBuild+MPITotTimebuild
    if (MPIVerbose.gt.0) write(iout,*)MPIStamp,' Rank:',MPIRank ,'Time to build jac:',MPITimeBuild

    call MPICollectBroadCastTime(real(MPITimeBuild,kind=8))
    !   end build jacobian ##############################################################
    !   collect jacobian ##############################################################

    MPITimeCollect=gettime(sec4)
    call MPICollectBroadCastJacobian(MPIiJacRow,MPIiJacCol,MPIrJacElem,nnz)
    MPITimeCollect=gettime(sec4)-MPITimeCollect
    if (istimingon .eq. 1) MPITotTimeCollect = MPITimeCollect+MPITotTimeCollect
    if (MPIVerbose.gt.0) write(iout,*)MPIStamp,' Time to collect/broadcast jac:',MPITimeCollect
    !   end collect jacobian ##############################################################


    !   for Debug purpose
    if (MPIWriteJacobian.eq.1) then
        write(filename,'(a,3i3,a)') "jac_MPI_",ijac(ig),".txt"
        call jac_write(filename,neq, rcsc, icsc, jcsc)
    endif

    if (MPICheckNaN.gt.0) then
        do i=1,nnzmx
            if (isnan(rcsc(i))) then
                write(iout,*) 'rcsc is NaN at i=',i,rcsc(i)
            endif
            if (icsc(i).ne.icsc(i)) then
                write(iout,*) 'icsc is NaN at i=',i,icsc(i)
            endif
        enddo
        do i=1,neq+1
            if (jcsc(i).ne.jcsc(i)) then
                write(iout,*) 'jcsc is NaN at i=',i
            endif
        enddo
    endif

    !   Convert Jacobian from cMPIressed sparse column to compressedsparse row format.
    time1=gettime(sec4)
    call csrcsc (neq, 1, 1, rcsc, icsc, jcsc, jac, ja, ia)

    !! ... If desired, calculate Jacobian elements for ion continuity
    !!     equation by using INEL routines, and combine elements with
    !!     those calculated above.
    !      if (iondenseqn .eq. "inel") then
    !         call iondens2 (neq, nnzmx, jac, ja, ia, rcsc, icsc, jcsc)
    !      endif

    !! ... If necessary, calculate Jacobian elements for rows corresponding
    !!     to impurity-density equations (for cells other than guard cells),
    !!     and combine elements with those calculated above.
    !      if (isimpon .eq. 3 .or. isimpon .eq. 4) then
    !         if (istimingon .eq. 1) tsimpjf = gettime(sec4)
    !         call impdens2 (neq, nnzmx, jac, ja, ia, rcsc, icsc, jcsc)
    !         if (istimingon .eq. 1) then
    !            dtimpjf = gettime(sec4) - tsimpjf
    !            ttimpjf = ttimpjf + dtimpjf
    !            ttotjf = ttotjf + dtimpjf
    !         endif
    !      endif

    if (istimingon .eq. 1) ttjstor = ttjstor + gettime(sec4) - tsjstor
    MPITimeJacCalc=gettime(sec4)-MPITimeJacCalc
    if (istimingon .eq. 1) MPITotJacCalc = MPITimeJacCalc+MPITotJacCalc
    if (ShowTime.gt.0) write(iout,*)'**** Time in jac_calc:',MPITimeJacCalc
    call MPI_barrier(MPI_COMM_WORLD,ierr)
    return
end subroutine jac_calc_mpi
!-------------------------------------------------------------------------------------------------
subroutine MPICollectBroadCastTime(TimeLocal)
    use Output
    use mpi
    use MpiOptions,only:Nprocs,MpiRank,MPIVerbose,MPIDebug,ioutmpi,MPIStamp
    use MpiJacobian,only: MPITimeLocalJac

    implicit none
    integer(kind=4)::iproc
    integer(kind=4) :: ierr,req0
    real,intent(in) :: TimeLocal
    real(kind=4) gettime
    real(kind=4) sec4, TimeCollect
    TimeCollect=gettime(sec4)
    if (MPIRank==0) MPITimeLocalJac(0)=Timelocal
    loopproc:do iproc=1,Nprocs-1
        if (MPIRank.eq.iproc) then
            ! send to the master proc
            if (MPIDebug.gt.0) then
                write(ioutmpi,*) MPIStamp,' Rank',MPIRAnk,' sending data to 0'
            endif
            CALL MPI_SEND(TimeLocal,1,MPI_REAL8,0,7,MPI_COMM_WORLD,ierr)
            if (MPIDebug.gt.0) then
                write(ioutmpi,*) MPIStamp,' Rank',iproc,'has sent data to 0'
            endif
            ! collect on the master proc
        elseif (MPIRank.eq.0) then
            if (MPIDebug.gt.0) then
                write(ioutmpi,*) MPIStamp,' Rank 0 getting data from:',iproc
            endif
            CALL MPI_RECV(Timelocal,1,MPI_REAL8,iproc,7,MPI_COMM_WORLD,req0,ierr)
            MPITimeLocalJac(iproc)=Timelocal
        endif
    enddo loopproc
    call MPI_barrier(MPI_COMM_WORLD,ierr)
    TimeCollect=gettime(sec4)-TimeCollect
    if (MPIVerbose.gt.0) write(iout,*)MPIStamp,' Time to collect time:',TimeCollect
    TimeCollect=gettime(sec4)
    ! And now broadcast it to all the procs
    call MPI_BCAST(MPITimeLocalJac(0:Nprocs-1), int(Nprocs,kind=4), MPI_Real8, 0, MPI_COMM_WORLD, IERR)
    call MPI_barrier(MPI_COMM_WORLD,ierr)

    TimeCollect=gettime(sec4)-TimeCollect
    if (MPIVerbose.gt.0) write(iout,*)MPIStamp,' Time to broadcast time:',TimeCollect
end subroutine MPICollectBroadCastTime
!-----------------------------------------------------------------------
subroutine MPICollectBroadCastJacobian(iJacRow,iJacCol,rJacElem,nnz)
    use Output
    use mpi
    use MpiOptions,only:Nprocs,MpiRank,mpineq,nnzmxperproc,MPIVerbose,MPIDebug,ioutmpi,MPIStamp
    use MpiJacobian,only: MPIivmin,MPIivmax
    use Jacobian_csc,only:rcsc,jcsc,icsc
    use LSode,only:neq
    Use Jacobian,only:nnzmx

    implicit none
    integer:: nnzcumlocal(Nprocs),ith
    integer(kind=4)::iproc,nnzmxperproclocal,mpineqlocal
    integer ,intent(in):: nnz
    integer,intent(in)::iJacCol(nnzmxperproc)
    integer,intent(in):: iJacRow(mpineq+1)
    real,intent(in):: rJacElem(nnzmxperproc)
    integer(kind=4) :: ierr,req0,req1,req2,req3
    real :: TimeCollect
    real(kind=4) gettime
    real(kind=4) sec4, TimeJacCalc
    TimeCollect=gettime(sec4)
    nnzmxperproclocal=int(nnzmxperproc,kind=4)
    mpineqlocal=int(mpineq,kind=4)
    ! First we collect data from the master proc (Mpirank=0)
    if (MPIDebug.gt.0) write(ioutmpi,*) 'Rank',MPIrank,'nnz-1:',nnz-1
    if (MPIRank.eq.0) then
        nnzcumlocal(1)=nnz-1
        jcsc(MPIivmin(0):MPIivmax(0))= iJacRow(MPIivmin(0):MPIivmax(0))
        rcsc(1:nnz-1)= rJacElem(1:nnz-1)
        icsc(1:nnz-1)= iJacCol(1:nnz-1)
    endif
    ! Then we collect from each non-master process (MPIrank>0)

    loopproc:do iproc=1,Nprocs-1
        call MPI_barrier(MPI_COMM_WORLD,ierr)
        ith=iproc+1
        if (MPIRank.eq.iproc) then
            ! send to the master proc
            if (MPIDebug.gt.0) then
                write(ioutmpi,*) MPIStamp,' Rank',MPIRAnk,' sending data to 0'
            endif
            CALL MPI_SEND(nnz,1,MPI_INTEGER8,0,7,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(iJacRow,mpineqlocal,MPI_INTEGER8,0,9,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(iJacCol,nnzmxperproclocal,MPI_INTEGER8,0,10,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(rJacElem,nnzmxperproclocal,MPI_REAL8,0,11,MPI_COMM_WORLD,ierr)
            if (MPIDebug.gt.0) then
                write(ioutmpi,*) MPIStamp,' Rank',iproc,'has sent data to 0'
            endif
            ! collect on the master proc
        elseif (MPIRank.eq.0) then
            if (MPIDebug.gt.0) then
                write(ioutmpi,*) MPIStamp,' Rank 0 getting data from:',iproc
            endif
            CALL MPI_RECV(nnz,1,MPI_INTEGER8,iproc,7,MPI_COMM_WORLD,req0,ierr)
            CALL MPI_RECV(iJacRow,mpineqlocal,MPI_INTEGER8,iproc,9,MPI_COMM_WORLD,req1,ierr)
            CALL MPI_RECV(iJacCol,nnzmxperproclocal,MPI_INTEGER8,iproc,10,MPI_COMM_WORLD,req2,ierr)
            CALL MPI_RECV(rJacElem,nnzmxperproclocal,MPI_REAL8,iproc,11,MPI_COMM_WORLD,req3,ierr)
            if (MPIDebug.gt.0) then
                write(ioutmpi,*) MPIStamp,' Rank 0 got data from:',iproc, '; ivmin:ivmax=',MPIivmin(iproc),MPIivmax(iproc)
            endif
            CALL MPI_WAIT(req0, MPI_STATUS_IGNORE, ierr)
            CALL MPI_WAIT(req1, MPI_STATUS_IGNORE, ierr)
            CALL MPI_WAIT(req2, MPI_STATUS_IGNORE, ierr)
            CALL MPI_WAIT(req3, MPI_STATUS_IGNORE, ierr)
            if (MPIDebug.gt.0) then
                write(ioutmpi,*) MPIStamp,' Rank 0 going next'
            endif

            ! Now build the jacobian vector iteratively
            nnzcumlocal(ith)=nnzcumlocal(ith-1)+nnz-1
            if (MPIDebug.gt.0) then
                write(ioutmpi,*) MPIStamp,' Rank:',iproc,' sent: nnz:',nnz, ';nnzcum:',nnzcumlocal(ith)
            endif
            if (nnzcumlocal(ith).gt.nnzmx) then
                write(iout,*) 'nnzcum=',nnzcumlocal(ith),'nnzmx=',nnzmx
                call xerrab(' Problem: nnzcum > nnzmx... Increase lenpfac')
            endif
            jcsc(MPIivmin(iproc):MPIivmax(iproc))= iJacRow(MPIivmin(iproc):MPIivmax(iproc))+nnzcumlocal(ith-1)
            rcsc(nnzcumlocal(ith-1)+1:nnzcumlocal(ith))=rJacElem(1:nnz-1)
            icsc(nnzcumlocal(ith-1)+1:nnzcumlocal(ith))=iJacCol(1:nnz-1)

        endif
        if (MPIDebug.gt.0) then
            write(ioutmpi,*) MPIStamp,' Rank:',MPIRank,' at end of loop'
        endif
        call MPI_barrier(MPI_COMM_WORLD,ierr)
    enddo loopproc

    if (MPIRank.eq.0) then
        jcsc(neq+1) = nnzcumlocal(Nprocs)+1
        if (MPIVerbose.gt.0) write(iout,'(a,i9)') '**** Number of non-zero Jacobian elems:',nnzcumlocal(Nprocs)
    endif

    call MPI_barrier(MPI_COMM_WORLD,ierr)
    TimeCollect=gettime(sec4)-TimeCollect
    if (MPIVerbose.gt.0) write(iout,*)MPIStamp,' Time to collect jac:',TimeCollect
    TimeCollect=gettime(sec4)
    ! And now broadcast it to all the procs
    call MPI_BCAST(rcsc, nnzmx, MPI_REAL8, 0, MPI_COMM_WORLD, IERR)
    call MPI_BCAST(jcsc, mpineq+1, MPI_INTEGER8, 0, MPI_COMM_WORLD, IERR)
    call MPI_BCAST(icsc, nnzmx, MPI_INTEGER8, 0, MPI_COMM_WORLD, IERR)

    call MPI_barrier(MPI_COMM_WORLD,ierr)
    TimeCollect=gettime(sec4)-TimeCollect
    if (MPIVerbose.gt.0) write(iout,*)MPIStamp,' Time to broadcast jac:',TimeCollect
end subroutine MPICollectBroadCastJacobian
!-----------------------------------------------------------------------
subroutine MPIJacBuilder(neq, t, yl,yldot00, ml,mu,wk,iJacCol,rJacElem,iJacRow,nnz)
    use MPIOptions,only:MPIDebug,Nprocs,nnzmxperproc,MPIRank,ioutmpi,MPIStamp
    use MPIJacobian, only:MPIivmin,MPIivmax
    use Output
    Use Jacobian,only:nnzmx


    implicit none
    integer,intent(inout)::nnz
    integer,intent(in):: neq      ! total number of equations (all grid points)
    integer,intent(in):: ml, mu   ! lower and upper bandwidths
    real,intent(in):: t           ! physical time
    real,intent(in) ::yl(*)       ! dependent variables
    real,intent(in) ::yldot00(neq+2) ! right-hand sides evaluated at yl
    real,intent(inout) :: wk(neq)
    integer,intent(out)::iJacCol(nnzmxperproc)
    integer,intent(out):: iJacRow(neq)
    real,intent(out):: rJacElem(nnzmxperproc)
    integer :: iproc

    if (MPIDebug.gt.0) write(ioutmpi,*) MPIStamp,' Building jacobian on proc:',MPIrank

    call LocalJacBuilder(MPIivmin(MPIrank),MPIivmax(MPIrank),neq, t, yl,yldot00,ml,mu,wk,&
        iJacCol,rJacElem,iJacRow,MPIrank,nnz,nnzmxperproc,Nprocs)


end subroutine MPIJacBuilder
!-----------------------------------------------------------------------
subroutine jac_calc_hybrid (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)

    !   Calculate Jacobian matrix (derivatives with respect to each
    !   dependent variable of the right-hand side of each rate equation).
    !   Lower and upper bandwidths are used to select for computation
    !   only those Jacobian elements that may be nonzero.
    !   Estimates of Jacobian elements are computed by finite differences.
    !   The Jacobian is stored in compressed sparse row format.
    use Output
    use Timing,only:istimingon,ttjstor,ttotjf,ttimpjf,ShowTime,MPITotJacCalc,MPITotTimeCollect,&
    OMPTotTimeCollect,MPITotTimeBuild,OMPTotTimeBuild,OMPTotJacCalc
    use Grid,only:ngrid,ig,ijac,ijactot
    use Jacobian_csc,only:rcsc,jcsc,icsc,yldot_pert,yldot_unpt
    use mpi
    use MPIOptions,only:Nprocs,nnzmxperproc,MPICheckNaN,MPIWriteJacobian,MPIRank,&
        MPILoadBalance,MPIAutoBalance,MPIBalanceStrength
    use OMPOptions,only:Nthreads,OMPLoadBalance,OMPAutoBalance,OMPBalanceStrength
    use MPIJacobian,only:MPIivmin,MPIivmax,MPIiJacCol,MPIrJacElem,MPIiJacRow,MPILoadWeight,MPITimeLocalJac
    use OmpJacobian,only:iJacCol,rJacElem,iJacRow,OMPivmin,OMPivmax,nnz,OMPLoadWeight,OMPTimeLocalJac
    use HybridOptions,only:HybridDebug,HybridVerbose,HybridCheckNaN,Hybridstamp
    use UEpar, only: svrpkg
    implicit none
    ! ... Input arguments:
    integer,intent(in):: neq      !      total number of equations (all grid points)
    real,intent(in)   :: t              ! physical time
    real,intent(in)   :: yl(*)          ! dependent variables
    real,intent(in)   :: yldot00(neq+2) ! right-hand sides evaluated at yl
    integer,intent(in):: ml, mu         ! lower and upper bandwidths
    integer,intent(in):: nnzmx          ! maximum number of nonzeros in Jacobian

    ! ... Output arguments:
    real,intent(out)   :: jac(nnzmx)     ! nonzero Jacobian elements
    integer,intent(out):: ja(nnzmx)   ! col indices of nonzero Jacobian elements
    integer,intent(out):: ia(neq+1)   ! pointers to beginning of each row in jac,ja

    ! ... Work-array argument:
    real wk(neq)     ! work space available to this subroutine

    ! ... Functions
    logical tstguardc
    real(kind=4) gettime

    ! ... Local variables:
    real factor
    real(kind=4)        :: sec4, tsjstor, tsimpjf, dtimpjf,time0,time1
    real(kind=4)        :: OMPTimeBuild,MPITimeJacCalc,OMPTimeCollect,MPITimeBuild,MPITimeCollect
    integer             :: i,iv,iproc,ithread,nnzcumout
    integer(kind=4)     :: ierr
    character(len = 80) ::  filename

    if (MPILoadBalance.ne.1 .and. MPIAutoBalance.ne.1) then
        MPILoadWeight(0:Nprocs-1)=1.0
    endif
    if (MPIAutoBalance.eq.1) then
        if (MPIBalanceStrength<=0) call xerrab('MPIBalanceStrength must be >0')
        !Check that time are not zero
        if (minval(MPITimeLocalJac(0:Nprocs-1)).gt.0.0) then
            do i=0,Nprocs-1
                factor=1/(MPITimeLocalJac(i)/sum(MPITimeLocalJac(0:Nprocs-1))*real(Nprocs))**MPIBalanceStrength
                MPILoadWeight(i)=MPILoadWeight(i)*factor
            enddo
        else
            MPILoadWeight(0:Nprocs-1)=1.0
        endif
    endif

    if (OMPLoadBalance.ne.1 .and. OMPAutoBalance.ne.1) then
        OMPLoadWeight(1:Nthreads)=1.0
    endif
    if (OMPAutoBalance.eq.1) then
    if (OMPBalanceStrength<=0) call xerrab('OMPBalanceStrength must be >0')
        !Check that time are not zero
        if (minval(OMPTimeLocalJac).gt.0.0) then
            do i=1,Nthreads
                factor=1/(OMPTimeLocalJac(i)/sum(OMPTimeLocalJac)*real(Nthreads))**OMPBalanceStrength
                OMPLoadWeight(i)=OMPLoadWeight(i)*factor
            enddo
        else
            OMPLoadWeight(1:Nthreads)=1.0
        endif
    endif


    if (HybridDebug.gt.0) write(*,*) Hybridstamp,' Starting jac_calc_hybrid'
    !   Get the range of the iv index for each thread
    call MPISplitIndex(neq,Nprocs,MPIivmin,MPIivmax,MPILoadWeight)

    call OMPSplitIndex(MPIivmin(MPIRank),MPIivmax(MPIRank),Nthreads,OMPivmin,OMPivmax,OMPLoadWeight)

    if (HybridVerbose.gt.1) then
        write(iout,*)HybridStamp, 'neq=',neq
        write(iout,*)HybridStamp,'| MPIivmin MPIivmax MPILoadWeight MPITimeJac | OMPivmin OMPivmax OMPLoadWeight '
        do iproc=0,Nprocs-1
            do ithread=1,Nthreads
                if (ithread==1) then
                    if (iproc==MPIRank) then
                    write(iout,'(a6,I3,a7,I3,a3,I10,I10,f10.1,f10.3,a3,I8,I8,f8.1)') 'rank', iproc,'thread', ithread,'|',&
                        MPIivmin(iproc),MPIivmax(iproc),MPILoadWeight(iproc),MPITimeLocalJac(iproc),&
                        '| ',OMPivmin(ithread),OMPivmax(ithread),OMPLoadWeight(ithread)
                    else
                    write(iout,'(a6,I3,a7,I3,a3,I10,I10,f10.1,f10.3,a3,a8,a8,a8)') 'rank', iproc,'thread', ithread,'|',&
                        MPIivmin(iproc),MPIivmax(iproc),MPILoadWeight(iproc),MPITimeLocalJac(iproc),&
                        '| ','-','-','-'
                        endif
                else
                if (HybridVerbose.gt.3) then
                    write(iout,'(a6,a3,a7,I3,a3,a10,a10,a10,a10,a3,I8,I8,f8.1)') ' ', ' ',  'thread', ithread,'|',&
                        ' ',' ',' ',' ',' | ',OMPivmin(ithread),OMPivmax(ithread),OMPLoadWeight(ithread)
                endif
                endif
            enddo
        enddo
    endif


    MPITimeJacCalc= gettime(sec4)

    !   Get initial value of system cpu timer.
    if (istimingon .eq. 1) tsjstor = gettime(sec4)

    !   Count Jacobian evaluations, both for total and for this case
    ijactot = ijactot + 1
    ijac(ig) = ijac(ig) + 1
    if (svrpkg.eq.'nksol') write(iout,*) ' Updating Jacobian, npe =  ',ijac(ig)

    !   Set up diagnostic arrays for debugging
    do iv = 1, neq
        yldot_unpt(iv) = yldot00(iv)  ! for diagnostic only
        yldot_pert(iv) = 0.
    enddo

    !   build jacobian #####################################################################
    OMPTimeBuild=gettime(sec4)
    call OMPJacBuilder(neq, t, yl,yldot00, ml,mu,wk,iJacCol,rJacElem,iJacRow,nnz)
    OMPTimeBuild=gettime(sec4)-OMPTimeBuild
    if (istimingon .eq. 1) OMPTotTimebuild = OMPTimeBuild+OMPTotTimebuild
    if (HybridVerbose.gt.2) write(iout,*)Hybridstamp,' Time to build OMP jac:',OMPTimeBuild

    ! Collect OMP part of the Jacobian ######################################################
    OMPTimeCollect=gettime(sec4)
    call OMPCollectJacobian(neq,nnzmxperproc,MPIrJacElem,MPIiJacCol,MPIiJacRow,nnzcumout)
    OMPTimeCollect=gettime(sec4)-OMPTimeCollect
    if (istimingon .eq. 1) OMPTotTimeCollect = OMPTimeCollect+OMPTotTimeCollect
    if (HybridVerbose.gt.2) write(iout,*)Hybridstamp,' Time to collect OMP jac:',OMPTimeCollect

    ! broadcast  MPItimebuild for load balancing #############################################
    MPITimeBuild=OMPTimeCollect+OMPTimeBuild
    if (istimingon .eq. 1) MPITotTimebuild = MPITimeBuild+MPITotTimebuild
    if (HybridVerbose.gt.2) write(iout,*)Hybridstamp,' Time to build MPI jac:',MPITimeBuild
    call MPICollectBroadCastTime(real(MPITimeBuild,kind=8))


    !   collect MPI jacobian    ##############################################################
    MPITimeCollect=gettime(sec4)
    !nnzcumout is the total amount of non-zero Jac elements stored in iJacCol and rJacElem
    ! but nnz=nnzcumout+1 is passed in the collector
    call MPICollectBroadCastJacobian(MPIiJacRow,MPIiJacCol,MPIrJacElem,nnzcumout+1)
    MPITimeCollect=gettime(sec4)-MPITimeCollect
    if (istimingon .eq. 1) MPITotTimeCollect = MPITimeBuild+MPITotTimeCollect
    if (HybridVerbose.gt.2) write(iout,*)Hybridstamp,' MPI: Time to collect/broadcast jac:',MPITimeCollect


    !   for Debug purpose
    if (MPIWriteJacobian.eq.1) then
        write(filename,'(a,3i3,a)') "jac_hybrid_",ijac(ig),".txt"
        call jac_write(filename,neq, rcsc, icsc, jcsc)
    endif

    if (HybridCheckNaN.gt.0) then
        do i=1,nnzmx
            if (isnan(rcsc(i))) then
                write(iout,*) 'rcsc is NaN at i=',i,rcsc(i)
            endif
            if (icsc(i).ne.icsc(i)) then
                write(iout,*) 'icsc is NaN at i=',i,icsc(i)
            endif
        enddo
        do i=1,neq+1
            if (jcsc(i).ne.jcsc(i)) then
                write(iout,*) 'jcsc is NaN at i=',i
            endif
        enddo
    endif

    !   Convert Jacobian from cMPIressed sparse column to compressedsparse row format.
    time1=gettime(sec4)
    call csrcsc (neq, 1, 1, rcsc, icsc, jcsc, jac, ja, ia)

    !! ... If desired, calculate Jacobian elements for ion continuity
    !!     equation by using INEL routines, and combine elements with
    !!     those calculated above.
    !      if (iondenseqn .eq. "inel") then
    !         call iondens2 (neq, nnzmx, jac, ja, ia, rcsc, icsc, jcsc)
    !      endif

    !! ... If necessary, calculate Jacobian elements for rows corresponding
    !!     to impurity-density equations (for cells other than guard cells),
    !!     and combine elements with those calculated above.
    !      if (isimpon .eq. 3 .or. isimpon .eq. 4) then
    !         if (istimingon .eq. 1) tsimpjf = gettime(sec4)
    !         call impdens2 (neq, nnzmx, jac, ja, ia, rcsc, icsc, jcsc)
    !         if (istimingon .eq. 1) then
    !            dtimpjf = gettime(sec4) - tsimpjf
    !            ttimpjf = ttimpjf + dtimpjf
    !            ttotjf = ttotjf + dtimpjf
    !         endif
    !      endif

    if (istimingon .eq. 1) ttjstor = ttjstor + gettime(sec4) - tsjstor
    MPITimeJacCalc=gettime(sec4)-MPITimeJacCalc
    if (istimingon .eq. 1) MPITotJacCalc = MPITimeJacCalc+MPITotJacCalc
    if (ShowTime.gt.0) write(iout,*)HybridStamp,'Time in jac_calc:',MPITimeJacCalc

    return
end subroutine jac_calc_hybrid
#endif
!-------------------------------------------------------------------------------------------------
subroutine MPISplitIndex(neq,Nprocs,ivmin,ivmax,weight)
    implicit none
    integer,intent(in) ::neq,Nprocs
    real,intent(inout) :: weight(0:Nprocs-1)
    integer,intent(out)::ivmin(0:Nprocs-1),ivmax(0:Nprocs-1)
    integer:: Nsize(0:Nprocs-1),imax,i

    if (Nprocs.gt.1) then
        do i=0,Nprocs-1
                if (weight(i)<=0) call xerrab('MPISplitIndex: weight <0')
        enddo
        ! Normalized weights
        weight(0:Nprocs-1)=weight(0:Nprocs-1)/sum(weight(0:Nprocs-1))*real(Nprocs)
        do i=0,Nprocs-1
            Nsize(i)=int(real(neq/Nprocs)*weight(i))
        enddo

        do i=0,Nprocs-1
            if (Nsize(i)<0) call xerrab('Nsize<0')
            if (Nsize(i)<2) Nsize(i)=Nsize(i)+1
        enddo

        if (neq.ne.sum(Nsize(0:Nprocs-1))) then
            imax=0
            do i=1,Nprocs-1
                if (Nsize(i)>Nsize(i-1)) then
                    imax=i
                endif
            enddo
            Nsize(imax) = Nsize(imax) + (neq-sum(Nsize(0:Nprocs-1)))
        endif
        if (neq.ne.sum(Nsize(0:Nprocs-1))) call xerrab('Nsize .ne. neq!!!')

        ivmin(0)=1
        ivmax(0)=1+Nsize(0)-1
        do i=1,Nprocs-1
            ivmin(i)=ivmax(i-1)+1
            ivmax(i)=ivmin(i)+Nsize(i)-1
        enddo
        if (ivmax(Nprocs-1)-ivmin(0)+1.ne.neq) call xerrab('ivmax(Nprocs-1)!=neq')
    else
        ivmin(0)=1
        ivmax(0)=neq
    endif

end subroutine MPISplitIndex
!-------------------------------------------------------------------------------------------------
subroutine OMPSplitIndex(ieqmin,ieqmax,Nthreads,ivmin,ivmax,weight)
    implicit none
    integer,intent(in) ::ieqmin,ieqmax,Nthreads
    real::weight(Nthreads)
    integer,intent(out)::ivmin(Nthreads),ivmax(Nthreads)
    integer:: Nsize(Nthreads),Msize,R,i,imax
    if (ieqmax-ieqmin+1<2) call xerrab('Number of equations to solve <2')
    if (Nthreads.gt.1) then
        !        if (OMPLoadWeight.eq.1) then

        do i=1,Nthreads
            if (weight(i)<=0) call xerrab('OMPSplitIndex: weight <0')
            !write(*,*) weight(i)
        enddo

        ! Normalized weights
        weight(1:Nthreads)=weight(1:Nthreads)/sum(weight(1:Nthreads))*real(Nthreads)
        do i=1,Nthreads
            Nsize(i)=int(real((ieqmax-ieqmin+1)/Nthreads)*weight(i))
            !write(*,*) Nsize(i),weight(i)
        enddo

        do i=1,Nthreads
            if (Nsize(i)<0) call xerrab('Nsize<0')
            if (Nsize(i)<2) Nsize(i)=Nsize(i)+1
        enddo
        if (ieqmax-ieqmin+1.ne.sum(Nsize)) then
            imax=1
            do i=2,Nthreads
                if (Nsize(i)>Nsize(i-1)) then
                    imax=i
                endif
            enddo
            Nsize(imax) = Nsize(imax) + ((ieqmax-ieqmin+1)-sum(Nsize))
        endif
        !write(*,*) Nsize,neq,sum(Nsize)
        if (ieqmax-ieqmin+1.ne.sum(Nsize)) call xerrab('Nsize .ne. neq!!!')
        ivmin(1)=ieqmin
        ivmax(1)=ieqmin+Nsize(1)-1
        do i=2,Nthreads
            ivmin(i)=ivmax(i-1)+1
            ivmax(i)=ivmin(i)+Nsize(i)-1
        enddo
        if (ivmax(Nthreads)-ivmin(1)+1.ne.(ieqmax-ieqmin+1)) call xerrab('ivmax(Nthreads)!=neq')
    else
        ivmin(Nthreads)=ieqmin
        ivmax(Nthreads)=ieqmax
    endif

end subroutine OMPSplitIndex
!-----------------------------------------------------------------------
! The routines below are for debugging of OMP implementation
subroutine WriteArrayReal(array,s,iu)
    implicit none
    real:: array(*)
    integer:: i,s,iu
    do i=1,s
        write(iu,*) array(i)
    enddo
end subroutine WriteArrayReal

subroutine WriteArrayInteger(array,s,iu)
    implicit none
    integer:: array(*)
    integer:: i,s,iu
    do i=1,s
        write(iu,*) array(i)
    enddo
end subroutine WriteArrayInteger
!JG 03/2020 Subroutine for debugging purpose: print out all the threadprivate variable.
!JG This routine was created with the OMPDebugger.py script, which also contains methods to analyze the outpout of this routine.
!JG This routine can be called within the jacobian constructor to determine which variables behave differenlty in serial vs openmp jacobian calculations
subroutine DebugHelper(FileName)
    Use Bcond
    Use Cfric
    Use Coefeq
    Use Comflo
    Use Comgeo
    Use Compla
    Use Comtra
    Use Conduc
    Use Dim
    Use Gradients
    Use Imprad
    Use Locflux
    Use MCN_sources
    Use PNC_params
    Use Rhsides
    Use Save_terms
    Use Selec
    Use Time_dep_nwt
    Use Timing
    Use UEpar
    Use Wkspace

    implicit none
    integer:: iunit
    character(len = *) ::  filename

    open (newunit = iunit, file = trim(filename))
    write(iunit,*) "alfneo"
    call WriteArrayReal(alfneo,size(alfneo),iunit)
    write(iunit,*) "betap"
    call WriteArrayReal(betap,size(betap),iunit)
    write(iunit,*) "cfneut"
    write(iunit,*) cfneut
    write(iunit,*) "cfneutdiv"
    write(iunit,*) cfneutdiv
    write(iunit,*) "cfvcsx"
    call WriteArrayReal(cfvcsx,size(cfvcsx),iunit)
    write(iunit,*) "cfvcsy"
    call WriteArrayReal(cfvcsy,size(cfvcsy),iunit)
    write(iunit,*) "cfvgpx"
    call WriteArrayReal(cfvgpx,size(cfvgpx),iunit)
    write(iunit,*) "cfvgpy"
    call WriteArrayReal(cfvgpy,size(cfvgpy),iunit)
    write(iunit,*) "cmneut"
    write(iunit,*) cmneut
    write(iunit,*) "cmneutdiv"
    write(iunit,*) cmneutdiv
    write(iunit,*) "coef1"
    write(iunit,*) coef1
    write(iunit,*) "coll_fe"
    call WriteArrayReal(coll_fe,size(coll_fe),iunit)
    write(iunit,*) "coll_fi"
    call WriteArrayReal(coll_fi,size(coll_fi),iunit)
    write(iunit,*) "conx"
    call WriteArrayReal(conx,size(conx),iunit)
    write(iunit,*) "conxe"
    call WriteArrayReal(conxe,size(conxe),iunit)
    write(iunit,*) "conxi"
    call WriteArrayReal(conxi,size(conxi),iunit)
    write(iunit,*) "cony"
    call WriteArrayReal(cony,size(cony),iunit)
    write(iunit,*) "conye"
    call WriteArrayReal(conye,size(conye),iunit)
    write(iunit,*) "conyi"
    call WriteArrayReal(conyi,size(conyi),iunit)
    write(iunit,*) "cs"
    write(iunit,*) cs
    write(iunit,*) "csh"
    write(iunit,*) csh
    write(iunit,*) "ctaue"
    write(iunit,*) ctaue
    write(iunit,*) "ctaui"
    write(iunit,*) ctaui
    write(iunit,*) "dclass_e"
    call WriteArrayReal(dclass_e,size(dclass_e),iunit)
    write(iunit,*) "dclass_i"
    call WriteArrayReal(dclass_i,size(dclass_i),iunit)
    write(iunit,*) "dif2_use"
    call WriteArrayReal(dif2_use,size(dif2_use),iunit)
    write(iunit,*) "dif_use"
    call WriteArrayReal(dif_use,size(dif_use),iunit)
    write(iunit,*) "diffusivwrk"
    call WriteArrayReal(diffusivwrk,size(diffusivwrk),iunit)
    write(iunit,*) "difp_use"
    call WriteArrayReal(difp_use,size(difp_use),iunit)
    write(iunit,*) "dp1"
    write(iunit,*) dp1
    write(iunit,*) "dtold"
    write(iunit,*) dtold
    write(iunit,*) "dtreal"
    write(iunit,*) dtreal
    write(iunit,*) "dutm_use"
    call WriteArrayReal(dutm_use,size(dutm_use),iunit)
    write(iunit,*) "eeli"
    call WriteArrayReal(eeli,size(eeli),iunit)
    write(iunit,*) "eqp"
    call WriteArrayReal(eqp,size(eqp),iunit)
    write(iunit,*) "eqpg"
    call WriteArrayReal(eqpg,size(eqpg),iunit)
    write(iunit,*) "erliz"
    call WriteArrayReal(erliz,size(erliz),iunit)
    write(iunit,*) "erlrc"
    call WriteArrayReal(erlrc,size(erlrc),iunit)
    write(iunit,*) "eta1"
    call WriteArrayReal(eta1,size(eta1),iunit)
    write(iunit,*) "ex"
    call WriteArrayReal(ex,size(ex),iunit)
    write(iunit,*) "ey"
    call WriteArrayReal(ey,size(ey),iunit)
    write(iunit,*) "fcdif"
    write(iunit,*) fcdif
    write(iunit,*) "fdiaxlb"
    call WriteArrayReal(fdiaxlb,size(fdiaxlb),iunit)
    write(iunit,*) "fdiaxrb"
    call WriteArrayReal(fdiaxrb,size(fdiaxrb),iunit)
    write(iunit,*) "feex"
    call WriteArrayReal(feex,size(feex),iunit)
    write(iunit,*) "feexy"
    call WriteArrayReal(feexy,size(feexy),iunit)
    write(iunit,*) "feey"
    call WriteArrayReal(feey,size(feey),iunit)
    write(iunit,*) "feey4ord"
    call WriteArrayReal(feey4ord,size(feey4ord),iunit)
    write(iunit,*) "feeycbo"
    call WriteArrayReal(feeycbo,size(feeycbo),iunit)
    write(iunit,*) "feix"
    call WriteArrayReal(feix,size(feix),iunit)
    write(iunit,*) "feixy"
    call WriteArrayReal(feixy,size(feixy),iunit)
    write(iunit,*) "feiy"
    call WriteArrayReal(feiy,size(feiy),iunit)
    write(iunit,*) "feiy4ord"
    call WriteArrayReal(feiy4ord,size(feiy4ord),iunit)
    write(iunit,*) "feiycbo"
    call WriteArrayReal(feiycbo,size(feiycbo),iunit)
    write(iunit,*) "flox"
    call WriteArrayReal(flox,size(flox),iunit)
    write(iunit,*) "floxe"
    call WriteArrayReal(floxe,size(floxe),iunit)
    write(iunit,*) "floxebgt"
    call WriteArrayReal(floxebgt,size(floxebgt),iunit)
    write(iunit,*) "floxi"
    call WriteArrayReal(floxi,size(floxi),iunit)
    write(iunit,*) "floxibgt"
    call WriteArrayReal(floxibgt,size(floxibgt),iunit)
    write(iunit,*) "floy"
    call WriteArrayReal(floy,size(floy),iunit)
    write(iunit,*) "floye"
    call WriteArrayReal(floye,size(floye),iunit)
    write(iunit,*) "floyi"
    call WriteArrayReal(floyi,size(floyi),iunit)
    write(iunit,*) "fmihxpt"
    call WriteArrayReal(fmihxpt,size(fmihxpt),iunit)
    write(iunit,*) "fmivxpt"
    call WriteArrayReal(fmivxpt,size(fmivxpt),iunit)
    write(iunit,*) "fmix"
    call WriteArrayReal(fmix,size(fmix),iunit)
    write(iunit,*) "fmixy"
    call WriteArrayReal(fmixy,size(fmixy),iunit)
    write(iunit,*) "fmiy"
    call WriteArrayReal(fmiy,size(fmiy),iunit)
    write(iunit,*) "fnix"
    call WriteArrayReal(fnix,size(fnix),iunit)
    write(iunit,*) "fnixcb"
    call WriteArrayReal(fnixcb,size(fnixcb),iunit)
    write(iunit,*) "fniy"
    call WriteArrayReal(fniy,size(fniy),iunit)
    write(iunit,*) "fniy4ord"
    call WriteArrayReal(fniy4ord,size(fniy4ord),iunit)
    write(iunit,*) "fniycb"
    call WriteArrayReal(fniycb,size(fniycb),iunit)
    write(iunit,*) "fniycbo"
    call WriteArrayReal(fniycbo,size(fniycbo),iunit)
    write(iunit,*) "fqp"
    call WriteArrayReal(fqp,size(fqp),iunit)
    write(iunit,*) "frice"
    call WriteArrayReal(frice,size(frice),iunit)
    write(iunit,*) "frici"
    call WriteArrayReal(frici,size(frici),iunit)
    write(iunit,*) "fricnrl"
    call WriteArrayReal(fricnrl,size(fricnrl),iunit)
    write(iunit,*) "fxe"
    write(iunit,*) fxe
    write(iunit,*) "fxi"
    write(iunit,*) fxi
    write(iunit,*) "ghxpt"
    write(iunit,*) ghxpt
    write(iunit,*) "gpex"
    call WriteArrayReal(gpex,size(gpex),iunit)
    write(iunit,*) "gpey"
    call WriteArrayReal(gpey,size(gpey),iunit)
    write(iunit,*) "gpix"
    call WriteArrayReal(gpix,size(gpix),iunit)
    write(iunit,*) "gpiy"
    call WriteArrayReal(gpiy,size(gpiy),iunit)
    write(iunit,*) "gprx"
    call WriteArrayReal(gprx,size(gprx),iunit)
    write(iunit,*) "gpry"
    call WriteArrayReal(gpry,size(gpry),iunit)
    write(iunit,*) "gtex"
    call WriteArrayReal(gtex,size(gtex),iunit)
    write(iunit,*) "gtey"
    call WriteArrayReal(gtey,size(gtey),iunit)
    write(iunit,*) "gtix"
    call WriteArrayReal(gtix,size(gtix),iunit)
    write(iunit,*) "gtiy"
    call WriteArrayReal(gtiy,size(gtiy),iunit)
    write(iunit,*) "gvxpt"
    write(iunit,*) gvxpt
    write(iunit,*) "hcxe"
    call WriteArrayReal(hcxe,size(hcxe),iunit)
    write(iunit,*) "hcxg"
    call WriteArrayReal(hcxg,size(hcxg),iunit)
    write(iunit,*) "hcxi"
    call WriteArrayReal(hcxi,size(hcxi),iunit)
    write(iunit,*) "hcxij"
    call WriteArrayReal(hcxij,size(hcxij),iunit)
    write(iunit,*) "hcxineo"
    call WriteArrayReal(hcxineo,size(hcxineo),iunit)
    write(iunit,*) "hcxn"
    call WriteArrayReal(hcxn,size(hcxn),iunit)
    write(iunit,*) "hcye"
    call WriteArrayReal(hcye,size(hcye),iunit)
    write(iunit,*) "hcyg"
    call WriteArrayReal(hcyg,size(hcyg),iunit)
    write(iunit,*) "hcyi"
    call WriteArrayReal(hcyi,size(hcyi),iunit)
    write(iunit,*) "hcyij"
    call WriteArrayReal(hcyij,size(hcyij),iunit)
    write(iunit,*) "hcyn"
    call WriteArrayReal(hcyn,size(hcyn),iunit)
    write(iunit,*) "i1"
    write(iunit,*) i1
    write(iunit,*) "i2"
    write(iunit,*) i2
    write(iunit,*) "i2p"
    write(iunit,*) i2p
    write(iunit,*) "i3"
    write(iunit,*) i3
    write(iunit,*) "i4"
    write(iunit,*) i4
    write(iunit,*) "i5"
    write(iunit,*) i5
    write(iunit,*) "i5m"
    write(iunit,*) i5m
    write(iunit,*) "i6"
    write(iunit,*) i6
    write(iunit,*) "i7"
    write(iunit,*) i7
    write(iunit,*) "i8"
    write(iunit,*) i8
    write(iunit,*) "ixf6"
    write(iunit,*) ixf6
    write(iunit,*) "ixs1"
    write(iunit,*) ixs1
    write(iunit,*) "iyf6"
    write(iunit,*) iyf6
    write(iunit,*) "iys1"
    write(iunit,*) iys1
    write(iunit,*) "j1"
    write(iunit,*) j1
    write(iunit,*) "j1p"
    write(iunit,*) j1p
    write(iunit,*) "j2"
    write(iunit,*) j2
    write(iunit,*) "j2p"
    write(iunit,*) j2p
    write(iunit,*) "j3"
    write(iunit,*) j3
    write(iunit,*) "j4"
    write(iunit,*) j4
    write(iunit,*) "j5"
    write(iunit,*) j5
    write(iunit,*) "j5m"
    write(iunit,*) j5m
    write(iunit,*) "j5p"
    write(iunit,*) j5p
    write(iunit,*) "j6"
    write(iunit,*) j6
    write(iunit,*) "j6p"
    write(iunit,*) j6p
    write(iunit,*) "j7"
    write(iunit,*) j7
    write(iunit,*) "j8"
    write(iunit,*) j8
    write(iunit,*) "k2neo"
    call WriteArrayReal(k2neo,size(k2neo),iunit)
    write(iunit,*) "ktneo"
    call WriteArrayReal(ktneo,size(ktneo),iunit)
    write(iunit,*) "kxbohm"
    call WriteArrayReal(kxbohm,size(kxbohm),iunit)
    write(iunit,*) "kxe_use"
    call WriteArrayReal(kxe_use,size(kxe_use),iunit)
    write(iunit,*) "kxi_use"
    call WriteArrayReal(kxi_use,size(kxi_use),iunit)
    write(iunit,*) "kybohm"
    call WriteArrayReal(kybohm,size(kybohm),iunit)
    write(iunit,*) "kye_use"
    call WriteArrayReal(kye_use,size(kye_use),iunit)
    write(iunit,*) "kyi_use"
    call WriteArrayReal(kyi_use,size(kyi_use),iunit)
    write(iunit,*) "lng"
    call WriteArrayReal(lng,size(lng),iunit)
    write(iunit,*) "mfl"
    write(iunit,*) mfl
    write(iunit,*) "msh"
    write(iunit,*) msh
    write(iunit,*) "msor"
    call WriteArrayReal(msor,size(msor),iunit)
    write(iunit,*) "msorold"
    call WriteArrayReal(msorold,size(msorold),iunit)
    write(iunit,*) "msorxr"
    call WriteArrayReal(msorxr,size(msorxr),iunit)
    write(iunit,*) "msorxrold"
    call WriteArrayReal(msorxrold,size(msorxrold),iunit)
    write(iunit,*) "na"
    call WriteArrayReal(na,size(na),iunit)
    write(iunit,*) "ne"
    call WriteArrayReal(ne,size(ne),iunit)
    write(iunit,*) "ney0"
    call WriteArrayReal(ney0,size(ney0),iunit)
    write(iunit,*) "ney1"
    call WriteArrayReal(ney1,size(ney1),iunit)
    write(iunit,*) "nfsp"
    write(iunit,*) nfsp
    write(iunit,*) "ng"
    call WriteArrayReal(ng,size(ng),iunit)
    write(iunit,*) "ngy0"
    call WriteArrayReal(ngy0,size(ngy0),iunit)
    write(iunit,*) "ngy1"
    call WriteArrayReal(ngy1,size(ngy1),iunit)
    write(iunit,*) "ni"
    call WriteArrayReal(ni,size(ni),iunit)
    write(iunit,*) "nit"
    call WriteArrayReal(nit,size(nit),iunit)
    write(iunit,*) "nity0"
    call WriteArrayReal(nity0,size(nity0),iunit)
    write(iunit,*) "nity1"
    call WriteArrayReal(nity1,size(nity1),iunit)
    write(iunit,*) "nixpt"
    call WriteArrayReal(nixpt,size(nixpt),iunit)
    write(iunit,*) "niy0"
    call WriteArrayReal(niy0,size(niy0),iunit)
    write(iunit,*) "niy0s"
    call WriteArrayReal(niy0s,size(niy0s),iunit)
    write(iunit,*) "niy1"
    call WriteArrayReal(niy1,size(niy1),iunit)
    write(iunit,*) "niy1s"
    call WriteArrayReal(niy1s,size(niy1s),iunit)
    write(iunit,*) "nm"
    call WriteArrayReal(nm,size(nm),iunit)
    write(iunit,*) "nratio"
    call WriteArrayReal(nratio,size(nratio),iunit)
    write(iunit,*) "ntau"
    call WriteArrayReal(ntau,size(ntau),iunit)
    write(iunit,*) "nucx"
    call WriteArrayReal(nucx,size(nucx),iunit)
    write(iunit,*) "nucxi"
    call WriteArrayReal(nucxi,size(nucxi),iunit)
    write(iunit,*) "nuelg"
    call WriteArrayReal(nuelg,size(nuelg),iunit)
    write(iunit,*) "nueli"
    call WriteArrayReal(nueli,size(nueli),iunit)
    write(iunit,*) "nuii"
    call WriteArrayReal(nuii,size(nuii),iunit)
    write(iunit,*) "nuiistar"
    call WriteArrayReal(nuiistar,size(nuiistar),iunit)
    write(iunit,*) "nuix"
    call WriteArrayReal(nuix,size(nuix),iunit)
    write(iunit,*) "nuiz"
    call WriteArrayReal(nuiz,size(nuiz),iunit)
    write(iunit,*) "nurc"
    call WriteArrayReal(nurc,size(nurc),iunit)
    write(iunit,*) "nuvl"
    call WriteArrayReal(nuvl,size(nuvl),iunit)
    write(iunit,*) "nz2"
    call WriteArrayReal(nz2,size(nz2),iunit)
    write(iunit,*) "nzloc"
    call WriteArrayReal(nzloc,size(nzloc),iunit)
    write(iunit,*) "openbox"
    write(iunit,*) openbox
    write(iunit,*) "parvis"
    call WriteArrayReal(parvis,size(parvis),iunit)
    write(iunit,*) "pg"
    call WriteArrayReal(pg,size(pg),iunit)
    write(iunit,*) "pgy0"
    call WriteArrayReal(pgy0,size(pgy0),iunit)
    write(iunit,*) "pgy1"
    call WriteArrayReal(pgy1,size(pgy1),iunit)
    write(iunit,*) "phi"
    call WriteArrayReal(phi,size(phi),iunit)
    write(iunit,*) "phiv"
    call WriteArrayReal(phiv,size(phiv),iunit)
    write(iunit,*) "phiy0"
    call WriteArrayReal(phiy0,size(phiy0),iunit)
    write(iunit,*) "phiy0s"
    call WriteArrayReal(phiy0s,size(phiy0s),iunit)
    write(iunit,*) "phiy1"
    call WriteArrayReal(phiy1,size(phiy1),iunit)
    write(iunit,*) "phiy1s"
    call WriteArrayReal(phiy1s,size(phiy1s),iunit)
    write(iunit,*) "pr"
    call WriteArrayReal(pr,size(pr),iunit)
    write(iunit,*) "prad"
    call WriteArrayReal(prad,size(prad),iunit)
    write(iunit,*) "pradc"
    call WriteArrayReal(pradc,size(pradc),iunit)
    write(iunit,*) "pradcff"
    call WriteArrayReal(pradcff,size(pradcff),iunit)
    write(iunit,*) "pradz"
    call WriteArrayReal(pradz,size(pradz),iunit)
    write(iunit,*) "pradzc"
    call WriteArrayReal(pradzc,size(pradzc),iunit)
    write(iunit,*) "pre"
    call WriteArrayReal(pre,size(pre),iunit)
    write(iunit,*) "prev"
    call WriteArrayReal(prev,size(prev),iunit)
    write(iunit,*) "pri"
    call WriteArrayReal(pri,size(pri),iunit)
    write(iunit,*) "priv"
    call WriteArrayReal(priv,size(priv),iunit)
    write(iunit,*) "priy0"
    call WriteArrayReal(priy0,size(priy0),iunit)
    write(iunit,*) "priy1"
    call WriteArrayReal(priy1,size(priy1),iunit)
    write(iunit,*) "prtv"
    call WriteArrayReal(prtv,size(prtv),iunit)
    write(iunit,*) "psor"
    call WriteArrayReal(psor,size(psor),iunit)
    write(iunit,*) "psor_tmpov"
    call WriteArrayReal(psor_tmpov,size(psor_tmpov),iunit)
    write(iunit,*) "psorbgg"
    call WriteArrayReal(psorbgg,size(psorbgg),iunit)
    write(iunit,*) "psorbgz"
    call WriteArrayReal(psorbgz,size(psorbgz),iunit)
    write(iunit,*) "psorc"
    call WriteArrayReal(psorc,size(psorc),iunit)
    write(iunit,*) "psorcxg"
    call WriteArrayReal(psorcxg,size(psorcxg),iunit)
    write(iunit,*) "psorcxgc"
    call WriteArrayReal(psorcxgc,size(psorcxgc),iunit)
    write(iunit,*) "psordis"
    call WriteArrayReal(psordis,size(psordis),iunit)
    write(iunit,*) "psorg"
    call WriteArrayReal(psorg,size(psorg),iunit)
    write(iunit,*) "psorgc"
    call WriteArrayReal(psorgc,size(psorgc),iunit)
    write(iunit,*) "psori"
    call WriteArrayReal(psori,size(psori),iunit)
    write(iunit,*) "psorold"
    call WriteArrayReal(psorold,size(psorold),iunit)
    write(iunit,*) "psorrg"
    call WriteArrayReal(psorrg,size(psorrg),iunit)
    write(iunit,*) "psorrgc"
    call WriteArrayReal(psorrgc,size(psorrgc),iunit)
    write(iunit,*) "psorxr"
    call WriteArrayReal(psorxr,size(psorxr),iunit)
    write(iunit,*) "psorxrc"
    call WriteArrayReal(psorxrc,size(psorxrc),iunit)
    write(iunit,*) "psorxrold"
    call WriteArrayReal(psorxrold,size(psorxrold),iunit)
    write(iunit,*) "pwrebkg"
    call WriteArrayReal(pwrebkg,size(pwrebkg),iunit)
    write(iunit,*) "pwribkg"
    call WriteArrayReal(pwribkg,size(pwribkg),iunit)
    write(iunit,*) "pwrze"
    call WriteArrayReal(pwrze,size(pwrze),iunit)
    write(iunit,*) "pwrzec"
    call WriteArrayReal(pwrzec,size(pwrzec),iunit)
    write(iunit,*) "q2cd"
    call WriteArrayReal(q2cd,size(q2cd),iunit)
    write(iunit,*) "qfl"
    write(iunit,*) qfl
    write(iunit,*) "qipar"
    call WriteArrayReal(qipar,size(qipar),iunit)
    write(iunit,*) "qsh"
    write(iunit,*) qsh
    write(iunit,*) "resco"
    call WriteArrayReal(resco,size(resco),iunit)
    write(iunit,*) "resee"
    call WriteArrayReal(resee,size(resee),iunit)
    write(iunit,*) "resei"
    call WriteArrayReal(resei,size(resei),iunit)
    write(iunit,*) "resmo"
    call WriteArrayReal(resmo,size(resmo),iunit)
    write(iunit,*) "rtau"
    call WriteArrayReal(rtau,size(rtau),iunit)
    write(iunit,*) "rtaue"
    call WriteArrayReal(rtaue,size(rtaue),iunit)
    write(iunit,*) "rtaux"
    call WriteArrayReal(rtaux,size(rtaux),iunit)
    write(iunit,*) "rtauy"
    call WriteArrayReal(rtauy,size(rtauy),iunit)
    write(iunit,*) "seec"
    call WriteArrayReal(seec,size(seec),iunit)
    write(iunit,*) "seev"
    call WriteArrayReal(seev,size(seev),iunit)
    write(iunit,*) "seg_ue"
    call WriteArrayReal(seg_ue,size(seg_ue),iunit)
    write(iunit,*) "seic"
    call WriteArrayReal(seic,size(seic),iunit)
    write(iunit,*) "seiv"
    call WriteArrayReal(seiv,size(seiv),iunit)
    write(iunit,*) "smoc"
    call WriteArrayReal(smoc,size(smoc),iunit)
    write(iunit,*) "smov"
    call WriteArrayReal(smov,size(smov),iunit)
    write(iunit,*) "sng_ue"
    call WriteArrayReal(sng_ue,size(sng_ue),iunit)
    write(iunit,*) "snic"
    call WriteArrayReal(snic,size(snic),iunit)
    write(iunit,*) "sniv"
    call WriteArrayReal(sniv,size(sniv),iunit)
    write(iunit,*) "sxyxpt"
    write(iunit,*) sxyxpt
    write(iunit,*) "te"
    call WriteArrayReal(te,size(te),iunit)
    write(iunit,*) "tev"
    call WriteArrayReal(tev,size(tev),iunit)
    write(iunit,*) "tey0"
    call WriteArrayReal(tey0,size(tey0),iunit)
    write(iunit,*) "tey1"
    call WriteArrayReal(tey1,size(tey1),iunit)
    write(iunit,*) "tg"
    call WriteArrayReal(tg,size(tg),iunit)
    write(iunit,*) "tgy0"
    call WriteArrayReal(tgy0,size(tgy0),iunit)
    write(iunit,*) "tgy1"
    call WriteArrayReal(tgy1,size(tgy1),iunit)
    write(iunit,*) "ti"
    call WriteArrayReal(ti,size(ti),iunit)
    write(iunit,*) "tiv"
    call WriteArrayReal(tiv,size(tiv),iunit)
    write(iunit,*) "tiy0"
    call WriteArrayReal(tiy0,size(tiy0),iunit)
    write(iunit,*) "tiy0s"
    call WriteArrayReal(tiy0s,size(tiy0s),iunit)
    write(iunit,*) "tiy1"
    call WriteArrayReal(tiy1,size(tiy1),iunit)
    write(iunit,*) "tiy1s"
    call WriteArrayReal(tiy1s,size(tiy1s),iunit)
    write(iunit,*) "travis"
    call WriteArrayReal(travis,size(travis),iunit)
    write(iunit,*) "trax_use"
    call WriteArrayReal(trax_use,size(trax_use),iunit)
    write(iunit,*) "tray_use"
    call WriteArrayReal(tray_use,size(tray_use),iunit)
    write(iunit,*) "ttnpg"
    write(iunit,*) ttnpg
    write(iunit,*) "ttotfe"
    write(iunit,*) ttotfe
    write(iunit,*) "ttotjf"
    write(iunit,*) ttotjf
    write(iunit,*) "up"
    call WriteArrayReal(up,size(up),iunit)
    write(iunit,*) "upe"
    call WriteArrayReal(upe,size(upe),iunit)
    write(iunit,*) "upi"
    call WriteArrayReal(upi,size(upi),iunit)
    write(iunit,*) "upxpt"
    call WriteArrayReal(upxpt,size(upxpt),iunit)
    write(iunit,*) "uu"
    call WriteArrayReal(uu,size(uu),iunit)
    write(iunit,*) "uup"
    call WriteArrayReal(uup,size(uup),iunit)
    write(iunit,*) "uz"
    call WriteArrayReal(uz,size(uz),iunit)
    write(iunit,*) "v2"
    call WriteArrayReal(v2,size(v2),iunit)
    write(iunit,*) "v2cb"
    call WriteArrayReal(v2cb,size(v2cb),iunit)
    write(iunit,*) "v2cd"
    call WriteArrayReal(v2cd,size(v2cd),iunit)
    write(iunit,*) "v2ce"
    call WriteArrayReal(v2ce,size(v2ce),iunit)
    write(iunit,*) "v2dd"
    call WriteArrayReal(v2dd,size(v2dd),iunit)
    write(iunit,*) "v2rd"
    call WriteArrayReal(v2rd,size(v2rd),iunit)
    write(iunit,*) "v2xgp"
    call WriteArrayReal(v2xgp,size(v2xgp),iunit)
    write(iunit,*) "ve2cb"
    call WriteArrayReal(ve2cb,size(ve2cb),iunit)
    write(iunit,*) "ve2cd"
    call WriteArrayReal(ve2cd,size(ve2cd),iunit)
    write(iunit,*) "vex"
    call WriteArrayReal(vex,size(vex),iunit)
    write(iunit,*) "vey"
    call WriteArrayReal(vey,size(vey),iunit)
    write(iunit,*) "veycb"
    call WriteArrayReal(veycb,size(veycb),iunit)
    write(iunit,*) "veycp"
    call WriteArrayReal(veycp,size(veycp),iunit)
    write(iunit,*) "visx"
    call WriteArrayReal(visx,size(visx),iunit)
    write(iunit,*) "visxneo"
    call WriteArrayReal(visxneo,size(visxneo),iunit)
    write(iunit,*) "visy"
    call WriteArrayReal(visy,size(visy),iunit)
    write(iunit,*) "visyxpt"
    call WriteArrayReal(visyxpt,size(visyxpt),iunit)
    write(iunit,*) "vsoree"
    call WriteArrayReal(vsoree,size(vsoree),iunit)
    write(iunit,*) "vsoreec"
    call WriteArrayReal(vsoreec,size(vsoreec),iunit)
    write(iunit,*) "vy"
    call WriteArrayReal(vy,size(vy),iunit)
    write(iunit,*) "vy_cft"
    call WriteArrayReal(vy_cft,size(vy_cft),iunit)
    write(iunit,*) "vy_use"
    call WriteArrayReal(vy_use,size(vy_use),iunit)
    write(iunit,*) "vyavis"
    call WriteArrayReal(vyavis,size(vyavis),iunit)
    write(iunit,*) "vycb"
    call WriteArrayReal(vycb,size(vycb),iunit)
    write(iunit,*) "vyce"
    call WriteArrayReal(vyce,size(vyce),iunit)
    write(iunit,*) "vycf"
    call WriteArrayReal(vycf,size(vycf),iunit)
    write(iunit,*) "vycp"
    call WriteArrayReal(vycp,size(vycp),iunit)
    write(iunit,*) "vycr"
    call WriteArrayReal(vycr,size(vycr),iunit)
    write(iunit,*) "vydd"
    call WriteArrayReal(vydd,size(vydd),iunit)
    write(iunit,*) "vygp"
    call WriteArrayReal(vygp,size(vygp),iunit)
    write(iunit,*) "vyhxpt"
    call WriteArrayReal(vyhxpt,size(vyhxpt),iunit)
    write(iunit,*) "vyrd"
    call WriteArrayReal(vyrd,size(vyrd),iunit)
    write(iunit,*) "vytan"
    call WriteArrayReal(vytan,size(vytan),iunit)
    write(iunit,*) "vyte_cft"
    call WriteArrayReal(vyte_cft,size(vyte_cft),iunit)
    write(iunit,*) "vyti_cft"
    call WriteArrayReal(vyti_cft,size(vyti_cft),iunit)
    write(iunit,*) "vyvxpt"
    call WriteArrayReal(vyvxpt,size(vyvxpt),iunit)
    write(iunit,*) "w"
    call WriteArrayReal(w,size(w),iunit)
    write(iunit,*) "w0"
    call WriteArrayReal(w0,size(w0),iunit)
    write(iunit,*) "w1"
    call WriteArrayReal(w1,size(w1),iunit)
    write(iunit,*) "w2"
    call WriteArrayReal(w2,size(w2),iunit)
    write(iunit,*) "w3"
    call WriteArrayReal(w3,size(w3),iunit)
    write(iunit,*) "wjdote"
    call WriteArrayReal(wjdote,size(wjdote),iunit)
    write(iunit,*) "xcnearlb"
    write(iunit,*) xcnearlb
    write(iunit,*) "xcnearrb"
    write(iunit,*) xcnearrb
    write(iunit,*) "zcoef"
    write(iunit,*) zcoef
    write(iunit,*) "zeff"
    call WriteArrayReal(zeff,size(zeff),iunit)
    write(iunit,*) "znot"
    call WriteArrayReal(znot,size(znot),iunit)
    close(iunit)
end subroutine DebugHelper

