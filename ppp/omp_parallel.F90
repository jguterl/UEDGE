subroutine InitOMPJac()

    Use OMPJacSettings,only:omplenpfac,OMPJacVerbose,OMPJacNchunks,OMPJacStamp
    use OMPJac,only: nnzmxperchunk,NchunksJac
    Use ParallelSettings,only: OMPParallelPandf1,Nthreads
    !Use OMPPandf1Settings,only: OMPPandf1FlagVerbose,OMPPandf1FirstRun
    Use MPIJacSettings,only:nnzmxperproc,MPIrank
    Use HybridSettings,only: HybridOMPMPI,Hybridstamp
    Use Jacobian,only:nnzmx
    Use Lsode, only:neq

    implicit none
    integer:: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM,OMP_GET_MAX_THREADS
    character*8 :: MPIRankTag
#ifndef _OPENMP
     call xerrab("UEDGE was not compiled with OpenMP. Cannot use OMP parallel features.")
#else
    ! prepare MPI/Hybrid stamps for output
    if (HybridOMPMPI>0) then
        write(MPIRankTag,'(I4)') MPIrank
        write(OMPJacStamp,'(a,a,a)') '[',trim(adjustl(trim(MPIRankTag))),'] OMPJac* '
    endif
    if (OMPJacVerbose.gt.1) write(*,*) OMPJacStamp,' Max number of threads available:',OMP_GET_MAX_THREADS()
    call OMP_SET_NUM_THREADS(OMP_GET_MAX_THREADS())

    !$omp parallel
    if (OMP_GET_THREAD_NUM().eq.0) then
        if (Nthreads.gt.OMP_GET_MAX_THREADS()) then
            if (OMPJacVerbose.gt.1) write(*,*) OMPJacStamp,' Warning: # threads requested > # threads available'

            Nthreads=OMP_GET_NUM_THREADS()


        endif
        if (Nthreads.le.0) then
            call xerrab('Nthread must be >0')
        endif
        if (OMPJacVerbose.gt.0) write(*,'(a,a,i3)') OMPJacStamp,' Number of threads for omp calculations:',Nthreads
    endif
    !$omp END parallel
    if (OMPJacNchunks.eq.0) then
    NchunksJac=neq
    elseif (OMPJacNchunks.lt.0) then
       NchunksJac=Nthreads
    else
       NchunksJac=OMPJacNchunks
     endif


    if (Nthreads.gt.1) then
        if (HybridOMPMPI>0) then
            nnzmxperchunk=ceiling(real(nnzmxperproc)/real(NchunksJac-1)*omplenpfac)
        else
            nnzmxperchunk=ceiling(real(nnzmx)/real(NchunksJac))*omplenpfac !nnzmx=neq*lenfac
        endif
    else
        if (HybridOMPMPI>0) then
            nnzmxperchunk=nnzmxperproc
        else
            nnzmxperchunk=ceiling(real(nnzmx)/real(NchunksJac))*omplenpfac !nnzmx=neq*lenfac
        endif
    endif

    write(*,*) OMPJacStamp,' Nthreads:', Nthreads, 'NchunksJac:',NchunksJac, 'nnzmxperchunk:',nnzmxperchunk,'neq:',neq
    call gchange('OMPJac',0)
!    if (OMPParallelPandf1.gt.0) then
!      OMPPandf1FlagVerbose=1
!      OMPPandf1FirstRun=1
!    endif

#endif
end subroutine InitOMPJAC

subroutine InitOMPPandf1
    use OMPPandf1Settings,only: OMPPandf1Nxchunks,OMPPandf1Nychunks,OMPPandf1Stamp,OMPPandf1Verbose
    use OMPPandf1, only: NchunksPandf1,Nxchunks,Nychunks
    use Dim,only:ny



    if (OMPPandf1Nychunks.lt.0) then
        call xerrab('Nychunks<0. Nxchunks must be >=0')
    endif

    if (OMPPandf1Nychunks.eq.0) then
        Nychunks=ny
    else
        Nychunks=OMPPandf1Nychunks
    endif

    if (OMPPandf1Nxchunks.ne.1) then
        call xerrab('OMPPandf1Nxchunks!=1. Only Nxchunks=1 is implemented for the moment...')
    else
        Nxchunks=1
    endif

    ! this is a placeholder for further parallelization but we need to implement handling of x-points discon. in ix indexing
    NchunksPandf1=Nychunks

    call gchange('OMPPandf1',0)

    if (OMPPandf1Verbose.gt.0) then
        write(*,*) OMPPandf1Stamp, ' NchunksPandf1 = ',NchunksPandf1
    endif

end subroutine InitOMPPandf1
!-------------------------------------------------------------------------------------------------

subroutine InitZeroOMP
#ifdef _OPENMP
call OmpInitZerobbb
call OmpInitZeroapi
call OmpInitZerocom
#endif
end subroutine InitZeroOMP
#ifdef _OPENMP
subroutine OMPCollectJacobian(neq,nnzmx,rcsc,icsc,jcsc,nnzcumout)
    use OMPJac,only:iJacCol,rJacElem,iJacRow,OMPivmin,OMPivmax,nnz,nnzcum,OMPTimeJacRow,NchunksJac
    use OMPJacSettings,only:OMPJacVerbose,OMPJacStamp,OMPTimingJacRow
    use ParallelDebug,only:OMPJacDebug
    integer,intent(in):: neq
    integer,intent(in):: nnzmx          ! maximum number of nonzeros in Jacobian
    real,intent(out)   :: rcsc(nnzmx)     ! nonzero Jacobian elements
    integer,intent(out):: icsc(nnzmx)   ! col indices of nonzero Jacobian elements
    integer,intent(out):: jcsc(neq+1)   ! pointers to beginning of each row in jac,ja
    integer,intent(out):: nnzcumout
    integer ichunk
    integer:: iunit,iv
    nnzcum(1:NchunksJac)=-1
    nnzcum(1)=nnz(1)-1

    do ichunk=2,NchunksJac
        nnzcum(ichunk)=nnzcum(ichunk-1)+nnz(ichunk)-1
    enddo
    if (OMPJacDebug.gt.0) then
        write(*,*) OMPJacStamp,' nnz:',nnz(1:NchunksJac)
        write(*,*) OMPJacStamp,' nnzcum:',nnzcum(1:NchunksJac)
    endif
    if (OMPJacVerbose.gt.0) write(*,'(a,i9)') '**** Number of non-zero Jacobian elems:',nnzcum(NchunksJac)

    if (nnzcum(NchunksJac).gt.nnzmx) then
        write(*,*) 'nnzcum=',nnzcum
        write(*,*) nnzmx
        call xerrab(' Problem: nnzcum > nnzmx...')
    endif



    jcsc(OMPivmin(1):OMPivmax(1))= iJacRow(OMPivmin(1):OMPivmax(1))
    do ichunk=2,NchunksJac
        jcsc(OMPivmin(ichunk):OMPivmax(ichunk))= iJacRow(OMPivmin(ichunk):OMPivmax(ichunk))+nnzcum(ichunk-1)
    enddo

    rcsc(1:nnz(1)-1)= rJacElem(1:nnz(1)-1,1)
    icsc(1:nnz(1)-1)= iJacCol(1:nnz(1)-1,1)
    do ichunk=2,NchunksJac
        rcsc(nnzcum(ichunk-1)+1:nnzcum(ichunk))=rJacElem(1:nnz(ichunk)-1,ichunk)
        icsc(nnzcum(ichunk-1)+1:nnzcum(ichunk))=iJacCol(1:nnz(ichunk)-1,ichunk)
    enddo
    nnzcumout=nnzcum(NchunksJac)
    if (OMPTimingJacRow.gt.0) then
    open (newunit = iunit, file = 'omptiming.dat')
    do iv=1,neq
    write(iunit,*) iv,OMPTimeJacRow(iv)
    enddo
    close(iunit)
    endif
end subroutine OMPCOllectJacobian
!-------------------------------------------------------------------------------------------------
subroutine jac_calc_omp (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)

    !   Calculate Jacobian matrix (derivatives with respect to each
    !   dependent variable of the right-hand side of each rate equation).
    !   Lower and upper bandwidths are used to select for computation
    !   only those Jacobian elements that may be nonzero.
    !   Estimates of Jacobian elements are computed by finite differences.
    !   The Jacobian is stored in compressed sparse row format.

    use Timing,only:istimingon,ttjstor,ttotjf,ttimpjf
    use PandfTiming
    use Grid,only:ngrid,ig,ijac,ijactot
    use Jacobian_csc,only:rcsc,jcsc,icsc,yldot_pert,yldot_unpt
    use OMPJac,only:NchunksJac,nnzmxperchunk
    use ParallelSettings,only: Nthreads
    use OMPJacSettings,only:OMPJacVerbose,OMPCheckNaN,&
    OMPLoadBalance,OMPAutoBalance,OMPJacStamp,OMPBalanceStrength
    use ParallelDebug,only: WriteJacobian,OMPJacDebug

    use OMPJac,only:iJacCol,rJacElem,iJacRow,OMPivmin,OMPivmax,nnz,nnzcum,OMPLoadWeight,OMPTimeLocalJac
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
    real tick,tock
    external tick,tock

    ! ... Local variables:
    real(kind=4) sec4, tsjstor, tsimpjf, dtimpjf,time0,time1
    integer:: i,thread,ichunk,iv,TID, OMP_GET_THREAD_NUM
    character(len = 80) ::  filename
    ! Calculate load distribution for threads
    if (OMPLoadBalance.ne.1 .and. OMPAutoBalance.ne.1) then
        OMPLoadWeight(1:NchunksJac)=1.0
    endif
    if (OMPAutoBalance.eq.1) then
        !Check that time are not zero
        if (OMPBalanceStrength<=0) call xerrab('OMPBalanceStrength must be >0')
        if (minval(OMPTimeLocalJac).gt.0.0) then
            do i=1,NchunksJac
                OMPLoadWeight(i)=OMPLoadWeight(i)*1/(OMPTimeLocalJac(i)/sum(OMPTimeLocalJac)*real(NchunksJac))**OMPBalanceStrength
            enddo
        else
            OMPLoadWeight(1:NchunksJac)=1.0
        endif
    endif
    !   Get the range of the iv index for each thread
    call OMPSplitIndex(1,neq,NchunksJac,OMPivmin,OMPivmax,OMPLoadWeight)

    if (OMPJacVerbose.gt.0) then
        write(*,*)' *OMPJac* neq=',neq,neqmx
        write(*,*)' *OMPJac* Ivmin(ichunk),Ivmax(ichunk), OMPLoadWeight(ichunk) : OMPTimeLocalJac(ichunk) ***'
        do ichunk=1,NchunksJac
write(*,'(a,I3,a,I7,I7,f6.2,a,f6.2)') '  * ichunk ', ichunk,':',OMPivmin(ichunk),OMPivmax(ichunk),OMPLoadWeight(ichunk)&
            ,' : ',OMPTimeLocalJac(ichunk)
        enddo
    endif


    OMPTimeJacCalc= tick()

    !   Get initial value of system cpu timer.
    if (istimingon .eq. 1) tsjstor = gettime(sec4)

    !   Count Jacobian evaluations, both for total and for this case
    ijactot = ijactot + 1
    ijac(ig) = ijac(ig) + 1
    if (svrpkg.eq.'nksol') write(*,*) ' Updating Jacobian, npe =  ',ijac(ig)

    !   Set up diagnostic arrays for debugging
    do iv = 1, neq
        yldot_unpt(iv) = yldot00(iv)  ! for diagnostic only
        yldot_pert(iv) = 0.
    enddo

    !   build jacobian ##############################################################
    OMPTimeBuild=tick()
    nnz(1:NchunksJac)=-1
    call OMPJacBuilder(neq, t, yl,yldot00, ml, mu,wk,iJacCol,rJacElem,iJacRow,nnz)
    OMPTotTimebuild = OMPTotTimeBuild+tock(OMPTimeBuild)
    if (OMPJacVerbose.gt.0) write(*,*)OMPJacStamp,' Time to build jac:',OMPTimeBuild
    !   end build jacobian ##############################################################

    !   collect jacobian ##############################################################
    OMPTimeCollect=tick()
    call OMPCollectJacobian(neq,nnzmx,rcsc,icsc,jcsc,nnzcumout)
    OMPTotTimeCollect = OMPTotTimeCollect+tock(OMPTimeCollect)
    if (OMPJacVerbose.gt.0) write(*,*)OMPJacStamp,' Time to collect jac:',OMPTimeCollect
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
                write(*,*) 'rcsc is NaN at i=',i,rcsc(i)
            endif
            if (icsc(i).ne.icsc(i)) then
                write(*,*) 'icsc is NaN at i=',i,icsc(i)
            endif
        enddo
        do i=1,neq+1
            if (jcsc(i).ne.jcsc(i)) then
                write(*,*) 'jcsc is NaN at i=',i
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

    
    OMPTotTimeJacCalc = OMPTotTimeJacCalc+tock(OMPTimeJacCalc)
    write(*,'(a,f8.2,a,I3,aI3,a)') 'Time in omp jac_calc:',OMPTimeJacCalc,'[',Nthreads,'|',NchunksJac,']'
    return
end subroutine jac_calc_omp
!-------------------------------------------------------------------------------------------------
subroutine OMPJacBuilder(neq, t, yl,yldot00, ml,mu,wk,iJacCol,rJacElem,iJacRow,nnz)
    use OMPJacSettings,only:OMPJacStamp,OMPJacVerbose,OMPLoopJacNchunk
    use ParallelDebug,only: OMPJacDebug,OMPCopyArray,OMPCopyScalar
    use OMPJac, only:OMPivmin,OMPivmax,OMPTimeLocalJac,OMPTimeJacRow,nnzmxperchunk,NchunksJac
    use Selec, only: yinc,xlinc,xrinc ! these variables are threadprivate because modify in pandf1 parallel loop. We copy them in ech thread with the copyin clause
    use OmpCopybbb
    use OmpCopycom
    use OmpCopyapi
    use omp_lib
    use PandfTiming

    implicit none
    integer,intent(inout)::nnz(NchunksJac)
    integer,intent(in):: neq      ! total number of equations (all grid points)
    integer,intent(in):: ml, mu   ! lower and upper bandwidths
    real,intent(in):: t           ! physical time
    real,intent(in) ::yl(*)       ! dependent variables
    real,intent(in) ::yldot00(neq+2) ! right-hand sides evaluated at yl
    real,intent(inout) :: wk(neq)
    integer,intent(out)::iJacCol(nnzmxperchunk,NchunksJac)
    integer,intent(out):: iJacRow(neq)
    real,intent(out):: rJacElem(nnzmxperchunk,NchunksJac)
    real ::wkcopy(neq)
    real::ylcopy(neq+2)
    real::tick,tock
    external tick,tock
    integer ::iJacColCopy(nnzmxperchunk),iJacRowCopy(neq)
    integer ::ivmincopy(NchunksJac),ivmaxcopy(NchunksJac)
    integer ::NchunksJaccopy,nnzmxperchunkcopy
    real :: rJacElemCopy(nnzmxperchunk),TimeJacRowcopy(neq)
    integer:: ichunk,tid,nnzlocal
    DOUBLE PRECISION :: TimeThread
    OMPTimeCopy=tick() 
    if (OMPJacDebug.gt.0)write(*,*) OMPJacStamp,' Copying data....'
    call pandf1 (-1, -1, 0, neq, 0, yl, ylcopy)
    if (OMPCopyArray.gt.0) then
        if (OMPJacDebug.gt.0)write(*,*) OMPJacStamp,' Copying array....'
        call OmpCopyPointerbbb
        call OmpCopyPointercom
        call OmpCopyPointerapi
    endif

    if (OMPCopyScalar.gt.0) then
        if (OMPJacDebug.gt.0)write(*,*) OMPJacStamp,' Copying scalar....'
        call OmpCopyScalarbbb
        call OmpCopyScalarcom
        call OmpCopyScalarapi
    endif

    
    !   We cannot use variables in the parallel construct declarations below when these variables are not in the scope of the subroutine
    NchunksJaccopy=NchunksJac
    nnzmxperchunkcopy=nnzmxperchunk
    ivmincopy(1:NchunksJac)=OMPivmin(1:NchunksJac)
    ivmaxcopy(1:NchunksJac)=OMPivmax(1:NchunksJac)
    iJacColCopy(1:nnzmxperchunk)=0
    rJacElemCopy(1:nnzmxperchunk)=0.0
    TimeJacRowcopy(1:neq)=0
    iJacRowCopy(1:neq)=0
    ylcopy(1:neq+2)=yl(1:neq+2) ! a very barbarian use of yl(neq+1) is implemented as a switch in pandf... Error-prone!
    wkcopy(1:neq)=wk(1:neq) ! Could be set equal to zero as well. The worker wk is not an output...

    if (OMPJacDebug.gt.0) then
        write(*,*) OMPJacStamp,' Starting parallel loop'
    endif
    tid=-1
    nnzlocal=-10000
    OMPTotTimeCopy=OMPTotTimeCopy+tock(OMPTimeCopy)
    OMPTimeLocal=tick()
    ! ivmincopy,ivmaxcopy,yldot00, neq an t  could be shared as well as well as
    !$omp parallel do schedule(dynamic,OMPLoopJacNchunk) default(shared)&
    !$omp& firstprivate(ivmincopy,ivmaxcopy,tid,nnzlocal,ylcopy)&
    !$omp& firstprivate(NchunksJaccopy,iJacRowCopy,iJacColCopy,rJacElemCopy,TimeJacRowcopy)&
    !$omp& private(TimeThread)  copyin(yinc,xlinc,xrinc)

    loopthread: do ichunk=1,NchunksJac !ichunk from 1 to Nthread, tid from 0 to Nthread-1
        Timethread = omp_get_wtime()
        tid=omp_get_thread_num()
        if (OMPJacDebug.gt.0) write(*,*) OMPJacStamp,' Thread id:',tid,' <-> ichunk:',ichunk
        ! we keep all these parameters as it is easier to debug LocalJacBuilder and deal wichunk private/shared attributes

        call LocalJacBuilder(ivmincopy(ichunk),ivmaxcopy(ichunk),neq, t, ylcopy,yldot00,ml,mu,&
            iJacColcopy,rJacElemcopy,iJacRowcopy,ichunk,nnzlocal,nnzmxperchunk,TimeJacRowcopy)
        if (OMPJacDebug.gt.0) write(*,*) OMPJacStamp,',',tid,' nzlocal:',nnzlocal

        !!!!$omp  critical
        iJacCol(1:nnzlocal,ichunk)=iJacColCopy(1:nnzlocal)
        rJacElem(1:nnzlocal,ichunk)=rJacElemCopy(1:nnzlocal)
        iJacRow(OMPivmin(ichunk):OMPivmax(ichunk))=iJacRowCopy(OMPivmin(ichunk):OMPivmax(ichunk))
        OMPTimeJacRow(ivmincopy(ichunk):ivmaxcopy(ichunk))=TimeJacRowcopy(ivmincopy(ichunk):ivmaxcopy(ichunk))
        nnz(ichunk)=nnzlocal
        OMPTimeLocalJac(tid)=omp_get_wtime() - Timethread
        !!!$omp  end critical

        if (OMPJacVerbose.gt.1) write(*,*) OMPJacStamp,' Time in thread #', tid,':',OMPTimeLocalJac(tid)
        if (OMPJacVerbose.gt.1) write(*,'(a,I3,a)') 'OMP thread ',tid,' exiting...'
    enddo loopthread
    !$omp  END PARALLEL DO
    OMPTotTimeLocal=OMPTotTimeLocal+tock(OMPTimeLocal)
    
    !nnz(1:NchunksJac)=nnzcopy(1:NchunksJac) !nnzcopy is not necssary as nnz would be shared anyway in the parallel construct

    if (OMPJacDebug.gt.0) then
        write(*,*) OMPJacStamp,' End of parallel loop....'
    endif


end subroutine OMPJacBuilder
#endif

subroutine OMPSplitIndex(ieqmin,ieqmax,NchunksJac,ivmin,ivmax,weight)
    implicit none
    integer,intent(in) ::ieqmin,ieqmax,NchunksJac
    real::weight(NchunksJac)
    integer,intent(out)::ivmin(NchunksJac),ivmax(NchunksJac)
    integer:: Nsize(NchunksJac),Msize,R,i,imax
    if (ieqmax-ieqmin+1<2) call xerrab('Number of equations to solve <2')
    if (NchunksJac.eq.ieqmax-ieqmin+1) then
    do i=1,NchunksJac
    ivmin(i)=i
    ivmax(i)=i
    enddo
    return
    endif

    if (NchunksJac.gt.1) then
        !        if (OMPLoadWeight.eq.1) then

        do i=1,NchunksJac
            if (weight(i)<=0) call xerrab('OMPSplitIndex: weight <0')
            !write(*,*) weight(i)
        enddo

        ! Normalized weights
        weight(1:NchunksJac)=weight(1:NchunksJac)/sum(weight(1:NchunksJac))*real(NchunksJac)
        do i=1,NchunksJac
            Nsize(i)=int(real((ieqmax-ieqmin+1)/NchunksJac)*weight(i))
            !write(*,*) Nsize(i),weight(i)
        enddo

        do i=1,NchunksJac
            if (Nsize(i)<0) call xerrab('Nsize<0')
            if (Nsize(i)<2) Nsize(i)=Nsize(i)+1
        enddo
        if (ieqmax-ieqmin+1.ne.sum(Nsize)) then
            imax=1
            do i=2,NchunksJac
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
        do i=2,NchunksJac
            ivmin(i)=ivmax(i-1)+1
            ivmax(i)=ivmin(i)+Nsize(i)-1
        enddo
        if (ivmax(NchunksJac)-ivmin(1)+1.ne.(ieqmax-ieqmin+1)) call xerrab('ivmax(NchunksJac)!=neq')
    else
        ivmin(NchunksJac)=ieqmin
        ivmax(NchunksJac)=ieqmax
    endif

end subroutine OMPSplitIndex


!!!! >>>>>> Routine for pandf1 parallelization <<<<<<<<<<<<<



!subroutine OMPSplitIndexyPandf(ieqmin,ieqmax,Nchunks,ic,inc,ivthread,neq)
!    Use Indexes,only: igyl
!    use OMPPandf1Settings,only: OMPPandf1Padyinc
!    implicit none
!    integer,intent(in) ::ieqmin,ieqmax,Nchunks,neq
!    integer,intent(out)::ic(Nchunks),inc(Nchunks),ivthread(1:neq)
!    integer :: ihigh(Nchunks),ilow(Nchunks),ixthread(ieqmin:ieqmax)
!    integer:: Nsize(Nchunks),Msize,R,i,imax,iv,imin,iend,istart,ichunk
!    if (ieqmax-ieqmin+1<2) call xerrab('Number of equations to solve <2')
!    i=ieqmin
!    if (ieqmax-ieqmin+1.eq.Nchunks) then
!    do ichunk=1,Nchunks
!    ic(ichunk)=i
!    ilow(ichunk)=i
!    ihigh(ichunk)=i
!    inc(ichunk)=OMPPandf1Padyinc
!    i=i+1
!    enddo
!    inc(1)=inc(1)+1
!    do iv=1,neq
!            do ichunk=1,Nchunks
!                if (igyl(iv,2).le.ihigh(ichunk) .and. igyl(iv,2).ge.ilow(ichunk)) then
!                    ivthread(iv)=ichunk
!                endif
!            enddo
!    enddo
!
!
!    return
!    endif
!
!        !        if (OMPLoadWeight.eq.1) then
!
!        !do i=1,Nchunks
!         !   if (weight(i)<=0) call xerrab('OMPSplitIndex: weight <0')
!            !write(*,*) weight(i)
!        !enddo
!
!        ! Normalized weights
!        !weight(1:Nchunks)=weight(1:Nchunks)/sum(weight(1:Nchunks))*real(Nchunks)
!        do i=1,Nchunks
!            !Nsize(i)=int(real((ieqmax-ieqmin+1)/Nchunks)*weight(i))
!            Nsize(i)=int(real((ieqmax-ieqmin+1)/Nchunks))
!        enddo
!
!        do i=1,Nchunks
!            if (Nsize(i)<1) call xerrab('Nsize<1')
!        enddo
!        if (Nchunks.ge.4) then
!        istart=2
!        iend=Nchunks-1
!        else
!        istart=1
!        iend=Nchunks
!        endif
!imin=istart
!        do while (ieqmax-ieqmin+1.gt.sum(Nsize))
!                Nsize(imin)=Nsize(imin)+1
!                if (imin.ge.iend) then
!                    imin=istart
!                else
!                    imin=imin+1
!                endif
!        enddo
!
!        if (ieqmax-ieqmin+1.ne.sum(Nsize)) then
!           write(*,*) Nsize,sum(Nsize)
!           call xerrab('OMPPandf1:Nsize .ne. ieqmax-ieqmin+1')
!       endif
!        ilow(1)=ieqmin
!        ihigh(1)=ieqmin+Nsize(1)-1
!        do i=2,Nchunks
!            ilow(i)=ihigh(i-1)+1
!            ihigh(i)=ilow(i)+Nsize(i)-1
!
!        enddo
!        if (ihigh(Nchunks).ne.ieqmax) call xerrab('ihigh(Nchunks)!=ieqmax')
!
!        ic=ilow+int((ihigh-ilow)/2)
!        do i=1,Nchunks
!           !write(*,*) max(ihigh(i)-ic(i),ic(i)-ilow(i)),OMPPandf1Padyinc
!           inc(i)=max(ihigh(i)-ic(i),ic(i)-ilow(i))+OMPPandf1Padyinc
!        enddo
!        do iv=1,neq
!            do i=1,Nchunks
!                if (igyl(iv,2).le.ihigh(i) .and. igyl(iv,2).ge.ilow(i)) then
!                    ivthread(iv)=i
!                endif
!            enddo
!        enddo
!
!end subroutine OMPSplitIndexyPandf

#ifdef _OPENMP
subroutine OMPPandf1Rhs(neq,time,yl,yldot)

    use omp_lib
    use OmpCopybbb
    use ParallelSettings,only: Nthreads,CheckPandf1
    use OMPPandf1Settings,only: OMPTimeParallelPandf1,OMPTimeSerialPandf1,OMPPandf1Stamp,OMPPandf1Verbose,OMPPandf1Debug
    use OMPPandf1,only:Nivchunk,ivchunk,yincchunk,xincchunk,iychunk,ixchunk,NchunksPandf1
    use OMPPandf1Settings,only:OMPPandf1loopNchunk
    use Dim,only:ny
    use Selec, only:yinc,xrinc,xlinc
    Use Grid,only:ijactot

    integer yinc_bkp,xrinc_bkp,xlinc_bkp,iv,tid
    integer,intent(in)::neq
    real,intent(in)::yl(*)
    real,intent(out)::yldot(*)
    real,intent(in)::time
    real::yldotcopy(1:neq)
    real yldotsave(1:neq),ylcopy(1:neq+2)
    character*80 ::FileName
    real time1,time2
    integer::ichunk
    !write(*,*)' ======================================== '
    ylcopy(1:neq+1)=yl(1:neq+1)

    if (ijactot.gt.0) then


    Time1=omp_get_wtime()
    call MakeChunksPandf1
    call OmpCopyPointerup
    !$omp parallel do default(shared) schedule(dynamic,OMPPandf1LoopNchunk) &
    !$omp& private(iv,tid) firstprivate(ylcopy) private(yldotcopy) copyin(yinc,xlinc,xrinc)
    loopthread: do ichunk=1,NchunksPandf1 !ichunk from 1 to Nthread, tid from 0 to Nthread-1
        !tid=omp_get_thread_num()
        !if (OMPPandf1Debug.gt.0) write(*,*) OMPPandf1Stamp,'Thread id:',tid,' <-> ichunk:',ichunk
        ! we keep all these parameters as it is easier to debug LocalJacBuilder and deal wichunk private/shared attributes
        yinc_bkp=yinc
        xlinc_bkp=xlinc
        xrinc_bkp=xrinc

        if (iychunk(ichunk).ne.-1) then
            yinc=yincchunk(ichunk)
        endif
        if (ixchunk(ichunk).ne.-1) then
            xrinc=xincchunk(ichunk)
            xlinc=xincchunk(ichunk)
        endif
        !xinc=OMPxinc(ichunk)
        !OMPTimeLocalPandf1(ichunk)=omp_get_wtime()

        call pandf1 (ixchunk(ichunk),iychunk(ichunk), 0, neq, time, ylcopy, yldotcopy)

        !OMPTimeLocalPandf1(ichunk)=omp_get_wtime()-OMPTimeLocalPandf1(ichunk)
        !OMPTimeCollectPandf1(ichunk)=omp_get_wtime()
        do iv=1,Nivchunk(ichunk)
            yldot(ivchunk(ichunk,iv))=yldotcopy(ivchunk(ichunk,iv))
        enddo


        !OMPTimeCollectPandf1(ichunk)=omp_get_wtime()-OMPTimeCollectPandf1(ichunk)
      yinc=yinc_bkp
      xlinc=xlinc_bkp
      xrinc=xrinc_bkp
    !if (OMPPandf1Verbose.gt.1) then
            !write(*,*) OMPPandf1Stamp,' Time in thread #', tid,':',OMPTimeLocalPandf1(ichunk),OMPTimeCollectPandf1(ichunk)
    !endif
    !OMPTimeLocalPandf1(ichunk)=OMPTimeCollectPandf1(ichunk)+OMPTimeLocalPandf1(ichunk)
    enddo loopthread
!$omp  END PARALLEL DO
Time1=omp_get_wtime()-Time1

      OMPTimeParallelPandf1=Time1+OMPTimeParallelPandf1
        !do iv=1,neq
        !yldot(iv)=yldotcopy(iv,OMPivthread(iv))
        !enddo




      if (CheckPandf1.gt.0) then
          Time2=omp_get_wtime()
         call pandf1 (-1, -1, 0, neq, time, ylcopy, yldotsave)
         Time2=omp_get_wtime()-Time2
          OMPTimeSerialPandf1=Time2+OMPTimeSerialPandf1
          write(*,*) "Timing Pandf1 serial:",OMPTimeSerialPandf1,"(",Time2,")/parallel:",OMPTimeParallelPandf1,'(',Time1,')'
          write(*,*) 'ijactot ',ijactot
          call Compare(yldot,yldotsave,neq)
          write(*,*) "serial and parallel pandf are identical"
        endif
  else
       call pandf1 (-1,-1, 0, neq, time, yl, yldot)
endif
end subroutine OMPPandf1Rhs

subroutine CreateBin(ieqmin,ieqmax,ichunkmin,ichunkmax,Padding,iCenterBin,iLeftBin,iRightBin,inc)
    integer,intent(in):: ieqmin, ieqmax,Padding,ichunkmin,ichunkmax
    integer,intent(out):: iCenterbin(Nchunk),inc(Nchunk),iLeftBin(Nchunk),iRightBin(Nchunk)
    integer ::N,SizeBin,i,Nchunk

    N=ieqmax-ieqmin+1
    Nchunk=ichunkmax-ichunkmin+1
    if (N>Nchunk) then
        SizeBin=int((N/Nchunk))
    else
        SizeBin=1
    endif

    iLeftBin(ichunkmin)=ieqmin
    iRightBin(ichunkmin)=iLeftBin(ichunkmin)+SizeBin-1
    iCenterBin(ichunkmin)=int((iLeftBin(ichunkmin)+iRightBin(ichunkmin))/2)
    inc(ichunkmin)=max(iCenterBin(ichunkmin)-iLeftBin(ichunkmin),iRightBin(ichunkmin)-iCenterBin(ichunkmin))+padding
    if (ichunkmax.gt.ichunkmin) then
    do i=ichunkmin+1,ichunkmax-1
        iLeftBin(i)=iRightBin(i-1)+1
        iRightBin(i)=iLeftBin(i)+SizeBin-1
        iCenterBin(i)=int((iLeftBin(i)+iRightBin(i))/2)
        inc(i)=max(iCenterBin(i)-iLeftBin(i),iRightBin(i)-iCenterBin(i))+padding
    enddo
    iLeftBin(ichunkmax)=iRightBin(ichunkmax-1)+1
    iRightBin(ichunkmax)=ieqmax
    iCenterBin(ichunkmax)=int((iLeftBin(ichunkmax)+iRightBin(ichunkmax))/2)
    inc(ichunkmax)=max(iCenterBin(ichunkmax)-iLeftBin(ichunkmax),iRightBin(ichunkmax)-iCenterBin(ichunkmax))+padding
    endif
    return
end subroutine CreateBin

subroutine MakeChunksPandf1()
    Use Indexes,only: igyl
    use OMPPandf1Settings,only: xpadding,ypadding,OMPPandf1Verbose
    use OMPPandf1, only: NchunksPandf1,yincchunk,xincchunk,iychunk,ixchunk,Nychunks_old,Nxchunks_old,neq_old,ivchunk,&
    Nivchunk,Nxchunks,Nychunks,iymaxchunk,ixmaxchunk,iyminchunk,ixminchunk
    use Lsode, only: neq
    use Dim, only:nx,ny
    implicit none
    integer:: remakechunk,i,ii,ichunk,iv,ix,iy
    integer:: iyCenterBin(Nychunks),iyRightBin(Nychunks),iyLeftBin(Nychunks),incy(Nychunks)
    integer::ixCenterBin(Nxchunks),ixRightBin(Nxchunks),ixLeftBin(Nxchunks),incx(Nxchunks)

       remakechunk=0
       if ((Nxchunks.ne.Nxchunks_old).or.(Nychunks.ne.Nychunks_old)) then


       if (Nychunks.gt.1) then
           if (Nychunks.eq.ny) then
              iyLeftBin(1)=0
              iyRightBin(1)=1
              iyCenterBin(1)=1
              incy(1)=ypadding
              call CreateBin(2,ny-1,2,Nychunks-1,ypadding,iyCenterBin,iyLeftBin,iyRightBin,incy)
              iyLeftBin(Nychunks)=ny
              iyRightBin(Nychunks)=ny+1
              iyCenterBin(Nychunks)=ny
              incy(Nychunks)=ypadding
              write(*,*) '----- Bins in y direction: ', Nychunks, ny+2
              do iy=1,Nychunks
                write(*,*) iyCenterBin(iy),iyLeftBin(iy),iyRightBin(iy),incy(iy)
              enddo
           else

            call CreateBin(0,ny+1,1,Nychunks,ypadding,iyCenterBin,iyLeftBin,iyRightBin,incy)

            write(*,*) '----- Bins in y direction: ', Nychunks, ny+2
            do iy=1,Nychunks
                    write(*,*) iyCenterBin(iy),iyLeftBin(iy),iyRightBin(iy),incy(iy)
            enddo

           ! now we check the first and last bins to check that ypadding is 3 if iyCenterBin=2
           if (iyCenterBin(1)==0) then
                incy(1)=incy(1)+1
           endif

           if (iyCenterBin(Nychunks)==ny+1) then
                incy(Nychunks)=incy(Nychunks)+1
           endif
           endif
       else
           iyCenterBin(1)=-1
           iyLeftBin(1)=0
           iyLeftBin(1)=ny+1
           incy(1)=0 !not used
       endif


       if (Nxchunks.gt.1) then
           call CreateBin(0,nx+1,1,Nxchunks,xpadding,ixCenterBin,ixLeftBin,ixRightBin,incx)
       else
           ixCenterBin(1)=-1
           ixLeftBin(1)=0
           ixRightBin(1)=nx+1
           incx(1)=0 !not used
       endif

     ichunk=1
    ! Build NchunksPandf1
    do iy=1,Nychunks
    if (iy==1) then
        iychunk(ichunk)=iyCenterBin(iy)
        iyminchunk(ichunk)=iyLeftBin(iy)
        iymaxchunk(ichunk)=iyRightBin(iy)
        yincchunk(ichunk)=incy(iy)
        ixchunk(ichunk)=-1
        ixminchunk(ichunk)=0
        ixmaxchunk(ichunk)=nx+1
        xincchunk(ichunk)=0
        ichunk=ichunk+1
        cycle
    endif
    if (iy==ny+1) then
        iychunk(ichunk)=iyCenterBin(iy)
        iyminchunk(ichunk)=iyLeftBin(iy)
        iymaxchunk(ichunk)=iyRightBin(iy)
        yincchunk(ichunk)=incy(iy)
        ixchunk(ichunk)=-1
        ixminchunk(ichunk)=0
        ixmaxchunk(ichunk)=nx+1
        xincchunk(ichunk)=0
        ichunk=ichunk+1
        cycle
    endif
        do ix=1,Nxchunks
            iychunk(ichunk)=iyCenterBin(iy)
            iyminchunk(ichunk)=iyLeftBin(iy)
            iymaxchunk(ichunk)=iyRightBin(iy)
            yincchunk(ichunk)=incy(iy)
            ixchunk(ichunk)=ixCenterBin(ix)
            ixminchunk(ichunk)=ixLeftBin(ix)
            ixmaxchunk(ichunk)=ixRightBin(ix)
            xincchunk(ichunk)=incx(ix)
            ichunk=ichunk+1
        enddo
    enddo

    Nychunks_old=Nychunks
    Nxchunks_old=Nxchunks
    remakechunk=1

    endif

    if ((neq.ne.neq_old) .or. remakechunk.gt.0) then

    do i=1,NchunksPandf1
            ii=1
            do iv=1,neq
                if ((igyl(iv,2).le.iymaxchunk(i) .and. igyl(iv,2).ge.iyminchunk(i)) &
                .and. (igyl(iv,1).le.ixmaxchunk(i) .and. igyl(iv,1).ge.ixminchunk(i))) then
                    ivchunk(i,ii)=iv
                    Nivchunk(i)=ii
                    ii=ii+1
                endif
            enddo
    enddo
    neq_old=neq
    if (OMPPandf1Verbose.gt.1) then
    write(*,"('Nychunks = ',I3,'; Nxchunks = ',I3)") Nychunks,Nxchunks
    write(*,"('Nchunks:',I3)") NchunksPandf1
    do i=1,NchunksPandf1
       write(*,'("ichunk: ",I3," | iyc = [",I3,";",I3,";",I3,"] ; ixc = ",I3, &
       " || xinc = ", I3,"; yinc = ",I3," | iv = [",I5,";",I5,"]")') &
       i,iyminchunk(i),iychunk(i),iymaxchunk(i), ixchunk(i),xincchunk(i),yincchunk(i), ivchunk(i,1),ivchunk(i,Nivchunk(i))
    enddo
    endif
    endif


end subroutine MakeChunksPandf1

!subroutine PrintChunks
!use OMPPandf1Settings,only: OMPPandf1Stamp
!use OMPPandf1, only: NchunksPandf1,yincchunk,xincchunk,iychunk,ixchunk,Nxchunks,Nychunks,iyminchunk,iymaxchunk
!integer ::i
!
!
!end subroutine PrintChunks

#endif
