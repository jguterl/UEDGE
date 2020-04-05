#include "bbb.h"
#include "../com/com.h"
#include "../mppl.h"
#include "../sptodp.h"

subroutine InitOMP()
    Use Output
    Use OmpOptions,only:Nthreads,nnzmxperthread,omplenpfac,ompneq,OMPVerbose
    Use Jacobian,only:nnzmx
    Use Lsode, only:neq
    Use ParallelOptions,only:OMPParallelJac
    
    implicit none
    integer:: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
    ! do the allocation for omp variables
#ifdef _OPENMP
    ifparalleljac: if (OMPParallelJac.eq.1) then
        write(iout,'(a,i3)') '*OMP* OMP actived'
        !$omp parallel
        if (OMP_GET_THREAD_NUM().eq.0) then
            if (Nthreads.gt.OMP_GET_NUM_THREADS()) then
                write(iout,*) '*OMP* Warning: Number of threads requested is larger the maximum number of threads available.&
                Nthreadsmax:'                ,OMP_GET_NUM_THREADS()
                write(iout,*) '*OMP* Resetting Nthreads to Nthreadsmax'
                Nthreads=OMP_GET_NUM_THREADS()
            endif
            if (Nthreads.le.0) then
                call xerrab('Nthread must be >0')
            endif
            if (OMPVerbose.gt.0) write(iout,'(a,i3)') '**** Number of threads for omp calculations:',Nthreads
        endif
        !$omp END parallel
#else
        Nthreads=1
        write(iout,'(a,i3)') '*OMP* OMP not actived'
#endif
        if (Nthreads.gt.1) then
            nnzmxperthread=ceiling(real(nnzmx)/real(Nthreads-1))*omplenpfac
        else
            nnzmxperthread=nnzmx
        endif
        ompneq=neq
        call gchange('OmpJacobian',0)
    else ifparalleljac
        write(iout,'(a)') '*OMP* OMP not actived'
    endif ifparalleljac

end subroutine InitOMP

subroutine jac_calc_omp (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)

    !   Calculate Jacobian matrix (derivatives with respect to each
    !   dependent variable of the right-hand side of each rate equation).
    !   Lower and upper bandwidths are used to select for computation
    !   only those Jacobian elements that may be nonzero.
    !   Estimates of Jacobian elements are computed by finite differences.
    !   The Jacobian is stored in compressed sparse row format.
    Use Output
    use Timing,only:istimingon,ttjstor,ttotjf,ttimpjf,ShowTime
    use Grid,only:ngrid,ig,ijac,ijactot
    use Jacobian_csc,only:rcsc,jcsc,icsc,yldot_pert,yldot_unpt
    use OmpOptions,only:OMPDebug,Nthreads,OMPVerbose,OMPDebug,nnzmxperthread,OMPCheckNaN,WriteJacobian
    use OmpJacobian,only:iJacCol,rJacElem,iJacRow,ivmin,ivmax,nnz,nnzcum
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
    integer,allocatable :: iJacConstructor(:,:)
    real,allocatable:: rJacConstructor(:,:)

    ! ... Functions
    logical tstguardc
    real(kind=4) gettime

    ! ... Local variables:
    real(kind=4) sec4, tsjstor, tsimpjf, dtimpjf,time0,time1,TimeBuild,TimeCollect,TimeJacCalc
    integer:: i,thread,ith,iv,TID, OMP_GET_THREAD_NUM
    character(len = 80) ::  filename

    !   Get the range of the iv index for each thread
    call SplitIndex(neq,Nthreads,ivmin,ivmax)

    if (OMPDebug.gt.0) then
        write(iout,*)' *OMP* neq=',neq
        write(iout,*)' *OMP* Ivmin(ithread),Ivmax(ithread) ***'
        do ith=1,Nthreads
            write(iout,'(a,I3,a,I7,I7)') 'ithread ', ith,':',ivmin(ith),ivmax(ith)
        enddo
    endif

    nnz(1:Nthreads)=-1  
    nnzcum(1:Nthreads)=-1
    iJacCol(1:nnzmxperthread,1:Nthreads)=0
    rJacElem(1:nnzmxperthread,1:Nthreads)=0.0
    iJacRow(1:neq,1:Nthreads)=0
    if (OMPDebug.gt.0) then
        write(iout,*) '*OMP* Jacobian arrays set to zero'
    endif
    TimeJacCalc= gettime(sec4)

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
    TimeBuild=gettime(sec4)
    
    call OMPJacBuilder(neq, t, yl,yldot00, ml, mu,wk,iJacCol,rJacElem,iJacRow,nnz)

    TimeBuild=gettime(sec4)-TimeBuild
    if (OMPVerbose.gt.0) write(iout,*)'*OMP* Time to build jac:',TimeBuild
    !   end build jacobian ##############################################################    

    !   collect jacobian ##############################################################
    TimeCollect=gettime(sec4)
    nnzcum(1)=nnz(1)-1
    
    do ith=2,Nthreads
        nnzcum(ith)=nnzcum(ith-1)+nnz(ith)-1
    enddo
    if (OMPDebug.gt.0) then
    write(iout,*) '*OMP* nnz:',nnz(1:Nthreads)
    write(iout,*) '*OMP* nnzcum:',nnzcum(1:Nthreads)
    endif
    if (OMPVerbose.gt.0) write(iout,'(a,i9)') '**** Number of non-zero Jacobian elems:',nnzcum(Nthreads)
    
    if (nnzcum(Nthreads).gt.nnzmx) then
        call xerrab(' Problem: nnzcum > nnzmx...')
    endif


    jcsc(ivmin(1):ivmax(1))= iJacRow(ivmin(1):ivmax(1),1)
    do ith=2,Nthreads
        jcsc(ivmin(ith):ivmax(ith))= iJacRow(ivmin(ith):ivmax(ith),ith)+nnzcum(ith-1)
    enddo

    rcsc(1:nnz(1)-1)= rJacElem(1:nnz(1)-1,1)
    icsc(1:nnz(1)-1)= iJacCol(1:nnz(1)-1,1)
    do ith=2,Nthreads
        rcsc(nnzcum(ith-1)+1:nnzcum(ith))=rJacElem(1:nnz(ith)-1,ith)
        icsc(nnzcum(ith-1)+1:nnzcum(ith))=iJacCol(1:nnz(ith)-1,ith)
    enddo
    jcsc(neq+1) = nnzcum(Nthreads)+1

    TimeCollect=gettime(sec4)-TimeCollect
    if (OMPVerbose.gt.0) write(iout,*)'*OMP* Time to collect jac:',TimeCollect
    !   end collect jacobian ##############################################################

 
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
    if (ShowTime.gt.0) write(iout,*)'**** Time in jac_calc:',gettime(sec4)-TimeJacCalc
    
    return
end subroutine jac_calc_omp
!-----------------------------------------------------------------------
subroutine SplitIndex(neq,Nthreads,ivmin,ivmax)
    implicit none
    integer,intent(in) ::neq,Nthreads
    integer,intent(out)::ivmin(Nthreads),ivmax(Nthreads)
    integer:: Nsize,R,i

    Nsize=neq/Nthreads
    R=MOD(neq, Nthreads)
    
    if (Nthreads.gt.1) then
        if (R.eq.0) then
            do i=1,Nthreads
                ivmin(i)=1+Nsize*(i-1)
                ivmax(i)=Nsize*i
            enddo
        else
            do i=1,Nthreads-1
                ivmin(i)=1+Nsize*(i-1)
                ivmax(i)=Nsize*i
            enddo
            ivmin(Nthreads)=ivmax(Nthreads-1)+1
            ivmax(Nthreads)=neq
        endif
    else
        ivmin(Nthreads)=1
        ivmax(Nthreads)=neq
    endif
    
end subroutine SplitIndex
!-----------------------------------------------------------------------
subroutine OMPJacBuilder(neq, t, yl,yldot00, ml,mu,wk,iJacCol,rJacElem,iJacRow,nnz)
    use OmpOptions,only:OMPDebug,OMPCopyArray,OMPCopyScalar,nthreads,nnzmxperthread
    use OMPJacobian, only:ivmin,ivmax
    use OmpCopybbb
    use OmpCopycom
    use OmpCopyapi
    
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
    integer:: ith,omp_get_thread_num,tid,nnzlocal,ithcopy


    if (OMPDebug.gt.0)write(iout,*) '*OMP* Copying data....'
    
    if (OMPCopyArray.gt.0) then
        if (OMPDebug.gt.0)write(iout,*) '*OMP* Copying array....'
        call OmpCopyPointerbbb
        call OmpCopyPointercom
        call OmpCopyPointerapi
    endif
    if (OMPCopyScalar.gt.0) then
        if (OMPDebug.gt.0)write(iout,*) '*OMP* Copying scalar....'
        call OmpCopyScalarbbb
        call OmpCopyScalarcom
        call OmpCopyScalarapi
    endif
    !   We cannot use variables in the parallel construct declarations below when these variables are not in the scope of the subroutine
    Nthreadscopy=Nthreads
    nnzmxperthreadcopy=nnzmxperthread
    ivmincopy(1:Nthreads)=ivmin(1:Nthreads)
    ivmaxcopy(1:Nthreads)=ivmax(1:Nthreads)
    iJacColCopy(1:nnzmxperthread)=0
    rJacElemCopy(1:nnzmxperthread)=0.0
    iJacRowCopy(1:neq)=0
    ylcopy(1:neq+2)=yl(1:neq+2) ! a very barbarian use of yl(neq+1) is implemented as a switch in pandf... Error-prone!
    wkcopy(1:neq)=wk(1:neq) ! Could be set equal to zero as well. The worker wk is not an output...

    if (OMPDebug.gt.0) then
        write(iout,*) '*OMP* Starting parallel loop'
    endif
    tid=-1
    nnzlocal=-10000
    ! ivmincopy,ivmaxcopy,yldot00, neq an t  could be shared as well as well as   
    !$omp parallel do default(shared)&
    !$omp& firstprivate(ithcopy,ivmincopy,ivmaxcopy,tid,nnzlocal,ylcopy,wkcopy,ml,mu,yldot00,t,neq)&
    !$omp& firstprivate(nnzmxperthreadcopy,nthreadscopy,iJacRowCopy,iJacColCopy,rJacElemCopy)
    loopthread: do ith=1,Nthreads !ith from 1 to Nthread, tid from 0 to Nthread-1
        tid=omp_get_thread_num()
        ithcopy=ith
        if (OMPDebug.gt.0) write(iout,*) '*OMP* Thread id:',tid,' <-> ith:',ithcopy
        ! we keep all these parameters as it is easier to debug LocalJacBuilder and deal with private/shared attributes
        call LocalJacBuilder(ivmincopy(ithcopy),ivmaxcopy(ithcopy),neq, t, ylcopy,yldot00,ml,mu,wkcopy,&
        iJacColcopy,rJacElemcopy,iJacRowcopy,ithcopy,nnzlocal,nnzmxperthreadcopy,nthreadscopy)
        if (OMPDebug.gt.0) write(iout,*) '*OMP*,',tid,' nzlocal:',nnzlocal
        !$omp  critical
        iJacCol(1:nnzlocal,ithcopy)=iJacColCopy(1:nnzlocal)
        rJacElem(1:nnzlocal,ithcopy)=rJacElemCopy(1:nnzlocal)
        iJacRow(1:neq,ithcopy)=iJacRowCopy(1:neq)
        nnzcopy(ithcopy)=nnzlocal
        !$omp  end critical

    enddo loopthread
    !$omp  END PARALLEL DO
    
    nnz(1:Nthreads)=nnzcopy(1:Nthreads) !nnzcopy is not necssary as nnz would be shared anyway in the parallel construct
    
    if (OMPDebug.gt.0) then
        write(iout,*) '*OMP* End of parallel loop....'
    endif

end subroutine OMPJacBuilder
! ----------------------------------------------------------------------
subroutine LocalJacBuilder(ivmin,ivmax,neq, t, yl,yldot00, ml, mu, wk,iJacCol,rJacElem,iJacRow,ith,nnz,nnzmxperthread,nthreads)

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
