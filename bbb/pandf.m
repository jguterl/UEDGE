c!include "bbb.h"
c!include "../com/com.h"
c!include "../mppl.h"
c!include "../sptodp.h"

c-----------------------------------------------------------------------
      subroutine pandf (xc, yc, neq, time, yl, yldot)

c... Calculates matrix A and the right-hand side depending on the values
c... of xc, yc.
c  Definitions for argument list
c
c  Input variables:
c    xc is poloidal index of perturbed variablefor Jacobian calc,
c       or =-1 for full RHS evaluation
c    yc is radial index for perturbed variable for Jacobian calc,
c       or =-1 for full RHS evaluation
c    neq is the total number of variables
c    time is the present physical time; useable by VODPK but not NKSOL
c    yl is the vector of unknowns
c  Output variables:
c    yldot is the RHS of ODE solver or RHS=0 for Newton solver (NKSOL)

      implicit none

*  -- input arguments
      integer xc, yc, neq
      real time, yl(*),yldot(*)

*  -- set local array dimension
      integer nigmx
      parameter (nigmx=100)

*  -- local variables
      integer ifld, jfld, zn, k, k1, k2, jx, ixt, ixt1, ixr, ixr1, iixt,
     .        ixt0
      real fxet, fxit, qr, vt0, vt1, vtn, vtn2, pradold, eeliold,
     .     erlizold, erlrcold, psorrgold(nigmx), psorcxgold(nigmx),
     .     nuizold(nigmx), nucxold(nigmx), nurcold(nigmx), nuixold(nigmx),
     .     psorgold(nigmx), tsfe, tsjf, niavex, niavey, teave, tiave, tgavex,
     .     zeffave, noavex, noavey, tiavey, tgavey, psordisold,
     .     nucxiold(nigmx), nueliold(nigmx), nuelgold(nigmx), rrfac, visxtmp,
     .     vttn, vttp, neavex, pwrebkgold, pwribkgold, feexflr, feixflr,
     .     naavex,naavey,nuelmolx,nuelmoly
      real fqpo, fqpom, friceo, friceom, upeo, upeom, fricio(100),
     .     friciom(100), upio(100), upiom(100), uupo(100), uupom(100)
      real nevol, ngvol, ngvolcap, kionz, krecz, kcxrz, kionm, krecm, kcxrm, nzbg,
     .     niz_floor, hflux, zflux, psorv, kionz0, pscx0, pxri, kcxrzig,
     .     nizm_floor, argx, massfac, ae, geyym, geyy0, geyyp, dgeyy0,
     .     dgeyy1, te_diss, wallfac, z1fac, bpolmin, rt2nus, epstmp, tv2
      real awoll,awll
      integer izch, ihyd, iimp, jg, jz, nsm1, ifld_fcs, ifld_lcs
      real uuv, ne_sgvi, nbarx, argth, fac_rad, ffyi, ffyo
      real grdnv, qflx, qfly, cshx, cshy, qshx, qshy, lxtec, lxtic
      real lmfpn, lmfppar, lmfpperp
      real temp1, temp2, temp3, temp4, cutlo3, lambd_ci, lambd_ce
      real upxavep1,upxave0,upxavem1,upf0,upfm1
      real teev
      logical xccuts, xcturb
      integer iy1, ixmp2, iyp1, iyp2, iym1, ixs, ixf, iys, iyf,
     .        methnx, methny, iy2, i2pwr, i5pwr, j2pwr, j5pwr,
     .        iysepu, ixpt1u, ixpt2u
      integer iy0, jy, jylo, jyhi, iy_min, iy_max
      real diffustotal
      real difnimix, kyemix, kyimix, v2dia, bscalf, bpfac
      real rdum,rdumaray(1),afqp,ltmax,lmfpe,lmfpi,flxlimf
      real rdumx,rdumy,dr1,dr2,qion
      real b_ctr,dbds_m,dbds_p,eta_h0,eta_hm,eta_hp,drag_1,drag_2,
     .     drag_3,mf_path,nu_ii,frac_col,fniy_recy
      real thetacc,dupdx,dupdy
      real dndym1,dndy0,dndyp1,d2ndy20,d2ndy2p1,d3ndy3
      real dtdym1,dtdy0,dtdyp1,d2tdy20,d2tdy2p1,d3tdy3,nhi_nha
      integer idum, idumaray(1)
      real(Size4) sec4, gettime, tsimpfe, tsimp, tsnpg, ueb
      integer impflag
      # former Aux module variables
      integer ix,iy,igsp,iv,iv1,iv2,iv3,ix1,ix2,ix3,ix4,ix5,ix6
      real tv,t0,t1,t2,a
      real length
cnxg      data igs/1/

      Use(Dim)      # nx,ny,nhsp,nusp,nzspt,nzsp,nisp,ngsp,nxpt
      Use(Xpoint_indices)      # ixlb,ixpt1,ixpt2,ixrb,iysptrx1,iysptrx2
                               # iysptrx
      Use(Math_problem_size)   # neqmx,numvar
      Use(Timing)   # istimingon,ttotfe,ttotjf,ttimpfe,ttimpjf,ttnpg
      Use(Share)    # nxpt,nxomit,geometry,nxc,cutlo,islimon,ix_lim,iy_lims
                    # istabon,isudsym
      Use(Phyvar)   # pi,me,mp,ev,qe,rt8opi
      Use(UEpar)    # methe,methu,methn,methi,methg,concap,lnlam,
                    # convis,icnuiz,icnucx,cnuiz,cnucx,isrecmon,
                    # ngbackg,ingb,eion,ediss,afix,coef,ce,ci,
                    # dp1,qfl,csh,qsh,mfl,msh,cs,fxe,ctaue,fxi,ctaui,
                    # zcoef,coef1,nurlxu,isphiofft,isintlog,isgxvon
                    # isgpye,frnnuiz,nzbackg,inzb,nlimix,nlimiy,
                    # isofric,isteaven
                    # isnionxy,isuponxy,isteonxy,istionxy,isngonxy,isphionxy,
                    # isupon,isteon,istion,isphion
      Use(Aux),only: ixmp
      Use(Coefeq)
      Use(Bcond)    # albedoo,albedoi,isfixlb,isfixrb
                    # xcnearlb,xcnearrb,openbox
      Use(Parallv)  # nxg,nyg
      Use(Rccoef)   # recylb,recyrb,recycw,recycz,sputtr
      Use(Fixsrc)
      Use(Selec)    # i1,i2,i3,i4,i5,i6,i7,i8,xlinc,xrinc,ixs1,
                    # j1,j1p,j2,j2p,j3,j4,j5,j6,j5p,j6p,j7,j8,j5m
      Use(Comgeo)   # isxptx,isxpty
      Use(Noggeo)   # fym,fy0,fyp,fymx,fypx,angfx,fxm,fx0,fxp,fxmy,fxpy
      Use(Compla)   # znucl,rtaux,rtauy,rtau,rt_scal
      Use(Comflo)
      Use(Conduc)   # lmfplim
      Use(Rhsides)
      Use(Save_terms)   # psorold,psorxrold
      Use(Indexes)
      Use(Ynorm)    # isflxvar,nnorm,ennorm,fnorm,n0,n0g
      Use(Poten)    # bcee,bcei,cthe,cthi
      Use(Comtra)   # parvis,travis,difni,difnit,difpr,difni2,
                    # difpr2,vcony,flalfe,flalfi,flgam,kxe,kxecore,
                    # kye,kyet,kxi,kxicore,kxn,kyi,kyit,kyn,feqp,
                    # flalfgx,flalfgy,flalfv,flalfvgx,flalfvgy,
                    # flalftgx,flalftgy,alfeqp,facbni,facbni2,facbee,
                    # facbei,isbohmcalc,diffusivity,diffusivwrk,
                    # diffusivloc,isdifxg_aug,isdifyg_aug,sigvi_floor,
                    # facbup
      Use(Bfield)   # btot,rbfbt,rbfbt2
      Use(Locflux)
      Use(Wkspace)
      Use(Gradients)
      Use(Imprad)   # isimpon,nzloc,impradloc,prad,pradz,na,ntau,nratio,
                    # afrac,atau,ismctab,rcxighg,pwrze

      Use(Volsrc)   # pwrsore,pwrsori,volpsor
      Use(Model_choice)   # iondenseqn
      Use(Time_dep_nwt)   # ylodt
      Use(Cfric)          # frice,frici
      Use(Turbulence)     # isturbnloc,isturbcons,diffuslimit,diffuswgts
      Use(Turbulence_diagnostics)   # chinorml,chinormh
      Use(MCN_dim)
      Use(MCN_sources)	  # uesor_ni,uesor_up,uesor_ti,uesor_te
      Use(Ext_neutrals)          # isextneuton, extneutopt
      Use(PNC_params)            # dtneut, dtold
      Use(PNC_data)              # ni_pnc, etc.
      Use(Reduced_ion_interface) # misotope,natomic
      Use(Indices_domain_dcl)    # ixmxbcl
      Use(Indices_domain_dcg)    # ndomain
      Use(Npes_mpi)              # mype
      Use(RZ_grid_info)  		 # bpol
      Use(Interp)				 # ngs, tgs
      Use(CapFloor)
      use ParallelEval,only: ParallelJac,ParallelPandf1
      use PandfTiming

*  -- procedures for atomic and molecular rates --
      integer zmax,znuc
      real dene,denz(0:1),radz(0:1)
      real rsa, rra, rqa, rcx, emissbs, erl1, erl2, radneq, radimpmc
      real radmc, svdiss, vyiy0, vyiym1, v2ix0, v2ixm1
      external rsa, rra, rqa, rcx, emissbs, erl1, erl2, radneq, radimpmc
      external radmc, svdiss
      real tick,tock
      external tick,tock

ccc      save

*  -- procedures --
      real ave, etaper, interp_yf, interp_xf
      ave(t0,t1) = 2*t0*t1 / (cutlo+t0+t1)
      etaper(ix,iy) = 3.234e-9*loglambda(ix,iy)/(max(te(ix,iy),temin*ev)
     .                                          /(1000.*ev))**(1.5)
      interp_yf(ix,iy,t0,t1) = (t0*gy(ix,iy) + t1*gy(ix,iy+1)) /
     .                                       (gy(ix,iy)+gy(ix,iy+1))
      interp_xf(ix,iy,t0,t1) =(t0*gx(ix,iy) + t1*gx(ixp1(ix,iy),iy)) /
     .                                (gx(ix,iy)+gx(ixp1(ix,iy),iy))
      cutlo3 = cutlo**0.3

c  Check array sizes
      if (ngsp > nigmx .or. nisp > nigmx) then
         call xerrab("***PANDF in oderhs.m: increase nigmx, recompile")
      endif
c... Timing of pandf components (added by J. Guterl)
        if (TimingPandf.gt.0
     . .and. yl(neq+1) .lt. 0 .and. ParallelPandf1.eq.0) then
        TimingPandfOn=1
        else
        TimingPandfOn=0
        endif
        if (TimingPandfOn.gt.0) TimePandf=tick()
c... Roadblockers for  call to pandf through openmp structures (added by J.Guterl)
      if ((isimpon.gt.0 .and. isimpon.ne.6) .and. (ParallelJac.gt.0 .or. ParallelPandf1.gt.0)) then
      call xerrab('Only isimpon=0 or 6 is validated with openmp.
     .Contact the UEDGE team to use other  options with openmp.')
      endif

      if ((ismcnon.gt.0) .and. (ParallelJac.gt.0 .or. ParallelPandf1.gt.0)) then
      call xerrab('Only ismcnon=0 is validated with openmp.
     .Contact the UEDGE team to use other options with openmp.')
      endif


************************************************************************
*  -- initialization --
************************************************************************
c     Set ix index for outer midplane turbulence
      if (isudsym==1) then
         ixmp2 = nxc + 2
      elseif (geometry=='dnull' .or. geometry(1:9)=="snowflake" .or.
     .        geometry=="dnXtarget" .or. geometry=="isoleg") then
         ixmp2 = ixmdp(2) + 1
      else
         ixmp2 = ixpt1(1) + nxomit + 3*(ixpt2(1)-ixpt1(1))/4
      endif

c     Set switches for neutrals-related source terms in plasma equations
c     (MER 1996/10/28)
c     (IJ  2015/04/06) add ismcnon>=3 for external call to run_neutrals
      if (ismcnon .eq. 1) then        # use MC sources only:
         if (svrpkg .eq. "cvode") then
           call xerrab('*** ismcnon=1 not allowed for cvode ***')
         endif
         cfneut=0.
         if (isupgon(1) .eq. 1) then
            cfvgpx(iigsp)=0.
            cfvgpy(iigsp)=0.
            cfvcsx(iigsp)=0.
            cfvcsy(iigsp)=0.
         endif
         cmneut=1.
      else if (ismcnon .eq. 2) then    # switch between two models:
         if (yl(neq+1) .gt. 0) then   # use fluid model for Jacobian
            cfneut=1.
            if (isupgon(1) .eq. 1) then
               cfvgpx(iigsp)=1.
               cfvgpy(iigsp)=1.
               cfvcsx(iigsp)=1.
               cfvcsy(iigsp)=1.
            endif
            cmneut=0.
         elseif (yl(neq+1) .lt. 0) then     # use MC model for evaluation
            cfneut=0.
            if (isupgon(1) .eq. 1) then
               cfvgpx(iigsp)=0.
               cfvgpy(iigsp)=0.
               cfvcsx(iigsp)=0.
               cfvcsy(iigsp)=0.
            endif
            cmneut=1.
         else
            call xerrab('*** PANDF: ismcnon=2 & yl(neq+1)=0 ???')
         endif
      else if (ismcnon .eq. 3) then         # switch between two neutral models internally
         if (yl(neq+1) .gt. 0) then         # use fluid model for preconditioner
            if (extneutmeth .eq. 1) then				#fluid source & implicit MC flux
               cfneut=1.     #turn on  fluid sources
               cfneutdiv=0.  #turn off fluid div fluxes
               cmneut=0.     #turn off MC sources
               cmneutdiv=1.  #turn on  MC div fluxes
               if (isupgon(1) .eq. 1) then
                  cfvgpx(iigsp)=1.
                  cfvgpy(iigsp)=1.
                  cfvcsx(iigsp)=1.
                  cfvcsy(iigsp)=1.
               endif
            else         								#fluid source & fluid flux
               cfneut=1.     #turn on fluid sources
               cfneutdiv=1.  #turn on fluid div fluxes
               cmneut=0.     #turn off MC sources
               cmneutdiv=0.  #turn off MC div fluxes
               if (isupgon(1) .eq. 1) then
                  cfvgpx(iigsp)=1.
                  cfvgpy(iigsp)=1.
                  cfvcsx(iigsp)=1.
                  cfvcsy(iigsp)=1.
               endif
            endif
         elseif (yl(neq+1) .lt. 0) then     # use MC model for RHS evaluation
            if (extneutmeth .eq. 1) then				#fluid source & implicit MC flux
               cfneut=1.     #turn on  fluid sources
               cfneutdiv=0.  #turn off fluid div fluxes
               cmneut=0.     #turn off MC sources
               cmneutdiv=1.  #turn on  MC div fluxes
               if (isupgon(1) .eq. 1) then
                  cfvgpx(iigsp)=1.
                  cfvgpy(iigsp)=1.
                  cfvcsx(iigsp)=1.
                  cfvcsy(iigsp)=1.
               endif
            else         								#MC source & fluid flux
               cfneut=0.     #turn off fluid sources
               cfneutdiv=1.  #turn on  fluid div fluxes
               cmneut=1.     #turn on  MC sources
               cmneutdiv=0.  #turn off MC div fluxes
               if (isupgon(1) .eq. 1) then
                  cfvgpx(iigsp)=0.
                  cfvgpy(iigsp)=0.
                  cfvcsx(iigsp)=0.
                  cfvcsy(iigsp)=0.
               endif
            endif

         else
            call xerrab('*** PANDF: ismcnon=3 & yl(neq+1)=0 ???')
         endif
         if (yl(neq+1) .lt. 0) then       # RHS eval (not precon eval)
			if (isextneuton .ne. 0) then  # implicit use of external neutrals inside exmain
                #Neutral step
                dtold=dtreal
                dtreal=dtneut
                call store_neutrals
                call run_neutrals		  # extneutopt sets choice of model
                call update_neutrals
                dtreal=dtold
            endif
        endif
      else if (ismcnon .eq. 4) then # test a different fluid model in the preconditioner
         if (yl(neq+1) .gt. 0) then   # Precon eval
            parvis=parvis*pnc_cfparvis
            travis=travis*pnc_cftravis
            do ifld=1,nisp
              ni(:,:,ifld)=ni(:,:,ifld)*pnc_cfni(ifld)
              up(:,:,ifld)=up(:,:,ifld)*pnc_cfup(ifld)
            enddo
c            write(*,*) 'ismcnon=4'
c            write(*,*) parvis
         endif
      end if #ismcnon

************************************************************************
*   This section is to use in the calculation of the jacobian locally.
************************************************************************

c ... Get initial value of system cpu timer.
      if(xc .lt. 0) then
         tsfe = gettime(sec4)
      else
         tsjf = gettime(sec4)
      endif

      if ( (xc .lt. 0) .or.
     .     ((0<=yc).and.(yc-yinc<=0)).and.isjaccorall==1 ) then
                                              # use full ix range near yc=0
                                              # with integrated core flux BC
         i1 = 0  # 1-ixmnbcl
         i2 = 1
         i2p = 1
         i3 = 0  # 1-ixmnbcl
         i4 = 0  # 1-ixmnbcl
         i5 = nx
         i5m = nx-1
         i6 = nx+1  # nx+ixmxbcl
         i7 = nx+1  # nx+ixmxbcl
         i8 = nx+1  # nx+ixmxbcl
      else
         i1 = max(0, xc-xlinc-1)
         i2 = max(1, xc-xlinc)
         i2p = max(1, xc-xrinc-1)
         i3 = xc-xlinc     # not used as loop indice (can be < 0)
         i4 = max(0, xc-xlinc)
         i5 = min(nx, xc+xrinc)
         i5m = min(nx-1, xc+xrinc)
         i6 = min(nx+1, xc+xrinc+1)
         i7 = xc+xrinc     # not used as loop indice (can be > nx)
         i8 = min(nx+1, xc+xrinc)
      endif
      if (yc .lt. 0) then
         j1 = 0
         j1p = 0
         j2 = 1
         j2p = 1
         j3 = 0
         j4 = 0
         j5 = ny
         j5m = ny-1
         j6 = ny+1
         j5p = ny
         j6p = ny+1
         j7 = ny+1
         j8 = ny+1
      else
         j1 = max(0, yc-yinc-1)
         j2 = max(1, yc-yinc)
         j1p = max(0, yc-yinc-2)
         j2p = max(1, yc-yinc-1)
         j3 = yc-yinc
         j4 = max(0, yc-yinc)
         j5 = min(ny, yc+yinc)
         j5m = min(ny-1, yc+yinc)
         j6 = min(ny+1, yc+yinc)
         j5p = min(ny, yc+yinc+1)
         j6p = min(ny+1, yc+yinc+1)
c         j6 = min(ny+1, yc+yinc+1)
         j7 = yc+yinc
         j8 = min(ny+1, yc+yinc)
      endif

c...  We will expand the range of possible responses when perturbing the
c...  plasma in a cell near one of the cuts.
      xccuts = .false.
      do jx = 1, nxpt
        if ( (xc-xlinc<=ixpt1(jx)+1) .and. (xc+xrinc+1>=ixpt1(jx)) .and.
     .       (yc-yinc<=iysptrx1(jx)) .and. (iysptrx1(jx)>0) ) xccuts=.true.
        if ( (xc-xlinc<=ixpt2(jx)+1) .and. (xc+xrinc+1>=ixpt2(jx)) .and.
     .       (yc-yinc<=iysptrx2(jx)) .and. (iysptrx2(jx)>0) ) xccuts=.true.
      enddo

c...  We must expand the range of ix in the vicinity of cells on
c...  which turbulent transport depends.
      xcturb = .false.
      do jx = 1, nxpt
         xcturb = xcturb .or. (xc.eq.ixlb(jx).and.ixmnbcl==1) .or.
     .                        (xc.eq.(ixrb(jx)+1).and.ixmxbcl==1)
      enddo
      xcturb = xcturb .or. (xc.eq.ixmp2)
      xcturb = xcturb .and. (isturbnloc.eq.1)
c...  NOTE: For a full double-null configuration, if there are 2 separatrices
c...  we use the innermost one (at iysptrx) to define the radial boundary
c...  of the turbulence.
      if (isturbcons .eq. 1) then
         xcturb = xcturb .and. yc .eq. iysptrx+1
      elseif (isturbcons .eq. 2) then
         xcturb = xcturb .and. yc .ge. iysptrx+1-diffuslimit
      else
         xcturb = xcturb .and. yc .ge. iysptrx+1
      endif

      if (xccuts .or. xcturb) then
         i1 = 0
         i2 = 1
         i3 = 0
         i4 = 0
         i5 = nx
         i6 = nx+1
         i7 = nx+1
         i8 = nx+1
      endif

c...  Define range for source terms to minimize calls to adpak-based routines
            ixs = i2
            ixf = i5
            iys = j2
            iyf = j5
            ixs1 = i1
            ixf6 = i6
            iys1 = j1
            iyf6 = j6
c...  Reset ioniz. and rad. indices to a point if this is a Jacobian calc.
         if (xc .ge.0 .and. yc .ge. 0) then
            ixs = xc
            ixf = xc
            iys = yc
            iyf = yc
            ixs1 = xc
            ixf6 = xc
            if (xrinc .ge. 20) then
               ixs1 = 0
               ixf6 = nx+1
            endif
            iys1 = yc
            iyf6 = yc
            if (yinc .ge. 20) then
               iys1 = 0
               iyf6 = ny+1
            endif
         endif

c...  Set flag that indicates wide open Jac'n "box" for subroutine bouncon.
      if (xc .lt. 0) then
         openbox = .true.
      elseif (xccuts .or. xcturb) then
         openbox = .true.
      elseif ( (0<=yc).and.(yc<=yinc) ) then # for integrated core flux BC
         openbox = .true.
      else
         openbox = .false.
      endif

c...  Set flags that indicate when the Jac'n "box" overlaps left or right
c...  boundary cells of a mesh region.  Used in subroutine bouncon.
      xcnearlb = .false.
      do jx = 1, nxpt
         xcnearlb = xcnearlb .or.
     .       ( (xc-xlinc.le.ixlb(jx)) .and. (xc+xrinc.ge.ixlb(jx)) )
      enddo
      if (xc .lt. 0) xcnearlb = .true.
      xcnearrb = .false.
      do jx = 1, nxpt
         xcnearrb = xcnearrb .or.
     .      ( (xc-xlinc.le.ixrb(jx)+1) .and. (xc+xrinc.ge.ixrb(jx)) )
      enddo
      if (xc .lt. 0) xcnearrb = .true.

************************************************************************
c... First, we convert from the 1-D vector yl to the plasma variables.
************************************************************************
         if (TimingPandfOn.gt.0) TimeConvert0=tick()
         call convsr_vo (xc, yc, yl)  # pre 9/30/97 was one call to convsr
         if (TimingPandfOn.gt.0) TotTimeConvert0=TotTimeConvert0+tock(TimeConvert0)
         if (TimingPandfOn.gt.0) TimeConvert1=tick()
         call convsr_aux (xc, yc)
         if (TimingPandfOn.gt.0) TotTimeConvert1=TotTimeConvert1+tock(TimeConvert1)

c ... Set variable controlling upper limit of species loops that
c     involve ion-density sources, fluxes, and/or velocities.

      nfsp = nisp
      if (isimpon .eq. 3 .or. isimpon .eq. 4) nfsp = nhsp

c ... Calculate the Bohm diffusion rates (units are m**2/s)
      do ifld = 1, nisp
       if (facbni+facbup+facbee+facbei>0 .and. isbohmcalc>0) then
         do iy = j1, j6
            iyp1 = min(ny+1, iy+1)
            do ix = i1, i6
               ix1 = ixp1(ix,iy)
               kybohm(ix,iy) = (te(ix,iy)+te(ix,iyp1)) /
     .                        (16*ev*(btot(ix,iy)+btot(ix,iyp1)))
               kxbohm(ix,iy) = (te(ix,iy)+te(ix1,iy)) /
     .                        (16*ev*(btot(ix,iy)+btot(ix1,iy)))
               dif_use(ix,iy,ifld)  = facbni*kybohm(ix,iy)
               dif2_use(ix,iy,ifld) = facbni2*kxbohm(ix,iy)
               tray_use(ix,iy,ifld)  = facbup*kybohm(ix,iy)
               kye_use(ix,iy)  = facbee*kybohm(ix,iy)
               kyi_use(ix,iy)  = facbei*kybohm(ix,iy)
	       dutm_use(ix,iy,ifld) = kybohm(ix,iy)
            enddo
         enddo
         if (isbohmcalc.eq.2) then  # calc. recip. average with const D
           fcdif = 0.   # used to zero constant diff. if recip. ave used
           do iy = j1, j6
             do ix = i1, i6
               dif_use(ix,iy,ifld)  = 0.5*ave(difni(ifld),  dif_use(ix,iy,ifld))
               dif2_use(ix,iy,ifld) = 0.5*ave(difni2(ifld), dif2_use(ix,iy,ifld))
               tray_use(ix,iy,ifld)  = 0.5*ave(travis(ifld), tray_use(ix,iy,ifld))
               kye_use(ix,iy)  = 0.5*ave(kye, kye_use(ix,iy))
               kyi_use(ix,iy)  = 0.5*ave(kyi, kyi_use(ix,iy))
               dutm_use(ix,iy,ifld)  = 0.5*ave(difutm(1), dutm_use(ix,iy,ifld))
             enddo
           enddo
         endif
       endif

c ... If isbohmcalc=3, then give (B0/B)**inbdif profile to diff
       if (isbohmcalc==3) then  # use inbtdif, inbpdif for btot, bpol scaling
         bpolmin = bpol(ixpt2(1)-ixbpmin,iysptrx,0)
         do iy = j1, j6
            do ix = i1, i6
              ix1 = ixp1(ix,iy)
	      bscalf=((.5*(btot(ixmp,iysptrx)/btot(ix,iy)) +
     .               .5*(btot(ixmp,iysptrx)/btot(ix1,iy)))**inbtdif)
     .       * ((bpol(ixmp,iysptrx,3)+bpol(ixmp,iysptrx,4))/
     .          (bpol(ix,iy,3)+bpol(ix,iy,4)+bpolmin))**inbpdif
	      dif_use(ix,iy,ifld)  = difniv(iy,ifld)*bscalf
	      difp_use(ix,iy,ifld) = difprv(iy,ifld)*bscalf
              dif2_use(ix,iy,ifld) = difniv2(iy,ifld)*bscalf
              tray_use(ix,iy,ifld)  = travisv(iy,ifld)*bscalf
              kye_use(ix,iy)  = kyev(iy)*bscalf
              kyi_use(ix,iy)  = kyiv(iy)*bscalf
              dutm_use(ix,iy,ifld) = difutmv(iy,ifld)*bscalf
	      vy_use(ix,iy,ifld) = vconyv(iy,ifld)*bscalf
            enddo
         enddo
       endif
      enddo  # loop over species lfld

c ,,, Add diffusion propto betap**iexpbp and (B0/B)**inbdif (as for isbohmcalc=3)
      if (isdifbetap == 1) then
       do ifld = 1, nisp
         if(zi(ifld) > 0.) then
           bpolmin = bpol(ixpt2(1)-ixbpmin,iysptrx,0)
           do iy = j1, j6
             do ix = i1, i6
               ix1 = ixp1(ix,iy)
               betap(ix,iy) = 8.*pi*1.e-7*pr(ix,iy)/bpol(ix,iy,0)**2
               bpfac = betap(ix,iy)**iexpbp
 	       bscalf = ((.5*(btot(ixmp,iysptrx)/btot(ix,iy)) +
     .                   .5*(btot(ixmp,iysptrx)/btot(ix1,iy)))**inbtdif)
     .                  * ((bpol(ixmp,iysptrx,3)+bpol(ixmp,iysptrx,4))/
     .                   (bpol(ix,iy,3)+bpol(ix,iy,4)+bpolmin))**inbpdif
	       dif_use(ix,iy,ifld)  = difniv(iy,ifld)*bscalf +
     .                                                    dfacbp*bpfac
	       difp_use(ix,iy,ifld) = difprv(iy,ifld)*bscalf
               dif2_use(ix,iy,ifld) = difniv2(iy,ifld)*bscalf +
     .                                                    dfacbp*bpfac
               tray_use(ix,iy,ifld)  = travisv(iy,ifld)*bscalf +
     .                                                   trfacbp*bpfac
               trax_use(ix,iy,ifld) = trfacbp*bpfac
               kye_use(ix,iy)  = kyev(iy)*bscalf +  kefacbp*bpfac
               kyi_use(ix,iy)  = kyiv(iy)*bscalf + kifacbp*bpfac
               kxe_use(ix,iy) = kefacbp*bpfac
               kxi_use(ix,iy) = kifacbp*bpfac
               dutm_use(ix,iy,ifld) = difutmv(iy,ifld)*bscalf
	       vy_use(ix,iy,ifld) = vconyv(iy,ifld)*bscalf
             enddo
           enddo
         endif
       enddo
      endif   # test on isdifbetap


************************************************************************
*  Transverse Drifts in y-direction and in 2-direction
*  (normal to B and y)
************************************************************************
*  ---------------------------------------------------------------------
*  compute drifts
*  ---------------------------------------------------------------------

c ... Compute log_lambda
      do iy = j1, j6
        do ix = i1, i6
          teev = te(ix,iy)/ev
          if (islnlamcon == 1) then
            loglambda(ix,iy) = lnlam  # set to constant
          elseif (teev < 50.) then    # Braginskii formula, teev<50
            loglambda(ix,iy) = 23.4-1.15*log10(1.e-6*ne(ix,iy))+
     .                              3.45*log10(teev)
          else                        #teev > 50
            loglambda(ix,iy) = 25.3-1.15*log10(1.e-6*ne(ix,iy))+
     .              2.33167537087122D+00*log10(teev)
          endif
          ctaui(ix,iy,1) = 2.1e13*sqrt(mi(1)/mp)/ loglambda(ix,iy)
          ctaue(ix,iy,1) = 3.5e11/loglambda(ix,iy) #both for zi=1
        enddo
      enddo

c ... Calculate collis. factors eta1 and rtaue for the simple Braginski model
      do iy = j1, j6
        do ix = i1, i6
           eta1(ix,iy) = cfeta1*0.3*nm(ix,iy,1)*ti(ix,iy)*
     .                   (1/(qe*btot(ix,iy))) / omgci_taui
           rtaue(ix,iy) = cfrtaue*(1/(qe*btot(ix,iy))) / omgce_taue
           dclass_i(ix,iy) = cfcl_i*eta1(ix,iy)/(0.3*nm(ix,iy,1))
           dclass_e(ix,iy) = cfcl_e*te(ix,iy)*rtaue(ix,iy)
        enddo
      enddo

*  -- loop over species number --

      do 100 ifld = 1, nfsp
c --- If this is the neutral species (zi(ifld).eq.0)) we dont want velocities
        if(zi(ifld) > 1.e-10) then  # if not, skip to end of 100 loop
         qion = zi(ifld)*qe
         do 18 iy = j1, j5
            iyp1 = min(iy+1,ny+1)
            iym1 = max(iy-1,0)
            do 17 ix = i1, i6
              ix3 = ixm1(ix,iy)
              ix4 = ixm1(ix,iy+1)
              temp1 =
     .                ( ex(ix ,iy) + ex(ix ,iy+1) +
     .                  ex(ix3,iy) + ex(ix4,iy+1) )
c... sknam: grad P from priv, prev
              temp1 = (-4.0)*(phiv(ix,iy) - phiv(ix3,iy))*gxc(ix,iy)
              temp2 = 4.0*(priv(ix,iy,ifld) - priv(ix3,iy,ifld))*gxc(ix,iy)
              temp3 = 4.0*(prev(ix,iy) - prev(ix3,iy))*gxc(ix,iy)

c...MER NOTE: For a full double-null configuration, the following test will
c...  use the radial index of the innermost separatrix (see iysptrx definition
c...  in subroutine nphygeo)
              if ( isxpty(ix,iy)==0 .and. iysptrx.gt.0 ) then
                temp1 = (-4.0)*(phiv(ix,iy) - phiv(ix3,iy))*gxc(ix,iy)
                temp2 = 4.0*(priv(ix,iy,ifld) - priv(ix3,iy,ifld))*gxc(ix,iy)
                temp3 = 4.0*(prev(ix,iy) - prev(ix3,iy))*gxc(ix,iy)
              endif
c...    Calc collisionality factors nu_s/(1 + nu_s) = 1/(1 + lambda_s)
              lambd_ci = 1e16*(ti(ix,iy)/ev)**2/nit(ix,iy)  # approx
              lambd_ce = 2e16*(te(ix,iy)/ev)**2/ne(ix,iy)   # approx
              coll_fi(ix,iy) = cfnus_i/(cfnus_i + (lambd_ci/(lconi(ix,iy))))
              coll_fe(ix,iy) = cfnus_e/(cfnus_e + (lambd_ce/(lcone(ix,iy))))
              vyce(ix,iy,ifld) = 0.125 * temp1
     .                          * ( rbfbt2(ix,iy) + rbfbt2(ix,iy+1) )
              vycb(ix,iy,ifld) = ( cfcurv*( 0.5*(ti(ix,iy)+ti(ix,iyp1)) +
     .                mi(ifld)*(0.25*(up(ix,iy,ifld)+up(ix,iyp1,ifld)+
     .                         up(ix3,iy,ifld)+up(ix4,iyp1,ifld)))**2 )*
     .                                       curvrby(ix,iy)/qion +
     .                   cfgradb*0.5*( ti(ix,iy)+ti(ix,iyp1) )*
     .                               gradby(ix,iy)/qion )*coll_fi(ix,iy)
              veycb(ix,iy) = ( -cfcurv*0.5*(te(ix,iy)+te(ix,iyp1))*
     .                                               curvrby(ix,iy)/qe -
     .                   cfgradb*0.5*( te(ix,iy)+te(ix,iyp1) )*
     .                                  gradby(ix,iy)/qe )*coll_fe(ix,iy)

              vycp(ix,iy,ifld) = -0.25 * temp2
     .               * (rbfbt2(ix,iy)+rbfbt2(ix,iy+1)) /
     .            (qion*(niy0(ix,iy,ifld)+niy1(ix,iy,ifld)))
              veycp(ix,iy) =  0.25 * temp3
     .               * (rbfbt2(ix,iy)+rbfbt2(ix,iy+1)) /
     .                    (qe*(ney0(ix,iy)+ney1(ix,iy)))
c...   zero the vy-diamagnetic velocity on the y guard-cell faces
              vycp(ix,0,ifld) = 0.
              vycp(ix,ny,ifld) = 0.
              veycp(ix,0) = 0.
              veycp(ix,ny) = 0.

c...  Precompute radial velocities from fixed BOUT turbulence fluxes
              vy_cft(ix,iy,ifld) = 2*fniyos_use(ix,iy,ifld)/
     .                            (niy0(ix,iy,ifld)+niy1(ix,iy,ifld))
              vyte_cft(ix,iy) = 2*feeyosn_use(ix,iy)/
     .                            (tey0(ix,iy)+tey1(ix,iy))
              vyti_cft(ix,iy) = 2*feiyosn_use(ix,iy)/
     .                            (tiy0(ix,iy)+tiy1(ix,iy))
              vyrd(ix,iy,ifld) = - 2. * gpry(ix,iy) /
     .                       (btot(ix,iy)**2/etaper(ix,iy) +
     .                        btot(ix,iy+1)**2/etaper(ix,iy+1) )
              vydd(ix,iy,ifld) = vcony(ifld) + vy_use(ix,iy,ifld) +
     .                                         vy_cft(ix,iy,ifld) -
     .                         (difpr(ifld) + difp_use(ix,iy,ifld)) *
     .                    ( 2*gpry(ix,iy)/(pr(ix,iy+1) + pr(ix,iy)) -
     .                      3.0*gtey(ix,iy)/(tey1(ix,iy)+tey0(ix,iy)) )
c ...   Note that the density grad. term for vydd added below
           if (cfrtaue.ne.0.) then  #special classical mom. transfer term
              vycr(ix,iy) = -0.5*(rtaue(ix,iy)+rtaue(ix,iyp1)) * (
     .                        (gpiy(ix,iy,1) + gpey(ix,iy))/
     .                        (0.5*(niy1(ix,iy,1)+niy0(ix,iy,1))) -
     .                          1.5*gtey(ix,iy) )
           endif
           if (cfeta1.ne.0. .and. iy.le.ny-1 .and. iy.gt.0) then
                                             #special classical vis. term
              geyym = 2*gpiy(ix,iym1,1)/(ney1(ix,iym1)+ney0(ix,iym1))-
     .                qe*ey(ix,iym1)
              geyy0 = 2*gpiy(ix,iy,1)/(ney1(ix,iy)+ney0(ix,iy)) -
     .                qe*ey(ix,iy)
              geyyp = 2*gpiy(ix,iyp1,1)/(ney1(ix,iyp1)+ney0(ix,iyp1))-
     .                qe*ey(ix,iyp1)
              dgeyy0 = (geyy0-geyym)*eta1(ix,iy)*gy(ix,iy)
              dgeyy1 = (geyyp-geyy0)*eta1(ix,iyp1)*gy(ix,iyp1)
              vycf(ix,iy) = 2*(dgeyy1-dgeyy0)*gy(ix,iy) /
     .                                 ( (ney1(ix,iy)+ney0(ix,iy))*
     .                     (qe*0.5*(btot(ix,iy)+btot(ix,iym1)))**2 )
           endif

              diffusivwrk(ix,iy)=fcdif*difni(ifld)+dif_use(ix,iy,ifld)
  17       continue
  18     continue

c
c ... Compute diffusive part of radial velocity.
c .. Needs further cleaning; no turbulence model used now TDR 9/1/15
         do iy = j1, j5
            do ix = i1, i6
              difnimix = diffusivwrk(ix,iy)

c ... Alter diffusivity in the SOL by mixing fixed diffusivity
c     with anomalous diffusivity computed in subroutine turb_diffus but
c     reduced by the factor difnit(ifld).  The mixing ratio is given by
c     cdifnit.  Diffusivity is unaltered if difnit(ifld) = 0.
c...MER NOTE: For a full double-null configuration, the SOL is defined to
c...  be the region outside the innermost separatrix (see iysptrx definition
c...  in subroutine nphygeo)
cc              if (difnit(ifld) .gt. 1.e-20 .and. zi(ifld) .eq. 1.
cc     .                                 .and. iy .gt. iysptrx) then
cc                 difnimix = (1. - cdifnit) *
cc     .                      (fcdif*difni(ifld) + dif_use(ix,iy,ifld)) +
cc     .                               cdifnit * difnit(ifld) * difnimix
cc              endif

              vydd(ix,iy,ifld) = vydd(ix,iy,ifld)
     .           -1. * difnimix * (
     .            2*(1-isvylog)*(niy1(ix,iy,ifld) - niy0(ix,iy,ifld)) *
     .              gyf(ix,iy) / (niy1(ix,iy,ifld) + niy0(ix,iy,ifld))+
     .              isvylog*(log(niy1(ix,iy,ifld)) -
     .                             log(niy0(ix,iy,ifld))) *gyf(ix,iy) )

c ... Compute total radial velocity.
              vy(ix,iy,ifld) = cfydd *bfacyrozh(ix,iy) *
     .                                 vycp(ix,iy,ifld) +
     .                         cfrd  * vyrd(ix,iy,ifld) +
     .                                 vydd(ix,iy,ifld) +
     .                         cfyef * vyce(ix,iy,ifld) +
     .                         cfybf * vycb(ix,iy,ifld) +
     .                        cfvycf * vycf(ix,iy) +
     .                        cfvycr * vycr(ix,iy)
c ... Compute radial vel v_grad_P eng eqn terms;cfydd+cfybf=1 or 0
              vygp(ix,iy,ifld) = (cfydd+cfybf)*bfacyrozh(ix,iy) *
     .                                         vycp(ix,iy,ifld) +
     .                                 cfrd  * vyrd(ix,iy,ifld) +
     .                                         vydd(ix,iy,ifld) +
     .                                 cfyef * vyce(ix,iy,ifld) +
     .                                cfvycf * vycf(ix,iy) +
     .                                cfvycr * vycr(ix,iy)
              if (isybdrywd == 1) then  #make vy diffusive in wall cells
                 if (iy==0 .and. matwalli(ix) > 0) then
                    vy(ix,iy,ifld) = vydd(ix,iy,ifld)
                 elseif (iy==ny .and. matwallo(ix) > 0) then
                    vy(ix,iy,ifld) = vydd(ix,iy,ifld)
                 endif
              endif
            enddo  #loop over iy
         enddo     #loop over ix

	 do 20 iy = j1, j6
	    do 19 ix = i1, i6
	      iy1 = max(0,iy-1)            # does iy=0 properly
              iy2 = min(ny+1,iy+1) # use ex*fqx since phi(0,) may be large
	      ix2 = ixp1(ix,iy)
	      ix4 = ixp1(ix,iy1)
              ix6 = ixp1(ix,iy2)
              do jx = 1, nxpt
                 if (ix==ixlb(jx) .and. ixmnbcl==1) then
                    temp1 =  (-4.)* (phiv(ix,iy) - phiv(ix,iy1))*gyc(ix,iy)
                    temp2 = 4.*(priv(ix,iy,ifld) - priv(ix,iy1,ifld))*gyc(ix,iy)
                    temp3 = 4.*(prev(ix,iy) - prev(ix,iy1))*gyc(ix,iy)
                    temp4 = 2.5*(tiv(ix,iy) - tiv(ix,iy1))*gyc(ix,iy)
                 elseif (ix==ixrb(jx) .and. ixmxbcl==1) then
                    temp1 =  (-4.)* (phiv(ix,iy) - phiv(ix,iy1))*gyc(ix,iy)
                    temp2 = 4.*(priv(ix,iy,ifld) - priv(ix,iy1,ifld))*gyc(ix,iy)
                    temp3 = 4.*(prev(ix,iy) - prev(ix,iy1))*gyc(ix,iy)
                    temp4 = 2.5*(tiv(ix,iy) - tiv(ix,iy1))*gyc(ix,iy)
                 else  # not a boundary
                    temp1 =  (-4.)* (phiv(ix,iy) - phiv(ix,iy1))*gyc(ix,iy)
                    temp2 = 4.*(priv(ix,iy,ifld) - priv(ix,iy1,ifld))*
     .                                                       gyc(ix,iy)
                    temp3 = 4.*(prev(ix,iy) - prev(ix,iy1))* gyc(ix,iy)
                    temp4 = 2.5*(tiv(ix,iy) - tiv(ix,iy1))*gyc(ix,iy)
                endif
              enddo  #vis end do-loop over nxpt mesh regions
c...MER NOTE: For a full double-null configuration, the following test will
c...  use the radial index of the innermost separatrix (see iysptrx definition
c...  in subroutine nphygeo)
              if ( isxptx(ix,iy)==0 .and. iysptrx.gt.0 ) then
                 temp1 =  (-4.)* (phiv(ix,iy) - phiv(ix,iy1))*gyc(ix,iy)
                 temp2 = 4.*(priv(ix,iy,ifld) - priv(ix,iy1,ifld))*
     .                                                       gyc(ix,iy)
              endif

              v2ce(ix,iy,ifld) = - 0.5 * temp1
     .             / ( btot(ix,iy) + btot(ix2,iy) )
              v2cb(ix,iy,ifld) =(cfcurv*( 0.5*(tiv(ix,iy)+tiv(ix,iy1)) +
     .                 mi(ifld)*up(ix,iy,ifld)**2 )*curvrb2(ix,iy) +
     .                     cfgradb*0.5*( tiv(ix,iy)+tiv(ix,iy1) )*
     .                            gradb2(ix,iy))/qion
              ve2cb(ix,iy) = -(cfcurv*0.5*(tev(ix,iy)+tev(ix,iy1))*
     .                                                 curvrb2(ix,iy) +
     .                          cfgradb*0.5*(tev(ix,iy)+tev(ix,iy1))*
     .                            gradb2(ix,iy))/qe
              v2cd(ix,iy,ifld) = temp2
     .                    / ((btot(ix,iy)+btot(ix2,iy))*qion*
     .                           (ni(ix,iy,ifld)+ni(ix2,iy,ifld)))
              ve2cd(ix,iy,1) = -temp3
     .                    / ((btot(ix,iy)+btot(ix2,iy))*qe*
     .                           (ni(ix,iy,ifld)+ni(ix2,iy,ifld)))
              q2cd(ix,iy,ifld) = (priv(ix,iy,ifld)+priv(ix,iy1,ifld))*temp4
     .                    / ( (btot(ix,iy)+btot(ix2,iy))*qion )

c...  Calculate plate electr diamag flux used to find sheath potential
              do jx = 1, nxpt
                 if (ix==ixlb(jx) .and. ixmnbcl==1) then
                 # use ix=ixlb+1 values to avoid BC variations
                    v2dia = -0.5*( gpey(ixlb(jx)+1,iy)+gpey(ixlb(jx)+1,iy1) ) /
     .                          ( btot(ixlb(jx)+1,iy)*qe*ne(ixlb(jx)+1,iy) )
                    fdiaxlb(iy,jx) = ne(ixlb(jx)+1,iy) * sx(ixlb(jx),iy) *
     .                              v2dia * rbfbt(ixlb(jx)+1,iy)
                 endif
                 if (ix==ixrb(jx) .and. ixmxbcl==1) then
                    v2dia = -0.5*( gpey(ixrb(jx),iy)+gpey(ixrb(jx),iy1) ) /
     .                          ( btot(ixrb(jx),iy)*qe*ne(ixrb(jx),iy) )
                    fdiaxrb(iy,jx) = ne(ixrb(jx),iy) * sx(ixrb(jx),iy) *
     .                              v2dia * rbfbt(ixrb(jx),iy)
                 endif
              enddo  # end do-loop over nxpt mesh regions

              v2rd(ix,iy,ifld) = - 2. * gprx(ix,iy) /
     .           ( btot(ix,iy)/(etaper(ix,iy)*rbfbt2(ix,iy)) +
     .             btot(ix2,iy)/(etaper(ix2,iy)*rbfbt2(ix2,iy)) )
              v2dd(ix,iy,ifld) = - 2. * difpr2(ifld) * gprx(ix,iy) /
     .                                ( pr(ix2,iy)/rbfbt(ix2,iy) +
     .                                  pr(ix,iy)/rbfbt(ix,iy) ) -
     .            2. * (fcdif*difni2(ifld) + dif2_use(ix,iy,ifld)) *
     .                               (ni(ix2,iy,ifld)-ni(ix,iy,ifld)) /
     .                      (ni(ix2,iy,ifld)/(rbfbt(ix2,iy)*gx(ix2,iy))+
     .                          ni(ix,iy,ifld)/(rbfbt(ix,iy)*gx(ix,iy)))
              v2(ix,iy,ifld) = cf2dd * bfacxrozh(ix,iy) *
     .                                 v2cd(ix,iy,ifld) +
     .                         cfrd  * v2rd(ix,iy,ifld) +
     .                                 v2dd(ix,iy,ifld) +
     .                         cf2ef * v2ce(ix,iy,ifld) +
     .                         cf2bf * v2cb(ix,iy,ifld)
c ...         Compute v2 for v2x_gradx_P eng terms; cf2dd+cf2bf=1 or 0
              v2xgp(ix,iy,ifld) =  0.5*(rbfbt(ix,iy)+rbfbt(ix2,iy)) * (
     .                 (cf2dd+cf2bf) * bfacxrozh(ix,iy) *
     .                                 v2cd(ix,iy,ifld) +
     .                         cfrd  * v2rd(ix,iy,ifld) +
     .                                 v2dd(ix,iy,ifld) +
     .                         cf2ef * v2ce(ix,iy,ifld) )
         if (isnonog.eq.1 .and. iy.le.ny) then
c            grdnv = ( 1/( fym (ix,iy,1)/ni(ix2,iy1,ifld) +
c     .                    fy0 (ix,iy,1)/ni(ix2,iy ,ifld) +
c     .                    fyp (ix,iy,1)/ni(ix2,iy2,ifld) +
c     .                    fymx(ix,iy,1)/ni(ix ,iy1,ifld) +
c     .                    fypx(ix,iy,1)/ni(ix, iy2,ifld) ) -
c     .                1/( fym (ix,iy,0)/ni(ix ,iy1,ifld) +
c     .                    fy0 (ix,iy,0)/ni(ix ,iy ,ifld) +
c     .                    fyp (ix,iy,0)/ni(ix ,iy2,ifld) +
c     .                    fymx(ix,iy,0)/ni(ix4,iy1,ifld) +
c     .                    fypx(ix,iy,0)/ni(ix6,iy2,ifld) ) )
c     .                                                 * gxfn(ix,iy)
cc            grdnv = ( exp( fym (ix,iy,1)*log(ni(ix2,iy1,ifld)) +
cc     .                     fy0 (ix,iy,1)*log(ni(ix2,iy ,ifld)) +
cc     .                     fyp (ix,iy,1)*log(ni(ix2,iy2,ifld)) +
cc     .                     fymx(ix,iy,1)*log(ni(ix ,iy1,ifld)) +
cc     .                     fypx(ix,iy,1)*log(ni(ix, iy2,ifld)) )
cc     .               -exp( fym (ix,iy,0)*log(ni(ix ,iy1,ifld)) +
cc     .                     fy0 (ix,iy,0)*log(ni(ix ,iy ,ifld)) +
cc     .                     fyp (ix,iy,0)*log(ni(ix ,iy2,ifld)) +
cc     .                     fymx(ix,iy,0)*log(ni(ix4,iy1,ifld)) +
cc     .                     fypx(ix,iy,0)*log(ni(ix6,iy2,ifld)) ) ) *
cc     .                                                      gxfn(ix,iy)
            grdnv = (    ( fym (ix,iy,1)*log(ni(ix2,iy1,ifld)) +
     .                     fy0 (ix,iy,1)*log(ni(ix2,iy ,ifld)) +
     .                     fyp (ix,iy,1)*log(ni(ix2,iy2,ifld)) +
     .                     fymx(ix,iy,1)*log(ni(ix ,iy1,ifld)) +
     .                     fypx(ix,iy,1)*log(ni(ix, iy2,ifld)) )
     .                  -( fym (ix,iy,0)*log(ni(ix ,iy1,ifld)) +
     .                     fy0 (ix,iy,0)*log(ni(ix ,iy ,ifld)) +
     .                     fyp (ix,iy,0)*log(ni(ix ,iy2,ifld)) +
     .                     fymx(ix,iy,0)*log(ni(ix4,iy1,ifld)) +
     .                     fypx(ix,iy,0)*log(ni(ix6,iy2,ifld)) ) ) *
     .                                                      gxfn(ix,iy)
            vytan(ix,iy,ifld)=(fcdif*difni(ifld) + dif_use(ix,iy,ifld)) *
     .                                      (grdnv/cos(angfx(ix,iy)) -
     .                       (log(ni(ix2,iy,ifld)) - log(ni(ix,iy,ifld)))
     .                                                 * gxf(ix,iy) )
            if (islimon.eq.1.and. ix.eq.ix_lim.and. iy.ge.iy_lims) then
              vytan(ix,iy,ifld) = 0.
            endif
            if (nxpt==2 .and. ix==ixrb(1)+1 .and. ixmxbcl==1) then
c             non-physical interface between upper target plates for dnull
              vytan(ix,iy,ifld) = 0.
            endif
         endif

 19        continue
 20     continue

         do 21 ix = i1, i6
            vy(ix,ny+1,ifld) = 0.0
 21      continue
        else    # test on zi > 1.e-10 to skip whole loop
c                  vy(:,:,ifld)=0
c                  vytan(:,:,ifld)=0
c                  v2cb(:,:,ifld)=0
c                  v2cd(:,:,ifld)=0
c                  v2ce(:,:,ifld)=0
c                  ve2cd(:,:,ifld)=0
c                  q2cd(:,:,ifld)=0
c                  v2rd(:,:,ifld)=0
c                  vygp(:,:,ifld)=0
        endif
  100 continue  # Giant loop over ifld (species)

c ... Save values returned by Hirshman mombal for Jacobian calc. to
c ... minimize calls - restore the "m" or ix-1 values at the end of pandf
c ... The Jacobian ix loop can then be reduced to only include ix-1 and ix
c ... Suffix "o" refers to "old" value at ix, and suffix "om" means "old"
c ... value at ix-1.

      if (xc.ge.0 .and. yc.ge.0) then
         ix1 = ixm1(xc,yc)
         fqpom = fqp(ix1,yc)
         friceom = frice(ix1,yc)
         upeom = upe(ix1,yc)
         fqpo = fqp(xc,yc)
         friceo = frice(xc,yc)
         upeo = upe(xc,yc)
         do ifld = 1, nfsp
            friciom(ifld) = frici(ix1,yc,ifld)    # dimension req. nfsp<101
            upiom(ifld) = upi(ix1,yc,ifld)
            uupom(ifld) = uup(ix1,yc,ifld)
            fricio(ifld) = frici(xc,yc,ifld)
            upio(ifld) = upi(xc,yc,ifld)
            uupo(ifld) = uup(xc,yc,ifld)
         enddo
      endif


c ... Need to calculate new currents (fqp) after saving old & before frice,i
      if(isphion+isphiofft .eq. 1)  call calc_currents

c ... Add anomalous perp vis vy using calc_currents result - awkward,change
      if (cfvyavis > 0.) then
        do ifld = 1, 1  # nfsp  # only good for ifld=1
          do iy = max(j1,2), min(j5,ny-1)
            do ix = max(i1,2), min(i6,nx-1)
              vyavis(ix,iy,ifld) = fqya(ix,iy)*2/(
     .                  qe*(niy1(ix,iy,1)+niy0(ix,iy,1))*sy(ix,iy) )
              vy(ix,iy,ifld) = vy(ix,iy,ifld) + cfvyavis*vyavis(ix,iy,ifld)
            enddo
          enddo
        enddo
      endif

c ... Calculate friction forces from Braginskii, if isimpon .ne. 5.

      if (isimpon.ne.5 .and. isimpon.ne.6 .and. isimpon.ne.7) then
         do iy = j1, j6    #iys1, iyf6
            do ix = i1, i6
               ix2 = ixp1(ix,iy)
               nbarx = 0.5*(ne(ix,iy)+ne(ix2,iy))
               ltmax = min( abs(te(ix,iy)/(rrv(ix,iy)*gtex(ix,iy) + cutlo)),
     .                   lcone(ix,iy) )
               lmfpe = 2e16*(te(ix,iy)/ev)**2/ne(ix,iy)
               flxlimf = flalftf*ltmax/(flalftf*ltmax + lmfpe)
               frice(ix,iy) = -cthe*flxlimf*nbarx*rrv(ix,iy)*gtex(ix,iy) +
     .                  cfnetap*qe*netap(ix,iy)*fqp(ix,iy)/sx(ix,iy)
               frici(ix,iy,1) = - frice(ix,iy)
               if (fac2sp .gt. 1.1 .and. nusp .eq. 2) then
                  frici(ix,iy,1) = - frice(ix,iy)/fac2sp
                  frici(ix,iy,2) = - frice(ix,iy)/fac2sp
               endif
            enddo
         enddo
      endif

c ... For use within subroutine mombal, the poloidal electric field is
c     calculated from || Ohms law if isphion = 0.  This field is not
c     intended for use in computing cross-field drifts, so a test of
c     cfyef is included. Both isphiofft=0 or 1 cases included in one loop

      if (isphion .eq. 0) then   # ex calc here assumes no parallel current
         do iy = iys1, iyf6
            do ix = i1, i6
               ix1 = ix
               do jx = 1, nxpt
                 if (ix==ixlb(jx) .and. ixmnbcl==1) then
                   ix1 = ixlb(jx) + 1
                 elseif (ix==ixrb(jx) .and. ixmxbcl==1) then
                   ix1 = ixrb(jx) - 1
                 endif
               enddo
               ix2 = ixp1(ix1,iy)
               ltmax = min( abs(te(ix,iy)/(rrv(ix,iy)*gtex(ix,iy) + cutlo)),
     .                   lcone(ix,iy) )
               lmfpe = 2e16*(te(ix,iy)/ev)**2/ne(ix,iy)
               flxlimf = flalftf*ltmax/(flalftf*ltmax + lmfpe)
               ex(ix,iy) = (1 - isphiofft) * (
     .                     -( gpex(ix1,iy)/(0.5*(ne(ix2,iy)+ne(ix1,iy)))+
     .                                cthe*flxlimf*gtex(ix1,iy) )/qe ) +
     .                      isphiofft * (
     .                            (phi(ix1,iy)-phi(ix2,iy))*gxf(ix1,iy) )
            enddo
         enddo
      endif

c ... Loop over cells (restricted to poloidal slice of box if doing
c     Jacobian), calling mombal if it is providing parallel flow, and
c     taking poloidal projection of parallel flow to get poloidal flow.
c     Unperturbed values of the parallel-flow contribution to uu are
c     saved here so they can be restored below.



      do iy = iys1, iyf6
         if (xc .gt. 0) then
            ix = xc
            ix1 = ixm1(ix,iy)
            if (isimpon .eq. 5) then   # Hirshmans reduced-ion approx.
               if (istimingon .eq. 1) tsimp = gettime(sec4)
               call mombal (ix1,ix,iy)
               if (istimingon .eq. 1) call timimpfj (tsimp, xc)
            elseif(isimpon .eq. 6 .or. isimpon .eq. 7) then # Force balance without inertia
               if (istimingon .eq. 1) tsimp = gettime(sec4)
               if (TimingPandfOn.gt.0) Timemombalni=tick()
                call mombalni (ix1,ix,iy)
          if (TimingPandfOn.gt.0) then
          TotTimemombalni=TotTimemombalni+tock(Timemombalni)
                endif
               if (istimingon .eq. 1) call timimpfj (tsimp, xc)
            endif
            do ifld = 1, nfsp
               if (ifld .le. nusp) then
                 upi(ix1,iy,ifld) = up(ix1,iy,ifld)
               else
                 do jx = 1, nxpt
                    if ( (ix1==ixlb(jx).and.ixmnbcl==1) .or.
     .                   (ix1==ixrb(jx).and.ixmxbcl==1) ) then
                       # constrain boundary velocity
                       if (zi(ifld) .gt. 1.e-10) then
                          argx = abs((2-2*upi(ix1,iy,ifld)/
     .                                     (upi(ix1,iy,1)+cutlo3))**3)
                          argx = min(20., argx)
                          upi(ix1,iy,ifld) = upi(ix1,iy,1) +
     .                       (upi(ix1,iy,ifld)-upi(ix1,iy,1))*exp(-argx)
                       endif  # end if-test on zi
                    endif  # end if-test on ix
                 enddo  # end do-loop over nxpt mesh regions
               endif
               uup(ix1,iy,ifld) = rrv(ix1,iy)*upi(ix1,iy,ifld)
            enddo
         endif
         do ix = ixs1, min(ixf6, nx+1-ixmxbcl)
            ix2 = ixp1(ix,iy)
            if (isimpon .eq. 5) then
               if (istimingon .eq. 1) tsimp = gettime(sec4)
               call mombal (ix,ix2,iy)
               if (istimingon .eq. 1) call timimpfj (tsimp, xc)
            elseif(isimpon .eq. 6 .or. isimpon .eq. 7) then # Force balance without inertia
               if (istimingon .eq. 1) tsimp = gettime(sec4)

               if (TimingPandfOn.gt.0) Timemombalni=tick()
                call mombalni (ix,ix2,iy)
          if (TimingPandfOn.gt.0) then
          TotTimemombalni=TotTimemombalni+tock(Timemombalni)
                endif
               if (istimingon .eq. 1) call timimpfj (tsimp, xc)
            endif
            do ifld = 1, nfsp
               if (ifld .le. nusp) then
                 upi(ix,iy,ifld) = up(ix,iy,ifld)
               else
                 do jx = 1, nxpt
                    if ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .                   (ix==ixrb(jx).and.ixmxbcl==1) ) then
                       # constrain boundary velocity
                       if (zi(ifld) .gt. 1.e-10) then
                          argx = abs((2-2*upi(ix,iy,ifld)/
     .                                     (upi(ix,iy,1)+cutlo3))**3)
                          argx = min(20., argx)
                          upi(ix,iy,ifld) = upi(ix,iy,1) +
     .                         (upi(ix,iy,ifld)-upi(ix,iy,1))*exp(-argx)
                       endif  # end if-test on zi
                    endif  # end if-test on ix
                 enddo  # end do-loop over nxpt mesh regions
               endif
               uup(ix,iy,ifld) = rrv(ix,iy)*upi(ix,iy,ifld)
            enddo
         enddo
      enddo

c ... Add contributions to poloidal velocity from cross-field drifts
c     to those from parallel flow.

      do ifld = 1, nfsp
         do iy = j1, j6
            if (i1 .gt. 0) then  # il is initial ix; here for uu(ixm1(i1-1,,)
               ix = i1
               ix1 = ixm1(ix,iy)
               uu(ix1,iy,ifld) = uup(ix1,iy,ifld) +
     .                           0.5 * (rbfbt(ix,iy) + rbfbt(ix1,iy)) *
     .                           v2(ix1,iy,ifld) - vytan(ix1,iy,ifld) -
     .                         difax(ifld) * 0.5 * ( ( 0.5*(
     .                         ni(ix1,iy,ifld)/ni(ix,iy,ifld) +
     .                         ni(ix,iy,ifld)/ni(ix1,iy,ifld)) -1)**2 ) *
     .                        (ni(ix,iy,ifld)-ni(ix1,iy,ifld))*gxf(ix1,iy)
     .                       /(ni(ix,iy,ifld)+ni(ix1,iy,ifld))
               uz(ix1,iy,ifld) = -uup(ix1,iy,ifld)/rrv(ix1,iy)*
     .          0.5*(rbfbt(ix,iy) + rbfbt(ix1,iy)) + sign(1.,b02d(ix,iy))*
     .               (cftef*v2ce(ix1,iy,ifld)+cftdd*v2cd(ix1,iy,ifld))*
     .                                                       rrv(ix1,iy)
            endif
            do ix = i1, i6    #now the remainder of the uu(ix,,)
               ix2 = ixp1(ix,iy)
               uu(ix,iy,ifld) = uup(ix,iy,ifld) +
     .                          0.5 * (rbfbt(ix,iy) + rbfbt(ix2,iy)) *
     .                          v2(ix,iy,ifld) - vytan(ix,iy,ifld) -
     .                         difax(ifld) * 0.5 * ( ( 0.5*(
     .                         ni(ix,iy,ifld)/ni(ix2,iy,ifld) +
     .                         ni(ix2,iy,ifld)/ni(ix,iy,ifld)) -1)**2 ) *
     .                        (ni(ix2,iy,ifld)-ni(ix,iy,ifld))*gxf(ix,iy)
     .                       /(ni(ix2,iy,ifld)+ni(ix,iy,ifld))
               uz(ix,iy,ifld) = -uup(ix,iy,ifld)/rrv(ix,iy)*
     .          0.5*(rbfbt(ix,iy) + rbfbt(ix2,iy)) + sign(1.,b02d(ix,iy))*
     .               (cftef*v2ce(ix,iy,ifld)+cftdd*v2cd(ix,iy,ifld))*
     .                                                       rrv(ix,iy)
            enddo
         enddo
      enddo

c...  If upi not from full ||mom eq (e.g.,isimpon=6), set impurity
c...  uu(ixrb,,) & upi(ixrb,,) via generalized Bohm cond.
      if(isimpon > 0) then
        do jx = 1, nxpt
	  ixt0 = ixlb(jx)
          ixt = ixrb(jx)+1
          ixt1 = ixrb(jx)
	  do ifld = nhsp+1, nfsp
            if(ifld > nusp) then  #species without full mom eqn
	      do iy = j1, j6
c ..          first left plate(s)
                if(isfixlb(jx) == 0) then # set upi for left plate
                  cs = csfacrb(ifld,jx)*sqrt( (te(ixt0,iy) +
     .                            csfacti*ti(ixt0,iy))/mi(ifld) )
                  ueb = cfueb*( cf2ef*v2ce(ixt0,iy,ifld)*rbfbt(ixt0,iy) -
     .                            vytan(ixt0,iy,ifld) )/rrv(ixt0,iy)
	          uu(ixt0,iy,ifld) = -rrv(ixt0,iy)*cs
                  upi(ixt0,iy,ifld) = -(cs - ueb)
                endif
c ..          switch to right plate(s)
                if(isfixrb(jx) == 0) then
                  cs = csfacrb(ifld,jx)*sqrt( (te(ixt1,iy) +
     .                            csfacti*ti(ixt1,iy))/mi(ifld) )
                  ueb = cfueb*( cf2ef*v2ce(ixt1,iy,ifld)*rbfbt(ixt,iy) -
     .                            vytan(ixt1,iy,ifld) )/rrv(ixt1,iy)
	          uu(ixt1,iy,ifld) = rrv(ixt1,iy)*cs
	          uu(ixt,iy,ifld) = uu(ixt1,iy,ifld)
                  upi(ixt1,iy,ifld) = cs - ueb
                  upi(ixt,iy,ifld) = upi(ixt1,iy,ifld)
                endif
              enddo
            endif   #checks if ifld > nusp
          enddo
        enddo
      endif         # checks if isimpon > 0

************************************************************************
*     Calculate the currents fqx, fqy, fq2 and fqp, if isphion = 1
*     or if isphiofft = 1.
************************************************************************
ccc      if(isphion+isphiofft .eq. 1)  call calc_currents

************************************************************************
*     Calculate the electron velocities, vex, upe, ve2, vey
************************************************************************

      do 25 iy = j1, j6
	  do 24 ix = i1, i6
	    vex(ix,iy) = 0.
	    vey(ix,iy) = 0.
   24    continue
   25 continue

      if (isimpon.eq.5) goto 29    # have upe from mombal

      do iy = j1, j6    #iys1, iyf6
         do ix = i1, i6
            upe(ix,iy) = 0.
         enddo
      enddo

      do 27 ifld = 1, nfsp
         do iy = j1, j6    #iys1, iyf6
	    do ix = i1, i6
               ix1 = ixp1(ix,iy)
	       upe(ix,iy) = upe(ix,iy) + upi(ix,iy,ifld)*zi(ifld)*0.5*
     .                      ( ni(ix,iy,ifld)+ni(ix1,iy,ifld) )
            enddo
         enddo
   27 continue
      afqp = 1.
      if (isimpon.eq.6 .or. isimpon.eq.7) afqp = fupe_cur  #allows gradual fix for old cases
      do iy = j1, j6    #iys1, iyf6
         do ix = i1, i6
            ix1 = ixp1(ix,iy)
	    upe(ix,iy) = (upe(ix,iy) -afqp*fqp(ix,iy)/
     .                               (rrv(ix,iy)*sx(ix,iy)*qe))/
     .                             (0.5*( ne(ix,iy)+ne(ix1,iy) ))
         enddo
      enddo

  29  continue

      do 731 iy = j1, j6   # ExB same all species;if cf2dd=1, no imp yet
	 do 730 ix = i1, i6
            ix1 = ixp1(ix,iy)
            vex(ix,iy) = upe(ix,iy)*rrv(ix,iy) +
     .                   (cf2ef*v2ce(ix,iy,1) + cf2bf*ve2cb(ix,iy) +
     .                         cf2dd*bfacxrozh(ix,iy)*ve2cd(ix,iy,1) ) *
     .                           0.5*(rbfbt(ix,iy) + rbfbt(ix1,iy)) -
     .                                               vytan(ix,iy,1)

  730    continue
  731 continue

      do 734 ifld = 1, nfsp
	 do 733 iy = j1, j5
	    do 732 ix = i1, i6   # grad_B will be ok as next fqy is subtr.
	       vey(ix,iy) = vey(ix,iy) + vy(ix,iy,ifld)*zi(ifld)*0.5*
     .                      ( niy0(ix,iy,ifld)+niy1(ix,iy,ifld) )
  732       continue
  733    continue
  734 continue

      do 36 iy = j1, j5
	 do 35 ix = i1, i6
	    vey(ix,iy) = (vey(ix,iy)-cfjve*fqy(ix,iy)/(sy(ix,iy)*qe))/
     .                    (0.5*( ney0(ix,iy)+ney1(ix,iy) ))
   35    continue
   36 continue

c ... if isnewpot=0, vey(,0) needs to be redone since fqy(,0)=0
      if (isnewpot==1) then
        do ix = i1, i6  # ExB vyce same all species
          vey(ix,0) = cfybf*veycb(ix,0) + vydd(ix,0,1) +
     .                cfyef*vyce(ix,0,1)
        enddo
      endif

c ... If isybdrywd = 1, make vey diffusive, just like vy
      if (isybdrywd == 1) then  #make vy diffusive in wall cells
        do ix = i1, i6
          if (matwalli(ix) > 0) vey(ix,0)  = vydd(ix,0,1)
          if (matwallo(ix) > 0) vey(ix,ny) = vydd(ix,ny,1)
        enddo
      endif
       if (TimingPandfOn.gt.0) TimeSource=tick()
************************************************************************
*   We Calculate the source terms now.
************************************************************************
*  ---------------------------------------------------------------------
*  Coefficients for the source terms.
*  ---------------------------------------------------------------------

      do 702 iy = j2, j5
         do 701 ix = i2, i5
            do ifld = 1, nfsp
               snic(ix,iy,ifld) = 0.0
               sniv(ix,iy,ifld) = 0.0
            enddo
            do ifld = 1, nusp
               smoc(ix,iy,ifld) = 0.0
               smov(ix,iy,ifld) = 0.0
            enddo
            seec(ix,iy) = 0.0
            seev(ix,iy) = 0.0
            seic(ix,iy) = 0.0
            seiv(ix,iy) = 0.0
	    psorbgz(ix,iy) = 0.    # diagnostic only
  701    continue
  702 continue

************************************************************************
*  -- steady sources
************************************************************************
*  ---------------------------------------------------------------------
*  volume sources. (old PHYSRC)
*  ---------------------------------------------------------------------

c ... Calculate effective Lyman-alpha optical depth factors used if
c ... istabon=14 or 15 for hydr. rate look-ups set rtauxfac<=0 to bypass
c ----------------------------------------------------------------------
      if (rtauxfac .gt. 0.) then

CC .. FIRST GET THE POLOIDAL OPTICAL DEPTH FACTOR
c ------------------------------------------------
         do iy = 0, ny+1

c ... get optical-depth to left (ix=0) boundary
            rdumx = 0.
            do ix = 0, nx+1
               dr1  = 0.5*dx(ix,iy)
               dr2  = 0.5*dx(ix,iy)
               if(ix==(nx+1)/2) rdumx=0.  #make two plate regions independent
               rdumx = rdumx + rtauxfac*dr1*ng(ix,iy,1)
               rtaux(ix,iy) = rdumx*rt_scal
               rdumx = rdumx + rtauxfac*dr2*ng(ix,iy,1)
            enddo

c ... get optical-depth to right (ix=nx+1) bdry; initial selection of min rtaux
            rdumx = 0.
            do ix = nx+1, 0, -1
               dr1  = 0.5*dx(ix,iy)
               dr2  = 0.5*dx(ix,iy)
               if(ix==(nx+1)/2) rdumx=0.  #make two plate regions independent
               rdumx = rdumx + rtauxfac*dr1*ng(ix,iy,1)
               rtaux(ix,iy) = min(rdumx*rt_scal, rtaux(ix,iy))
               rdumx = rdumx + rtauxfac*dr2*ng(ix,iy,1)
            enddo
         enddo   # iy loop

CC .. NOW GET THE RADIAL OPTICAL DEPTH FACTOR
c --------------------------------------------
         do ix = 0, nx+1

c ... get optical-depth to inner (iy=0) bdry
            rdumy = 0.
            do iy = 0, ny+1
               dr1  = 0.5*dy(ix,iy)
               dr2  = 0.5*dy(ix,iy)
               rdumy = rdumy + rtauyfac*dr1*ng(ix,iy,1)
               rtauy(ix,iy) = rdumy*rt_scal
               rdumy = rdumy + rtauyfac*dr2*ng(ix,iy,1)
            enddo

c ... get optical-depth to outer (iy=ny+1) bdry; selection of min rtau
            rdumy = 0.
            do iy = ny+1, 0, -1
               dr1  = 0.5*dy(ix,iy)
               dr2  = 0.5*dy(ix,iy)
               rdumy = rdumy + rtauyfac*dr1*ng(ix,iy,1)
               rtauy(ix,iy) = min(rdumy*rt_scal, rtauy(ix,iy))
               rtau(ix,iy) = min(rtaux(ix,iy), rtauy(ix,iy))
               rdumy = rdumy + rtauyfac*dr2*ng(ix,iy,1)
            enddo
         enddo   # ix loop

      endif     # test on rtauxfac, skip if rtauxfac is negative

*     The following is a temporary recycling model.

*  -- recalculate particle source psor if ifixpsor=0 --

c...  Initialize save-variables if this is a Jacobian (xc,yc > -1)
         if (xc .ge. 0 .and. yc .ge. 0..and.Saveold.gt.0) then
            psordisold = psordis(xc,yc)
cc            write(*,*) 'Just after psordisold; xc,yc=',xc,yc
            do ifld = 1, nfsp
               psorold(ifld) = psorc(xc,yc,ifld)
               psorxrold(ifld) = psorxr(xc,yc,ifld)
               msorold(ifld) = msor(xc,yc,ifld)
               msorxrold(ifld) = msorxr(xc,yc,ifld)
               nucxiold(ifld) = nucxi(xc,yc,ifld)
               nueliold(ifld) = nueli(xc,yc,ifld)
            enddo
            do igsp = 1, ngsp
               nucxold(igsp) = nucx(xc,yc,igsp)
               nurcold(igsp) = nurc(xc,yc,igsp)
               nuizold(igsp) = nuiz(xc,yc,igsp)
               nuixold(igsp) = nuix(xc,yc,igsp)
               nuelgold(igsp) = nuelg(xc,yc,igsp)
               psorgold(igsp) = psorgc(xc,yc,igsp)
               psorrgold(igsp) = psorrgc(xc,yc,igsp)
               psorcxgold(igsp) = psorcxgc(xc,yc,igsp)
            enddo
         endif

c...  The particle source can be frozen if ifixpsor.ne.0
      if(ifixpsor .eq. 0) then

        igsp = 0
        do ifld = 1, nhsp  # Hydrogen-only loop
         if (zi(ifld) > 0.) then  #calc only for hydrogen ions
          igsp = igsp + 1
          do iy = iys1, iyf6
            do ix = ixs1, ixf6

c     Ionization of neutral hydrogen by electrons and recombination--
               if (icnuiz .eq. 0) then
                  ne_sgvi = ne(ix,iy)
                  if (ifxnsgi.eq.1) ne_sgvi = cne_sgvi  # fix density dependence
                  nuiz(ix,iy,igsp) = chioniz *  ne(ix,iy) * (
     .                           rsa(te(ix,iy),ne_sgvi,rtau(ix,iy),0)
     .                         + sigvi_floor )
                  if (xc .ge. 0) then        # limit Jacobian element
                     nuiz(ix,iy,igsp) = fnnuiz*nuiz(ix,iy,igsp) +
     .                               (1-fnnuiz)*nuizold(igsp)
                  endif
               elseif (icnuiz .eq. 1) then
                  nuiz(ix,iy,igsp) = cnuiz
               endif
               if (isrecmon == 1) then
                  nurc(ix,iy,igsp) = cfrecom * ne(ix,iy)
     .                         * rra(te(ix,iy),ne(ix,iy),rtau(ix,iy),1)
                  if (xc .ge. 0) then        # limit Jacobian element
                     nurc(ix,iy,igsp) = fnnuiz*nurc(ix,iy,igsp) +
     .                             (1-fnnuiz)*nurcold(igsp)
                  endif
               else
                   nurc(ix,iy,igsp) = 0.
               endif
               psorbgg(ix,iy,igsp) = ngbackg(igsp)*( (0.9 + 0.1*
     .                            (ngbackg(igsp)/ng(ix,iy,igsp))**ingb) ) *
     .                             nuiz(ix,iy,igsp) * vol(ix,iy)
               psorgc(ix,iy,igsp) = -ngcap(ix,iy,igsp)*nuiz(ix,iy,igsp)*vol(ix,iy) +
     .                              psorbgg(ix,iy,igsp)
               psorc(ix,iy,ifld) = - psorgc(ix,iy,igsp)
               psordis(ix,iy) = psorc(ix,iy,1)  # changed below if ishymol=1
               psorxrc(ix,iy,ifld) = -nicap(ix,iy,ifld)*nurc(ix,iy,igsp)*vol(ix,iy)
               psorrgc(ix,iy,igsp) = -psorxrc(ix,iy,ifld)
               msor(ix,iy,ifld) = 0.
               msorxr(ix,iy,ifld) = 0.


c     Charge exchange on neutral hydrogen --
              if (icnucx .eq. 0) then
	         t0 = max(ti(ix,iy),temin*ev)
ccc   we omit the weak velocity dependence as it brings in ni(ix+1) in Jac
                 t1 = t0/(mi(ifld)/mp)
                 nucx(ix,iy,igsp) = ni(ix,iy,ifld) * rcx(t1,ni(ix,iy,ifld),1)
              elseif (icnucx .eq. 1) then
                 nucx(ix,iy,igsp) = cnucx
              elseif (icnucx == 2) then
	         t0 = max(ti(ix,iy),temin*ev)
                 nucx(ix,iy,igsp) = sqrt(t0/mi(ifld))*
     .                         sigcx*(ni(ix,iy,ifld)+rnn2cx*ng(ix,iy,igsp))
              endif
              nuix(ix,iy,igsp) = fnuizx*nuiz(ix,iy,igsp) +
     .                           fnucxx*nucx(ix,iy,igsp)
                              #dont use neutral-neutral collisions here
c
c   neutral particle source/sink for isupgon=1 (reset below if multispecies
c   models are on [isimpon = 5 or 6 or 7])
              if(isupgon(igsp) .eq. 1)then # inertia gas species is ifld+1
                 psorc(ix,iy,ifld+1)= -psorc(ix,iy,ifld)
                 psorxrc(ix,iy,ifld+1)= -psorxrc(ix,iy,ifld)
                 msor(ix,iy,ifld+1)= 0.
                 msorxr(ix,iy,ifld+1)= 0.
              endif
c
            enddo   #end loop over ix
          enddo     #end loop over iy
         endif      #omit whole loop if zi(ifld) = 0. (neutrals)
        enddo       #end loop over hydrogen species (ifld)

c*****************************************************************
c ... Average psorgc and psorc over cell vol with simple 5pt ave
c*****************************************************************
        igsp = 0
        do ifld = 1, nhsp  # Hydrogen-only loop
         if (zi(ifld) > 0.) then  #calc only for hydr ions, not neuts
           igsp = igsp + 1
           if (ispsorave.eq.0.) then  #use only single-cell value
             do iy = iys1, iyf6
               do ix = ixs1, ixf6
                 psorg(ix,iy,igsp) = psorgc(ix,iy,igsp)
                 psor(ix,iy,ifld) =  psorc(ix,iy,ifld)
                 psorxr(ix,iy,ifld) = psorxrc(ix,iy,ifld)
                 psorrg(ix,iy,igsp) = psorrgc(ix,iy,igsp)
                 if(isupgon(igsp) .eq. 1) then # inertia gas is ifld+1
                   psor(ix,iy,ifld+1)= -psor(ix,iy,ifld)
                   psorxr(ix,iy,ifld+1)= -psorxr(ix,iy,ifld)
                 endif
               enddo
             enddo

           elseif (ispsorave > 0.) # use 5pt ave; first divide by vol

             if (xc < 0) then  #full RHS eval
               j2pwr = j2
               j5pwr = j5
             else  # Jacobian eval
               j2pwr = max(1, yc-1)
               j5pwr = min(ny, yc+1)
             endif
             do iy = j2pwr, j5pwr
               if (xc < 0) then #full RHS eval
                 i2pwr = i2
                 i5pwr = i5
               else  #Jacobian eval
                 i2pwr = max(1,ixm1(xc,yc))
                 i5pwr = min(nx, ixp1(xc,yc))
               endif
               do ix = i2pwr, i5pwr
                 ix1 = ixm1(ix,iy)
                 ix2 = ixp1(ix,iy)
                 psorg(ix,iy,igsp) = (1.-ispsorave*0.5)*
     .                                  psorgc(ix,iy,igsp)+
     .                               0.125*ispsorave*vol(ix,iy)*
     .                          ( psorgc(ix,iy-1,igsp)/vol(ix,iy-1) +
     .                            psorgc(ix,iy+1,igsp)/vol(ix,iy+1) +
     .                            psorgc(ix1,iy,igsp)/vol(ix1,iy)   +
     .                            psorgc(ix2,iy,igsp)/vol(ix2,iy) )
                 psorxr(ix,iy,ifld) = (1.-ispsorave*0.5)*
     .                                 psorxrc(ix,iy,ifld) +
     .                                  0.125*ispsorave*vol(ix,iy)*
     .                           ( psorxrc(ix,iy-1,ifld)/vol(ix,iy-1) +
     .                             psorxrc(ix,iy+1,ifld)/vol(ix,iy+1) +
     .                             psorxrc(ix1,iy,ifld)/vol(ix1,iy)   +
     .                             psorxrc(ix2,iy,ifld)/vol(ix2,iy) )
                 psor(ix,iy,ifld) = -psorg(ix,iy,igsp)
                 psorrg(ix,iy,igsp) = -psorxr(ix,iy,ifld)
                 if(isupgon(igsp) .eq. 1) then # inertia gas is ifld+1
                   psor(ix,iy,ifld+1)= -psor(ix,iy,ifld)
                   psorxr(ix,iy,ifld+1)= -psorxr(ix,iy,ifld)
                 endif

               enddo   #end loop over ix
             enddo     #end loop over iy
           endif       #if-loop on ipsorave
         endif         #omit whole loop if zi(ifld) = 0. (neutrals)
        enddo          #end loop over hydrogen species (ifld)

c ... Can now calc current from nucx since it is updated
      if (cfqyn .gt. 0.) call calc_curr_cx

c ... Ionization and recombination of impurities.
c     (Note that charge-exchange recombination is implemented for
c     impurities but the corresponding terms need to be added to
c     equations for hydrogenic-ion and gas continuity and electron
c     and ion energy.) ???
c     The total source is written as psor+psorxr, where
c     psor = n_(z-1) ne K^i_(z-1) - n_z ne K^i_z     # ionization gain/loss
c     psorxr = -n_z[ne K^r_z+ng K^cx_z]              # cx/r loss from z to z-1
c              +n_(z+1)[ne K^r_(z+1)+ng K^cx_(z+1)]  # cx/r gain to z from z+1

        if (isimpon .ge. 5 .and. nzspt .gt. 0) then
          do 604 iy = iys1, iyf6
             do 603 ix = ixs1, ixf6

                  if (istimingon .eq. 1) tsimp = gettime(sec4)
                  nevol = ne(ix,iy) * vol(ix,iy)
                  ngvol = ng(ix,iy,1) * vol(ix,iy)
                  ngvolcap = ngcap(ix,iy,1) * vol(ix,iy)
                  jg = nhgsp
                  ifld_lcs = nhsp
                  do jz = 1, ngspmx-1           # for all impurity isotopes
                     if (nzsp(jz)==0) break
                     ifld_fcs = ifld_lcs + 1
                     ifld_lcs = ifld_fcs + nzsp(jz) - 1
                     if (ngsp .gt. nhgsp) then  # impurity gas is present
                         jg = jg + 1
                         if (ismctab .eq. 1) then  # rates for Z=0 gas
                            call imprates(te(ix,iy), 0, nzsp(jz), kionz0,
     .                                    krecz, kcxrz)
                         elseif (ismctab .eq. 2) then
                            call mcrates(ne(ix,iy),te(ix,iy),
     .                                   ti(ix,iy)*mp/mi(1),
     .                                   0, nzsp(jz), znucl(ifld_fcs),
     .                                   kionz0, krecz, kcxrz)
                         endif
                         kionz0 = kionz0 + sigvi_floor
			 psorbgg(ix,iy,jg)= ngbackg(jg)*
     .                     (0.9+0.1*(ngbackg(jg)/ng(ix,iy,jg))**ingb) *
     .                                                      nevol*kionz0
                         psorg(ix,iy,jg) = -ngcap(ix,iy,jg)*nevol*kionz0 +
     .                                      psorbgg(ix,iy,jg)
                         psor(ix,iy,ifld_fcs) = - psorg(ix,iy,jg)
                         msor(ix,iy,ifld_fcs)= 0.  # zero gas mom. assumed
                         if (ismctab .eq. 1) then  # rates for Z=1 ions
                            call imprates(te(ix,iy), 1, nzsp(jz), kionz,
     .                                    krecz, kcxrz)
                         elseif (ismctab .eq. 2) then
                            call mcrates(ne(ix,iy), te(ix,iy),
     .                                   ti(ix,iy)*mp/mi(1),
     .                                   1, nzsp(jz), znucl(ifld_fcs),
     .                                   kionz, krecz, kcxrz)
                            kcxrz = cfcximp1*kcxrz  # rescale if desired
                         endif
                         kionz = kionz + sigvi_floor # only to set kionm below
                         kcxrzig = rcxighg(jg)*kcxrz  # K_cx of ng(jg)+ni(1)->
                         niz_floor = nzbackg(ifld_fcs) * (0.9 + 0.1*
     .                          (nzbackg(ifld_fcs)/ni(ix,iy,ifld_fcs))**inzb)
                         pscx0 = ngvol*(nicap(ix,iy,ifld_fcs)-niz_floor)*kcxrz -
     .                           ngcap(ix,iy,jg)*ni(ix,iy,1)*vol(ix,iy)*
     .                                                        kcxrzig
                         psorcxg(ix,iy,jg) = pscx0
                         psorcxg(ix,iy,1) = -pscx0
                         psorrg(ix,iy,jg) = nevol*(nicap(ix,iy,ifld_fcs)-
     .                                                   niz_floor)*krecz
                         psorxr(ix,iy,ifld_fcs)= -psorrg(ix,iy,jg) - pscx0
                         psorxr(ix,iy,1) = psorxr(ix,iy,1) + pscx0
cc                    Note: summed over ion/neutrals here backgrd source=0
                         massfac = cfmassfac*16*mi(1)/(3*(mg(jg)+mi(1)))
                         nuiz(ix,iy,jg) = kionz0*ne(ix,iy)
                         nuix(ix,iy,jg) = fnuizx*nuiz(ix,iy,jg) +
     .                                   kcxrzig*ni(ix,iy,1) +
     .                         massfac*( kelighi(jg)*ni(ix,iy,1) +
     .                                   kelighg(jg)*ng(ix,iy,1) )
                         nucxi(ix,iy,ifld_fcs) = sigcxms(ifld_fcs,jg)*
     .                          sqrt(ti(ix,iy)/mi(ifld_fcs))*ng(ix,iy,jg)
                         nucx(ix,iy,jg) = sigcxms(ifld_fcs,jg)*
     .                        sqrt(ti(ix,iy)/mi(ifld_fcs))*ni(ix,iy,ifld_fcs)
                         massfac = cfmassfac*16*mi(ifld_fcs)/
     .                                               (3*(mg(jg)+mi(ifld_fcs)))
                         nueli(ix,iy,ifld_fcs) = massfac*( keligii(jg)*
     .                                                    ng(ix,iy,jg) )
                         nuelg(ix,iy,jg) = massfac*( keligii(jg)*
     .                                                    ni(ix,iy,ifld_fcs) )

                     else                       # no impurity gas present
                         izch = nint(zi(ifld_fcs))
                         if (ismctab .eq. 1) then
                            call imprates(te(ix,iy), izch, nzsp(jz),
     .                                   kionz, krecz, kcxrz)
                         elseif (ismctab .eq. 2) then
                            call mcrates(ne(ix,iy), te(ix,iy),
     .                                   ti(ix,iy)*mp/mi(1),
     .                                   izch, nzsp(jz), znucl(ifld_fcs),
     .                                   kionz, krecz, kcxrz)
                         endif
                         kionz = kionz + sigvi_floor
                         psor(ix,iy,ifld_fcs) = 0.   # in case nzspt=1
                         psorxr(ix,iy,ifld_fcs) = 0.
                         msor(ix,iy,ifld_fcs) = 0.   # in case nzspt=1
                         msorxr(ix,iy,ifld_fcs) = 0.
                         krecz = 0.
                         kcxrz = 0.
                     endif   # end if-branches for Z=1 with/wo impurity gas
                     if (znucl(ifld_fcs)==1) then
                         # hydrogenic impurity: include cx on ifld_fcs in nuix
                         nuix(ix,iy,jg) = nuix(ix,iy,jg) +
     .                                   kcxrz*ni(ix,iy,ifld_fcs)
                         nuix(ix,iy,1) = nuix(ix,iy,1) +
     .                                   kcxrzig*ni(ix,iy,ifld_fcs)
                     endif

                     kionm = kionz     # set values as previous charge state
                     krecm = krecz
                     kcxrm = kcxrz

                     do ifld = ifld_fcs + 1, ifld_lcs  # for charge states Z > 1
                        izch = nint(zi(ifld))
                        if (ismctab .eq. 1) then
                           call imprates(te(ix,iy), izch, nzsp(jz),
     .                                  kionz, krecz, kcxrz)
                        elseif (ismctab .eq. 2) then
                           call mcrates(ne(ix,iy), te(ix,iy),
     .                                  ti(ix,iy)*mp/mi(1),
     .                                  izch, nzsp(jz), znucl(ifld),
     .                                  kionz, krecz, kcxrz)
                            kcxrz = cfcximp2*kcxrz   #rescale if desired
                        endif
                        kionz = kionz + sigvi_floor
			if (ifld==ifld_lcs) kionz = 0. #ensure no lcs ioniz
                        pxri = 0.    # gets reset if ifld.eq.ifld_fcs+1
                        z1fac = 1.   # gets reset = 0 if ifld.eq.ifld_fcs+1

                        if (ifld .eq. ifld_fcs+1) then #for 2nd charge-state
                           nizm_floor = nzbackg(ifld_fcs) * (0.9 + 0.1*
     .                         (nzbackg(ifld_fcs)/ni(ix,iy,ifld_fcs))**inzb)
                           psor(ix,iy,ifld_fcs) = psor(ix,iy,ifld_fcs)-
     .                            nevol*(nicap(ix,iy,ifld_fcs)-nizm_floor)*
     .                                                            kionm
			   psorbgz(ix,iy) = psorbgz(ix,iy)-nevol*nizm_floor*
     .                                                            kionm
	                   msor(ix,iy,ifld_fcs) = msor(ix,iy,ifld_fcs)-
     .                            nevol*nicap(ix,iy,ifld_fcs)*
     .                            kionm*mi(ifld_fcs)*up(ix,iy,ifld_fcs)
                           pxri = psorxr(ix,iy,ifld_fcs) #set in Z=1 loop
                           z1fac = 0.
                        endif

                        niz_floor = nzbackg(ifld) * (0.9 + 0.1*
     .                            (nzbackg(ifld)/ni(ix,iy,ifld))**inzb)
                        psor(ix,iy,ifld) = nevol *
     .                                     ( nicap(ix,iy,ifld-1)* kionm -
     .                             (nicap(ix,iy,ifld)-niz_floor) * kionz )
			psorbgz(ix,iy) = psorbgz(ix,iy)-nevol*niz_floor*kionz
                        msor(ix,iy,ifld) = nevol *
     .                                  ( nicap(ix,iy,ifld-1)* kionm *
     .                                   mi(ifld-1) * up(ix,iy,ifld-1) -
     .                            (nicap(ix,iy,ifld)) * kionz *
     .                                     mi(ifld) * up(ix,iy,ifld) )
                        psorxr(ix,iy,ifld-1) = pxri-z1fac*(nevol*krecm +
     .                                                     ngvol*kcxrm) *
     .                                    (nicap(ix,iy,ifld-1)-nizm_floor) +
     .                                   (nevol*krecz + ngvol*kcxrz) *
     .                                                   nicap(ix,iy,ifld)
			psorbgz(ix,iy) = psorbgz(ix,iy) + z1fac*
     .                                   (nevol*krecm + ngvol*kcxrm) *
     .                                   nizm_floor
                        msorxr(ix,iy,ifld-1) = 0. - (nevol*krecm +
     .                                               ngvol*kcxrm) *
     .                                    nicap(ix,iy,ifld-1)*
     .                                      mi(ifld-1)*up(ix,iy,ifld-1) +
     .                        (nevol*krecz + ngvol*kcxrz)*nicap(ix,iy,ifld)*
     .                                            mi(ifld)*up(ix,iy,ifld)
                        psorxr(ix,iy,1) = psorxr(ix,iy,1) + ngvolcap*
     .                                               ni(ix,iy,ifld)*kcxrz
                        psorcxg(ix,iy,1) = psorcxg(ix,iy,1) - ngvolcap*
     .                                               ni(ix,iy,ifld)*kcxrz
                        nucxi(ix,iy,ifld) = sigcxms(ifld,jg)*
     .                              sqrt(ti(ix,iy)/mi(ifld))*ng(ix,iy,jg)
                        nucx(ix,iy,jg) = nucx(ix,iy,jg) + sigcxms(ifld,jg)*
     .                          sqrt(ti(ix,iy)/mi(ifld))*ni(ix,iy,ifld)
                        massfac = cfmassfac*16*mi(ifld)/
     .                                              (3*(mg(jg)+mi(ifld)))
                        nueli(ix,iy,ifld) = massfac*( keligii(jg)*
     .                                                   ng(ix,iy,jg) )
                        nuelg(ix,iy,jg) = nuelg(ix,iy,jg) + massfac*
     .                                  ( keligii(jg)*ni(ix,iy,ifld) )

                        kionm = kionz
                        krecm = krecz
                        kcxrm = kcxrz
                        nizm_floor = niz_floor
                        if (ifld .eq. ifld_lcs) then  # last charge-state
                          psorxr(ix,iy,ifld) = -(nevol * krecz +
     .                                           ngvol * kcxrz) *
     .                                          (nicap(ix,iy,ifld)-niz_floor)
			  psorbgz(ix,iy) = psorbgz(ix,iy) + niz_floor *
     .                                      (nevol*krecz + ngvol*kcxrz)
                          msorxr(ix,iy,ifld) = -(nevol * krecz +
     .                                           ngvol * kcxrz) *
     .                                          (nicap(ix,iy,ifld))*
     .                                             mi(ifld)*up(ix,iy,ifld)
                          nuix(ix,iy,jg) = nuix(ix,iy,jg) + nucx(ix,iy,jg) +
     .                                     nuelg(ix,iy,jg)
                        endif


                     enddo   # end do-loop on charge states for isotope jz
                  enddo  # end do-loop on impurity isotopes

                  if (istimingon .eq. 1) call timimpfj (tsimp, xc)
c
c   neutral particle source/sink for isupgon=1; make consistent with impurity
c   contributions just calculated for multispecies
              if (isupgon(1) .eq. 1) then #should be generalized to D & T
                 psor(ix,iy,iigsp)= -psor(ix,iy,1)
                 psorxr(ix,iy,iigsp)= -psorxr(ix,iy,1)
              endif
c
  603        continue
  604     continue
        endif            # end of if (isimpon .ge. 5 .and. nzspt .gt. 0)

c ... Add volume particle loss terms for quasi 1-D radial model
      if (l_parloss .le. 1e9) then
        do iy = iys1, iyf6  # core region has no loss
          do ix = ixs1, ixf6
            do ifld = 1, nfsp
              if (iy .le. iysptrx) then # inside the LCFS when nxpt>1
                                        # (see definition of iysptrx in nphygeo)
                nuvl(ix,iy,ifld) = 0.
              else
                nuvl(ix,iy,ifld) = ( cfvlh*
     .                             sqrt((te(ix,iy)+ti(ix,iy))/mi(1))+
     .                             cfvli(ifld)*
     .              sqrt((zi(ifld)*te(ix,iy)+ti(ix,iy))/mi(ifld)) ) /
     .                                                      l_parloss
              endif
            enddo
          enddo
        enddo
      endif

c ... Set up nuiz & sources for hydrogen molecular gas
      if (ishymol .eq. 1) then
         if (nhgsp .eq. 1) then
          call xerrab('*** nhgsp must exceed 1 for ishymol=1 ***')
        endif
        do iy = iys1, iyf6
         do ix = ixs1, ixf6
           nuiz(ix,iy,2) = ne(ix,iy) * (
     .                          svdiss( te(ix,iy) ) + sigvi_floor )
           massfac = 16*mi(1)/(3*(mg(2)+mi(1)))
           nuix(ix,iy,2)= fnuizx*nuiz(ix,iy,2) +
     .                           massfac*( kelighi(2)*ni(ix,iy,1)+
     .                                     kelighg(2)*ng(ix,iy,1) )
c ...  molecule-molecule collisions would enter viscosity, not nuix
           psorbgg(ix,iy,2) = ngbackg(2)*
     .                     (0.9+0.1*(ngbackg(2)/ng(ix,iy,2))**ingb ) *
     .                                        nuiz(ix,iy,2) * vol(ix,iy)
           psorgc(ix,iy,2) = - ngcap(ix,iy,2)*nuiz(ix,iy,2)*vol(ix,iy)+
     .                        psorbgg(ix,iy,2)
           psorg(ix,iy,2) = psorgc(ix,iy,2)  # no mol sor averaging
           psordis(ix,iy) = -2*psorgc(ix,iy,2)  # 2 atoms per molecule
           if(isupgon(1) .eq. 1) then
             psor(ix,iy,iigsp) = psor(ix,iy,iigsp) + psordis(ix,iy)
           endif
         enddo
        enddo
      endif  # end of loop for ishymol=1 (hydrogen molecules on)


c  *** Now integrate sources over cell volume if ishosor=1 & yl(neq+1)=-1,
c  *** where the last condition means this is only a full RHS eval, not
c  *** a Jacobian calculation

         if (ishosor.eq.1) then  #full RHS eval

           if (svrpkg.eq."cvode") then    # cannot access yl(neq+1)
            call xerrab('*** svrpkg=cvode not allowed for ishosor=1 **')
           endif

          if (yl(neq+1).lt.0) then  #full RHS eval

c ...    integ. sources over cells (but not for Jac) for higher-order accuracy

             do ifld = 1, nfsp  # loop over ions
                call volave(nx, ny, j2, j5, i2, i5, ixp1(0,0), ixm1(0,0),
     .                 fsprd, vol(0,0), psor_tmpov(0,0), psor(0,0,ifld))
                call volave(nx, ny, j2, j5, i2, i5, ixp1(0,0), ixm1(0,0),
     .                 fsprd, vol(0,0), psor_tmpov(0,0), psorxr(0,0,ifld))
             enddo

c *** Now do the gas
             do igsp = 1, ngsp  # now loop over gas species
                call volave(nx, ny, j2, j5, i2, i5, ixp1(0,0), ixm1(0,0),
     .                fsprd, vol(0,0), psor_tmpov(0,0), psorg(0,0,igsp))
                call volave(nx, ny, j2, j5, i2, i5, ixp1(0,0), ixm1(0,0),
     .                fsprd, vol(0,0), psor_tmpov(0,0), psorrg(0,0,igsp))
                call volave(nx, ny, j2, j5, i2, i5, ixp1(0,0), ixm1(0,0),
     .                fsprd, vol(0,0), psor_tmpov(0,0), psorcxg(0,0,igsp))
             enddo

          endif   # end of if (yl(neq+1).lt.0) test
         endif    # end of integrating over sources and ishosor test

      endif              # end of big loop starting if (ifixpsor .eq. 0)

*-----------------------------------------------------------------------
*  -- Calculates the fixed source if it is on
*-----------------------------------------------------------------------

      if (ifixsrc .ne. 0) then
         do 920 iy = j2, j5
            do 910 ix = i2, i5
               snic(ix,iy,1) = snic(ix,iy,1) + vol(ix,iy) * a1n *
     .                          exp(-b1n*(xcs(ix)-xxsrc)**2) *
     .                          exp(-c1n*(yyc(iy)-yysrc)**2)
               seic(ix,iy) = seic(ix,iy) + vol(ix,iy) * a1i * ev *
     .                          exp(-b1i*(xcs(ix)-xxsrc)**2) *
     .                          exp(-c1i*(yyc(iy)-yysrc)**2)
               seec(ix,iy) = seec(ix,iy) + vol(ix,iy) * a1e * ev *
     .                          exp(-b1e*(xcs(ix)-xxsrc)**2) *
     .                          exp(-c1e*(yyc(iy)-yysrc)**2)
 910        continue
 920     continue
      endif

      if (TimingPandfOn.gt.0) TotTimeSource=TotTimeSource+tock(TimeSource)

*****************************************************************
c In the case of neutral parallel mom, call neudif to get
c flux fngy, vy and uu, now that we have evaluated nuix etc.
*****************************************************************
      if (TimingPandfOn.gt.0) TimeNeudif=tick()
ccc      if (isupgon .eq. 1 .and. zi(ifld) .eq. 0.0) call neudif
      if (ineudif .eq. 1) then
         call neudif
      elseif (ineudif .eq. 2) then
c ..Timing
      if(istimingon==1) tsnpg=gettime(sec4)
         call neudifpg
c ..Timing
      if(istimingon==1) ttnpg=ttnpg+(gettime(sec4)-tsnpg)
      elseif (ineudif .eq. 3) then
         call neudifl
      else
         call neudifo
      endif
      if (TimingPandfOn.gt.0) TotTimeNeudif=TotTimeNeudif+tock(TimeNeudif)
*****************************************************************
*  Other volume sources calculated in old SRCMOD
*****************************************************************
*  ---------------------------------------------------------------------
*  electron-ion transfer terms and an
*  approximation to the ion-ion thermal force.
*  cfw: For the neutral momentum eq there is also a v_gas*grad(p_gas)
*       term which is evaluated using ng, ti and gpiy
*  ---------------------------------------------------------------------

c...  Force fluxes and gradients on cuts to be zero for half-space problems
      if (isfixlb(1).eq.2.or. isfixrb(1).eq.2) then
         if (isfixlb(1).eq.2) then
            ix = ixpt2(1)
         else
            ix = ixpt1(1)
         endif
         if (ix.ge.i2 .and. ix.le.i5+1 .and. iysptrx1(1) > 0) then
            do iy = 0, iysptrx1(1)
               gpex(ix,iy) = 0.
               frice(ix,iy) = 0.
               ex(ix,iy) = 0.
               upe(ix,iy) = 0.
               do ifld = 1, nfsp
                  gpix(ix,iy,ifld) = 0.
                  frici(ix,iy,ifld) = 0.
                  uu(ix,iy,ifld) = 0.
                  upi(ix,iy,ifld) = 0.
               enddo
            enddo
         endif
      endif

*  -- Set up electron parallel contribution to seec & smoc
      do iy = j2, j5
         do ix = i2, i5
            ix1 = ixm1(ix,iy)
            ix2 = ixp1(ix,iy)
            t1 =.5*cvgp*(upe(ix,iy)*rrv(ix,iy)*
     .              ave(gx(ix,iy),gx(ix2,iy))*gpex(ix,iy)/gxf(ix,iy) +
     .                   upe(ix1,iy)*rrv(ix1,iy)*
     .              ave(gx(ix,iy),gx(ix1,iy))*gpex(ix1,iy)/gxf(ix1,iy) )
            t2 = 1.e-20* 0.25*(fqp(ix,iy)+fqp(ix1,iy))*
     .                (ex(ix,iy)+ex(ix1,iy))/gx(ix,iy)
            seec(ix,iy) = seec(ix,iy) + t1*vol(ix,iy) - t2
            if (nusp-isupgon(1).eq.1) smoc(ix,iy,1) = -cpgx*gpex(ix,iy)*
     .                                   rrv(ix,iy)*sx(ix,iy)/gxf(ix,iy)
         enddo
      enddo

*  -- Now loop over all ion species for seec, seic, and smoc --

      do 101 ifld = 1, nusp  #not nfsp; up only for ifld<=nusp
*     -- coupling in the x-direction --

         do 31 iy = j2, j5
            do 30 ix = i2, i5
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               tv =  gpix(ix ,iy,ifld)/gxf(ix,iy)
               t1 = gpix(ix1,iy,ifld)/gxf(ix1,iy)
               t1 = .5*cvgp*( up(ix,iy,ifld)*rrv(ix,iy)*
     .                                 ave(gx(ix2,iy),gx(ix,iy))*tv
     .                      + up(ix1,iy,ifld)*rrv(ix1,iy)*
     .                                 ave(gx(ix,iy),gx(ix1,iy))*t1 )
               seic(ix,iy) = seic(ix,iy) + cfvgpx(ifld)*t1*vol(ix,iy)
               if (zi(ifld) .ne. 0) then
                 t0 = - cpiup(ifld)*gpix(ix,iy,ifld)*rrv(ix,iy)*
     .                                              sx(ix,iy)/gxf(ix,iy)
                 if (nusp-isupgon(1) .eq. 1) then  # single ion mom. eq.
                    smoc(ix,iy,1) = smoc(ix,iy,1) + cpgx*t0
                 else                # multiple mom. eq., so nusp=nisp
                    t0 = t0 +( qe*zi(ifld)*0.5*( ni(ix2,iy,ifld)+
     .                         ni(ix,iy,ifld) )*ex(ix,iy)*rrv(ix,iy) +
     .                         frici(ix,iy,ifld) )* sx(ix,iy)/gxf(ix,iy)
                    if (ifld <= nusp) smoc(ix,iy,ifld) =
     .                                        smoc(ix,iy,ifld) + cpgx*t0
                 endif
c...  Add friction part of Q_e here
                 tv = 0.25*(frice(ix,iy)+frice(ix1,iy))*
     .                ( upe(ix,iy)     + upe(ix1,iy) -
     .                upi(ix,iy,ifld) - upi(ix1,iy,ifld) )
                 seec(ix,iy) = seec(ix,iy) - zi(ifld)**2*ni(ix,iy,ifld)*
     .                                          tv*vol(ix,iy)/nz2(ix,iy)
c... REMEMBER TO ADD CONTRIBUTION TO SEEC FROM V2 1/26/95
               endif
   30       continue
   31    continue

*     -- coupling in the x & y-directions --
         do 34 iy = j2, j5
           do 33 ix = i2, i5
             if (isgpye == 0) then
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               vyiy0 = fracvgpgp*vygp(ix,iy,ifld)
     .                    +(1.-fracvgpgp)*vycb(ix,iy,ifld)
               vyiym1 = fracvgpgp*vygp(ix,iy-1,ifld)
     .                    +(1.-fracvgpgp)*vycb(ix,iy-1,ifld)
               v2ix0 = fracvgpgp*v2xgp(ix,iy,ifld)
     .                    +(1.-fracvgpgp)*v2cb(ix,iy,ifld)
               v2ixm1 = fracvgpgp*v2xgp(ix1,iy,ifld)
     .                    +(1.-fracvgpgp)*v2cb(ix1,iy,ifld)
               t1 =.5*cvgp*( vygp(ix,iy  ,ifld)*gpiy(ix,iy  ,ifld) +
     .                       vygp(ix,iy-1,ifld)*gpiy(ix,iy-1,ifld) +
     .                    v2xgp(ix ,iy,ifld)*ave(gx(ix,iy),gx(ix2,iy))*
     .                               gpix(ix ,iy,ifld)/gxf(ix ,iy) +
     .                    v2xgp(ix1,iy,ifld)*ave(gx(ix,iy),gx(ix1,iy))*
     .                               gpix(ix1,iy,ifld)/gxf(ix1,iy) )
               t2 = t1
             elseif (isgpye == 1) then    # Old B2 model with Jperp=0
               t1 = -0.5*( vy(ix,iy  ,ifld)*gpey(ix,iy  ) +
     .                     vy(ix,iy-1,ifld)*gpey(ix,iy-1) )
               t2 = t1
             elseif (isgpye == 2) then    # Knoll expression
               t1 = -0.5*( vy(ix,iy  ,ifld)*gpey(ix,iy  ) +
     .                     vy(ix,iy-1,ifld)*gpey(ix,iy-1) )
             endif
             if (zi(ifld).gt.1e-10) then
               seec(ix,iy) = seec(ix,iy) - fluxfacy*t1 * vol(ix,iy)
               seic(ix,iy) = seic(ix,iy) + fluxfacy*cfvgpy(ifld)*t2*
     .                                                     vol(ix,iy)
             endif
   33      continue
   34    continue

  101 continue

*****************************************************************
*  Other physics coefficients. (old PHYVIS)
*****************************************************************

* -- loop over species number --

      do 102 ifld = 1, nfsp

c
c     neutral viscosity for isupgon=1
c
         if(isupgon(1) .eq. 1 .and. zi(ifld) .eq. 0)then
            do 936 iy = j1,j6
               iyp1 = min(iy+1,ny+1)
               do 937 ix = i1,i6
c
                  ix1 = ixm1(ix,iy)
                  vtn = sqrt(max(tg(ix,iy,1),temin*ev)/mi(ifld))
 		  qfl = flalfvgxa(ix)*nm(ix,iy,ifld)*vtn**2
                  if(isvisxn_old == 1) then
                    lmfpn = 1./(sigcx *
     .                          (ni(ix,iy,1) + rnn2cx*ni(ix,iy,ifld)))
                  elseif(isvisxn_old==0 .and. ishymol==0) then
                    lmfppar = vtn/(kelhihg*ni(ix,iy,1) +
     .                                         kelhghg*ni(ix,iy,ifld))
                    lmfpperp = vtn/( vtn*sigcx*ni(ix,iy,1) +
     .                      kelhihg*ni(ix,iy,1)+kelhghg*ni(ix,iy,ifld) )
                    rrfac = rr(ix,iy)*rr(ix,iy)
                    lmfpn = lmfppar*rrfac + lmfpperp*(1-rrfac)
                  else   # (isvisxn_old=0 .and. ishymol=1) then #with mols
                    lmfppar = vtn/(kelhihg*ni(ix,iy,1) +
     .                     kelhghg*ni(ix,iy,ifld) + kelhmhg*ng(ix,iy,2))
                    lmfpperp = vtn/( vtn*sigcx*ni(ix,iy,1) +
     .                     kelhihg*ni(ix,iy,1) +kelhghg*ni(ix,iy,ifld) +
     .                     kelhmhg*ng(ix,iy,2) )
                    rrfac = rr(ix,iy)*rr(ix,iy)
                    lmfpn = lmfppar*rrfac + lmfpperp*(1-rrfac)
                  endif
                  csh = lmfpn*nm(ix,iy,ifld)*
     .                            sqrt(max(ti(ix,iy),temin*ev)/mi(ifld))*
     .                                      lgvmax/(lgvmax + lmfpn)
                  qsh = csh * (up(ix1,iy,ifld)-up(ix,iy,ifld)) *
     .                                                       gx(ix,iy)
                  visx(ix,iy,ifld)= cfvisxn*csh/
     .               (1 + (abs(qsh/(qfl+cutlo))**flgamvg))**(1./flgamvg)
     .               + cfanomvisxg*travis(ifld)*nm(ix,iy,ifld)

c    Now do y-direction; csh is the same for neutrals, but should use den face
 		  qfl = flalfvgya(iy)*nm(ix,iy,ifld)*vtn**2
                  qsh = csh * (up(ix,iy,ifld)-up(ix,iyp1,ifld)) *
     .                                                        gyf(ix,iy)
                  visy(ix,iy,ifld)= cfvisyn*csh /
     .               (1 + (abs(qsh/(qfl+cutlo))**flgamvg))**(1./flgamvg)
     .               + cfanomvisyg*travis(ifld)*nm(ix,iy,ifld)
c
 937           continue
 936        continue
         endif
c
c
       if(zi(ifld) > 1.e-20) then
         do 39 iy = j1, j6
            do 38 ix = i1, i6
               w(ix,iy) = 0.0e0
 38         continue
 39      continue

         do 42 jfld = 1, nisp
            tv = zi(jfld)**2 / sqrt((mi(ifld)+mi(jfld))/(2*mp))
            do 41 iy = j1, j6
               do 40 ix = i1, i6
                  w(ix,iy) = w(ix,iy) + tv*ni(ix,iy,jfld)
   40          continue
   41       continue
   42    continue

         do 44 iy = j1, j6
            do 43 ix = i1, i6
	       ctaui(ix,iy,ifld) = 2.1e13/(loglambda(ix,iy)*zi(ifld)**2) # mass fac?
               tv2 = ctaui(ix,iy,ifld)/(ev*sqrt(ev))
               if (convis .eq. 0) then
                  a = max (ti(ix,iy), temin*ev)
               else
                  a = afix*ev
               endif
	       epstmp = max(epsneo(ix,iy), 1.e-50)
               visxtmp = tv2 * coef * rr(ix,iy) * rr(ix,iy) *
     .                            a*a*sqrt(a) * ni(ix,iy,ifld)/w(ix,iy)
               visx(ix,iy,ifld) = parvis(ifld)*visxtmp +
     .                            trax_use(ix,iy,ifld)*nm(ix,iy,ifld)
               nuii(ix,iy,ifld) = w(ix,iy)/(tv2*a*sqrt(a))
               nuiistar(ix,iy,ifld) = ( lconneo(ix,iy)*nuii(ix,iy,ifld)/
     .                                 (epstmp**1.5*(2*a/mi(ifld))**0.5)
     .                                           + 1.e-50 )
               visxneo(ix,iy,ifld) = visxtmp*
     .                 (1./(1.+epstmp**(-1.5)/nuiistar(ix,iy,ifld)))*
     .                               (1./(1.+1./nuiistar(ix,iy,ifld)))
               rt2nus = 1.414*nuiistar(ix,iy,ifld)
               ktneo(ix,iy,ifld) = (-0.17 + 1.05*rt2nus**.5 +
     .                     2.7*rt2nus**2*epstmp**3) / ( 1.+
     .                0.7*rt2nus**.5 + rt2nus**2*epstmp**3 )
               alfneo(ix,iy,ifld) = (8./15.)*(ktneo(ix,iy,ifld) - 1.)*
     .                 (1./(1.+epstmp**(-1.5)/nuiistar(ix,iy,ifld)))*
     .                            ( 1./(1.+1./nuiistar(ix,iy,ifld)) )
               k2neo(ix,iy,ifld) =(.66 + 1.88*epstmp**.5 - 1.54*epstmp)/
     .                            (1. + 1.03*rt2nus**.5 + 0.31*rt2nus) +
     .             1.17*epstmp**3*rt2nus/(1. + 0.74*epstmp**1.5*rt2nus)
c...  flux limit the viscosity; beware of using visx(0,iy) and
c...  visx(nx+1,iy) as they are meaningless when flux limited
               ix1 = ixm1(ix,iy)
               t0 = max (ti(ix,iy), temin*ev)
               vtn = sqrt(t0/mi(ifld))
               mfl = flalfv * nm(ix,iy,ifld) * rr(ix,iy) *
     .               vol(ix,iy) * gx(ix,iy) * (t0/mi(ifld))
ccc  Distance between veloc. cell centers:
               if (isgxvon .eq. 0) then     # dx(ix)=1/gx(ix)
                 csh = visx(ix,iy,ifld) * vol(ix,iy) * gx(ix,iy)
     .                            * gx(ix,iy)
               elseif (isgxvon .eq. 1) then # dx(ix)=.5/gxf(ix-1)+.5/gxf(ix)
                 csh = visx(ix,iy,ifld) * vol(ix,iy) * gx(ix,iy)
     .               * 2*gxf(ix,iy)*gxf(ix1,iy)/(gxf(ix,iy)+gxf(ix1,iy))
               endif
ccc
               msh = abs( csh*(upi(ix1,iy,ifld) - upi(ix,iy,ifld)) )
               visx(ix,iy,ifld) = visx(ix,iy,ifld)
     .               / (1 + (msh/(mfl+1.e-20*msh))**flgamv )**(1/flgamv)
               visy(ix,iy,ifld)=(fcdif*travis(ifld)+ tray_use(ix,iy,ifld))*
     .                              nm(ix,iy,ifld) +  4*eta1(ix,iy)
   43       continue
   44    continue
       endif      # test if zi(ifld) > 1.e-20
  102 continue    # large loop for ifld = 1, nfsp

*****************************************************************
*****************************************************************
*  Heat Conduction. (old PHYTHC)
*****************************************************************
*  ---------------------------------------------------------------------
*  compute conductivities on cell faces
*  ---------------------------------------------------------------------

*  -- initialize to 0 --

      do 706 iy = j1, j6
         do 705 ix = i1, i6
            hcxe(ix,iy) = 0.0e0
            hcxi(ix,iy) = 0.0e0
            hcxineo(ix,iy) = 0.0e0
            hcye(ix,iy) = 0.0e0
            hcyi(ix,iy) = 0.0e0
            do ifld = 1, nisp
               hcxij(ix,iy,ifld) = 0.0e0
               hcyij(ix,iy,ifld) = 0.0e0
            enddo
  705    continue
  706 continue

*  -- loop over species number --

      do 103 ifld = 1, nisp
c -- Skip this if these are the neutrals (zi(ifld).eq.0)
         if (zi(ifld) .eq. 0.0e0) goto 103

c...  Initialize w1 and w2 for each species
         do 49 iy = j1, j6
            do 48 ix = i1, i6
               w1(ix,iy) = 0.0e0
               w2(ix,iy) = 0.0e0
 48         continue
 49      continue

*     -- conductivities --
*        The poloidal conductivities  are initially computed without
*        the factor rr**2 * tv**2.5)

         do 52 jfld = 1, nisp
            tv = zi(jfld)**2
            a = zi(jfld)**2 *
     .            sqrt(2*mi(ifld)*mi(jfld)/(mi(ifld)+mi(jfld)))
            do 51 iy = j1, j6
               do 50 ix = i1, i6
                  ix1 = ixp1(ix,iy)
                  w1(ix,iy) = w1(ix,iy) + tv*(ni(ix,iy,jfld)*gx(ix,iy) +
     .                                     ni(ix1,iy,jfld)*gx(ix1,iy)) /
     .                                        (gx(ix,iy) + gx(ix1,iy))
                  w2(ix,iy) = w2(ix,iy) + a*(ni(ix,iy,jfld)*gx(ix,iy) +
     .                                     ni(ix1,iy,jfld)*gx(ix1,iy)) /
     .                                        (gx(ix,iy) + gx(ix1,iy))
   50          continue
   51       continue
   52    continue

         do 59 iy = j1, j6
            do 58 ix = i1, i6
               ix1 = ixp1(ix,iy)
               iyp1 = min(ny+1, iy+1)
               ctaue(ix,iy,ifld) = 3.5e11*zi(ifld)/loglambda(ix,iy)
               ctaui(ix,iy,ifld) =2.1e13/(loglambda(ix,iy)*zi(ifld)**2)
               fxe = kxe * ce * ctaue(ix,iy,ifld) / (me*ev*sqrt(ev))
               fxi = kxi * ci * ctaui(ix,iy,ifld) / (ev*sqrt(ev*mp))
               fxet = fxe
               fxit = fxi
               do jx = 1, nxpt  #reduce kxe inside sep by rkxecore fac
                  if ( (iy.le.iysptrx) .and.
     .                    ix.gt.ixpt1(jx) .and. ix.le.ixpt2(jx) ) then
                     fxet = fxe/( 1. + (rkxecore-1.)*
     .                              (yyf(iy)/(yyf(0)+4.e-50))**inkxc )
                     fxit = kxicore * fxi
                  endif
               enddo
               niavex = ( ni(ix ,iy,ifld)*gx(ix ,iy) +
     .                    ni(ix1,iy,ifld)*gx(ix1,iy)) /
     .                     (gx(ix,iy) + gx(ix1,iy))
               niavey = ( niy0(ix,iy,ifld)*gy(ix,iy) +
     .                    niy1(ix,iy,ifld)*gy(ix,iyp1)) /
     .                     (gy(ix,iy) + gy(ix,iyp1))
               hcxe(ix,iy) = hcxe(ix,iy)+fxet*niavex/w1(ix,iy)

c ... Use fixed diffusivity inside the separatrix, anomalous outside,
c     if anomalous-diffusivity multiplier is nonzero.
               kyemix = fcdif*kye + kye_use(ix,iy)
cccMER NOTE: when there are multiple x-points, as in 'dnull' configuration,
ccc          iysptrx is the last closed flux surface (see S.R. nphygeo)
               if(kyet .gt. 1.e-20 .and. iy .gt. iysptrx) then
                  kyemix = (1. - ckyet) * kyemix +
     .               ckyet * kyet * diffusivwrk(ix,iy)
               endif
               hcye(ix,iy) = hcye(ix,iy) + ( kyemix +
     .                       2.33*(dclass_e(ix,iy)+dclass_e(ix,iyp1)) )*
     .                                               zi(ifld) * niavey

               hcxij(ix,iy,ifld) = fxit*niavex/w2(ix,iy)
               kyimix = fcdif*kyi + kyi_use(ix,iy)
cccMER NOTE: when there are multiple x-points, as in 'dnull' configuration,
ccc          iysptrx is the last closed flux surface (see S.R. nphygeo)
               if(kyit .gt. 1.e-20 .and. iy .gt. iysptrx) then
                  kyimix = (1. - ckyit) * kyimix +
     .               ckyit * kyit * diffusivwrk(ix,iy)
               endif
               hcyij(ix,iy,ifld) = hcyij(ix,iy,ifld) + ( kyimix +
     .                           (dclass_i(ix,iy)+dclass_i(ix,iyp1)) )*
     .                                                         niavey
   58       continue
   59    continue

  103 continue

c ... Add ion temp. dep. for pol. terms, flux limit, & build total ion hcx,yi
      do 595 ifld = 1, nisp
         if (zi(ifld) .eq. 0.e0) goto 595
         do iy = j1, j6
            do ix = i1, i6
               ix1 = ixp1(ix,iy)
               if (concap .eq. 0) then
                  tiave = (ti(ix,iy)*gx(ix,iy) + ti(ix1,iy)*gx(ix1,iy)) /
     .                                         (gx(ix,iy) + gx(ix1,iy))
                  do jx = 1, nxpt
                    if (ix==ixlb(jx).and.ixmnbcl==1) tiave=ti(ixlb(jx)+1,iy)
                    if (ix==ixrb(jx).and.ixmxbcl==1) tiave=ti(ixrb(jx),iy)
                  enddo
                  a = max (tiave, temin*ev)
               else
                  a = afix*ev
               endif
               hcxij(ix,iy,ifld) = hcxij(ix,iy,ifld)*rrv(ix,iy)*
     .                                            rrv(ix,iy)*a*a*sqrt(a)
c...  reduce hcxij if ti very flat; prevents large conduction for high ti
c...  or if lmfpi exceeds a mean-free path limit, lmfplim
               lmfpi = 1.e16*(tiave/ev)**2/ni(ix,iy,1) # i-mean-free-path
               niavex = ( ni(ix ,iy,ifld)*gx(ix ,iy) +
     .                    ni(ix1,iy,ifld)*gx(ix1,iy)) /
     .                                     (gx(ix,iy) + gx(ix1,iy))
               hcxij(ix,iy,ifld) = hcxij(ix,iy,ifld)/(1.+lmfpi/lmfplim)
               hcxij(ix,iy,ifld) = hcxij(ix,iy,ifld) *
     .                         ( cutlo + (ti(ix,iy)-ti(ix1,iy))**2 ) /
     .                         ( cutlo + (ti(ix,iy)-ti(ix1,iy))**2 +
     .                      (0.5*alfkxi*(ti(ix,iy)+ti(ix1,iy)))**2 ) +
     .                         kxi_use(ix,iy)*niavex

c ... Flux limit individ. hcxij in poloidal direction if isflxldi=2
               if (isflxldi .eq. 2) then
                  niavex = ( ni(ix ,iy,ifld)*gx(ix ,iy) +
     .                       ni(ix1,iy,ifld)*gx(ix1,iy)) /
     .                                        (gx(ix,iy) + gx(ix1,iy))
                  wallfac = 1.
                  do jx = 1, nxpt
                     if ( ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .                      (ix==ixrb(jx).and.ixmxbcl==1) )
     .                    .and. (isplflxl==0) ) wallfac = flalfipl/flalfi
                  enddo
                  qflx = wallfac*flalfi * rrv(ix,iy) *
     .                                    sqrt(a/mi(ifld)) * niavex * a
                  cshx = hcxij(ix,iy,ifld)
	          lxtic = 0.5*(ti(ix,iy)+ti(ix1,iy)) /
     .               (abs(ti(ix,iy)-ti(ix1,iy))*gxf(ix,iy) + 100.*cutlo)
	          qshx = cshx * (ti(ix,iy)-ti(ix1,iy)) * gxf(ix,iy) *
     .                                              (1. + lxtic/lxtimax)
                  hcxij(ix,iy,ifld) = cshx  / (1 + abs(qshx/qflx))
               endif
               hcxi(ix,iy) = hcxi(ix,iy) + hcxij(ix,iy,ifld)
               hcxineo(ix,iy) = hcxineo(ix,iy) + hcxij(ix,iy,ifld)*
     .                      1.5676*epsneo(ix,iy)**1.5/k2neo(ix,iy,ifld)
               hcyi(ix,iy) = hcyi(ix,iy) + hcyij(ix,iy,ifld)
               qipar(ix,iy,ifld) = hcxij(ix,iy,ifld)*gxf(ix,iy)*
     .                                (ti(ix,iy)-ti(ix1,iy))/rrv(ix,iy)
            enddo
         enddo
 595  continue

c...  Now include elec. temp and other dep. in poloidal terms + diff. neut.
      do 61 iy = j1, j6
         do 60 ix = i1, i6
            ix1 = ixp1(ix,iy)
            iyp1 = min(ny+1, iy+1)
            if (concap .eq. 0) then
               teave = (te(ix,iy)*gx(ix,iy) + te(ix1,iy)*gx(ix1,iy)) /
     .                                       (gx(ix,iy) + gx(ix1,iy))
               do jx = 1, nxpt
                 if(ix==ixlb(jx).and.ixmnbcl==1) teave=te(ixlb(jx)+1,iy)
                 if(ix==ixrb(jx).and.ixmxbcl==1) teave=te(ixrb(jx),iy)
               enddo
               a = max (teave, temin*ev)
            else
               a = afix*ev
            endif
            zeffave = (zeff(ix,iy)*gx(ix,iy) + zeff(ix1,iy)*gx(ix1,iy)) /
     .                                         (gx(ix,iy) + gx(ix1,iy))
            zcoef = 0.308 + 0.767*zeffave - 0.075*zeffave**2
            hcxe(ix,iy) = hcxe(ix,iy)*rrv(ix,iy)*rrv(ix,iy)*a*a*sqrt(a)
     .                    *zcoef

c...  reduce hcxe if te very flat; prevents very large conduction for high te
c...  or if the mean-free path exceeds lmfplim
            lmfpe = 2e16*(te(ix,iy)/ev)**2/ne(ix,iy)  #mfp for elec [m]
            neavex = ( ne(ix ,iy)*gx(ix ,iy) +
     .                 ne(ix1,iy)*gx(ix1,iy))/(gx(ix,iy) + gx(ix1,iy))
            hcxe(ix,iy) = hcxe(ix,iy) *
     .                         ( cutlo + (te(ix,iy)-te(ix1,iy))**2 ) /
     .                         ( cutlo + (te(ix,iy)-te(ix1,iy))**2 +
     .                     (0.5*alfkxe*(te(ix,iy)+te(ix1,iy)))**2 ) +
     .                        kxe_use(ix,iy)*neavex
            hcxe(ix,iy) = hcxe(ix,iy) /( (1. + lmfpe/lmfplim) *
     .                  (1+hcxe(ix,iy)*gx(ix,iy)**2*tdiflim/ne(ix,iy)) )
            if (isupgon(1).eq.0) then   # add diff. gas cx contrib. to hcxi
               hcxn(ix,iy) = 0.
               hcyn(ix,iy) = 0.
c..1dn0802
c IJ 2016/10/10	add cfneutsor_ei multiplier to control fraction of neutral energy to add
               hcxi(ix,iy) = hcxi(ix,iy)
     .                + cfneut*cfneutsor_ei*kxn*( ng(ix ,iy,1)*ti(ix ,iy)
     .                              +ng(ix1,iy,1)*ti(ix1,iy) ) /
     .                  (mi(1)*(nucx(ix,iy,1) + nucx(ix1,iy,1)))
               hcyi(ix,iy) = hcyi(ix,iy)
     .                + cfneut*cfneutsor_ei*kyn*( ngy0(ix,iy,1)*tiy0(ix,iy)
     .                              +ngy1(ix,iy,1)*tiy1(ix,iy) ) /
     .                  (mi(1)*(nucx(ix,iy,1) + nucx(ix,iyp1,1)))
            endif

   60    continue
 61   continue
c
c
      if (isupgon(1).eq.1) then
c
c ----- Section for the inertial neutral fluid; we need to do different
c ----- things than for the ions. Note third index=iigsp is neutral species
c ----- The inertial neutrals coeff. are flux-limited and add to total here
         do 62 iy = j1, j6
            iy1 = min(iy,ny)   #dont use j5 because hcx also in loop (not imp.)
            do 63 ix = i1, i6
               ix1 = ixp1(ix,iy)
               tgavex = max(0.5*(tg(ix,iy,1) + tg(ix1,iy,1)), temin*ev)
               tgavey= max(0.5*(tgy0(ix,iy,1)+tgy1(ix,iy,1)), temin*ev)
               niavex = 0.5*(ni(ix,iy,1) + ni(ix1,iy,1)) #only for coll. term
               niavey = 0.5*(niy0(ix,iy1,1) + niy1(ix,iy1,1)) #only coll. term
               noavex = ( ni(ix ,iy,iigsp)*gx(ix ,iy) +
     .                    ni(ix1,iy,iigsp)*gx(ix1,iy)) /
     .                     (gx(ix,iy) + gx(ix1,iy))
               noavey = 0.5*(niy0(ix,iy1,iigsp) + niy1(ix,iy1,iigsp))

c          Set up flux-limit variables (no rrv here)
c          First limit the poloidal coeff, then radial
c IJ 2016/10/10	add cfneutsor_ei multiplier to control fraction of neutral energy to add
               qflx = flalftgxa(ix) * sqrt(tgavex/mi(iigsp)) * noavex *
     .                                                     tgavex
               lmfpn = 1./(sigcx * (niavex + rnn2cx*noavex))
               cshx = lmfpn*sqrt(tgavex/mi(iigsp))*noavex *
     .                         lgtmax(iigsp)/(lmfpn + lgtmax(iigsp))
               qshx = cshx * (tg(ix,iy,1)-tg(ix1,iy,1)) * gxf(ix,iy)
	       hcxn(ix,iy) = cshx  /
     .                      (1 + (abs(qshx/qflx))**flgamtg)**(1./flgamtg)
               hcxi(ix,iy) = hcxi(ix,iy) + cfneut*cfneutsor_ei*hcxn(ix,iy)
c          Now for the radial flux limit - good for nonorthog grid too
               qfly = flalftgya(iy) * sqrt(tgavey/mi(iigsp)) * noavey *
     .                                                     tgavey
               lmfpn = 1./(sigcx * (niavey + rnn2cx*noavey))
               cshy = lmfpn*sqrt(tgavey/mi(iigsp))*noavey *
     .                         lgtmax(iigsp)/(lmfpn + lgtmax(iigsp))
               qshy = cshy * (tgy0(ix,iy1,1)-tgy1(ix,iy1,1)) * gyf(ix,iy)
               hcyn(ix,iy) = cshy /
     .                      (1 + (abs(qshy/qfly))**flgamtg)**(1./flgamtg)
               hcyi(ix,iy) = hcyi(ix,iy) + cfneut*cfneutsor_ei*hcyn(ix,iy)
c
  63        continue
  62     continue
      endif

*  Equipartition (old PHYEQP)
*****************************************************************
*  ---------------------------------------------------------------------
*  compute equipartition.
*  ---------------------------------------------------------------------

*     -- initialize w3 --
      do 69 iy = j1, j6
         do 68 ix = i1, i6
            w3(ix,iy) = 0.0e0
 68      continue
 69   continue

      do 72 ifld = 1, nisp
         tv = zi(ifld)**2/mi(ifld)
         do 71 iy = j2, j5
            do 70 ix = i2, i5
               w3(ix,iy) = w3(ix,iy) + tv*ni(ix,iy,ifld)
   70       continue
   71    continue
   72 continue

*  -- compute equipartition --
ccc In detail, coef1 = qe**4*sqrt(me)*lnlam / ((2*pi)**1.5*eps0**2)
      do 74 iy = j2, j5
         do 73 ix = i2, i5
            a = max (te(ix,iy), temin*ev)
            coef1 = feqp*4.8e-15*loglambda(ix,iy)*sqrt(ev)*ev*mp
            eqp(ix,iy) = coef1 * w3(ix,iy) * ne(ix,iy) / (a*sqrt(a))
c...       reduce eqp when (te-ti)/(te+ti) << 1
            eqp(ix,iy) = eqp(ix,iy) * (a-ti(ix,iy))**2 / ( cutlo +
     .                (a-ti(ix,iy))**2 + (alfeqp*(a+ti(ix,iy)))**2 )
   73    continue
   74 continue

*********************************************************
c ... Gas thermal coefficients, initially for molecules *
*********************************************************
*
c ... Gas thermal conductivity coeffs - from self-collisions
*****************************************************************

      if (nisp >= 2) then  # uses ni(,,2), so must have atoms
       do igsp = 1, ngsp
        do iy = j1, j6
        iy1 = min(iy,ny)
          do ix = i1, i6
            ix1 = ixp1(ix,iy)
            tgavex = max( (tg(ix,iy,igsp)*gx(ix,iy) +
     .                              tg(ix1,iy,igsp)*gx(ix1,iy)) /
     .                             (gx(ix,iy) + gx(ix1,iy)), temin*ev )
            tgavey=max(0.5*(tgy0(ix,iy,igsp)+tgy1(ix,iy,igsp)),temin*ev)
            noavex = ( ng(ix,iy,igsp)*gx(ix,iy) +
     .                                   ng(ix1,iy,igsp)*gx(ix1,iy)) /
     .                                   (gx(ix,iy) + gx(ix1,iy))
            niavex = ( ni(ix,iy,1)*gx(ix,iy) +
     .                                   ni(ix1,iy,1)*gx(ix1,iy)) /
     .                                   (gx(ix,iy) + gx(ix1,iy))
            naavex = ( ni(ix,iy,2)*gx(ix,iy) +
     .                                   ni(ix1,iy,2)*gx(ix1,iy)) /
     .                                   (gx(ix,iy) + gx(ix1,iy))
            noavey = 0.5*(ngy0(ix,iy1,igsp) + ngy1(ix,iy1,igsp))
            niavey = 0.5*(niy0(ix,iy1,1) + niy1(ix,iy1,1))
            naavey = 0.5*(niy0(ix,iy1,2) + niy1(ix,iy1,2))
            nuelmolx = noavex*kelhmhm + niavex*kelhmhg +
     .                 naavex*kelhmhg
            qflx = flalftmx*sqrt(tgavex/mg(igsp))*noavex*tgavex
            cshx = cftgcond*noavex*tgavex/(mg(igsp)*nuelmolx)  #assume K not fcn Tg
            qshx = cshx * (tg(ix,iy,igsp)-tg(ix1,iy,igsp)) * gxf(ix,iy)
            hcxg(ix,iy,igsp) = cshx /
     .                     (1.+ (abs(qshx/qflx))**flgamtg)**(1./flgamtg)
            hcxg(ix,iy,igsp)=(1.-cfhcxgc(igsp))*hcxg(ix,iy,igsp)+
     .                     cfhcxgc(igsp)*noavex*kxg_use(ix,iy,igsp)
c..   Now radial direction
            nuelmoly = noavey*kelhmhm + niavey*kelhmhg +
     .                 naavey*kelhmhg
            qfly = flalftmy*sqrt(tgavey/mg(igsp))*noavey*tgavey
            cshy = cftgcond*noavey*tgavey/(mg(igsp)*nuelmoly)  #assume Kel_s not fcn Tg
            qshy = cshy*(tgy0(ix,iy1,igsp)-tgy1(ix,iy1,igsp))*gyf(ix,iy)
           hcyg(ix,iy,igsp) = cshy /
     .                     (1 + (abs(qshy/qfly))**flgamtg)**(1./flgamtg)
            hcyg(ix,iy,igsp)=(1-cfhcygc(igsp))*hcyg(ix,iy,igsp)+
     .                     cfhcygc(igsp)*noavey*kyg_use(ix,iy,igsp)
          enddo
        enddo
       enddo
      endif

c ... Gas molecule thermal equipartition with hydrogen ions and atoms
*****************************************************************
      if (nisp >= 2) then   # uses ni(,,2), so must have atoms
       do igsp = 1, ngsp
        do iy = j1, j6
          do ix = i1, i6
	    nhi_nha = ni(ix,iy,1)+ni(ix,iy,2)
            eqpg(ix,iy,igsp) = cftgeqp*ng(ix,iy,igsp)*nhi_nha*
     .                                            keligig(igsp)
          enddo
        enddo
       enddo
      endif

c ... Call routine to evaluate gas energy fluxes
****************************************************************
      call engbalg


*****************************************************************
*****************************************************************
*  Here starts the old routine PARBAL
*****************************************************************
      do 104 ifld = 1, nfsp
*  ---------------------------------------------------------------------
*     compute flux, residual
*     The residual is: res := snic + sniv * ni - outflow(ni).
*  ---------------------------------------------------------------------

*  -- compute fnix --

         methnx = mod(methn, 10)
         methny = methn/10
         do 81 iy = j4, j8
            do 80 ix = i1, i5
              if ( zi(ifld).eq.0. .and. ineudif.ne.0 .and.
     .                                   1.-rrv(ix,iy) > 1.e-4 ) then
                 fnix(ix,iy,ifld) = fngx(ix,iy,1)
              else
               ix2 = ixp1(ix,iy)

               if (methnx .eq. 2) then   # central differencing
                  t2 = ( ni(ix, iy,ifld) + ni(ix2,iy,ifld) ) / 2

               elseif (methnx .eq. 3) then   # upwind differencing

                  if( uu(ix,iy,ifld) .ge. 0.) then
                     t2 = ni(ix,iy,ifld)
                  else
                     t2 = ni(ix2,iy,ifld)
                  endif

               elseif (methnx .eq. 6) then   # log central differencing
                  t2 = exp(0.5*
     .                ( log(ni(ix,iy,ifld)) + log(ni(ix2,iy,ifld)) ))

               else   # interp. ave or harmonic ave depending on wind*grad

                  t0 = ( ni(ix, iy,ifld)*gx(ix, iy) +
     .                   ni(ix2,iy,ifld)*gx(ix2,iy) ) /
     .                                      ( gx(ix,iy)+gx(ix2,iy) )
                  t1 = ( gx(ix,iy)+gx(ix2,iy) ) * ni(ix,iy,ifld) *
     .                   ni(ix2,iy,ifld) / ( cutlo + ni(ix,iy,ifld)*
     .                   gx(ix2,iy) + ni(ix2,iy,ifld)*gx(ix,iy) )
                  if( uu(ix,iy,ifld)*(ni(ix,iy,ifld)-ni(ix2,iy,ifld))
     .                                                     .ge. 0.) then
                     t2 = t0
                  else
                     t2 = t1
                  endif

               endif

               fnix(ix,iy,ifld) = cnfx*uu(ix,iy,ifld) * sx(ix,iy) * t2
               fnixcb(ix,iy,ifld)=cnfx*sx(ix,iy) * t2 * 0.5*
     .                 (rbfbt(ix,iy) + rbfbt(ix2,iy))*v2cb(ix,iy,ifld)
                  fnix(ix,iy,ifld) = fnix(ix,iy,ifld)/sqrt( 1 +
     .              (nlimix(ifld)*ni(ix ,iy,ifld)/ni(ix2,iy,ifld))**2 +
     .              (nlimix(ifld)*ni(ix2,iy,ifld)/ni(ix ,iy,ifld))**2 )
              endif
   80      continue
           if ((isudsym==1.or.geometry.eq.'dnXtarget') .and. nxc > 1) then
              fnix(nxc-1,iy,ifld)=0.
              fnix(nxc,  iy,ifld)=0.
              fnix(nxc+1,iy,ifld)=0.
              uu(nxc-1,iy,ifld) = 0.
              uu(nxc  ,iy,ifld) = 0.
              uu(nxc+1,iy,ifld) = 0.
              vytan(nxc-1,iy,ifld) = 0.
              vytan(nxc  ,iy,ifld) = 0.
              vytan(nxc+1,iy,ifld) = 0.

           endif
           if (islimon.ne.0 .and. iy.ge.iy_lims) fnix(ix_lim,iy,ifld)=0.
           if (nxpt==2 .and. ixmxbcl==1) fnix(ixrb(1)+1,iy,ifld)=0.
   81    continue


*  -- compute fniy  --

         do 83 iy = j1, j5
            do 82 ix = i4, i8
               if (zi(ifld).eq.0.) then #inertial gas must follow ion index
                  fniy(ix,iy,ifld) = fngy(ix,iy,ifld-1)
               else
                  if (methny .eq. 2) then   # central differencing
                     t2 = ( niy0(ix,iy,ifld) + niy1(ix,iy,ifld) ) / 2

                  elseif (methny .eq. 3) then   # upwind differencing

                     if( vy(ix,iy,ifld) .ge. 0.) then
                        t2 = niy0(ix,iy,ifld)
                     else
                        t2 = niy1(ix,iy,ifld)
                     endif

                  elseif (methny .eq. 6) then   # log central differencing
                     t2 = exp( 0.5*
     .                   (log(niy0(ix,iy,ifld))+log(niy1(ix,iy,ifld))) )

                  else    # interp. ave or harmonic ave depending on wind*grad

                     t0 = ( niy0(ix,iy,ifld)*gy(ix,iy  ) +
     .                    niy1(ix,iy,ifld)*gy(ix,iy+1) ) /
     .                    ( gy(ix,iy)+gy(ix,iy+1) )
                     t1 = ( gy(ix,iy)+gy(ix,iy+1) ) * niy0(ix,iy,ifld)*
     .                    niy1(ix,iy,ifld) / ( cutlo + niy0(ix,iy,ifld)*
     .                    gy(ix,iy+1) + niy1(ix,iy,ifld)*gy(ix,iy) )
                     if( (niy0(ix,iy,ifld)-niy1(ix,iy,ifld))*
     .                    vy(ix,iy,ifld) .ge. 0.) then
                        t2 = t0
                     else
                        t2 = t1
                     endif

                  endif

                  fniy(ix,iy,ifld) = cnfy*vy(ix,iy,ifld)*sy(ix,iy)*t2
                  fniycb(ix,iy,ifld) = cnfy*vycb(ix,iy,ifld)*sy(ix,iy)*t2
                  if (vy(ix,iy,ifld)*(ni(ix,iy,ifld)-ni(ix,iy+1,ifld))
     .                                                      .lt. 0.) then
                     fniy(ix,iy,ifld) = fniy(ix,iy,ifld)/( 1 +
     .                               (nlimiy(ifld)/ni(ix,iy+1,ifld))**2 +
     .                               (nlimiy(ifld)/ni(ix,iy  ,ifld))**2 )
                  endif
c...  Note: nonorthogonality comes in through calc. of vy
               endif
 82         continue
 83      continue

c ... cosmetic setting of fniy - not used
         do ix = i4, i8
            fniy(ix,ny+1,ifld) = 0.0e0
         enddo

 104  continue

c ... Add rad flux of 4th order diff operator; damp grid-scale oscillations
      do ifld = 1, nfsp
        if (abs(dif4order(ifld)) > 1.e-50) then
          do iy = j2p, j5m   #limits to range iy=1:ny-1 for fniy4ord
            iym1 = max(iy-1,0)
            iyp1 = min(iy+1,ny+1)
            iyp2 = min(iy+2,ny+1)
            do ix = i4, i8
              dndym1 = (ni(ix,iy,ifld)-ni(ix,iym1,ifld))*gyf(ix,iym1)
              dndy0 = (ni(ix,iyp1,ifld)-ni(ix,iy,ifld))*gyf(ix,iy)
              dndyp1 = (ni(ix,iyp2,ifld)-ni(ix,iyp1,ifld))*gyf(ix,iyp1)
              d2ndy20 = (dndy0 - dndym1)*gy(ix,iy)
              d2ndy2p1 = (dndyp1 - dndy0)*gy(ix,iyp1)
              d3ndy3 = (d2ndy2p1 - d2ndy20)*gyf(ix,iy)
              fniy4ord(ix,iy,ifld) = dif4order(ifld)*d3ndy3*sy(ix,iy)/
     .                                                  gyf(ix,iy)**2
              fniy(ix,iy,ifld) = fniy(ix,iy,ifld) + fniy4ord(ix,iy,ifld)
            enddo
          enddo
        endif
      enddo

c ... Setup a correction to surface-flux for grad_B and grad_P effects at iy=0
      do ifld = 1, nfsp
        do ix = i4, i8
           fniycbo(ix,ifld) = 0.0
        enddo
         do ix = i4, i8
            fniycbo(ix,ifld) = ( ni(ix,0,ifld)*sy(ix,0) ) *
     .                         ( (1-cfniybbo)*cfybf*vycb(ix,0,ifld) -
     .                            cfniydbo*(1-cfydd)*vycp(ix,0,ifld) )
         enddo
      enddo


c----------------------------------------------------------------------c
c          SCALE SOURCE TERMS FROM MONTE-CARLO-NEUTRALS MODEL
c
c     These sources are used in the residuals (resco,resmo,resee,resei)
c     so the call to scale_mcn must occur BEFORE these residuals are
c     evaluated.  Since they scale with fnix at the divertor plates,
c     the call to scale_mcn must occur AFTER fnix has been calculated.


      if (ismcnon .ne. 0) then
c 	     write(*,*) 'TEST ISMCNON START: ismcnon=',ismcnon
c         call scale_mcn
         call scale_mcnsor
      endif
c----------------------------------------------------------------------c

*  -- compute the residual if isnion = 1 --

      do ifld = 1, nfsp
       do 86 iy = j2, j5
         do 85 ix = i2, i5
	   if(isnionxy(ix,iy,ifld) == 1) then
              resco(ix,iy,ifld) =
     .           snic(ix,iy,ifld)+sniv(ix,iy,ifld)*ni(ix,iy,ifld) +
     .           volpsor(ix,iy,ifld) +
     .           cfneut * cfneutsor_ni * cnsor * psor(ix,iy,ifld) +
     .           cfneut * cfneutsor_ni * cnsor * psorxr(ix,iy,ifld)
     .           -nuvl(ix,iy,ifld)*vol(ix,iy)*ni(ix,iy,ifld) +
     .           voljcsor(ix,iy)/qe
           endif
c           if (ifld .ne. iigsp) then
	       if(zi(ifld) .ne. 0) then # IJ 2016 skip if neutral zi(ifld)=0
              resco(ix,iy,ifld) = resco(ix,iy,ifld) + cmneut * uesor_ni(ix,iy,ifld)
           else # IJ 2016 zi==0, assume neutral->ifld and ion->ifld-1
              resco(ix,iy,ifld) = resco(ix,iy,ifld) - cmneut * uesor_ni(ix,iy,ifld-1)
           endif
   85    continue
   86  continue

       do 302 iy = j2, j5
         do 301 ix = i2, i5
	       if(isnionxy(ix,iy,ifld) == 1) then
              ix1 = ixm1(ix,iy)
	        if(zi(ifld) .ne. 0) then # IJ 2016 skip if neutral  zi(ifld)=0
                 resco(ix,iy,ifld) = resco(ix,iy,ifld)
     .                                   - ( (fnix(ix,iy,ifld) - fnix(ix1,iy, ifld))
     .                            + fluxfacy*(fniy(ix,iy,ifld) - fniy(ix,iy-1,ifld)) )
              else ## zi==0
                 resco(ix,iy,ifld) = resco(ix,iy,ifld)
     .              - cfneutdiv*cfneutdiv_fng*( (fnix(ix,iy,ifld)-fnix(ix1,iy, ifld))
     .                               + fluxfacy*(fniy(ix,iy,ifld)-fniy(ix,iy-1,ifld)) )

c ... IJ 2016/10/19 add MC neutral flux if flags set
                 if (get_neutral_moments .and. cmneutdiv_fng .ne. 0.0) then
                    jfld=1  ## assume main ions in ifld=1
                    sng_ue(ix,iy,jfld) = - ( (fngx_ue(ix,iy,jfld) - fngx_ue(ix1,iy, jfld))
     .                        +   fluxfacy*(fngy_ue(ix,iy,jfld) - fngy_ue(ix,iy-1,jfld)) )
     .                        *( (ng(ix,iy,jfld)*ti(ix,iy))/(ng(ix,iy,jfld)*ti(ix,iy)) )
c                   if (ix .eq. 1 .and. iy .eq. 1) write(*,*) 'sng_ue', ifld, jfld
                    resco(ix,iy,ifld) = resco(ix,iy,ifld) +
     .                                  cmneutdiv*cmneutdiv_fng*sng_ue(ix,iy,jfld)
                 endif
              endif
           endif
 301     continue
 302   continue
       enddo       # end of ifld loop

if (TimingPandfOn.gt.0) TimeMomBal=tick()
*********************************************************************
c  Here we do the neutral gas diffusion model
c  The diffusion is flux limited using the thermal flux
**********************************************************************

ccc         if(isngon .eq. 1) call neudif

*****************************************************************
*****************************************************************
*  Here starts the old MOMBAL_B2
*****************************************************************


*  ---------------------------------------------------------------------
*  loop over all species.
*  ---------------------------------------------------------------------

      do 105 ifld = 1, nusp
      if(isupon(ifld) .eq. 0) goto 105
*     ------------------------------------------------------------------
*     compute the residual.
*     ------------------------------------------------------------------

*  -- evaluate flox and conx --

         do 91 iy = j4, j8
            flox(0,iy) = 0.0e0
            conx(0,iy) = 0.0e0
            do 90 ix = i2, i6
               ix1 = ixm1(ix,iy)
               if (isimpon.ge.5 .and. ifld.eq.1) then
                   #up(,,1) is total mass vel, whereas uu(,,i) for each ion
                  uuv =0.5*( (up(ix1,iy,ifld)*rrv(ix1,iy)+
     .                        up(ix,iy,ifld)*rrv(ix,iy)) +
     .                     (v2(ix1,iy,ifld)+v2(ix,iy,ifld))*rbfbt(ix,iy)-
     .                     (vytan(ix1,iy,ifld)+vytan(ix,iy,ifld)) )
               else
                  uuv = 0.5 * (uu(ix1,iy,ifld)+uu(ix,iy,ifld))
               endif
               flox(ix,iy) = cmfx * nm(ix,iy,ifld) * uuv *
     .                          vol(ix,iy) * gx(ix,iy)
ccc Distance between veloc. cell centers:
               if (isgxvon .eq. 0) then     # dx(ix)=1/gx(ix)
                 conx(ix,iy) = visx(ix,iy,ifld) * vol(ix,iy) * gx(ix,iy)
     .                            * gx(ix,iy)
               elseif (isgxvon .eq. 1) then # dx(ix)=.5/gxf(ix-1) + .5/gxf(ix)
                 conx(ix,iy) = visx(ix,iy,ifld) * vol(ix,iy) * gx(ix,iy)
     .               * 2*gxf(ix,iy)*gxf(ix1,iy)/(gxf(ix,iy)+gxf(ix1,iy))
               endif
   90       continue
   91    continue

*  -- evaluate floy and cony without averaging over two ix cells --

         do 93 iy = j1, j5
            if (nxpt == 1 .or. iy <= iysptrx1(1)) then
              iysepu = iysptrx1(1)
              if (ndomain > 1) iysepu = iysptrxg(mype+1)  # and ixpt1,2u??
              ixpt1u = ixpt1(1)
              ixpt2u = ixpt2(1)
            else  # nxpt=2 and iy > iysptrx1(1), use second separatrix
              iysepu = iysptrx1(2)
              ixpt1u = ixpt1(2)
              ixpt2u = ixpt2(2)
            endif
            do 92 ix = i4, i8
               ix2 = ixp1(ix,iy)
               ix4 = ixp1(ix,iy+1)
               cony(ix,iy) = .5 * syv(ix,iy) *
     .                       (ave(visy(ix,iy,ifld)*gy(ix,iy),
     .                           visy(ix,iy+1,ifld)*gy(ix,iy+1)) +
     .                        ave(visy(ix2,iy,ifld)*gy(ix2,iy),
     .                            visy(ix4,iy+1,ifld) * gy(ix4,iy+1)))
               if (iy==iysepu .and. (ix==ixpt1u .or. ix==ixpt2u)) then
                 cony(ix,iy) = syv(ix,iy) *
     .                       ( ave(visy(ix,iy,ifld)*gy(ix,iy),
     .                             visy(ix,iy+1,ifld)*gy(ix,iy+1)) )
                 floy(ix,iy) = (cmfy/2) * syv(ix,iy) *(
     .                    ave(nm(ix,iy,ifld),nm(ix,iy+1,ifld)) )*
     .                                         vy(ix,iy,ifld)
                 if (ifld==1) then  # add user-specified convection
                   floy(ix,iy) = floy(ix,iy)+(cmfy/2)*syv(ix,iy) *(
     .                    ave(nm(ix,iy,ifld),nm(ix,iy+1,ifld)) )*
     .                                                 vyup_use(ix,iy)
                 endif
               elseif(isugfm1side == 1 .and. zi(ifld) == 0.) then
                 floy(ix,iy) = (cmfy/4) * syv(ix,iy) *(
     .                    ave(nm(ix,iy,ifld),nm(ix,iy+1,ifld)) +
     .                    ave(nm(ix2,iy,ifld),nm(ix4,iy+1,ifld))) *
     .                       ( vy(ix,iy,ifld) + vy(ix,iy,ifld) )
               else
                 floy(ix,iy) = (cmfy/4) * syv(ix,iy) *(
     .                    ave(nm(ix,iy,ifld),nm(ix,iy+1,ifld)) +
     .                    ave(nm(ix2,iy,ifld),nm(ix4,iy+1,ifld))) *
     .                       ( vy(ix,iy,ifld) + vy(ix2,iy,ifld) )
                 if (ifld==1) then  # add user-specified convection
                    floy(ix,iy) = floy(ix,iy)+(cmfy/4)*syv(ix,iy) *(
     .                    ave(nm(ix,iy,ifld),nm(ix,iy+1,ifld)) +
     .                    ave(nm(ix2,iy,ifld),nm(ix4,iy+1,ifld))) *
     .                       ( vyup_use(ix,iy) + vyup_use(ix2,iy) )
                 endif
               endif
	     if(ishavisy == 1) then
               cony(ix,iy) = .5 * syv(ix,iy) *
     .                       (ave(visy(ix,iy,ifld)*gy(ix,iy),
     .                           visy(ix,iy+1,ifld)*gy(ix,iy+1)) +
     .                        ave(visy(ix2,iy,ifld)*gy(ix2,iy),
     .                            visy(ix4,iy+1,ifld) * gy(ix4,iy+1)))
             else
               cony(ix,iy) = .25 * cfaccony*syv(ix,iy) *
     .                       ( visy(ix,iy,ifld)*gy(ix,iy) +
     .                         visy(ix,iy+1,ifld)*gy(ix,iy+1) +
     .                         visy(ix2,iy,ifld)*gy(ix2,iy) +
     .                         visy(ix4,iy+1,ifld)*gy(ix4,iy+1) )
             endif

   92       continue
   93    continue

*  -- compute the momentum transport --

         call fd2tra (nx,ny,flox,floy,conx,cony,
     .                up(0,0,ifld),fmix(0,0,ifld),
     .                fmiy(0,0,ifld),1, methu)

      if (isnonog .eq. 1) then

c     Compute y-component fmixy of nonorthogonal diffusive momentum flux.
c     The convective component is already already added through uu(ix,iy).
c     Average fym, etc in ix to get staggered velocity-grid values fymv, etc.
c     The density-stencil dxnog has to be averaged as well.
         do 96 iy = j2, j5
            iy1 = max(iy-1,0)
            do 97 ix = i2, i5+1    # ixp1(i5,iy)
               ix1 = ixm1(ix,iy)
               ix3 = ixm1(ix,iy1)
               ix5 = ixm1(ix,iy+1)
               grdnv = (
     .                  fymv (ix,iy,1)*up(ix ,iy1 ,ifld)+
     .                  fy0v (ix,iy,1)*up(ix ,iy  ,ifld)+
     .                  fypv (ix,iy,1)*up(ix ,iy+1,ifld)+
     .                  fymxv(ix,iy,1)*up(ix3,iy1 ,ifld)+
     .                  fypxv(ix,iy,1)*up(ix5,iy+1,ifld)-
     .                  fymv (ix,iy,0)*up(ix3,iy1 ,ifld)-
     .                  fy0v (ix,iy,0)*up(ix1,iy  ,ifld)-
     .                  fypv (ix,iy,0)*up(ix5,iy+1,ifld)-
     .                  fymxv(ix,iy,0)*up(ix ,iy1 ,ifld)-
     .                  fypxv(ix,iy,0)*up(ix ,iy+1,ifld) ) *
     .                     (2*gxfn(ix,iy)*gxfn(ix1,iy) /
     .                      (gxfn(ix,iy)+gxfn(ix1,iy)))
               if (isgxvon .eq. 0) then
                  fmixy(ix,iy,ifld) = cfvisxy(ifld)*visy(ix,iy,ifld) *
     .              ( grdnv/cos(0.5*(angfx(ix1,iy)+angfx(ix,iy))) -
     .               (up(ix,iy,ifld) - up(ix1,iy,ifld))*gx(ix,iy) ) *
     .              0.5*(sx(ix1,iy)+sx(ix,iy))
               elseif (isgxvon .eq. 1) then
                  fmixy(ix,iy,ifld) = cfvisxy(ifld)*visy(ix,iy,ifld) *
     .              ( grdnv/cos(0.5*(angfx(ix1,iy)+angfx(ix,iy))) -
     .               (up(ix,iy,ifld) - up(ix1,iy,ifld))*
     .                     ( 2*gxf(ix,iy)*gxf(ix1,iy) /
     .                        (gxf(ix,iy)+gxf(ix1,iy)) ) ) *
     .                     0.5*(sx(ix1,iy)+sx(ix,iy))
               endif
c...  Now flux limit with flalfvgxy if ifld=2
               if (ifld==2) then
                 t0 = max(tg(ix,iy,1),temin*ev)
                 t1 = max(tg(ix2,iy,1),temin*ev)
                 vtn = sqrt(t0/mg(1))
                 qfl = flalfvgxya(ix)*0.5*(sx(ix,iy)+sx(ix1,iy))*vtn**2*
     .                                        nm(ix,iy,ifld) + cutlo
                 fmixy(ix,iy,ifld) = fmixy(ix,iy,ifld) /
     .                             sqrt(1+(fmixy(ix,iy,ifld)/qfl)**2)
               endif

 97         continue
 96      continue
      endif

c...  Compute viscous drag from nonuniform B-field, then add to smoc
      if (isupdrag .eq. 1 .and. ifld .eq. 1) then
        do iy = j2, j5
          do ix = i2, i5
            ix1 = ixm1(ix,iy)
            ix2 = ixp1(ix,iy)
c ...   First, the short mfp drag
            b_ctr = 0.5*(btot(ix,iy)+btot(ix2,iy))
                    # derviatives dbds_m and dbds_p are one-sided
            dbds_m = (btot(ix,iy) - btot(ix1,iy))*
     .                                      gxf(ix1,iy)*rrv(ix1,iy)
            dbds_p = (btot(ix2,iy) - btot(ix,iy))*
     .                                        gxf(ix,iy)*rrv(ix,iy)
            eta_h0 = visx(ix,iy,1)/b_ctr**2.5
            eta_hm = 0.5*(visx(ix,iy,1)+visx(ix1,iy,1))/
     .                                              btot(ix,iy)**2.5
            eta_hp = 0.5*(visx(ix2,iy,1)+visx(ix,iy,1))/
     .                                              btot(ix2,iy)**2.5
            drag_1 = -2*eta_h0*0.5*(up(ix2,iy,1)-up(ix1,iy,1))*
     .                            (btot(ix2,iy)-btot(ix,iy))*
     .                                   (gxf(ix,iy)*rrv(ix,iy))**2
            drag_2 = up(ix,iy,1)*(eta_hp*dbds_p - eta_hm*dbds_m)*
     .                                        gxf(ix,iy)*rrv(ix,iy)
c ...   now for the trapped particle drag (sloppy)
            nu_ii = ni(ix,iy,1)*(2*mp/mi(1))**0.5/
     .                                   (3e12*(ti(ix,iy)*ev)**1.5)
            drag_3 = -mi(1)*ni(ix,iy,1)*up(ix,iy,1)*nu_ii*frac_pt
            mf_path = (2*ti(ix,iy)/mi(1))**0.5 / nu_ii
            frac_col = 1 / (1 + (mf_path/con_leng)**2)
            smoc(ix,iy,1) = smoc(ix,iy,1) + (0.6666667*frac_col*
     .                      b_ctr**2.5*(drag_1 + drag_2) +
     .                      (1 - frac_col)*drag_3) * volv(ix,iy)
          enddo
        enddo
      endif

c...  Compute total viscosity for nonuniform B-field; put in visvol_v,q
      if (cfvisxneov+cfvisxneoq > 0.) call upvisneo

c...  Now fix the fluxes touching the x-point(s):


         do k = 1, nxpt   # loop over all x-points
           k1 = k      # region argument of ixpt1 that touches this x-point
           k2 = k-1    # region argument of ixpt2 that touches this x-point
           if (k==1) k2 = nxpt

           if (nxpt==2) then      # set ghxpt,gvxpt,sxyxpt for full double null
              if (k==1) then      # this is the lower x-point
                 ghxpt = ghxpt_lower
                 gvxpt = gvxpt_lower
                 sxyxpt = sxyxpt_lower
              elseif (k==2) then  # this is the upper x-point
                 ghxpt = ghxpt_upper
                 gvxpt = gvxpt_upper
                 sxyxpt = sxyxpt_upper
              endif
           endif

	   if( ((2*(yc-iysptrx1(k1))-1)/4 .le. 1) .or. j1 == 0 ) then
           if( ((2*(xc-ixpt1(k1))-1)/4)*((2*(xc-ixpt2(k2))-1)/4).eq.0 .or.
     .                                                        i1.eq.0 ) then
           if(isnfmiy .eq. 1) then

           fmiy(ixpt1(k1),iysptrx1(k1),ifld) = 0.
           fmiy(ixpt2(k2),iysptrx2(k2),ifld) = 0.
           nixpt(ifld,k1) = 0.125 * (
     .           ni(ixpt1(k1),iysptrx1(k1)  ,ifld) + ni(ixpt1(k1)+1,iysptrx1(k1)  ,ifld)
     .         + ni(ixpt1(k1),iysptrx1(k1)+1,ifld) + ni(ixpt1(k1)+1,iysptrx1(k1)+1,ifld)
     .         + ni(ixpt2(k2),iysptrx2(k2)  ,ifld) + ni(ixpt2(k2)+1,iysptrx2(k2)  ,ifld)
     .         + ni(ixpt2(k2),iysptrx2(k2)+1,ifld) + ni(ixpt2(k2)+1,iysptrx2(k2)+1,ifld) )
           visyxpt(ifld,k1) = 0.125 * (
     .        visy(ixpt1(k1),iysptrx1(k1)  ,ifld) + visy(ixpt1(k1)+1,iysptrx1(k1)  ,ifld)
     .      + visy(ixpt1(k1),iysptrx1(k1)+1,ifld) + visy(ixpt1(k1)+1,iysptrx1(k1)+1,ifld)
     .      + visy(ixpt2(k2),iysptrx2(k2)  ,ifld) + visy(ixpt2(k2)+1,iysptrx2(k2)  ,ifld)
     .      + visy(ixpt2(k2),iysptrx2(k2)+1,ifld) + visy(ixpt2(k2)+1,iysptrx2(k2)+1,ifld) )
           upxpt(ifld,k1) = 0.25 * (
     .           up(ixpt1(k1),iysptrx1(k1)  ,ifld) + up(ixpt2(k2),iysptrx2(k2)  ,ifld)
     .         + up(ixpt1(k1),iysptrx1(k1)+1,ifld) + up(ixpt2(k2),iysptrx2(k2)+1,ifld) )
           vyvxpt(ifld,k1) = (0.707*0.25) * (
     .           vy(ixpt1(k1)  ,iysptrx1(k1),ifld) - vy(ixpt2(k2)  ,iysptrx2(k2),ifld)
     .         - vy(ixpt1(k1)+1,iysptrx1(k1),ifld) + vy(ixpt2(k2)+1,iysptrx2(k2),ifld) )
           vyhxpt(ifld,k1) = (0.707*0.25) * (
     .         - vy(ixpt1(k1)  ,iysptrx1(k1),ifld) + vy(ixpt2(k2)  ,iysptrx2(k2),ifld)
     .         - vy(ixpt1(k1)+1,iysptrx1(k1),ifld) + vy(ixpt2(k2)+1,iysptrx2(k2),ifld) )
cccMER The convective contributions to fmihxpt and fmivxpt seem to have an
cccMER erroneous multiplicative factor -1/2 (from original code) ???
           fmihxpt(ifld,k1) = cfnfmiy*( - cmfy*(mi(ifld)/2)*sxyxpt*nixpt(ifld,k1)*
     .                     vyhxpt(ifld,k1)*upxpt(ifld,k1)
     .                   - sxyxpt*visyxpt(ifld,k1)*(up(ixpt2(k2),iysptrx2(k2)+1,ifld)
     .                     - up(ixpt1(k1),iysptrx1(k1)+1,ifld))*ghxpt )
           fmivxpt(ifld,k1) = cfnfmiy*(- cmfy*(mi(ifld)/2)*sxyxpt*nixpt(ifld,k1)*
     .                     vyvxpt(ifld,k1)*upxpt(ifld,k1)
     .                   - sxyxpt*visyxpt(ifld,k1)*(up(ixpt2(k2),iysptrx2(k2),ifld)
     .                     - up(ixpt1(k1),iysptrx1(k1),ifld))*gvxpt )
           smoc(ixpt1(k1),iysptrx1(k1)+1,ifld) = smoc(ixpt1(k1),iysptrx1(k1)+1,ifld)
     .                                - fmihxpt(ifld,k1)
           smoc(ixpt2(k2),iysptrx2(k2)+1,ifld) = smoc(ixpt2(k2),iysptrx2(k2)+1,ifld)
     .                                + fmihxpt(ifld,k1)
           smoc(ixpt1(k1),iysptrx1(k1)  ,ifld) = smoc(ixpt1(k1),iysptrx1(k1)  ,ifld)
     .                                - fmivxpt(ifld,k1)
           smoc(ixpt2(k2),iysptrx2(k2)  ,ifld) = smoc(ixpt2(k2),iysptrx2(k2)  ,ifld)
     .                                + fmivxpt(ifld,k1)

           endif # end if-test on isnfmiy
           endif # end if-test on xc
           endif # end if-test on yc

         enddo # end do-loop over nxpt x-points

*  -- source term and pressure gradient --

         do 99 iy = j2, j5
            do 98 ix = i2, i5
               ix2 = ixp1(ix,iy)
               if (zi(ifld) .ne. 0) then  # additions only for charged ions
                  dp1 =  cngmom(ifld)*(1/fac2sp)*
     .                      ( ng(ix2,iy,1)*tg(ix2,iy,1)-
     .                        ng(ix ,iy,1)*tg(ix ,iy,1) )
                  resmo(ix,iy,ifld) = 0.
                  resmo(ix,iy,ifld) =
     .                  smoc(ix,iy,ifld)
     .                + smov(ix,iy,ifld) * up(ix,iy,ifld)
     .                - cfneut * cfneutsor_mi * sx(ix,iy) * rrv(ix,iy) * dp1
     .                - cfneut * cfneutsor_mi * cmwall(ifld)*0.5*(ng(ix,iy,1)+ng(ix2,iy,1))
     .                      * mi(ifld)*up(ix,iy,ifld)*0.5
     .                      *(nucx(ix,iy,1)+nucx(ix2,iy,1))*volv(ix,iy)
     .                + cmneut * cmneutsor_mi * uesor_up(ix,iy,ifld)
     .                + cfmsor*(msor(ix,iy,ifld) + msorxr(ix,iy,ifld)) #### IJ 2017: needs *cfneut for multi-charge state ions & MC neutrals?
     .                + volmsor(ix,iy,ifld)
     .                + cfvisxneov*visvol_v(ix,iy,ifld)
     .                + cfvisxneoq*visvol_q(ix,iy,ifld)
c  Add drag with cold, stationary impurity neutrals
                  if (ifld > nhsp) then
                    resmo(ix,iy,ifld) = resmo(ix,iy,ifld) - cfupimpg*
     .                    0.25*mi(ifld)*(ni(ix,iy,ifld)+ni(ix2,iy,ifld))*
     .                     ( nucxi(ix,iy,ifld)+nucxi(ix2,iy,ifld)+
     .                       nueli(ix,iy,ifld)+nueli(ix2,iy,ifld) )*
     .                       up(ix,iy,ifld)*volv(ix,iy)
                  endif
               endif

               if (isupgon(1) .eq. 1) then

c     If we are solving the parallel neutral mom eq. we need different/addtnl
c     source terms. Beware that cngmom and cmwall should be zero so that the
c     main ions do not get coupled twice to the neutrals!
c     The CX, ionization, recomb. friction for the parallel momentum eqs
c     for the neutrals and main ions are included here.
c     Assumes the neutrals have index iigsp and corresponding ions index 1

                  if (ifld .eq. 1) then
c     The main ions, momentum coupling:
                     resmo(ix,iy,ifld) = resmo(ix,iy,ifld)
     .                   + cfneut * cfneutsor_mi * cfupcx*0.25*volv(ix,iy)*
     .                       (nucx(ix,iy,1)+nucx(ix2,iy,1))*
     .                       (nm(ix,iy,iigsp)+nm(ix2,iy,iigsp))*
     .                       (up(ix,iy,iigsp)-up(ix,iy,1))
     .                   + cfneut * cfneutsor_mi * 0.25*volv(ix,iy)*
     .                       (  (nuiz(ix,iy,1)+nuiz(ix2,iy,1))*
     .                          (nm(ix,iy,iigsp)+nm(ix2,iy,iigsp))*
     .                          up(ix,iy,iigsp)
     .                        - (nurc(ix,iy,1)+nurc(ix2,iy,1))*
     .                          (nm(ix,iy,1)+nm(ix2,iy,1))*up(ix,iy,1)
     .                       )
                  elseif ((isupgon(1) .eq. 1) .and. ifld .eq. iigsp) then
c     The neutral species, momentum coupling AND other source terms:
                      resmo(ix,iy,iigsp) =   # TR resmo(ix,iy,ifld) #IJ 2016
     .                    - cmneut * cmneutsor_mi * uesor_up(ix,iy,1)
     .                    -sx(ix,iy) * rrv(ix,iy) *
     .                         cpgx*( ni(ix2,iy,iigsp)*ti(ix2,iy)-
     .                                ni(ix,iy,iigsp)*ti(ix,iy) )
     .                    -cfupcx*0.25*volv(ix,iy)*
     .                       (nucx(ix,iy,1)+nucx(ix2,iy,1))*
     .                       (nm(ix,iy,iigsp)+nm(ix2,iy,iigsp))*
     .                       (up(ix,iy,iigsp)-up(ix,iy,1))
     .                    -0.25*volv(ix,iy)*(
     .                       (nuiz(ix,iy,1)+nuiz(ix2,iy,1))*
     .                       (nm(ix,iy,iigsp)+nm(ix2,iy,iigsp))*
     .                                            up(ix,iy,iigsp)
     .                       -(nurc(ix,iy,1)+nurc(ix2,iy,1))*
     .                       (nm(ix,iy,1)+nm(ix2,iy,1))*up(ix,iy,1) )
                  endif
               endif
 98         continue
 99      continue

*  -- divergence of momentum flow --

         if (isnonog.eq.1) then
            do 3051 iy = j2, j5
               do 3061 ix = i2, i5
                  ix2 = ixp1(ix,iy)
c ... IJ 2016/10/10 use cfneutdiv_fmg multiplier for neutrals
c                 if (ifld .ne. iigsp) then
	          if(zi(ifld) .ne. 0) then # IJ 2016; depends if ion or neut
                    resmo(ix,iy,ifld) = resmo(ix,iy,ifld)
     .                 + (fmixy(ix2,iy,ifld) - fmixy(ix,iy,ifld))
                  else
                    resmo(ix,iy,ifld) = resmo(ix,iy,ifld)
     .                 + cfneutdiv*cfneutdiv_fmg*(fmixy(ix2,iy,ifld) - fmixy(ix,iy,ifld))
c***	IJ 2017/09/21: Need to add similar fmgxy calculation for MC neutrals on nonorthogonal mesh ***
                  endif
 3061          continue
 3051       continue
         endif

         do 305 iy = j2, j5
            do 306 ix = i2, i5
               ix2 = ixp1(ix,iy)
c IJ 2016/10/10 add cfneutdiv_fmg multiplier for neutrals to control fraction of momentum to add
c               if (ifld .ne. iigsp) then
c ... IJ 2016 resmo contrib changes if ion or neut
	       if(zi(ifld) .ne. 0) then
                 resmo(ix,iy,ifld) = resmo(ix,iy,ifld)
     .                                 - (fmix(ix2,iy,ifld) - fmix(ix,iy  ,ifld)
     .                        + fluxfacy*(fmiy(ix ,iy,ifld) - fmiy(ix,iy-1,ifld)) )
               else
                 resmo(ix,iy,ifld) = resmo(ix,iy,ifld)
     .                      - cfneutdiv*cfneutdiv_fmg*(fmix(ix2,iy,ifld) - fmix(ix,iy  ,ifld)
     .                        + fluxfacy*(fmiy(ix ,iy,ifld) - fmiy(ix,iy-1,ifld)) )
                 if(cmneutdiv_fmg .ne. 0.0) then
                    jfld=1
                    resmo(ix,iy,ifld) = resmo(ix,iy,ifld)
     .                    - cmneutdiv*cmneutdiv_fmg*( (fmgx_ue(ix2,iy,jfld) - fmgx_ue(ix,iy  ,jfld))
     .                        + fluxfacy*(fmgy_ue(ix ,iy,jfld) - fmgy_ue(ix,iy-1,jfld)) )
     .                        * (ni(ix,iy,ifld)*ti(ix,iy))/(ni(ix,iy,ifld)*ti(ix,iy))
                 endif
              endif
  306       continue
  305    continue

c  -- Include frictional drag in parallel flow here if isofric=1; otherwise
c  -- it is included in frici from mombal or mombalni

        if (isofric.eq.1 .and. nusp .gt.1) then
*  -- w0 now accumulates friction coefficient --
*     -- set w2 = vol*ti**(-1.5) --
         do iy = j1, j6
           do ix = i1, i6
             fricnrl(ix,iy,ifld) = 0.  #diagnostic ~ ni*mi*nu*(up1-up2)
             w0(ix,iy) = 0.0e0
             w2(ix,iy) = vol(ix,iy) / (ti(ix,iy)*sqrt(ti(ix,iy)))
           enddo
         enddo

*  -- consider all other species --

         do jfld = 1, nusp
           if (jfld .ne. ifld) then
*     -- common factor in collision frequency --
             awoll = zi(ifld)**2 * zi(jfld)**2 *
     .            (qe**4/(12*pi**2*eps0**2)) *
     .            sqrt (2*pi*mi(ifld)*mi(jfld)/(mi(ifld)+mi(jfld)))

*     -- frictional coupling --
             do iy = j1, j6
               do ix = i1, i5
                 ix2 = ixp1(ix,iy)
                 t0 = ni(ix,iy,ifld) * ni(ix,iy,jfld) * w2(ix,iy)
                 t1 = ni(ix2,iy,ifld)*ni(ix2,iy,jfld)*w2(ix2,iy)
                 awll = awoll*loglambda(ix,iy)
                 tv  = awll*(t0+t1)/2
                 resmo(ix,iy,ifld) = resmo(ix,iy,ifld) +
     .                        tv * (up(ix,iy,jfld)-up(ix,iy,ifld))
                 fricnrl(ix,iy,ifld) = fricnrl(ix,iy,ifld) +
     .                tv*(up(ix,iy,jfld)-up(ix,iy,ifld))/vol(ix,iy)
                 w0(ix,iy) = w0(ix,iy) - tv
               enddo
             enddo
           endif
         enddo
        endif   # if test on isofric.eq.1

 105  continue
if (TimingPandfOn.gt.0) TotTimeMomBal=TotTimeMomBal+tock(TimeMomBal)
if (TimingPandfOn.gt.0) TimeEngBal=tick()
*****************************************************************
*****************************************************************
*  Here starts the old ENEBAL
*****************************************************************
*  ---------------------------------------------------------------------
*  compute temperature conductances.
*  ---------------------------------------------------------------------
*  -- initialize to 0 --

      do 708 iy = j1, j6
         do 707 ix = i1, i6
            floxe(ix,iy) = 0.0e0
            floxi(ix,iy) = 0.0e0
            floye(ix,iy) = 0.0e0
            floyi(ix,iy) = 0.0e0
            feiycbo(ix) = 0.0e0
            feeycbo(ix) = 0.0e0
            w0(ix,iy) = 0.0e0
            w1(ix,iy) = 0.0e0
            wvh(ix,iy,1) = 0.0e0
	    wvh(ix,iy,2) = 0.0e0
  707    continue
  708 continue

*  -- compute conxe and conxi --

*     (The computation of conxe involves a flux limit)

      do 121 iy = j4, j8
         do 120 ix = i1, i5
            ix2 = ixp1(ix,iy)
            t0 = max (te(ix,iy), temin*ev)
            t1 = max (te(ix2,iy), temin*ev)
            vt0 = sqrt(t0/me)
            vt1 = sqrt(t1/me)
            wallfac = 1.
            do jx = 1, nxpt
               if ( ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .                (ix==ixrb(jx).and.ixmxbcl==1) )
     .              .and. (isplflxl==0) ) wallfac = flalfepl/flalfe
            enddo
            qfl = wallfac*flalfe * sx(ix,iy) * rrv(ix,iy) *
     .          (ne(ix,iy)*vt0*t0 + ne(ix2,iy)*vt1*t1) / 2
            csh = sx(ix,iy) * hcxe(ix,iy) * gxf(ix,iy)
	    lxtec = 0.5*(te(ix,iy)+te(ix2,iy)) /
     .             ( abs(te(ix,iy)-te(ix2,iy))*gxf(ix,iy) + 100.*cutlo )
	    qsh = csh * (te(ix,iy)-te(ix2,iy)) * (1. + lxtec/lxtemax)
            qr = (1-isflxlde)*abs(qsh/qfl)
            conxe(ix,iy) = (1-isflxlde)*csh / (1 + qr)**2 +
     .                isflxlde*csh /(1 + abs(qsh/qfl)**flgam)**(1/flgam)
            floxe(ix,iy) = floxe(ix,iy) + (sign(qr*qr,qsh)/(1 + qr)**2)
     .                  *flalfea(ix) * sx(ix,iy) * ( ne(ix,iy)*
     .                   rr(ix,iy)*vt0 + ne(ix2,iy)*rr(ix2,iy)*vt1 ) / 2
c.... Now do the ions (hcxi is flux-limited previously when it is built)
          if (isflxldi .ne. 2) then    # Else flux limit done on hcxij
            t0 = max (ti(ix,iy), temin*ev)
            t1 = max (ti(ix2,iy), temin*ev)
            vt0 = sqrt(t0/mi(1))
            vt1 = sqrt(t1/mi(1))
            wallfac = 1.
            do jx = 1, nxpt
               if ( ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .                (ix==ixrb(jx).and.ixmxbcl==1) )
     .              .and. (isplflxl==0) ) wallfac = flalfipl/flalfi
            enddo
            qfl = wallfac*flalfia(ix) * sx(ix,iy) * rrv(ix,iy) *
     .                    (ne(ix,iy)*vt0*t0 + ne(ix2,iy)*vt1*t1) / 2
            csh = sx(ix,iy) * hcxi(ix,iy) * gxf(ix,iy)
            lxtic = 0.5*(ti(ix,iy)+ti(ix2,iy)) /
     .             ( abs(ti(ix,iy)-ti(ix2,iy))*gxf(ix,iy) + 100.*cutlo)
	    qsh = csh * (ti(ix,iy)-ti(ix2,iy)) * (1. + lxtic/lxtimax)
            qipar(ix,iy,1) = qsh/(rrv(ix,iy)*sx(ix,iy))
            qr = (1-isflxldi)*abs(qsh/qfl)
            conxi(ix,iy) = (1-isflxldi)*csh / (1 + qr)**2 +
     .                isflxldi*csh / (1 + abs(qsh/qfl)**flgam)**(1/flgam)
            floxi(ix,iy) = floxi(ix,iy) + (sign(qr*qr,qsh)/(1 + qr)**2)
     .                  *flalfia(ix) * sx(ix,iy) *( ne(ix,iy)*rr(ix,iy)*
     .                   vt0 + ne(ix2,iy)*rr(ix2,iy)*vt1 ) / 2
          else
            conxi(ix,iy) = sx(ix,iy) * hcxi(ix,iy) * gxf(ix,iy)
          endif
  120    continue
         conxe(nx+1,iy) = 0
         conxi(nx+1,iy) = 0
  121 continue

*  -- compute conye and conyi --

      do 123 iy = j1, j5
         do 122 ix = i4, i8
            conye(ix,iy) = sy(ix,iy) * hcye(ix,iy) * gyf(ix,iy)
            conyi(ix,iy) = sy(ix,iy) * hcyi(ix,iy) * gyf(ix,iy)
  122    continue
  123 continue

      do 124 ix = i1, i6
         conye(ix,ny+1) = 0.0e0
         conyi(ix,ny+1) = 0.0e0
  124 continue

*  ---------------------------------------------------------------------
*  compute the strength of convection for te and ti.
*  ---------------------------------------------------------------------
*  includes a correction to account for e-velocity .ne. i-velocity
*  (ccn term), and also a conduction contribution to electron heat flux
*  due to friction (cthe term).  (M.E. Rensink, 07/20/90)
*  ---------------------------------------------------------------------
*     floxe, floxi  contain the cross-derivative terms now
*                        JLM      5/3/90
*  ---------------------------------------------------------------------

      do 126 iy = j4, j8
         do 125 ix = i1, i5
            ix1 = ixp1(ix,iy)
            ltmax = min( abs(te(ix,iy)/(rrv(ix,iy)*gtex(ix,iy) + cutlo)),
     .                   lcone(ix,iy) )
            lmfpe = 2e16*(te(ix,iy)/ev)**2/ne(ix,iy)
            flxlimf = flalftf*ltmax/(flalftf*ltmax + lmfpe)
            floxe(ix,iy) = floxe(ix,iy) + cfcvte*1.25*
     .                  (ne(ix,iy)+ne(ix1,iy))*vex(ix,iy)*sx(ix,iy)
     .                   - cthe*flxlimf*cfjhf*fqp(ix,iy)/ev
  125    continue
         floxe(nx+1,iy) = 0.0e0
  126 continue

c IJ 2016/10/10	add cfneutsor_ei multiplier to control fraction of neutral energy to add
      do 729 ifld = 1, nfsp
         if ((isupgon(1) .eq. 1) .and. (ifld .eq. iigsp)) then  #neutrals
            do 726 iy = j4, j8
               do 725 ix = i1, i5
                  floxi(ix,iy) = floxi(ix,iy) +
     .                           cfcvti*2.5*cfneut*cfneutsor_ei*fnix(ix,iy,ifld)
 725           continue   # next correct for incoming neut pwr = 0
               do jx = 1, nxpt  #if at plate, sub (1-cfloxiplt)*neut-contrib
                 if(ixmnbcl==1) then  #real plate-need for parallel UEDGE
                   iixt = ixlb(jx) #left plate
                   if(fnix(iixt,iy,ifld) > 0.) then
                     floxi(iixt,iy) = floxi(iixt,iy) - (1.-cfloxiplt)*
     .                 cfcvti*2.5*cfneut*cfneutsor_ei*fnix(iixt,iy,ifld)
                   endif
                 endif
                 if(ixmxbcl==1) then #real plate-need for parallel UEDGE
                   iixt = ixrb(jx) # right plate
                   if(fnix(iixt,iy,ifld) < 0.) then
                     floxi(iixt,iy) = floxi(iixt,iy) - (1.-cfloxiplt)*
     .                 cfcvti*2.5*cfneut*cfneutsor_ei*fnix(iixt,iy,ifld)
                   endif
                   floxi(ixrb(jx)+1,iy) = 0.0e0  #cosmetic
                 endif
               enddo
 726        continue
         else  #ions
            do 728 iy = j4, j8
               do 727 ix = i1, i5
                  floxi(ix,iy) = floxi(ix,iy) +
     .                           cfcvti*2.5*fnix(ix,iy,ifld)
 727           continue
               floxi(nx+1,iy) = 0.0e0
 728        continue

         endif
 729  continue

*  -- compute floye and floyi --

      do 129 iy = j1, j5    # note: cfloye usually = 2.5 or 1.5 (ExB turb)
         do 128 ix = i4, i8
            floye(ix,iy) = floye(ix,iy) + (cfloye/2.)*
     .                    (ney0(ix,iy)+ney1(ix,iy))*vey(ix,iy)*sy(ix,iy)
     .                + (vyte_use(ix,iy)+vyte_cft(ix,iy))*0.5*sy(ix,iy)*
     .                     (ney0(ix,iy)+ney1(ix,iy))
            if (iy == 0) then
               feeycbo(ix) =  cfloye*
     .                          ( ne(ix,0)*te(ix,0)*sy(ix,0) ) *
     .                         ( (1-cfeeybbo)*cfybf*veycb(ix,0) -
     .                             cfeeydbo*(1-cfydd)*veycp(ix,0) )
            endif
  128    continue
  129 continue

      do 629 ifld = 1, nfsp
         if ((isupgon(1) .eq. 1) .and. (ifld .eq. iigsp)) then
            do iy = j1, j5
               do ix = i4, i8
                  floyi(ix,iy) = floyi(ix,iy)
     .                           + cfneut * cfneutsor_ei * 2.5 * fniy(ix,iy,ifld)
               enddo
            enddo
c ...       Make correction at walls to prevent recyc neutrals injecting pwr
            do ix = i4, i8
              if (matwallo(ix) > 0 .and. recycwot(ix,1)>0.) then
                fniy_recy = max(recycwot(ix,1)*fac2sp*fniy(ix,ny,1), 0.)
                floyi(ix,ny) = floyi(ix,ny) +
     .                           cfneut*cfneutsor_ei*2.5*(1.-cfloygwall)*fniy_recy
              endif
              if (matwalli(ix) > 0 .and. recycwit(ix,1,1)>0.) then
                fniy_recy = min(recycwit(ix,1,1)*fac2sp*fniy(ix,0,1), 0.)
                floyi(ix,0) = floyi(ix,0) +
     .                          cfneut*cfneutsor_ei*2.5*(1.-cfloygwall)*fniy_recy
              endif
            enddo

         else
            do 628 iy = j1, j5 # note: cfloyi usually = 2.5 or 1.5 (ExB turb)
               do 627 ix = i4, i8
                  floyi(ix,iy) = floyi(ix,iy)
     .                            + cfloyi * fniy(ix,iy,ifld)
     .                            + (vyti_use(ix,iy)+vyti_cft(ix,iy))*
     .                                                  0.5*sy(ix,iy)*
     .                              (niy0(ix,iy,ifld)+niy1(ix,iy,ifld))
                  if (iy == 0) then
                     feiycbo(ix) = feiycbo(ix) + cfloyi*fniycbo(ix,ifld)*
     .                                                  ti(ix,0)
                  endif
 627           continue
 628        continue
         endif
 629  continue

c...  Next B x grad(T), first for the ions
      if(abs(cfbgt) .gt. 0 .or. cfeexdbo+cfeixdbo > 0.) then

      do 133 ifld = 1, nfsp
        do 132 iy = j4, j8
           do 131 ix = i1, i5
	     iy1 = max(0,iy-1)
             ix1 = ixp1(ix,iy)
             ix5 = ixp1(ix,iy1)
	     if (iy.eq.0 .or. iy.eq.ny+1) goto 131
c... sknam: grad T from tiv
             temp1 = 4.0*(tiv(ix,iy) - tiv(ix,iy1))*gyc(ix,iy)
cccMER For full double-null configuration, iysptrx is last closed flux surface.
             if ( isxptx(ix,iy)==0 .and. iysptrx.gt.0 )
     .                 temp1 = 4.0*(tiv(ix,iy) - tiv(ix,iy1))*gyc(ix,iy)

             if (zi(ifld) > 1.e-10) then
                floxibgt(ix,iy,ifld)=(5*sx(ix,iy)/(32*qe*zi(ifld) )) *
     .                           ( ni(ix,iy,ifld) + ni(ix1,iy,ifld) ) *
     .                           ( rbfbt2(ix,iy) + rbfbt2(ix1,iy) ) *
     .                            temp1
c             else
c                 floxibgt(ix,iy,ifld)=0.0
             endif
             floxi(ix,iy) = floxi(ix,iy) + cfbgt*floxibgt(ix,iy,ifld)
  131      continue
  132  continue
  133 continue

      do 136 ifld = 1, nfsp
        do 135 iy = j1, j5
           do 134 ix = i4, i8
             ix3 = ixm1(ix,iy)
             ix4 = ixm1(ix,iy+1)
             do jx = 1, nxpt
                if ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .               (ix==ixrb(jx)+1.and.ixmxbcl==1) ) goto 134
             enddo
c... sknam: grad T from tiv
             temp1 = 4.0*(tiv(ix,iy) - tiv(ix3,iy))*gxc(ix,iy)
cccMER For full double-null configuration, iysptrx is last closed flux surface.
             if ( isxpty(ix,iy)==0 .and. iysptrx.gt.0 )
     .              temp1 = 4.0*(tiv(ix,iy) - tiv(ix3,iy))*gxc(ix,iy)
             if (zi(ifld) > 1.e-10) then
                floyi(ix,iy) = floyi(ix,iy) -
     .                      cfbgt*( 5*sy(ix,iy) / (32*qe*zi(ifld) )) *
     .                          ( ni(ix,iy,ifld) + ni(ix,iy+1,ifld) ) *
     .                          ( rbfbt2(ix,iy) + rbfbt2(ix,iy+1) ) *
     .                           temp1
             endif
 134      continue
 135   continue
 136  continue

c...  Now B x grad(T) for the electrons

      do 138 iy = j4, j8
         do 137 ix = i1, i5
	     iy1 = max(0,iy-1)
             ix1 = ixp1(ix,iy)
             ix5 = ixp1(ix,iy1)
	     if (iy.eq.0 .or. iy.eq.ny+1) goto 137
c             temp1 = (gtey(ix,iy ) + gtey(ix1,iy ) +
c     .                gtey(ix,iy1) + gtey(ix5,iy1))
c... sknam: grad T from tev
             temp1 = 4.0*(tev(ix,iy) - tev(ix,iy1))*gyc(ix,iy)
cccMER For full double-null configuration, iysptrx is last closed flux surface.
             if ( isxptx(ix,iy)==0 .and. iysptrx.gt.0 )
     .                    temp1 = 4.0*(tev(ix,iy) - tev(ix,iy1))*gyc(ix,iy)
             floxebgt(ix,iy) = ( 5*sx(ix,iy) / (32*qe) ) *
     .                           ( ne(ix,iy) + ne(ix1,iy) ) *
     .                          ( rbfbt2(ix,iy) + rbfbt2(ix1,iy) ) *
     .                           temp1
             floxe(ix,iy) = floxe(ix,iy) - cfbgt*floxebgt(ix,iy)
  137    continue
  138 continue

      do 140 iy = j1, j5
	 do 139 ix = i4, i8
	    ix3 = ixm1(ix,iy)
	    ix4 = ixm1(ix,iy+1)
            do jx = 1, nxpt
               if ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .              (ix==ixrb(jx)+1.and.ixmxbcl==1) ) goto 139
            enddo
c... sknam: grad T from tev
            temp1 = 4.0*(tev(ix,iy) - tev(ix3,iy))*gxc(ix,iy)
cccMER For full double-null configuration, iysptrx is last closed flux surface.
            if ( isxpty(ix,iy)==0 .and. iysptrx.gt.0 )
     .               temp1 = 4.0*(tev(ix,iy) - tev(ix3,iy))*gxc(ix,iy)
            floye(ix,iy) = floye(ix,iy) +
     .                       cfbgt*( 5*sy(ix,iy) / (32*qe) ) *
     .                         ( ne(ix,iy) + ne(ix,iy+1) ) *
     .                         ( rbfbt2(ix,iy) + rbfbt2(ix,iy+1) ) *
     .                          temp1

  139    continue
  140 continue

      endif

c...Add the charge-exhange neutral contributions to ion+neutral temp eq.


         do 141 iy = j4, j8
            do 142 ix = i1, i5
               floxi(ix,iy) = floxi(ix,iy) +
     .                    cfneut*cfneutsor_ei*cngtgx(1)*cfcvti*2.5*fngx(ix,iy,1)
 142        continue
            floxi(nx+1,iy) = 0.0e0
 141        continue

*  --Adds to floyi --

         do 145 iy = j1, j5
            do 144 ix = i4, i8
               floyi(ix,iy) = floyi(ix,iy)
     .                       + cfneut*cfneutsor_ei*cngtgy(1)*2.5*fngy(ix,iy,1)
 144        continue
 145     continue

*  ---------------------------------------------------------------------
*  compute the electron and the ion energy flow.
*  ---------------------------------------------------------------------

*  -- compute the electron energy flow --
      if(isteon .eq. 1) call fd2tra (nx,ny,floxe,floye,conxe,conye,
     .                               te,feex,feey,0,methe)

*  -- compute the ion thermal energy flow --
      if(istion .eq. 1) call fd2tra (nx,ny,floxi,floyi,conxi,conyi,
     .                               ti,feix,feiy,0,methi)

c  -- Add rad flux of 4th order diff operator; damp grid-scale oscillations
      if (abs(kye4order)>1.e-50 .or. abs(kyi4order)>1.e-50) then
        do iy = j2p, j5m   #range to iy=1:ny-1 for feey4ord,feiy4ord
          iym1 = max(iy-1,0)
          iyp1 = min(iy+1,ny+1)
          iyp2 = min(iy+2,ny+1)
          do ix = i4, i8
            dtdym1 = (te(ix,iy)-te(ix,iym1))*gyf(ix,iym1)
            dtdy0 = (te(ix,iyp1)-te(ix,iy))*gyf(ix,iy)
            dtdyp1 = (te(ix,iyp2)-te(ix,iyp1))*gyf(ix,iyp1)
            d2tdy20 = (dtdy0 - dtdym1)*gy(ix,iy)
            d2tdy2p1 = (dtdyp1 - dtdy0)*gy(ix,iyp1)
            d3tdy3 = (d2tdy2p1 - d2tdy20)*gyf(ix,iy)
            feey4ord(ix,iy) = kye4order*d3tdy3*ney1(ix,iy)*
     .                                     sy(ix,iy)/gyf(ix,iy)**2
            feey(ix,iy) = feey(ix,iy) + feey4ord(ix,iy)
          enddo
          do ix = i4, i8
            dtdym1 = (ti(ix,iy)-ti(ix,iym1))*gyf(ix,iym1)
            dtdy0 = (ti(ix,iyp1)-ti(ix,iy))*gyf(ix,iy)
            dtdyp1 = (ti(ix,iyp2)-ti(ix,iyp1))*gyf(ix,iyp1)
            d2tdy20 = (dtdy0 - dtdym1)*gy(ix,iy)
            d2tdy2p1 = (dtdyp1 - dtdy0)*gy(ix,iyp1)
            d3tdy3 = (d2tdy2p1 - d2tdy20)*gyf(ix,iy)
            feiy4ord(ix,iy) = kyi4order*d3tdy3*niy1(ix,iy,1)*
     .                                     sy(ix,iy)/gyf(ix,iy)**2
            feiy(ix,iy) = feiy(ix,iy) + feiy4ord(ix,iy)
          enddo
        enddo
      endif

*  ---------------------------------------------------------------------
*  compute the energy residuals.
*  ---------------------------------------------------------------------

*  -- source terms --

      do 150 iy = j2, j5
         do 149 ix = i2, i5
            resee(ix,iy) =
     .             seec(ix,iy) + seev(ix,iy) * te(ix,iy)
     .           + pwrsore(ix,iy)
     .           + cmneut * cmneutsor_ee * uesor_te(ix,iy)
     .           - nuvl(ix,iy,1)*vol(ix,iy)*bcee*ne(ix,iy)*te(ix,iy)
            resei(ix,iy) =
     .             seic(ix,iy) + seiv(ix,iy) * ti(ix,iy)
     .           + pwrsori(ix,iy)
     .           + cmneut * cmneutsor_ei * uesor_ti(ix,iy)
     .           - nuvl(ix,iy,1)*vol(ix,iy)*bcei*ne(ix,iy)*ti(ix,iy)
  149    continue
  150 continue

*  -- divergence of electron and ion energy flows --

c...  Add y-component of nonorthogonal diffusive flux; convective component
c...  already added to uu(ix,iy)
      if (isnonog .eq. 1) then

         do 3094 iy = j1, j6
            if (iy .gt. ny) goto 3094
            iy1 = max(iy-1,0)
            do 3093 ix = i1, i6
c...  First do the Te equation
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               ix3 = ixm1(ix,iy1)
               ix4 = ixp1(ix,iy1)
               ix5 = ixm1(ix,iy+1)
               ix6 = ixp1(ix,iy+1)
                 grdnv=(    ( fym (ix,iy,1)*log(te(ix2,iy1 )) +
     .                        fy0 (ix,iy,1)*log(te(ix2,iy  )) +
     .                        fyp (ix,iy,1)*log(te(ix2,iy+1)) +
     .                        fymx(ix,iy,1)*log(te(ix ,iy1 )) +
     .                        fypx(ix,iy,1)*log(te(ix ,iy+1)) )
     .                     -( fym (ix,iy,0)*log(te(ix ,iy1 )) +
     .                        fy0 (ix,iy,0)*log(te(ix ,iy  )) +
     .                        fyp (ix,iy,0)*log(te(ix ,iy+1)) +
     .                        fymx(ix,iy,0)*log(te(ix4,iy1 )) +
     .                        fypx(ix,iy,0)*log(te(ix6,iy+1)) ) ) *
     .                                                   gxfn(ix,iy)
               feexy(ix,iy) = exp( 0.5*
     .                         (log(te(ix2,iy)) + log(te(ix,iy))) )*
     .                               (fcdif*kye+kye_use(ix,iy))*0.5*
     .                                       (ne(ix2,iy)+ne(ix,iy))*
     .                                     (grdnv/cos(angfx(ix,iy)) -
     .                         (log(te(ix2,iy)) - log(te(ix,iy)))*
     .                                         gxf(ix,iy))*sx(ix,iy)

c...  Now do the Ti equation.
c --- If we are using the parallel neutral momentum equation, we automatically
c --- change to a combined neutral+ion energy equation. We thus need to
c --- include the neutral heat conductivity. Since it is is isotropic
c --- we could use hcxn though we take the radial derivative; but this is
c --- only true if we dont flux limit.  Thus, we use 4-pt average of hcyn.
c --- Note: this four-point average results in not getting the full Jac. for
c --- a nonorthogonal mesh because of niy1,0 - see def. of hcyn

                 grdnv =(    ( fym (ix,iy,1)*log(ti(ix2,iy1 )) +
     .                         fy0 (ix,iy,1)*log(ti(ix2,iy  )) +
     .                         fyp (ix,iy,1)*log(ti(ix2,iy+1)) +
     .                         fymx(ix,iy,1)*log(ti(ix ,iy1 )) +
     .                         fypx(ix,iy,1)*log(ti(ix ,iy+1)) )
     .                      -( fym (ix,iy,0)*log(ti(ix ,iy1 )) +
     .                         fy0 (ix,iy,0)*log(ti(ix ,iy  )) +
     .                         fyp (ix,iy,0)*log(ti(ix ,iy+1)) +
     .                         fymx(ix,iy,0)*log(ti(ix4,iy1 )) +
     .                         fypx(ix,iy,0)*log(ti(ix6,iy+1)) ) ) *
     .                                                   gxfn(ix,iy)
               feixy(ix,iy) = exp( 0.5*
     .                       (log(ti(ix2,iy)) + log(ti(ix,iy))) )*
     .                           ( (fcdif*kyi+kyi_use(ix,iy))*0.5*
     .                                     (nit(ix2,iy)+nit(ix,iy))
     .                   + cfneut*cfneutsor_ei*0.25*(hcyn(ix ,iy)+hcyn(ix ,iy1)
     .                              +hcyn(ix2,iy)+hcyn(ix4,iy1)) ) *
     .                                 (  grdnv/cos(angfx(ix,iy))
     .                         - (log(ti(ix2,iy)) - log(ti(ix,iy)))*
     .                                        gxf(ix,iy) )*sx(ix,iy)
c...  Flux limit with flalftxt even though hcys have parallel FL built in
               t0 = max(ti(ix,iy),temin*ev)
               t1 = max(ti(ix2,iy),temin*ev)
               vttn = t0*sqrt( t0/mi(1) )
               vttp = t1*sqrt( t1/mi(1) )
               qfl = flalftxy * 0.125 * sx(ix,iy) * (vttn + vttp) *
     .               (ni(ix,iy,1)+ng(ix,iy,1)+ni(ix2,iy,1)+ng(ix2,iy,1))
               feixy(ix,iy) = feixy(ix,iy) /
     .                              sqrt(1. + (feixy(ix,iy)/qfl)**2)

 3093      continue
 3094    continue

c...  Fix the fluxes with the same indice range as in fd2tra
         do iy = j4, j8
            do ix = i1, i5
               feex(ix,iy) = feex(ix,iy) - feexy(ix,iy)
               feix(ix,iy) = feix(ix,iy) - feixy(ix,iy)
            enddo
         enddo

      endif
c...  Finished with nonorthogonal mesh part

c ... Demand that net feex cannot be out of the plates
      if (isfeexpl0 == 1) then
        do iy = j4, j8
          do jx = 1, nxpt
            if(feex(ixlb(jx),iy) > 0. .and. ixmnbcl==1) then
              feexflr = ni(ixlb(jx),iy,1)*1.e4*ev*sx(ixlb(jx),iy)
              feex(ixlb(jx),iy) = feex(ixlb(jx),iy)/
     .                (1.+ (feex(ixlb(jx),iy)/feexflr)**4)
            endif
            if(feex(ixrb(jx),iy) < 0. .and. ixmxbcl==1) then
              feexflr = ni(ixrb(jx),iy,1)*1.e4*ev*sx(ixrb(jx),iy)
              feex(ixrb(jx),iy) = feex(ixrb(jx),iy)/
     .                (1.+ (feex(ixrb(jx),iy)/feexflr)**4)
            endif
          enddo
        enddo
      endif

      if (isfeixpl0 == 1) then
        do iy = j4, j8
          do jx = 1, nxpt
            if(feix(ixlb(jx),iy) > 0.) then
              feixflr = ni(ixlb(jx),iy,1)*1.e3*ev*sx(ixlb(jx),iy)
              feix(ixlb(jx),iy) = feix(ixlb(jx),iy)/
     .                (1.+ (feix(ixlb(jx),iy)/feixflr)**4)
            endif
            if(feix(ixrb(jx),iy) < 0.) then
              feixflr = ni(ixrb(jx),iy,1)*1.e3*ev*sx(ixrb(jx),iy)
              feix(ixrb(jx),iy) = feix(ixrb(jx),iy)/
     .                (1.+ (feix(ixrb(jx),iy)/feixflr)**4)
            endif
          enddo
        enddo
      endif

      do 310 iy = j2, j5
	    if((isudsym==1.or.geometry.eq.'dnXtarget') .and. nxc > 1) then
            feex(nxc-1,iy) = 0.
            feix(nxc-1,iy) = 0.
            feex(nxc,iy) = 0.
            feix(nxc,iy) = 0.
            feex(nxc+1,iy) = 0.
            feix(nxc+1,iy) = 0.
         endif
         if (islimon .ne. 0 .and. iy .ge. iy_lims) then
            feex(ix_lim,iy) = 0.
            feix(ix_lim,iy) = 0.
         endif
         if (nxpt==2 .and. ixmxbcl==1) then
            feex(ixrb(1)+1,iy) = 0.
            feix(ixrb(1)+1,iy) = 0.
         endif
         do 309 ix = i2, i5
            ix1 = ixm1(ix,iy)
            resee(ix,iy) = resee(ix,iy)
     .                  - ( feex(ix,iy) - feex(ix1,iy)
     .          + fluxfacy*(feey(ix,iy) - feey(ix,iy-1)) )
c ... ## IJ 2017 cfneutsor_ei flags above control neutral contrib.
            resei(ix,iy) = resei(ix,iy)
     .                  - ( feix(ix,iy) - feix(ix1,iy)
     .          + fluxfacy*(feiy(ix,iy) - feiy(ix,iy-1)) )

c ... ## IJ 2016/10/19 add MC neutral flux
           if(get_neutral_moments .and. cmneutdiv_feg .ne. 0.0) then
              jfld=1
              seg_ue(ix,iy,jfld)=-( (fegx_ue(ix,iy,jfld)-fegx_ue(ix1,iy,  jfld))
     .                   + fluxfacy*(fegy_ue(ix,iy,jfld)-fegy_ue(ix, iy-1,jfld)) )
     .                  *( (ni(ix,iy,jfld)*ti(ix,iy))/(ni(ix,iy,jfld)*ti(ix,iy)) )
              resei(ix,iy) = resei(ix,iy) +
     .                              cmneutdiv*cmneutdiv_feg*seg_ue(ix,iy,jfld)
            endif
  309    continue
  310 continue


*  -- total energy residual and equipartition --

c...  Electron radiation loss -- ionization and recombination
            do 316 iy = iys1, iyf6  #iys,iyf
               do 315 ix = ixs1, ixf6  #iys, iyf
                  erlizold = erliz(ix,iy)
                  erlrcold = erlrc(ix,iy)
                  ne_sgvi = ne(ix,iy)
                  if (ifxnsgi.eq.1) ne_sgvi = cne_sgvi  # fix density dependence
                  if (istabon==16) then      # compute from b2frates data file
                     zmax=1
                     znuc=1
                     denz(0)=ng(ix,iy,1)     # use ngbackg as below ?
                     denz(1)=ni(ix,iy,1)     # use fac2sp  as below ?
                     dene=ne_sgvi
                     rdum=radmc(zmax,znuc,te(ix,iy),dene,denz,radz)
                     erliz(ix,iy)=chradi*radz(0)*vol(ix,iy)
                     if (isrecmon .ne. 0) erlrc(ix,iy)=chradr*radz(1)*vol(ix,iy)
                  else                       # compute from other data files
                     erliz(ix,iy) = chradi *
     .                           erl1(te(ix,iy),ne_sgvi,rtau(ix,iy))
     .                                  * (ng(ix,iy,1)-ngbackg(1)*
     .                    (0.9+0.1*(ngbackg(1)/ng(ix,iy,1))**ingb) ) *
     .                                                       vol(ix,iy)
                     if (isrecmon .ne. 0) erlrc(ix,iy) = chradr *
     .                               erl2(te(ix,iy),ne_sgvi,rtau(ix,iy))
     .                             * fac2sp*ni(ix,iy,1) * vol(ix,iy)
                  endif
                  eeliold = eeli(ix,iy)
                  if (icnuiz.le.1 .and. psor(ix,iy,1).ne.0.)
     .                                           eeli(ix,iy) = 13.6*ev +
     .                               erliz(ix,iy)/(fac2sp*psor(ix,iy,1))

                  pradhyd(ix,iy)= ( (eeli(ix,iy)-ebind*ev)*psor(ix,iy,1)+
     .                                         erlrc(ix,iy) )/vol(ix,iy)
 315           continue
 316        continue

      do iy = iys1, iyf6  #j2, j5
        do ix = ixs1, ixf6  #i2, i5
          vsoreec(ix,iy) =
     .          - cfneut*cfneutsor_ee*cnsor*13.6*ev*fac2sp*psorc(ix,iy,1)
     .          + cfneut*cfneutsor_ee*cnsor*13.6*ev*fac2sp*psorrgc(ix,iy,1)
     .          - cfneut*cfneutsor_ee*cnsor*erliz(ix,iy)
     .          - cfneut*cfneutsor_ee*cnsor*erlrc(ix,iy)
     .          - cfneut*cfneutsor_ee*cnsor*ediss*ev*(0.5*psordis(ix,iy))
        enddo
      enddo

ccc         if (ishosor.eq.1) then  #full RHS eval
ccc
ccc           if (svrpkg.eq."cvode") then    # cannot access yl(neq+1)
ccc            call xerrab('*** svrpkg=cvode not allowed for ishosor=1 **')
ccc           endif
ccc
ccc           if (yl(neq+1).lt.0) then  #full RHS eval
ccc
cccc ...    integ. source over cells (but not for Jac) for higher-order accuracy
ccc
ccc             call volave(nx, ny, j2, j5, i2, i5, ixp1(0,0), ixm1(0,0),
ccc     .                         fsprd, vol(0,0), psor_tmpov(0,0), vsoree)
ccc
ccc           endif   # end of if (yl(neq+1).lt.0) test
ccc         endif    # end of integrating over sources and ishosor test

c*************************************************************
c   Perform 5pt average of source terms as volume integral
c*************************************************************
         if (iseesorave == 0.) then  #use only single-cell value
           do iy = iys1, iyf6
             do ix = ixs1, ixf6
               vsoree(ix,iy) = vsoreec(ix,iy)
             enddo
           enddo

         elseif (iseesorave > 0.)

            if (xc < 0) then  #full RHS eval
              j2pwr = j2
              j5pwr = j5
            else  # Jacobian
              j2pwr = max(1, yc-1)
              j5pwr = min(ny, yc+1)
            endif
            do iy = j2pwr, j5pwr
              if (xc < 0) then #full RHS eval
                i2pwr = i2
                i5pwr = i5
              else  #Jacobian eval
                i2pwr = max(1,ixm1(xc,yc))
                i5pwr = min(nx, ixp1(xc,yc))
              endif
              do ix = i2pwr, i5pwr
                ix1 = ixm1(ix,iy)
                ix2 = ixp1(ix,iy)
                vsoree(ix,iy) = (1.-iseesorave*0.5)*
     .                                  vsoreec(ix,iy)+
     .                               0.125*iseesorave*vol(ix,iy)*
     .                           ( vsoreec(ix,iy-1)/vol(ix,iy-1) +
     .                             vsoreec(ix,iy+1)/vol(ix,iy+1) +
     .                             vsoreec(ix1,iy)/vol(ix1,iy)   +
     .                             vsoreec(ix2,iy)/vol(ix2,iy) )
              enddo
            enddo

         endif    # end of integrating over sources and iseesorave test


      do 152 iy = j2, j5
         do 151 ix = i2, i5
            ix1 = ixm1(ix,iy)
            w0(ix,iy) = vol(ix,iy) * eqp(ix,iy) * (te(ix,iy)-ti(ix,iy))
            resee(ix,iy) = resee(ix,iy) - w0(ix,iy) + vsoree(ix,iy)
            if (isupgon(1).eq.1) then
c These terms include electron-ion equipartition as well as terms due
c to the friction force between neutrals and ions
               t1 = 0.5*(up(ix,iy,1)+up(ix1,iy,1))
               t2 = 0.5*(up(ix,iy,iigsp)+up(ix1,iy,iigsp))
               resei(ix,iy) = resei(ix,iy) + w0(ix,iy)
     .             + cfneut * cfneutsor_ei * cfnidh * 0.5*mi(1) * (t1-t2)*(t1-t2) *
     .                    (  psor(ix,iy,1) + psorrg(ix,iy,1)
     .              + 2*cfticx*nucx(ix,iy,1)*ng(ix,iy,1)*vol(ix,iy)  )
     .             + cfneut * cfneutsor_ei * cnsor*eion*ev*psordis(ix,iy)
            else
               resei(ix,iy) = resei(ix,iy) + w0(ix,iy)
     .             + cfneut * cfneutsor_ei * ctsor*1.25e-1*mi(1)*
     .                    (upi(ix,iy,1)+upi(ix1,iy,1))**2*
     .                    fac2sp*psor(ix,iy,1)
     .             + cfneut * cfneutsor_ei * ceisor*cnsor* eion*ev*psordis(ix,iy)
     .             - cfneut * cfneutsor_ei * ccoldsor*ng(ix,iy,1)*nucx(ix,iy,1)*
     .                    (  1.5*ti(ix,iy)
     .                     - 0.125*mi(1)*(upi(ix,iy,1)+upi(ix1,iy,1))**2
     .                     - eion*ev  ) * vol(ix,iy)
            endif
  151    continue
  152 continue


c ... If molecules are present as gas species 2, add ion/atom cooling
      if(ishymol == 1) then
        do iy = j2, j5
          do ix = i2, i5
            resei(ix,iy) = resei(ix,iy) - vol(ix,iy)*eqpg(ix,iy,2)*
     .                                     (ti(ix,iy)-tg(ix,iy,2))
          enddo
        enddo
      endif

*  -- Energy transfer to impurity neutrals at tg(,,igsp)
      if (ngsp >= 2) then   # for now, specialized to igsp=2 only
        do ifld = nhsp+1, nisp
          do iy = j2, j5    # iys,iyf limits dont seem to work(?)
            do ix = i2, i5
              resei(ix,iy) =resei(ix,iy) -cftiimpg*1.5*ni(ix,iy,ifld)*
     .                      (nucxi(ix,iy,ifld)+nueli(ix,iy,ifld))*
     .                      (ti(ix,iy) - tg(ix,iy,2))*vol(ix,iy)
            enddo
          enddo
        enddo
      endif

*  -- impurity radiation --

      if (isimpon .ge. 2) then
         if (istimingon .eq. 1) tsimp = gettime(sec4)

         do 533 iy = iys1, iyf6  #if Jacobian, only 1 cell done - local sor
            do 532 ix = ixs1, ixf6
               ntau(ix,iy) = atau(ix,iy) * ne(ix,iy)
               nratio(ix,iy) = ng(ix,iy,1)/ne(ix,iy)
               pradold = pwrzec(ix,iy)
               if (isimpon .eq. 2) then   # fixed-fraction model
                  na(ix,iy) = afrac(ix,iy) * ne(ix,iy)
                  pradcff(ix,iy) = na(ix,iy) * ne(ix,iy) *
     .               emissbs (te(ix,iy), nratio(ix,iy), ntau(ix,iy))
                  pradc(ix,iy) = pradcff(ix,iy)
                  pwrzec(ix,iy) = pradc(ix,iy)
               elseif (isimpon .eq. 3) then   # average-ion model
                  na(ix,iy) = ni(ix,iy,nhsp+1)
                  pradc(ix,iy) = na(ix,iy) * ne(ix,iy) *
     .               radneq (te(ix,iy), nratio(ix,iy))
                  pwrzec(ix,iy) = pradc(ix,iy)
               elseif (nzspt .eq. 0) then
                  # Multicharge group not allocated, so avoid radimpmc
                  pradc(ix,iy) = 0.
                  pwrzec(ix,iy) = pradc(ix,iy)
               elseif (isimpon .ge. 4) then  # multi-charge model
                  pradc(ix,iy) = 0.
                  pwrzec(ix,iy) = 0.
                  nsm1 = nhsp
                  do igsp = nhgsp+1, ngsp  # loop over diff imp species
                     jz = igsp - nhgsp
                     zn = znucl(nsm1+nzsp(jz))
                     if (ngsp .ge. nhgsp) nzloc(0) = ng(ix,iy,igsp)
                     do ifld = 1, nzsp(jz)
                        nzloc(ifld) = ni(ix,iy,nsm1+ifld)
                     enddo
                     nsm1 = nsm1 + nzsp(jz)   # setup for next igsp
                     argth = (te(ix,iy)-1.*ev)/(del_te_ro*ev)
                     fac_rad = 1.
                     if(del_te_ro.lt. 100.) fac_rad=0.5*(1+tanh(argth))
                     if (ismctab .eq. 1) then
                        pwrzec(ix,iy)= pwrzec(ix,iy) + fac_rad*
     .                                radimpmc (nzsp(jz), te(ix,iy),
     .                                    ne(ix,iy), nzloc, impradloc)
                     elseif (ismctab .eq. 2) then
                        pwrzec(ix,iy)= pwrzec(ix,iy) + fac_rad*
     .                                radmc(nzsp(jz), zn, te(ix,iy),
     .                                   ne(ix,iy), nzloc, impradloc)
                     endif

                     do ifld = 0, nzsp(jz)
                        pradzc(ix,iy,ifld,jz) = impradloc(ifld)
                        pradc(ix,iy) = pradc(ix,iy)+impradloc(ifld)
                     enddo

                  enddo

		  if (isimpon .eq. 7) then  # add fixed-fraction contrib
                     na(ix,iy) = afrac(ix,iy) * ne(ix,iy)
                     pradcff(ix,iy) = na(ix,iy)* ne(ix,iy)*
     .                     emissbs(te(ix,iy), nratio(ix,iy), ntau(ix,iy))
                     pradc(ix,iy) = pradc(ix,iy) + pradcff(ix,iy)
                     pwrzec(ix,iy) = pwrzec(ix,iy) + pradcff(ix,iy)
                  endif

               endif
 532        continue
 533     continue

c*************************************************************
c   Perform 5pt average of source terms as volume integral
c*************************************************************
cc         if (ishosor.eq.0) then  #use only single-cell value
         if (iseesorave.eq.0.) then  #use only single-cell value
           do iy = iys1, iyf6
             do ix = ixs1, ixf6
               pwrze(ix,iy) = pwrzec(ix,iy)
               prad(ix,iy) = pradc(ix,iy)
               do igsp = nhgsp+1, ngsp
                 jz = igsp - nhgsp
                 do ifld = 0, nzsp(jz)
                   pradz(ix,iy,ifld,jz) = pradzc(ix,iy,ifld,jz)
                 enddo
               enddo
             enddo
           enddo

cc         elseif (ishosor .ne. 0)
         elseif (iseesorave > 0.)

            if (xc < 0) then  #full RHS eval
              j2pwr = j2
              j5pwr = j5
            else  # Jacobian
              j2pwr = max(1, yc-1)
              j5pwr = min(ny, yc+1)
            endif
            do iy = j2pwr, j5pwr
              if (xc < 0) then #full RHS eval
                i2pwr = i2
                i5pwr = i5
              else  #Jacobian eval
                i2pwr = max(1,ixm1(xc,yc))
                i5pwr = min(nx, ixp1(xc,yc))
              endif
              do ix = i2pwr, i5pwr
                ix1 = ixm1(ix,iy)
                ix2 = ixp1(ix,iy)
                pwrze(ix,iy) = (1.-iseesorave*0.5)*pwrzec(ix,iy) +
     .                                            0.125*iseesorave*
     .                         ( pwrzec(ix,iy-1)+ pwrzec(ix,iy+1)+
     .                           pwrzec(ix1,iy) + pwrzec(ix2,iy) )
                if (isimpon < 4) prad(ix,iy) = pwrze(ix,iy)
                if (isimpon >= 4) then  #prad, pradz only diagnostic
                  prad(ix,iy) = (1.-iseesorave*0.5)*pradc(ix,iy) +
     .                                         0.125*iseesorave*
     .                          ( pradc(ix,iy-1)+ pradc(ix,iy+1)+
     .                            pradc(ix1,iy) + pradc(ix2,iy) )
                  do igsp = nhgsp+1, ngsp
                    jz = igsp - nhgsp
                    do ifld = 0, nzsp(jz)
                      pradz(ix,iy,ifld,jz) = (1.-iseesorave*0.5)*
     .                                        pradzc(ix,iy,ifld,jz) +
     .                                            0.125*iseesorave*
     .              ( pradzc(ix,iy-1,ifld,jz)+ pradzc(ix,iy+1,ifld,jz)+
     .                pradzc(ix1,iy,ifld,jz) + pradzc(ix2,iy,ifld,jz) )
                    enddo
                  enddo
                endif
              enddo
            enddo

         endif    # end of integrating over sources and iseesorave test

c*******************************************************************
c ... Define a background elec energy source to prevent very low Te
c******************************************************************
      do iy = iys, iyf  #j2, j5
        do ix = ixs, ixf  #i2, i5
          pwrebkgold = pwrebkg(ix,iy)
          if (isimpon == 0) then
            pwrebkg(ix,iy) = (tebg*ev/te(ix,iy))**iteb*pwrbkg_c
          else  #add impurity rad loss
            pwrebkg(ix,iy) = (tebg*ev/te(ix,iy))**iteb*pwrbkg_c
          endif
        enddo
      enddo

c******************************************************************
c...  Update resee over whole "box" because initially set to zero
c******************************************************************
         do 536 iy = j2, j5
            do 535 ix = i2, i5
               resee(ix,iy) = resee(ix,iy) -
     .                            cnimp*pwrze(ix,iy)*vol(ix,iy) +
     .                                pwrebkg(ix,iy)*vol(ix,iy)
 535        continue
 536     continue

         if (istimingon .eq. 1) call timimpfj (tsimp, xc)
      endif  #loop for isimpon==2

*  -- joule heating --

      if (jhswitch > 0) then  # relies on div(J)=0, so omit iy=1 & ny
         if (isnewpot .eq. 1) then
            iy_min = 2
            iy_max = ny-1
         else
            iy_min = 1
            iy_max = ny
         endif
         if (jhswitch == 1) then   # div(J)=0 gives -grad(phi).J=-div(phi.J)
           do iy = max(iy_min, j2), min(iy_max, j5)
             do ix = i2, i5
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               wjdote(ix,iy) =
     .                          - 0.5*(fqp(ix,iy)+fq2(ix,iy))*
     .                                (phi(ix2,iy)+phi(ix,iy))
     .                          + 0.5*(fqp(ix1,iy)+fq2(ix1,iy))*
     .                                (phi(ix,iy)+phi(ix1,iy))
     .                          - 0.5*fqygp(ix,iy)*
     .                                (phi(ix,iy+1)+phi(ix,iy))
     .                          + 0.5*fqygp(ix,iy-1)*
     .                                (phi(ix,iy)+phi(ix,iy-1))
               resee(ix,iy) = resee(ix,iy) + wjdote(ix,iy) / ( 1. +
     .                             cfwjdotelim*(tebg*ev/te(ix,iy))**iteb )
             enddo
           enddo
         else  # for jhswitch > 1
           do iy = max(iy_min, j2), min(iy_max, j5)
             do ix = i2, i5    # use ex*fqx since phi(0,) may be large
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               wjdote(ix,iy) =
     .                       0.5*( ex(ix1,iy) *fqx(ix1,iy) +
     .                             ex(ix, iy) *fqx(ix, iy) )/gx(ix,iy)
     .                     + 0.5*( ey(ix, iy) *fqy(ix, iy) +
     .                             ey(ix,iy-1)*fqy(ix,iy-1) )/gy(ix,iy)
               resee(ix,iy) = resee(ix,iy) + wjdote(ix,iy)
             enddo
           enddo
         endif
      endif

*  -- Now we introduce the viscous heating; one-side derviatives are used
*  -- on either side of the x-point where isxpty = 0

      do 157 iy = j2, j5
         do 156 ix = i2, i5
            do 155 ifld = 1, nusp  # if nusp --> nfsp, problems from y-term
               ix1 = ixm1(ix,iy)
               ix2 = ixm1(ix,iy+1)
               ix3 = ixm1(ix,iy-1)
	       thetacc = 0.5*(angfx(ix1,iy) + angfx(ix,iy))
	       dupdx = gx(ix,iy)*(upi(ix,iy,ifld)-upi(ix1,iy,ifld))
               wvh(ix,iy,ifld) = cfvcsx(ifld)*cfvisx*cos(thetacc)*
     .                                    visx(ix,iy,ifld)*dupdx**2
               if ( isxpty(ix,iy)==0 ) then  #1-sided deriv down in y
                 dupdy = 0.5*( upi(ix,iy,  ifld)+upi(ix1,iy  ,ifld) -
     .                         upi(ix,iy-1,ifld)-upi(ix3,iy-1,ifld) )*
     .                                                    gyf(ix,iy-1)
               elseif (isxpty(ix,iy)== -1) then #1-sided up in y
                 dupdy = 0.5*( upi(ix,iy+1,ifld)+upi(ix2,iy+1,ifld) -
     .                         upi(ix,iy  ,ifld)-upi(ix1,iy  ,ifld) )*
     .                                                    gyf(ix,iy)
               elseif (isxpty(ix,iy)==1.and.isvhyha==1) then
                                 #use harm y-ave for up face-values
                                 #take abs() to avoid near-zero denomin;
                                 #small err in wvh because up then small
                 upxavep1 = 0.5*(upi(ix,iy+1,ifld)+upi(ix2,iy+1,ifld))
                 upxave0 =  0.5*(upi(ix,iy  ,ifld)+upi(ix1,iy  ,ifld))
                 upxavem1 = 0.5*(upi(ix,iy-1,ifld)+upi(ix3,iy-1,ifld))
                 upf0  = 2.*upxavep1*upxave0*(upxavep1+upxave0) /
     .                           ( (upxavep1+upxave0)**2 + upvhflr**2 )
                 upfm1 = 2.*upxave0*upxavem1*(upxave0+upxavem1) /
     .                           ( (upxave0+upxavem1)**2 + upvhflr**2 )
                 dupdy = (upf0 - upfm1)*gy(ix,iy)
               else	#V7.08.04 option - linear ave in y-direction
		 dupdy = 0.25*( (upi(ix,iy+1,ifld)+upi(ix2,iy+1,ifld) -
     .                           upi(ix,iy  ,ifld)-upi(ix1,iy  ,ifld))*
     .                                                     gyf(ix,iy) +
     .                          (upi(ix,iy  ,ifld)+upi(ix1,iy  ,ifld) -
     .                           upi(ix,iy-1,ifld)-upi(ix3,iy-1,ifld))*
     .                                                     gyf(ix,iy-1) )
               endif
               wvh(ix,iy,ifld) = wvh(ix,iy,ifld) + cfvcsy(ifld)*cfvisy*
     .                                   visy(ix,iy,ifld)*dupdy**2
	       wvh(ix,iy,ifld) = wvh(ix,iy,ifld) -
     .                             sin(thetacc)*cfvcsy(ifld)*cfvisy*
     .                                   visy(ix,iy,ifld)*dupdx*dupdy
            resei(ix,iy) = resei(ix,iy) + wvh(ix,iy,ifld)*vol(ix,iy)
  155       continue   # loop over up species ifld
  156    continue
 157  continue


c*******************************************************************
c ... Define a background ion energy source to prevent very low Ti
c******************************************************************
      do iy = iys, iyf  #j2, j5
        do ix = ixs, ixf  #i2, i5
          pwribkgold = pwribkg(ix,iy)
          pwribkg(ix,iy) = (tibg*ev/ti(ix,iy))**iteb*pwribkg_c
        enddo
      enddo

      do iy = j2, j5
        do ix = i2, i5
          resei(ix,iy) = resei(ix,iy) + pwribkg(ix,iy)*vol(ix,iy)
        enddo
      enddo
if (TimingPandfOn.gt.0) TotTimeEngBal=TotTimeEngBal+tock(TimeEngBal)

**********************************************************************
*  --  Equations to be solved --
**********************************************************************
      do 270 iy = j2, j5
         do 260 ix = i2, i5
            do 254 ifld = 1, nisp
	       if(isnionxy(ix,iy,ifld) .eq. 1) then
                  iv = idxn(ix,iy,ifld)
                  yldot(iv) = (1-iseqalg(iv)) *
     .                        resco(ix,iy,ifld)/(vol(ix,iy)*n0(ifld))
               endif
 254        continue
            do 255 ifld = 1, nusp
	       if(isuponxy(ix,iy,ifld) .eq. 1) then
                  iv = idxu(ix,iy,ifld)
                  yldot(iv) = (1-iseqalg(iv)) *
     .                        resmo(ix,iy,ifld)/(volv(ix,iy)*fnorm(ifld))
                  do jx = 1, nxpt
                     if (ix.eq.ixrb(jx) .and. ixmxbcl.eq.1) yldot(iv) =
     .                        resmo(ix,iy,ifld)/(volv(ix,iy)*fnorm(ifld))
                  enddo
               endif
 255        continue
            iv =  idxte(ix,iy)
            iv1 = idxti(ix,iy)
	    if(isteonxy(ix,iy).eq.1) yldot(iv) = (1-iseqalg(iv)) *
     .                                  resee(ix,iy)/(vol(ix,iy)*ennorm)
	    if(istionxy(ix,iy).eq.1) yldot(iv1) = (1-iseqalg(iv1)) *
     .                                 resei(ix,iy)/(vol(ix,iy)*ennorm)
            do 256 igsp = 1, ngsp
	      if(isngonxy(ix,iy,igsp).eq.1) then
                iv2 = idxg(ix,iy,igsp)
                yldot(iv2) = (1-iseqalg(iv2)) *
     .                      resng(ix,iy,igsp)/(vol(ix,iy)*n0g(igsp))
              endif
	      if(istgonxy(ix,iy,igsp).eq.1) then
                iv2 = idxtg(ix,iy,igsp)
                yldot(iv2) = (1-iseqalg(iv2)) *
     .                      reseg(ix,iy,igsp)/(vol(ix,iy)*ennorm)
              endif
 256        continue
 260     continue
 270  continue
c ... The factor (1-iseqalg(iv)) above forces yldot=0 for algebraic
c ... equations, except up(nx,,); these yldot are subsequently set in
c ... subroutine bouncon.


c  POTEN calculates the electrostatic potential, and BOUNCON calculates the
c  equations for the boundaries. For the vodpk solver, the B.C. are ODEs
c  in time (rate equations).  Both bouncon and poten must be called before
c  the perturbed variables are reset below to get Jacobian correct

      if (isphion.eq.1) call poteneq (neq, yl, yldot)

      call bouncon (neq, yl, yldot)

c...  Finally, reset some source terms if this is a Jacobian evaluation
         if (xc .ge. 0 .and. yc .ge. 0.and.SaveOld.gt.0) then
            ix1 = ixm1(xc,yc)
            if(isimpon.gt.0) pwrzec(xc,yc) = pradold
            pwrebkg(xc,yc) = pwrebkgold
            pwribkg(xc,yc) = pwribkgold
            erliz(xc,yc) = erlizold
            erlrc(xc,yc) = erlrcold
            eeli(xc,yc) = eeliold
            fqp(ix1,yc) = fqpom
            fqp(xc,yc) = fqpo
            frice(ix1,yc) = friceom
            frice(xc,yc) = friceo
            upe(ix1,yc) = upeom
            upe(xc,yc) = upeo
            psordis(xc,yc) = psordisold
            do ifld = 1, nfsp
               psorc(xc,yc,ifld) = psorold(ifld)
               psorxr(xc,yc,ifld) = psorxrold(ifld)
               frici(ix1,yc,ifld) = friciom(ifld)
               frici(xc,yc,ifld) = fricio(ifld)
               upi(ix1,yc,ifld) = upiom(ifld)
               upi(xc,yc,ifld) = upio(ifld)
               uup(ix1,yc,ifld) = uupom(ifld)
               uup(xc,yc,ifld) = uupo(ifld)
               nucxi(xc,yc,ifld) = nucxiold(ifld)
               nueli(xc,yc,ifld) = nueliold(ifld)
            enddo
            do igsp = 1, ngsp
               nucx(xc,yc,igsp) = nucxold(igsp)
               nurc(xc,yc,igsp) = nurcold(igsp)
               nuiz(xc,yc,igsp) = nuizold(igsp)
               nuelg(xc,yc,igsp) = nuelgold(igsp)
               nuix(xc,yc,igsp) = nuixold(igsp)
               psorgc(xc,yc,igsp) = psorgold(igsp)
               psorrgc(xc,yc,igsp) = psorrgold(igsp)
               psorcxgc(xc,yc,igsp) = psorcxgold(igsp)
            enddo
         endif

      if (ismcnon .eq. 4) then # test a different fluid model in the preconditioner
         if (yl(neq+1) .gt. 0) then   # Precon eval
            parvis=parvis/pnc_cfparvis
            travis=travis/pnc_cftravis
            do ifld=1,nisp
              ni(:,:,ifld)=ni(:,:,ifld)/pnc_cfni(ifld)
              up(:,:,ifld)=up(:,:,ifld)/pnc_cfup(ifld)
            enddo
         endif
      end if #ismcnon

c ... Accumulate cpu time spent here.
      if(xc .lt. 0) then
         ttotfe = ttotfe + gettime(sec4) - tsfe
      else
         ttotjf = ttotjf + gettime(sec4) - tsjf
      endif
      if (TimingPandfOn.gt.0) TotTimePandf=TotTimePandf+tock(TimePandf)
      return
      end
c****** end of subroutine pandf ************



c-----------------------------------------------------------------------
      subroutine pandf1(xc, yc, ieq, neq, time, yl, yldot)

c ... Calculates matrix A and the right-hand side depending on the
c     values of xc, yc.
c  Definitions for argument list
c
c  Input variables:
c    xc is poloidal index of perturbed variablefor Jacobian calc,
c       or =-1 for full RHS evaluation
c    yc is radial index for perturbed variable for Jacobian calc,
c       or =-1 for full RHS evaluation
c    ieq is the eqn number for Jacobian eval; not presently used
c    neq is the total number of variables
c    time is the present physical time; useable by VODPK but not NKSOL
c    yl is the vector of unknowns
c  Output variables:
c    yldot is the RHS of ODE solver or RHS=0 for Newton solver (NKSOL)

      implicit none
      Use(Dim)     # nusp,nisp,ngsp
      Use(Math_problem_size)   # neqmx(for arrays not used here)
      Use(UEpar)   # svrpkg,isbcwdt,isnionxy,isuponxy,isteonxy,istionxy,
                   # isngonxy,isphionxy
cc      Use(Selec)   # i2,i5,j2,j5
      Use(Time_dep_nwt)   # nufak,dtreal,ylodt,dtuse
      Use(Indexes) # idxn,idxg,idxu,dxti,idxte,idxphi
      Use(Ynorm)   # isflxvar,isrscalf
      Use(Share)    # geometry,nxc,isnonog,cutlo
      Use(Indices_domain_dcl) # ixmnbcl,ixmxbcl,iymnbcl,iymxbcl
      Use(Compla)  # zi
      Use(Xpoint_indices)      # ixpt1,ixpt2,iysptrx

*  -- arguments
      integer xc, yc, ieq, neq     # ieq is the equation index for Jac. calc
      real time, yl(neqmx),yldot(neq)

*  -- local variables
      integer ix,iy,igsp,iv,iv1,ifld,j2l,j5l,i2l,i5l
      character*80 msgjm
      integer nrcv, ierrjm, ijmgetmr

ccc      save

c
c     Check if "k" or "kaboom" has been typed to jump back to the parser
c
      if (((svrpkg.eq.'nksol') .or. (svrpkg.eq.'petsc')) .and. iskaboom.eq.1) then
                              #can only call once - preserves 's' in vodpk
        ierrjm = ijmgetmr(msgjm,80,1,nrcv)
        if (ierrjm .eq. 0) then
          if (msgjm(1:nrcv).eq.'kaboom' .or. msgjm(1:nrcv).eq.'k')then
            call xerrab("")
          endif
        endif
      endif

c     check if a "ctrl-c" has been type to interrupt - from basis
      call ruthere

c
c  PANDF calculates the equations in the interior of the grid, plus calls
c  bouncon for B.C. and poten for potential
c
      call pandf (xc, yc, neq, time, yl, yldot)
c
c...  If isflxvar=0, we use ni,v,Te,Ti,ng as variables, and the ODEs need
c...  to be modified as original equations are for d(nv)/dt, etc
c...  If isflxvar=2, variables are ni,v,nTe,nTi,ng. Boundary equations and
c...  potential equations are not reordered.

      if(isflxvar.ne.1 .and. isrscalf.eq.1) call rscalf(yl,yldot)
c
c ... Now add psuedo or real timestep for nksol method, but not both
      if (nufak.gt.1.e5 .and. dtreal.lt.1.e-5) then
         call xerrab('***Both 1/nufak and dtreal < 1.e5 - illegal***')
      endif

c...  Add a real timestep, dtreal, to the nksol equations
c...  NOTE!! condition yl(neq+1).lt.0 means a call from nksol, not jac_calc

      if(dtreal < 1.e15) then
       if((svrpkg=='nksol' .and. yl(neq+1)<0) .or. svrpkg == 'petsc') then
         if (isbcwdt .eq. 0) then  # omit b.c. eqns
cccMER   NOTE: what about internal guard cells (for dnbot,dnull,limiter) ???
            j2l = 1
            j5l = ny
            i2l = 1
            i5l = nx
         else                      # include b.c. eqns
            j2l = (1-iymnbcl)
            j5l = ny+1-(1-iymxbcl)
            i2l = (1-ixmnbcl)
            i5l = nx+1-(1-ixmxbcl)
         endif
         do iy = j2l, j5l    # if j2l=j2, etc., omit the boundary equations
            do ix = i2l, i5l
              do ifld = 1, nisp
                if(isnionxy(ix,iy,ifld) .eq. 1) then
                  iv = idxn(ix,iy,ifld)
                  yldot(iv) = (1.-fdtnixy(ix,iy,ifld))*yldot(iv)
                  if(zi(ifld).eq.0. .and. ineudif.eq.3) then
                    yldot(iv) = yldot(iv) - (1/n0(ifld))*
     .                          (exp(yl(iv))-exp(ylodt(iv)))/dtuse(iv)
                  else
                    yldot(iv) =yldot(iv)-(yl(iv)-ylodt(iv))/dtuse(iv)
                  endif
                endif
              enddo
               if(ix.ne.nx+2*isbcwdt) then
                              # nx test - for algebr. eq. unless isbcwdt=1
                  do ifld = 1, nusp
                    if(isuponxy(ix,iy,ifld).eq.1) then
                      iv = idxu(ix,iy,ifld)
                      yldot(iv) = (1.-fdtupxy(ix,iy,ifld))*yldot(iv)
                      yldot(iv) = yldot(iv)-(yl(iv)-ylodt(iv))/dtuse(iv)
                    endif
                  enddo
               endif
               if (isteonxy(ix,iy) == 1) then
                 iv =  idxte(ix,iy)
                 yldot(iv) = (1.-fdttexy(ix,iy))*yldot(iv)
                 yldot(iv) = yldot(iv) - (yl(iv)-ylodt(iv))/dtuse(iv)
               endif
               if (istionxy(ix,iy) == 1) then
                 iv1 = idxti(ix,iy)
                 yldot(iv1) = (1.-fdttixy(ix,iy))*yldot(iv1)
                 yldot(iv1)=yldot(iv1) - (yl(iv1)-ylodt(iv1))/dtuse(iv1)
               endif
               do igsp = 1, ngsp
                  if(isngonxy(ix,iy,igsp).eq.1) then
                     iv = idxg(ix,iy,igsp)
                     yldot(iv) = (1.-fdtngxy(ix,iy,igsp))*yldot(iv)
                     if(ineudif.eq.3) then
                       yldot(iv) = yldot(iv) - (1/n0g(igsp))*
     .                            (exp(yl(iv))-exp(ylodt(iv)))/dtuse(iv)
                     else
                       yldot(iv) =yldot(iv)-(yl(iv)-ylodt(iv))/dtuse(iv)
                     endif
                  endif
               enddo
               do igsp = 1, ngsp
                  if(istgonxy(ix,iy,igsp).eq.1) then
                     iv = idxtg(ix,iy,igsp)
                     yldot(iv) = (1.-fdttgxy(ix,iy,igsp))*yldot(iv)
                     yldot(iv) =yldot(iv)-(yl(iv)-ylodt(iv))/dtuse(iv)
                  endif
               enddo
               if (isphionxy(ix,iy).eq.1 .and. isbcwdt.eq.1) then
                  iv = idxphi(ix,iy)
                  yldot(iv) = (1.-fdtphixy(ix,iy))*yldot(iv)
                  yldot(iv) = yldot(iv) - (yl(iv)-ylodt(iv))/dtuse(iv)
               endif

            enddo
         enddo

C...  Now do an additional relaxation of the potential equations with
c...  timestep dtphi
        if (dtphi < 1e10) then
          do iy = 0, ny+1
            do ix = 0, nx+1
              if (isphionxy(ix,iy) == 1) then
                iv = idxphi(ix,iy)
                yldot(iv) = yldot(iv) - (yl(iv)-ylodt(iv))/dtphi
              endif
            enddo
          enddo
        endif

       endif   #if-test on svrpkg and yl(neq+1)
      endif    #if-test on dtreal

      return
      end
c****** end of subroutine pandf1 ************
