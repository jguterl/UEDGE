#include "bbb.h"
#include "../com/com.h"
#include "../mppl.h"
#include "../sptodp.h"

c-----------------------------------------------------------------------
      subroutine mombal0 (nisp, nhsp, nzsp, minu, ziin,
     .                                      misotope, natomic, nchstate)
c ... Compute 'misotope', 'nchstate', and 'natomic', and allocate memory
c     for arrays used in subroutine mombal.

      implicit none

c ... Input arguments:
      integer nisp      ! total number of ion species
      integer nhsp      ! total number of hydrogenic ion species
      integer nzsp(ngspmx-1)   ! number of charge states for each imp isotope
      real minu(nisp)   ! mass (in amu) of ion species
      real ziin(nisp)   ! charge (in units of e) of ion species

c ... Output arguments:
      integer misotope  ! total number of isotopes (including electrons)
      integer natomic(*)   ! maximum charge state of each isotope
      integer nchstate     ! maximum charge state among all isotopes

c ... Local variables:
      integer misa, ifld, jz

c ... Loop over ion species, looking for change to a new isotope, and
c     finding maximum charge state.
      natomic(1) = 1   ! electrons are "isotope 1"
      nchstate = 0
      misa = 2
      do ifld = 1, nhsp
         natomic(misa) = max(nint(ziin(ifld)), 1)   ! must be .ge. 1
         nchstate = max(nchstate, natomic(misa))
         if (ifld .eq. nhsp) go to 50
         if (minu(ifld+1) .ne. minu(ifld)) misa = misa + 1
      enddo
 50   misotope = misa
      do jz = 1, ngspmx-1
         if (nzsp(jz)==0) exit
         misotope = misotope + 1
	 if (misotope .gt. MXMISO) then
	   write(*,*) 'misotope:',misotope,MXMISO
           call remark("subroutine mombal0 error: ")
           call remark("To avoid write out-of-bounds for array natomic")
           call remark("increase the value of MXMISO and recompile.")
	       call xerrab("")
         endif
         natomic(misotope) = nzsp(jz)
         nchstate = max(nchstate, natomic(misotope))
      enddo

c ... Allocate memory for arrays used in subroutine mombal.
      call gallot("Reduced_ion_interface", 0)

      return
      end
c****** end of subroutine mombal0 ************
c-----------------------------------------------------------------------
      subroutine mombal (ix,ix1,iy)
c ... Prepare information needed to call Steve Hirshman reduced-ion
c     momentum-balance routine, and distribute results from it.
c     Note that results are computed at a poloidal density-cell face.
c     The parallel flow velocities are returned in arrays up and upe.

c ... Common blocks:
            use Dim           ! nx,ny,nisp
            use Selec         ! isupgoon
            use Comgeo        ! sx,rrv
            use UEpar         ! lnlam,isupgon
            use Cfric         ! frice,frici
            use Compla        ! ni,up,upe,upi,te,ti,ng,ne
            use Gradients     ! gpix,gtix,ex,gtex,gpex
            use Comflo        ! fqp
            use UEint         ! minu,ziin
            use Reduced_ion_interface  ! misotope,nchstate,natomic,
                                 ! amu,tempa,qneut,uneut,den,gradp,
                                 ! gradt,friction,nuion,nurec,qcond,ucond
                                 ! friccomp
            use Imprad        ! ismctab
            use Phyvar        ! ev
            use Comtra        ! fricflf,sigvi_floor
            use Share         ! cutlo

c ... External functions:
      implicit none

c ... Input arguments:
      integer ix    ! poloidal index of density cell to left of face
      integer ix1   ! poloidal index of density cell to right of face
      integer iy    ! radial index of density cell

      real rra, rsa

c ... Local variables:
      integer ifld, misa, nz
      integer ldir
      real rdum, dloglam, epar, parcurrent, umass, lmfpe, lmfpi,
     .     ltmax, tif, flxlimf, umassni, massni

c ... Determine a flux-limit factor for all input terms by finding the min
c ... scale length. First consider the electrons
      den(1,1) = 0.5 * (ne(ix,iy) + ne(ix1,iy))
      tempa(1) = 0.5 * (te(ix,iy) + te(ix1,iy))
      tif = 0.5 * (ti(ix,iy) + ti(ix1,iy))
      ltmax = min( abs(tempa(1)/(rrv(ix,iy)*gtex(ix,iy) + cutlo)),
     .             abs(tif/(rrv(ix,iy)*gtix(ix,iy) + cutlo)),
     .             abs(den(1,1)*tempa(1)/
     .                          (rrv(ix,iy)*gpex(ix,iy) + cutlo)) )
      lmfpe = 1.e16*((tempa(1)/ev)**2/den(1,1)) !Approx Coulomb e-mean-free-path
      lmfpi = 1.e16*(tif/ev)**2/den(1,1)  ! Approx. Coulomb i-mean-free-path
c ... Now check ion pressure gradient scale-lengths
      do nz = 1, 1   ! previously nisp - only do majority species now
         if (zi(nz) .gt. 1.e-10) then  ! omits neutrals species
            ltmax = min( ltmax, abs( 0.5*(pri(ix,iy,nz)+pri(ix1,iy,nz))/
     .                           (rrv(ix,iy)*gpix(ix,iy,nz) + cutlo) ) )
         endif
      enddo
      flxlimf = 1 / (1 + fricflf*((lmfpe+lmfpi)/ltmax)**2)

c ... Store electron density and temperature, gradients of electron
c     pressure and temperature, (hydrogenic) gas density, ionization rate
c     of gas, and a few other (physically meaningless) quantities for which
c     there are array locations.
      amu(1) = 5.45e-4
      gradp(1,1) = flxlimf*rrv(ix,iy) * gpex(ix,iy)
      gradt(1,1) = flxlimf*rrv(ix,iy) * gtex(ix,iy)
      den(2,0) = 0.5 * (ng(ix,iy,1) + ng(ix1,iy,1))
      nuion(2,0) = den(1,1) * (rsa(tempa(1), den(1,1), 0., 0)
     .                         + sigvi_floor)
      den(1,0) = 0.
      qneut(1) = 0.
      uneut(1) = up(ix,iy,1)  ! netural velocity. Cant be zero?; was umass
      nuion(1,0) = 0.
      nuion(1,1) = 0.
      nurec(1,1) = 0.

c ... Loop over isotopes.
      ifld = 1
      do misa = 2, misotope
         amu(misa) = minu(ifld)   ! Store mass of this isotope

c ... Store ionization rate of neutral if this isotope is an impurity.
         if (natomic(misa) .gt. 1) then
            den(misa,0) = 0.              ! impurity neutral density
            if (ismctab .eq. 1) then
               call imprates(tempa(1), 0, natomic(misa),
     .                       nuion(misa,0), rdum, rdum)
               nuion(misa,0) = den(1,1)*nuion(misa,0) !sigv-->ne*sigv
            elseif (ismctab .eq. 2) then

               call mcrates(den(1,1), tempa(1), 0., 0, natomic(misa),
     .                      znucl(ifld), nuion(misa,0), rdum, rdum)
               nuion(misa,0) = den(1,1)*nuion(misa,0) !sigv-->ne*sigv
            endif
         endif

c ...    Loop over charged states, storing ni and parallel gradients
c        of pressure and Ti at poloidal density-cell faces.
         do nz = 1, natomic(misa)
            den(misa,nz) = 0.5 * (ni(ix,iy,ifld) + ni(ix1,iy,ifld))
            gradp(misa,nz) = flxlimf*rrv(ix,iy) * gpix(ix,iy,ifld)
            tempa(misa) = 0.5 * (ti(ix,iy) + ti(ix1,iy))
            gradt(misa,nz) = flxlimf*rrv(ix,iy) * gtix(ix,iy)
c .......   Get ionization and recombination rates.
c           Note that nuion has no meaning for the fully-stripped state,
c           but space is available to store zero returned by imprates.
            if (isupimpap .eq. 0) then  ! omit atomic physics coupling in up
               nuion(misa,nz) = 0.
               nurec(misa,nz) = 0.
            elseif (natomic(misa) .eq. 1) then   ! hydrogenic isotope
               nuion(misa,1) = 0.
               nurec(misa,1) = den(1,1) *
     .                         rra(tempa(1), den(1,1), 0., 1)
            else                             ! impurity isotope
               if (ismctab .eq. 1) then
                  call imprates(tempa(1), nz, natomic(misa),
     .                          nuion(misa,nz), nurec(misa,nz), rdum)
                  nuion(misa,nz) = den(1,1)*nuion(misa,nz) !sigv-->ne*sigv
                  nurec(misa,nz) = den(1,1)*nurec(misa,nz)
               elseif (ismctab .eq. 2) then
               write(*,*) 'mcrates with misa,misotope,natomic',misa,misotope,natomic(misa)
                  call mcrates(den(1,1), tempa(1), 0., nz, natomic(misa),
     .                znucl(ifld), nuion(misa,nz), nurec(misa,nz), rdum)
                  nuion(misa,nz) = den(1,1)*nuion(misa,nz) !sigv-->ne*sigv
                  nurec(misa,nz) = den(1,1)*nurec(misa,nz)
               endif
            endif
            ifld = ifld + 1
            if (natomic(misa).eq.1.and.isupgon(ifld-1).eq.1) ifld=ifld+1
         enddo    ! end of loop over charge states with index nz
      enddo    ! end of loop over isotopes with index misa

c ... Set up other inputs for fmombal, including flow velocities.
 50   epar = flxlimf*rrv(ix,iy) * ex(ix,iy)
      umassni = 0.
      massni = 0.
      do ifld = 1, nusp
        if(zi(ifld) .gt. 1.e-20) then
          umassni = umassni + mi(ifld)*ni(ix,iy,ifld)*up(ix,iy,ifld)
          massni = massni + mi(ifld)*ni(ix,iy,ifld)
        endif
      enddo
      umass = umassni/massni
      do misa = 2, misotope
         qneut(misa) = 0.
         uneut(misa) = up(ix,iy,1)  ! use hydr ion as default (orig)
      enddo
      ldir = 2
      dloglam = lnlam

c ... Call Steve Hirshman reduced-ion momentum-balance routine.
      call fmombal(amu,den,dloglam,epar,friction,gradp,gradt,
     >	 nuion,nurec,qcond,qneut,ucond,uneut,umass,
     >   parcurrent,tempa,natomic,misotope,nchstate,ldir,friccomp)

c ... Distribute results into arrays used in pandf.  Note that we do
c     nothing with qcond.
      fqp(ix,iy) = cfparcur*parcurrent * rrv(ix,iy)*sx(ix,iy)
      frice(ix,iy) = friction(1,1)
      upe(ix,iy) = ucond(1,1)
      ifld = 0
      do misa = 2, misotope
         do nz = 1, natomic(misa)
 60         ifld = ifld + 1
            if (ziin(ifld) .lt. 1.e-10) goto 60    !omits neutral species
            frici(ix,iy,ifld) = friction(misa,nz)
            upi(ix,iy,ifld) = ucond(misa,nz)
            upifmb(ix,iy,ifld) = ucond(misa,nz)    !diagnostic only
         enddo
      enddo

      return
      end
c****** end of subroutine mombal ************
c-----------------------------------------------------------------------
      subroutine mombalni (ix,ix1,iy)
c ... Use force balance for the impurity momentum equation, neglecting
c     inertia, viscosity, and atomic-physics coupling ala Knoll, Campbell.
c     We also use Keilhacker, et al., Nucl. Fusion., Vol. 31, 537 (1991)
c     which differs somewhat from Campbell: alfi gets divided by zeffv
c     but we retain the extra term 0.6*... in betai from Campbell.
c     Note that results are computed at a poloidal density-cell face.
c     The parallel flow velocities are returned in arrays up and upe.

c ... Common blocks:
            use Dim           ! nx,ny,nisp,nhsp,nusp
            use Selec         ! isupgon
            use Comgeo        ! sx,rrv,vol
            use UEpar         ! lnlam,isupgon,isofric,is_z0_imp_const,z0_imp_const
            use Cfric         ! frice,frici,cfgte,cfgti,cftaud,alfe,betai
            use Compla        ! ni,up,upe,upi,te,ti,ng,ne,zeff,netap
            use Gradients     ! gpix,gtix,ex,gtex,gpex
            use Comflo        ! fqp
            use UEint         ! minu,ziin
            use Reduced_ion_interface  ! misotope,nchstate,natomic,dztot,
                                 ! amu,tempa,qneut,uneut,den,gradp,
                                 ! gradt,friction,nuion,nurec,qcond,ucond
            use Imprad        ! ismctab
            use Phyvar        ! ev,qe
            use Comtra        ! fricflf,fupe_cur
            use Share         ! cutlo
            use Coefeq        ! cfnetap
            use Volsrc        ! volmsor
            use Npes_mpi 


c ... Local variables:
      implicit none

c ... Input arguments:
      integer,intent(in):: ix    ! poloidal index of density cell to left of face
      integer,intent(in):: ix1   ! poloidal index of density cell to right of face
      integer,intent(in):: iy    ! radial index of density cell

      integer ifld, misa, nz
      integer ldir
      real rdum, tdum, dloglam, epar, parcurrent, umass, lmfpe, lmfpi,
     .     zeffv, z0, taud, taudeff, ltmax, tif, flxlimf,
     .     dzj, dzz2tot


c ... Determine a flux-limit factor for all input terms by finding the min
c ... scale length. First consider the electrons
      den(1,1) = 0.5 * (ne(ix,iy) + ne(ix1,iy))
      tempa(1) = 0.5 * (te(ix,iy) + te(ix1,iy))
      tif = 0.5 * (ti(ix,iy) + ti(ix1,iy))
      ltmax = min( abs(tempa(1)/(rrv(ix,iy)*gtex(ix,iy) + cutlo)),
     .             abs(tif/(rrv(ix,iy)*gtix(ix,iy) + cutlo)),
     .             abs(den(1,1)*tempa(1)/
     .                          (rrv(ix,iy)*gpex(ix,iy) + cutlo)) )
      lmfpe = 1.e16*((tempa(1)/ev)**2/den(1,1)) !Approx Coulomb e-mean-free-path
      lmfpi = 1.e16*((tif/ev)**2/den(1,1))  ! Approx. Coulomb i-mean-free-path
c ... Now check ion pressure gradient scale-lengths
      do nz = 1, 1      ! previously nisp - only do majority species now
         if (zi(nz) .gt. 1.e-10) then  ! omits neutrals species
            ltmax = min( ltmax, abs( 0.5*(pri(ix,iy,nz)+pri(ix1,iy,nz))/
     .                           (rrv(ix,iy)*gpix(ix,iy,nz) + cutlo) ) )
         endif
      enddo
      flxlimf = 1 / (1 + fricflf*((lmfpe+lmfpi)/ltmax)**2)

c ... Store electron density and temperature, gradients of electron
c     pressure and temperature, (hydrogenic) gas density, ionization rate
c     of gas, and a few other (physically meaningless) quantities for which
c     there are array locations.
      gradp(1,1) = rrv(ix,iy) * gpex(ix,iy)
      gradt(1,1) = rrv(ix,iy) * gtex(ix,iy)
      frice(ix,iy) = -0.71*flxlimf*den(1,1)*gradt(1,1) +
     .                cfnetap*qe*netap(ix,iy)*fqp(ix,iy)/sx(ix,iy)

c ... Loop over charge states to get total impurity density
      ifld = nhsp
      dztot = 0.
      dzz2tot = 0.
      do misa = 3, misotope
         do nz = 1, natomic(misa)
            ifld = ifld + 1
            dzj = 0.5 * (ni(ix,iy,ifld) + ni(ix1,iy,ifld))
            dztot = dztot + dzj
            dzz2tot = dzz2tot + dzj*zi(ifld)**2
         enddo
      enddo

c ... Set the hydrogen values based on electrons (trace limit)
ccc      frici(ix,iy,1) = - frice(ix,iy)   ! needed for hydrogen
ccc   For arbitrary impurity concentration use the following frici:
      frici(ix,iy,1) = - frice(ix,iy) *
     .                   ni(ix,iy,1)*zi(1)**2/(ne(ix,iy)*zeff(ix,iy))
      upe(ix,iy) = up(ix,iy,1) - fupe_cur*fqp(ix,iy)/( sx(ix,iy)*qe*
     .                         rrv(ix,iy)*0.5*(ne(ix,iy)+ne(ix1,iy)) )
      upi(ix,iy,1) = up(ix,iy,1)
      den(2,1) = 0.5 * (ni(ix,iy,1) + ni(ix1,iy,1))

c ... Loop over isotopes for friction coefficients
      ifld = 1
      do misa = 3, misotope      ! only executed if impurities are present

c ... Loop over charged states, storing ni and parallel gradients
c ... of pressure and Ti at poloidal density-cell faces.
         do nz = 1, natomic(misa)  ! note: hydrogen friction set above
            ifld = ifld + 1
            if (ziin(ifld) .lt. 1.e-10) ifld = ifld+1 !skip gas index
            den(misa,nz) = 0.5 * (ni(ix,iy,ifld) + ni(ix1,iy,ifld))
            gradp(misa,nz) = rrv(ix,iy) * gpix(ix,iy,ifld)
            tempa(misa) = 0.5 * (ti(ix,iy) + ti(ix1,iy))
            gradt(misa,nz) = rrv(ix,iy) * gtix(ix,iy)
            zeffv = 0.5*(zeff(ix,iy)+zeff(ix1,iy))
            if (is_z0_imp_const == 0) then
              z0 = den(1,1)*zeffv/den(2,1) - 1.
            else ! fixed user input
              z0 = z0_imp_const
            endif
            if (isbetaicalc(ifld) == 1) then
              betai(ifld)=cfgti*1.56*zi(ifld)**2*(1+1.414*z0)*(1+.52*z0)/
     .                       ( (1+2.65*z0)*(1+.285*z0)*( z0 +
     .                         sqrt( 0.5*(mi(1)+mi(ifld))/mi(ifld)) ) )
     .                       + 0.6*(zi(ifld)**2*den(misa,nz)/dzz2tot - 1.)
            endif
            if (isalfecalc(ifld) == 1) then
              alfe(ifld) = cfgte*2.2*zi(ifld)**2*(1+.52*zeffv) /
     .                      ( (1+2.65*zeffv)*(1+.285*zeffv)*zeffv )
            endif
c... NOTE:next coefficient 12*pi*sqrt(pi/2)*epsilon**2/e**4 = 5.624e54 in mks
            taud = cftaud*5.624e54*sqrt(mi(1))*mi(ifld)*tempa(misa)**1.5 /
     .             ( lnlam*den(misa,nz)*zi(ifld)**2*(mi(1)+mi(ifld)) )
            taudeff = flxlimf*taud*den(misa,nz)*(1+2.65*z0)*(1+.285*z0) /
     .                         ( den(1,1)*(1+.24*z0)*(1+.93*z0) )
            upi(ix,iy,ifld) = up(ix,iy,1) + (taudeff/mi(1)) * (
     .                         - gradp(misa,nz)/den(misa,nz)
     .                         + alfe(ifld)*gradt(1,1)
     .                         + betai(ifld)*gradt(misa,nz)
     .                         + qe*zi(ifld)*rrv(ix,iy)*ex(ix,iy)
     .                         + volmsor(ix,iy,ifld)/
     .                                       (den(misa,nz)*vol(ix,iy)) )
c ...       For force balance, frici just balances E-field and pressure
c ...       No flxlimf for 1st option; it only enhances (1/taudeff)
            if (nusp-isupgon(1) .eq. 1) then  !only frici(,,1) used here
              frici(ix,iy,ifld) =-qe*zi(ifld)*den(misa,nz)*
     .                           rrv(ix,iy)*ex(ix,iy) + gradp(misa,nz)
            else ! multi ion mom eqns; drag calc elsewhere (w0) if isofric=1
              frici(ix,iy,ifld) =flxlimf*den(misa,nz)*(
     .                                          alfe(ifld)*gradt(1,1) +
     .                                     betai(ifld)*gradt(misa,nz) +
     .                                               (1-isofric)*mi(1)*
     .                           (up(ix,iy,1)-up(ix,iy,ifld))/ taudeff )
            endif
         enddo
      enddo

      return
      end
c****** end of subroutine mombalni ************
c-----------------------------------------------------------------------
      subroutine timimpfj (tsimp, xc)
            use Timing    ! ttimpfe,ttimpjf
      real(kind=4):: tsimp
      integer xc
      real(kind=4):: sec4, gettime, dtimp

      dtimp = gettime(sec4) - tsimp
      if (xc .lt. 0) then
         ttimpfe = ttimpfe + dtimp
      else
         ttimpjf = ttimpjf + dtimp
      endif

      return
      end
c---- end of subroutine timimpfj ---------------------------------------
c-----------------------------------------------------------------------


c --------------------------------------------------------------------------
c SUBROUTINE TO SET UP ENERGY EQUATION FOR NEUTRAL GAS
c --------------------------------------------------------------------------

      subroutine engbalg

            use Dim       ! nx,ny,nhsp,nisp,ngsp,nxpt
            use Xpoint_indices       ! ixlb,ixpt1,ixpt2,ixrb,iysptrx1
            use Share     ! geometry,nxc,isnonog,cutlo,islimon,ix_lim,iy_lims
            use Phyvar    ! pi,me,mp,ev,qe,rt8opi
            use UEpar     ! methg,
                    ! qfl,csh,qsh,cs,
                    ! isupgon,iigsp,nlimgx,nlimgy,nlimiy,rld2dx,rld2dy
            use Coefeq    ! cngfx,cngfy,cngmom,cmwall,cdifg,rld2dxg,rld2dyg
            use Bcond     ! albedoo,albedoi
            use Parallv   ! nxg,nyg
            use Rccoef    ! recylb,recyrb,recycw,recycz,sputtr
            use Selec     ! i1,i2,i3,i4,i5,i6,i7,i8,j1,j2,j3,j4,j5,j6,j7,j8
                    ! xlinc,xrinc,yinc,ixm1,ixp1,stretcx
            use Comgeo    ! vol, gx, gy, sx ,xy
            use Noggeo    ! fym,fy0,fyp,fymx,fypx,angfx
            use Compla    ! mi, zi, ni, uu, up, v2, v2ce, vygtan, mg
            use Comflo    ! fngx,fngy,fngxy,fnix,fniy
            use Conduc    ! nuiz, nucx, nuix
            use Rhsides   ! resng,psor,psorg,psorrg,sniv
            use Comtra    ! flalfgx,flalfgy
            use Locflux   ! floxg,floyg,conxg,conyg
            use Indices_domain_dcl     ! iymnbcl,iymxbcl
            use Volsrc    ! volpsorg
            use Wkspace   ! w0,w1,etc
            use MCN_dim   !
            use MCN_sources    ! cfneutsor_ei


      implicit none

*  -- local variables
      real vtn, vtnp, qr, qtgf, nconv, grdnv, difgx2
      real tnuiz,ngnot,lmfp,ty0,ty1,nlmt,nu1,nu2,ffyi,ffyo
      real vt0,vt1,wallfac,lxtgc,dupdx,dupdy,fniy_recy,thetacc
      real vttn,vttp
      integer ifld,iixt,iy1, methgx, methgy, iy2, jx
      integer ix,iy,igsp,iv,iv1,iv2,iv3,ix1,ix2,ix3,ix4,ix5,ix6
      real t0,t1,t,a

*  -- initialize some arrays to 0
      do igsp = 1, ngsp
        do iy = j1, j6
          do ix = i1, i6
            floxge(ix,iy,igsp) = 0.0e0
            floyge(ix,iy,igsp) = 0.0e0
            conxge(ix,iy,igsp) = 0.0e0
            conyge(ix,iy,igsp) = 0.0e0
          enddo
        enddo
      enddo

*  ---------------------------------------------------------------------
*  Compute thermal conductances
*  ---------------------------------------------------------------------
c ... Compute poloidal conduction
      do igsp = 1,ngsp
        do iy = j4, j8
          do ix = i1, i5
            ix2 = ixp1(ix,iy)

            t0 = max (tg(ix,iy,igsp), temin*ev)
            t1 = max (tg(ix2,iy,igsp), temin*ev)
            vt0 = sqrt(t0/mg(igsp))
            vt1 = sqrt(t1/mg(igsp))
c... flux-limit occurs in building hcxg - do not flux-limit 2nd time
            conxge(ix,iy,igsp) = sx(ix,iy) * hcxg(ix,iy,igsp) * gxf(ix,iy)
          enddo
          conxge(nx+1,iy,igsp) = 0
        enddo
      enddo

*  -- compute radial conduction conyge
      do igsp = 1, ngsp
        do iy = j1, j5
          do ix = i4, i8
            conyge(ix,iy,igsp) = sy(ix,iy)*hcyg(ix,iy,igsp)*gyf(ix,iy)
          enddo
        enddo
      enddo

      do igsp = 1, ngsp
        do ix = i1, i6
          conyge(ix,ny+1,igsp) = 0.0e0
        enddo
      enddo
*  ---------------------------------------------------------------------
*  compute convective flow of tg
*  ---------------------------------------------------------------------

*  -- compute floxge --

      do igsp = 1, ngsp
        do iy = j4, j8
          do ix = i1, i5
            floxge(ix,iy,igsp) = cfcvtg*2.5*fngx(ix,iy,igsp)
          enddo
          floxge(nx+1,iy,igsp) = 0.
        enddo
      enddo

*  -- Correct bdry:remove any inward power from plates; ok in parallel
      do igsp = 1, ngsp
        do iy = j4, j8
          do jx = 1, nxpt  !if at plate, sub (1-cfloxiplt)*neut-contrib
            if(ixmnbcl==1) then  !real plate-need for parallel UEDGE
              iixt = ixlb(jx) !left plate
              if(fngx(iixt,iy,igsp) > 0.) then
                floxge(iixt,iy,igsp) = floxge(iixt,iy,igsp) -
     .                    (1.-cfloxiplt)*cfcvti*2.5*fngx(iixt,iy,igsp)
              endif
            endif
            if(ixmxbcl==1) then !real plate-need for parallel UEDGE
              iixt = ixrb(jx) ! right plate
              if(fngx(iixt,iy,igsp) < 0.) then
                floxge(iixt,iy,igsp) = floxge(iixt,iy,igsp) -
     .                   (1.-cfloxiplt)*cfcvti*2.5*fngx(iixt,iy,igsp)
              endif
              floxge(ixrb(jx)+1,iy,igsp) = 0.0e0 !cosmetic
            endif
          enddo
        enddo
      enddo

*  -- compute floyge --

      do igsp = 1, ngsp
        do iy = j1, j5
          do ix = i4, i8
            floyge(ix,iy,igsp) = cfcvtg*2.5*fngy(ix,iy,igsp)
          enddo
        enddo
      enddo

*  -- Combine conduction/convection to compute thermal energy flow --
      do igsp = 1,ngsp
        if(istgon(igsp) == 1) then
          call fd2tra (nx,ny,floxge(0,0,igsp),floyge(0,0,igsp),
     .                 conxge(0,0,igsp),conyge(0,0,igsp),tg(0,0,igsp),
     .                 fegx(0,0,igsp),fegy(0,0,igsp),0,methi)
        endif
      enddo

*  ---------------------------------------------------------------------
*  compute the energy residuals.
*  ---------------------------------------------------------------------

*  -- total energy residual and equipartition --

      do igsp = 1, ngsp
        do iy = j2, j5
          do ix = i2, i5
            ix1 = ixm1(ix,iy)
            reseg(ix,iy,igsp)= -( fegx(ix,iy,igsp)-fegx(ix1,iy,  igsp)+
     .                            fegy(ix,iy,igsp)-fegy(ix, iy-1,igsp) )
            reseg(ix,iy,igsp)= reseg(ix,iy,igsp) + vol(ix,iy)*
     .                      eqpg(ix,iy,igsp)*(ti(ix,iy)-tg(ix,iy,igsp))+
     .                   cftgdiss(igsp)*psorg(ix,iy,igsp)*tg(ix,iy,igsp)
          enddo
        enddo
      enddo

*  -- Energy transfer to impurity neutrals at tg(,,igsp)
      if (ngsp >= 2) then   ! for now, specialized to igsp=2 only
        do ifld = nhsp+1, nisp
          do iy = j2, j5    ! iys,iyf limits dont seem to work(?)
            do ix = i2, i5
              resei(ix,iy) =resei(ix,iy) -cftiimpg*1.5*ni(ix,iy,ifld)*
     .                      (nucxi(ix,iy,ifld)+nueli(ix,iy,ifld))*
     .                      (ti(ix,iy) - tg(ix,iy,2))*vol(ix,iy)
            enddo
          enddo
        enddo
      endif

      return
      end

c-----------------------------------------------------------------------
c  END subroutine engbalg - THE NEUTRAL GAS ENERGY EQUATION
c-----------------------------------------------------------------------



c-----------------------------------------------------------------------
      subroutine rscalf(yl, yldot)

c...  This routine reorders the equation for yldot depending on what
c...  variables are used, i.e., n,nv,nT, or n,v,T, or n,v,nT

            use Dim        ! nx,ny,nusp,nisp
            use Math_problem_size    ! neqmx for arrays not used here 
            use Ynorm      ! isflxvar,nnorm,n0
            use Compla 
            use Indexes 
            use Selec      ! i2,i5,j2,j5,ixp1
            use Coefeq     ! cngtgx
            use UEpar      ! isnion,isupon,isteon,istion,isngon,isnionxy,isuponxy,
                     ! isteonxy,istionxy,isngonxy,isphionxy
            use Rhsides    ! resco
            use Comgeo     ! vol

      implicit none
*  -- Input parameters
      real,intent(in):: yl(*)
      real::yldot(*)

*  -- Local variables
      integer ifld,ix,iy,iv,iv1,iv2,ix1,igsp
      real nbv, nbvdot, nbidot, nbedot, nbgdot, yldot_np1, nbg2dot(ngsp)

c...  If isflxvar=0, we use ni,v,Te,Ti,ng as variables, and the ODEs need
c...  to be modified as original equations are for d(nv)/dt, etc
c...  If isflxvar=2, variables are ni,v,nTe,nTi,ng.  Boundary eqns and
c...  potential are not reordered (loops from i2-i5 and j2-j5).

*****************************************************************
*  --  ODE Equations to be solved - rescaling
*****************************************************************

      do 270 iy = j2, j5
         do 260 ix = i2, i5
            nbedot = 0.
            nbidot = 0.
            nbgdot = 0.
            nbvdot = 0.
ccc            if (isngonxy(ix,iy,1) .eq. 1) nbidot = cngtgx(1)*yldot(idxg(ix,iy,1))
            do 255 ifld = 1, nisp
	       if (isnionxy(ix,iy,ifld) .eq. 1) then
                  iv = idxn(ix,iy,ifld)
                  nbidot = nbidot + yldot(iv)*n0(ifld)
                  if (isupgon(1)==1 .and. zi(ifld)==0) then  !neutral hyd
                    nbgdot = yldot(iv)*n0(ifld)
                  endif
                  nbedot = nbedot + zi(ifld)*yldot(iv)*n0(ifld)
               endif
  255       continue
            do igsp = 1, ngsp
               nbg2dot(igsp) = 0.

c              if(isngonxy(ix,iy,ifld) .eq. 1) then !JG: ifld out of bound !!!!!!
               if(isngonxy(ix,iy,igsp) .eq. 1) then
                 iv = idxg(ix,iy,igsp)
                 nbg2dot(igsp) = yldot(iv)*n0g(igsp)
               endif
            enddo
c...    Omit cases where iseqalg=1 indicating algebraic b.c. and ix=nx for up.
c...    Could cause trouble for nisp.ne.nusp, but majority species should
c...    be the first ones in the nisp series, so probably ok.
            do 257 ifld = 1, nusp
	       if (isuponxy(ix,iy,ifld) .eq. 1) then
                  ix1 = ixp1(ix,iy)
                  iv2 = idxu(ix,iy,ifld)
		if (iseqalg(iv2)== 0.and.isnionxy(ix,iy,ifld)==1) then
                  iv = idxn(ix,iy,ifld)
                  iv1 = idxn(ix1,iy,ifld)
                     yldot_np1 = resco(ix1,iy,ifld)/(vol(ix1,iy)*n0(ifld))
! need to use resco rather than yldot if dtreal is added; recursive prob.
c ....            Fix limiter case with algebraic eqns, not ODEs
                     if (iseqalg(iv).eq.1) then
                        if (isnupdot1sd == 0) then
                           nbvdot = yldot_np1*n0(ifld)
                        else
                           nbvdot = yldot(iv1)*n0(ifld)
                        endif
                        nbv = ni(ix1,iy,ifld)
                     elseif (iseqalg(iv1).eq.1) then
                        nbvdot = yldot(iv)*n0(ifld)
                        nbv = ni(ix,iy,ifld)
                     else
                        if (isnupdot1sd == 0) then
                          nbvdot = 0.5*(yldot(iv) + yldot_np1)*n0(ifld)
                        else
                          nbvdot = yldot(iv)*n0(ifld)
                        endif
                         nbv = 0.5*(ni(ix,iy,ifld) + ni(ix1,iy,ifld))
                     endif
                     yldot(iv2) = (yldot(iv2)*n0(ifld) - yl(iv2)*nbvdot)
     .                                                              /nbv
                endif
               endif
  257       continue
            if (isflxvar .eq. 0) then
               iv =  idxte(ix,iy)
               iv1 = idxti(ix,iy)
	       if(isteonxy(ix,iy).eq.1 .and. iseqalg(iv).eq.0) then
                 yldot(iv) = ( yldot(iv)*nnorm -
     .                                       yl(iv)*nbedot ) / ne(ix,iy)
               endif
	       if(istionxy(ix,iy).eq.1 .and. iseqalg(iv1).eq.0) then
                 if(isupgon(1)==1) then
                   yldot(iv1) = ( yldot(iv1)*nnorm -
     .                              yl(iv1)*(nbidot + nbgdot) ) /
     .                               (nit(ix,iy) + ni(ix,iy,2) )
                 else
                   yldot(iv1) = ( yldot(iv1)*nnorm - yl(iv1)*nbidot ) /
     .                            (nit(ix,iy) + cngtgx(1)*ng(ix,iy,1))
                 endif
               endif
               do igsp = 1, ngsp
                 if(istgonxy(ix,iy,igsp) == 1) then
                   iv = idxtg(ix,iy,igsp)
                   yldot(iv) = ( yldot(iv)*n0g(igsp) -
     .                            yl(iv)*nbg2dot(igsp) )/ng(ix,iy,igsp)
                 endif
              enddo
            endif
  260    continue
 270  continue

      return
      end

c   -------------------------------------------------------------------------
      subroutine volsor

*     VOLSOR defines volume particle and power sources for the ions and
*     electrons. Input is total current, power, and shape of 2-D Gaussians

            use Dim       ! nx,nxm,ny,nisp
            use Math_problem_size    ! neqmx for arrays not used here 
            use Comgeo    ! vol
            use RZ_grid_info    ! rm,zm
            use Volsrc    ! pwrsore,pwrsori,volpsor,pvole,pvoli,ivolcur,
                    ! z0pe,z0pi,r0pe,r0pi,zwpe,zwpi,rwpe,rwpi,
                    ! z0ni,r0ni,zwni,rwni,voljcsor,jcvsor,
                    ! ix_sjcsor, ix_ejcsor, iy_sjcsor, iy_ejcsor,
                    ! thetarot,rcutmin,zcutmin,effvng
            use Phyvar    ! ev
            use Bcond     ! islimsor,rlimiter
            use Parallv   ! nxg,nyg
            use Xpoint_indices   ! ixpt1,ixpt2,iysptrx
            use Share     ! nxomit

      implicit none

*  -- local scalars --
      real effvni, effvup, effvpe, effvpi, effvjel, zc, rc, ivolcurt,
     &     ivolcurgt, mvolcurt
      real argr, argz
      integer isxjcsor,iexjcsor,isyjcsor,ieyjcsor, ifld, nj,ix,iy,igsp


c...  Initialize values and arrays
      nj = nxomit
      effvni = 0.
      effvup = 0.
      effvpe = 0.
      effvpi = 0.
      effvjel = 0.
      ivolcurt = 0.
      mvolcurt = 0.
      ivolcurgt = 0.
      do igsp = 1, ngsp
        effvng(igsp) = 0.
      enddo
      call s2fill (nx+2, ny+2, 0., pwrsori, 1, nx+2)
      call s2fill (nx+2, ny+2, 0., pwrsore, 1, nx+2)
      call s2fill (nx+2, ny+2, 0., voljcsor, 1, nx+2)
      do ifld = 1, nisp
        call s2fill (nx+2, ny+2, 0., volpsor(0,0,ifld), 1, nx+2)
        call s2fill (nx+2, ny+2, 0., volmsor(0,0,ifld), 1, nx+2)
        ivolcurt = ivolcurt + ivolcur(ifld)
        mvolcurt = mvolcurt + mvolcur(ifld)
      enddo
      do igsp = 1, ngsp
        call s2fill (nx+2, ny+2, 0., volpsor(0,0,igsp), 1, nx+2)
        call s2fill (nx+2, ny+2, 0., volmsor(0,0,igsp), 1, nx+2)
        ivolcurgt = ivolcurgt + ivolcurg(igsp)
      enddo
cccMER NOTE: generalize the following for multiple x-points
cc Define index ranges for a localized ion-loss sink; crude & temporary
      if (ix_sjcsor .gt. 0) then
        isxjcsor = ix_sjcsor
      else
        isxjcsor = (ixpt1(1) + ixpt2(1))/2
      endif
      if (ix_ejcsor .gt. 0) then
        iexjcsor = ix_ejcsor
      else
        iexjcsor = ixpt2(1)
      endif
      if (iy_sjcsor > 0) then
        isyjcsor = iy_sjcsor
      else
        isyjcsor = 1
      endif
      if (iy_ejcsor > 0) then
        ieyjcsor = iy_ejcsor
      else
        ieyjcsor = iysptrx
      endif

      if ( abs(pvoli+pvole+ivolcurt+mvolcurt+ivolcurgt+jcvsor)
     .                      > 0. ) then   ! skip, big jump to end

      do 20 iy = 0, ny+1
         if (rm(0+nj,iy,0).lt.rlimiter .or. islimsor.eq.1) then
         do 10 ix = 0, nx+1
            zc = z0ni - (rm(ix+nj,iy,0)-r0ni)*sin(thetarot) +
     .                  (zm(ix+nj,iy,0)-z0ni)*cos(thetarot)
            rc = r0ni + (rm(ix+nj,iy,0)-r0ni)*cos(thetarot) +
     .                  (zm(ix+nj,iy,0)-z0ni)*sin(thetarot)
            if (zc.lt.zcutmin .or. rc.lt.rcutmin) goto 10
             argz = min(25., ((zc-z0ni)/zwni)**2)
             argr = min(25., ((rc-r0ni)/rwni)**2)
            effvni = effvni + vol(ix,iy) * exp(-argz -argr)
             argz = min(25., ((zc-z0up)/zwup)**2)
             argr = min(25., ((rc-r0up)/rwup)**2)
            effvup = effvup + vol(ix,iy) * exp(-argz -argr)
             argz = min(25., ((zc-z0pe)/zwpe)**2)
             argr = min(25., ((rc-r0pe)/rwpe)**2)
            effvpe = effvpe + vol(ix,iy) * exp(-argz -argr)
             argz = min(25., ((zc-z0pi)/zwpi)**2)
             argr = min(25., ((rc-r0pi)/rwpi)**2)
            effvpi = effvpi + vol(ix,iy) * exp(-argz -argr)
            do igsp = 1, ngsp
               argz = min(25., ((zc-z0ng(igsp))/zwng(igsp))**2)
               argr = min(25., ((rc-r0ng(igsp))/rwng(igsp))**2)
              effvng(igsp) = effvng(igsp) + vol(ix,iy)*exp(-argz-argr)
            enddo
cccMER For full double-null configuration, iysptrx is last closed flux surface.
cc  Temporary localized current source (for prompt ion loss)
            if (iy >= isyjcsor .and. iy <= ieyjcsor) then
              if (ix .ge. isxjcsor .and. ix .le. iexjcsor) then
                effvjel = effvjel + vol(ix,iy)
              endif
            endif
 10      continue
         endif
 20   continue

      do 40 iy = 0, ny+1
         if (rm(0+nj,iy,0).lt.rlimiter .or. islimsor.eq.1) then
         do 30 ix = 0, nx+1
            zc = z0ni - (rm(ix+nj,iy,0)-r0ni)*sin(thetarot) +
     .                  (zm(ix+nj,iy,0)-z0ni)*cos(thetarot)
            rc = r0ni + (rm(ix+nj,iy,0)-r0ni)*cos(thetarot) +
     .                  (zm(ix+nj,iy,0)-z0ni)*sin(thetarot)
            if (zc.lt.zcutmin .or. rc.lt.rcutmin) goto 30
             argz = min(25., ((zc-z0pe)/zwpe)**2)
             argr = min(25., ((rc-r0pe)/rwpe)**2)
            pwrsore(ix,iy) =pvole*vol(ix,iy)*exp(-argz-argr)/effvpe
             argz = min(25., ((zc-z0pi)/zwpi)**2)
             argr = min(25., ((rc-r0pi)/rwpi)**2)
            pwrsori(ix,iy) =pvoli*vol(ix,iy)*exp(-argz-argr)/effvpi
            do ifld = 1, nisp
               argz = min(25., ((zc-z0ni)/zwni)**2)
               argr = min(25., ((rc-r0ni)/rwni)**2)
              volpsor(ix,iy,ifld) = ivolcur(ifld) * vol(ix,iy) *
     .                                   exp(-argz-argr)/(effvni*ev)
               argz = min(25., ((zc-z0up)/zwup)**2)
               argr = min(25., ((rc-r0up)/rwup)**2)
              volmsor(ix,iy,ifld) = mvolcur(ifld) * vol(ix,iy) *
     .                                   exp(-argz-argr)/(effvup*ev)
            enddo
            do igsp = 1, ngsp
               argz = min(25., ((zc-z0ng(igsp))/zwng(igsp))**2)
               argr = min(25., ((rc-r0ng(igsp))/rwng(igsp))**2)
              volpsorg(ix,iy,igsp) = ivolcurg(igsp) * vol(ix,iy) *
     .                                 exp(-argz-argr)/(effvng(igsp)*ev)
            enddo
cccMER For full double-null configuration, iysptrx is last closed flux surface
cc  Temporary localized current source (for prompt ion loss)
            if (iy >= isyjcsor .and. iy <= ieyjcsor) then
              if (ix .ge. isxjcsor .and. ix .le. iexjcsor) then
                voljcsor(ix,iy) = jcvsor*vol(ix,iy) / effvjel
              endif
            endif
 30      continue
         endif
 40   continue

      endif  ! end of large if beginning if (abs(pvoli+pvole+...

      return
      end
c-----------------------------------------------------------------------
      subroutine resid (t, yl, yldot, cj, delta, ires, rpar, ipar)

c ... Calculate residuals for differential-algebraic solution of the
c     UEDGE equations.
c
c     Note that the residuals can be expressed in terms of f(i) =
c     right-hand sides of rate equations used in differential
c     solution of the UEDGE equations.
c
c     The residuals have the form
c        delta(i) = f(i)
c     for the algebraic equations, and
c        delta(i) = f(i) - yldot(i)
c     for the differential equations.

      implicit none

c ... Input arguments:
      real t          ! physical time
      real yl(*)      ! most recent iterate of solution vector
      real yldot(*)   ! most recent iterate of time-derivative of yl
      real cj         ! proportional to 1/del_t; can rescale alg. constraint eq.
      real rpar(*)    ! real parameter communication (not used)
      integer ipar(*) ! integer parameter communication

c ... In-out argument:
      integer ires    ! 0 on entry, set to -1 or -2 if errors occur

c ... Output argument:
      real delta(*)   ! residuals

c ... Function:
      integer res_algeb

c ... Local variables:
      integer neq, i, ifail

c ... Get total number of equations (all grid points).
      neq = ipar(1)

c ... Calculate f(i), storing them in delta(i).
      call rhsdpk (neq, t, yl, delta, ifail)
      if (ifail .ne. 0) then
         ires = -1
         return
      endif

c ... Loop through all i, skipping those that correspond to algebraic
c     equations, and subtract yldot(i) from delta(i).
      do i = 1, neq
         if (res_algeb (i) .eq. 0) delta(i) = delta(i) - yldot(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ffun (neq, t, yl, yldot)

c ... Calculate the right-hand sides of the UEDGE rate equations.

            use Math_problem_size    ! neqmx
            use Constraints          ! icflag,rlx,icnstr,ylprevc,ylchng

c ... Input arguments used for all entry points:
      implicit none

      integer neq
      real yl(*)   ! most recent iterate of solution vector

c ... Additional input argument required for ffun and rhsvd:
      real t         ! physical time

c ... Additional input arguments required for rhsvd:
      real rpar(*)
      integer ipar(*)

c ... Output argument used for all entry points:
      real yldot(neq)   ! right-hand sides

c ... Output argument for rhsvd and rhsdpk:
      integer ifail

c ... Local variables:
      integer i, ivar
      real tloc, tau, rlxl
      data tau /1.e0/, rlxl /1.e20/   ! dummy argument here for cnstrt

c ... Beginning of execution for call ffun (by newton or lsode).
      goto 8

c ... Beginning of execution for call rhsvd (by vodpk), check constraints
      entry rhsvd (neq, t, yl, yldot, rpar, ipar, ifail)

      if (icflag .gt. 0) then
         if (icflag .eq. 2) rlxl = rlx
         do 5 i = 1, neq
            ylchng(i) = yl(i) - ylprevc(i)
 5       continue
         call cnstrt (neq,ylprevc,ylchng,icnstr,tau,rlxl,ifail,ivar)
         if (ifail .ne. 0) then
            call remark ('***Constraint failure in VODPK, dt reduced***')
            write (*,*) 'variable index = ',ivar,'   time = ',t
            goto 20
         endif
cJG         ylprevc(1:neq)=yl(1:neq)
         call scopy (neq, yl, 1, ylprevc, 1)  !put yl into ylprevc
      else
         ifail = 0
      endif

      go to 8

c ... Beginning of execution for call rhsdpk (by daspk), check constraints
      entry rhsdpk (neq, t, yl, yldot, ifail)

      if (icflag .gt. 0 .and. t .gt. 0.) then
         if (icflag .eq. 2) rlxl = rlx
         do 6 i = 1, neq
            ylchng(i) = yl(i) - ylprevc(i)
 6       continue
         call cnstrt (neq,ylprevc,ylchng,icnstr,tau,rlxl,ifail,ivar)
         if (ifail .ne. 0) then
            call remark ('***Constraint failure in DASPK, dt reduced***')
            write (*,*) 'variable index = ',ivar,'   time = ',t
            goto 20
         endif
      else
         ifail = 0
      endif
cJG      ylprevc(1:neq)=yl(1:neq)
      call scopy (neq, yl, 1, ylprevc, 1)  !put yl into ylprevc

 8    tloc = t
      go to 10

c ... Beginning of execution for call rhsnk (by nksol).
      entry rhsnk (neq, yl, yldot)
      tloc = 0.

c ... Calculate right-hand sides for interior and boundary points.
ccc 10   call convsr_vo (-1,-1, yl)  ! test new convsr placement
ccc      call convsr_aux (-1,-1, yl) ! test new convsr placement
 10   call pandf1 (-1, -1, 0, neq, tloc, yl, yldot)

 20   continue
      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine jac_calc (neq, t, yl, yldot00, ml, mu, wk,
     .                     nnzmx, jac, ja, ia)
c ... Calculate Jacobian matrix (derivatives with respect to each
c     dependent variable of the right-hand side of each rate equation).
c     Lower and upper bandwidths are used to select for computation
c     only those Jacobian elements that may be nonzero.
c     Estimates of Jacobian elements are computed by finite differences.
c     The Jacobian is stored in compressed sparse row format.
c ... Common blocks:
            use Dim                      ! nx,ny,
                                   ! nusp[for fnorm not used here]
            use Timing                   ! istimingon,ttjstor,ttotjf,ttimpjf
            use Math_problem_size        ! neqmx,numvar
            use Grid                     ! ngrid,ig,ijac,ijactot
            use Indexes                  ! igyl,iseqalg
            use Variable_perturbation    ! del,dylconst
            use Jacobian_clipping        ! jaccliplim,istopjac,irstop,icstop
            use Jacobian_csc             ! rcsc,jcsc,icsc,yldot_pert
            use Ynorm                    ! suscal,sfscal
            use UEpar                    ! isphion,isnewpot,svrpkg,isbcwdt
            use Model_choice             ! iondenseqn
            use Imprad                   ! isimpon
            use Bcond                    ! isextrnpf,isextrtpf,isextrngc,
                                   ! isextrnw,isextrtw
            use Parallv                  ! nxg,nyg
            use Time_dep_nwt             ! nufak,dtreal,ylodt,dtuse
            use Selec                    ! yinc
c ... Functions:
      implicit none
c ... Input arguments:
      integer neq      ! total number of equations (all grid points)
      real t           ! physical time
      real yl(*)       ! dependent variables
      real yldot00(neq+2) ! right-hand sides evaluated at yl
      integer ml, mu   ! lower and upper bandwidths
      integer nnzmx    ! maximum number of nonzeros in Jacobian
c ... Work-array argument:
      real wk(neq)     ! work space available to this subroutine
c ... Output arguments:
      real jac(nnzmx)     ! nonzero Jacobian elements
      integer ja(nnzmx)   ! col indices of nonzero Jacobian elements
      integer ia(neq+1)   ! pointers to beginning of each row in jac,ja
      intrinsic isnan
      logical tstguardc
      real(kind=4):: gettime
cc      real(kind=4):: ranf
c ... Local variables:
      integer nnz, ii, iv, ii1, ii2, xc, yc, ix, iy
      real yold, dyl, jacelem
      real(kind=4):: sec4, tsjstor, tsimpjf, dtimpjf,TimeBuild
ccc      save
c ... Get initial value of system cpu timer.
      if (istimingon .eq. 1) tsjstor = gettime(sec4)
c ... Pause from BASIS if a ctrl_c is typed
c     call ruthere
c ... Count Jacobian evaluations, both for total and for this case
      ijactot = ijactot + 1
      ijac(ig) = ijac(ig) + 1
      if (svrpkg.eq.'nksol') write(*,*) ' Updating Jacobian, npe =  ',
     .                                                          ijac(ig)
       write(*,*) 'ml,mu,neq',ml, mu,neq
c ... Set up diagnostic arrays for debugging
      do iv = 1, neq
        yldot_unpt(iv) = yldot00(iv)  ! for diagnostic only
        yldot_pert(iv) = 0.
      enddo
      wk(1:neq)=0
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c ... Begin loop over dependent variables.
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      TimeBuild=gettime(sec4)
      nnz = 1
      do iv = 1, neq
ccc ... Only perturb variables that are being solved for (for Daspk option)
ccc      if (iseqon(iv) .eq. 0) goto 18
c ... Set beginning and ending indices of right-hand sides that might be
c     perturbed.
         ii1 = max(iv-mu, 1)
         ii2 = min(iv+ml, neq)
c ... Reset range if this is a potential perturbation with isnewpot=1
ccc         if (isphion*isnewpot.eq.1 .and. mod(iv,numvar).eq.0) then
         if (isphion*isnewpot.eq.1) then
            ii1 = max(iv-4*numvar*nx, 1)      ! 3*nx may be excessive
            ii2 = min(iv+4*numvar*nx, neq)    ! 3*nx may be excessive
         endif
c ... Reset range if extrapolation boundary conditions are used
         if (isextrnpf+isextrtpf+isextrngc+isextrnw+isextrtw.gt.0) then
            ii1 = max(iv-2*numvar*(nx+3), 1)      ! guess to include extrap. bc
            ii2 = min(iv+2*numvar*(nx+3), neq)    ! guess to include extrap. bc
         endif
c ... Initialize all of those right-hand sides to their unperturbed
c     values.
         do ii = ii1, ii2   ! below wk is reset, but only over limited range
            wk(ii) = yldot00(ii)
         enddo
c ... Set spatial-location indices for this dependent variable.
         xc = igyl(iv,1)
         yc = igyl(iv,2)
c ... Save value of dependent variable, then perturb it.
c     The perturbation to the variable is proportional to parameter
c     del and to a measure of the size of the variable.  That measure
c     increases with the absolute value of the variable if it exceeds
c     the typical size given by dylconst/suscal but can never be less
c     than that typical size.
         yold = yl(iv)
         dyl = delperturb * (abs(yold) + dylconst / suscal(iv))
         yl(iv) = yold + dyl
c ... Calculate right-hand sides near location of perturbed variable.
ccc         call convsr_vo (xc, yc, yl)  ! test new convsr placement
ccc         call convsr_aux (xc, yc, yl) ! test new convsr placement
         call pandf1 (xc, yc, iv, neq, t, yl, wk)
ccc         yl(iv) = yold - dyl    ! for 2nd order Jac
ccc         call pandf1 (xc, yc, iv, neq, t, yl, yldot0) ! for 2nd order Jac
c ... Calculate possibly nonzero Jacobian elements for this variable,
c     and store nonzero elements in compressed sparse column format.
         jcsc(iv) = nnz      ! sets index for first Jac. elem. of var. iv
         do ii = ii1, ii2
ccc            if (iondenseqn .eq. "inel") then   ! skip ion cont. eqn
ccc               ix = igyl(ii,1)               ! (but not for guard cells)
ccc               iy = igyl(ii,2)
ccc               if (isnion(1) .eq. 1 .and. mod(ii,numvar) .eq. 1 .and.
ccc     .             .not. tstguardc (ix, iy)) go to 10  !end of iv loop removed
ccc            endif
            jacelem = (wk(ii) - yldot00(ii)) / dyl
ccc            jacelem = (wk(ii) - yldot0(ii)) / (2*dyl)  ! for 2nd order Jac
c ...  Add diagonal 1/dt for nksol
            if (((svrpkg.eq."nksol") .or. (svrpkg.eq."petsc")) .and. iv.eq.ii) then
              if (iseqalg(iv)*(1-isbcwdt).eq.0) then
                jacelem = jacelem - 1/dtuse(iv)
              endif
              ix = igyl(iv,1)
              iy = igyl(iv,2)
              if (idxphi(ix,iy)==iv .and. dtphi<1e10) then !selects phi eqn
                jacelem = jacelem - 1/dtphi
              endif
            endif
c ...  Add a pseudo timestep to the diagonal !! if eqn is not algebraic
            if (svrpkg .ne. "cvode" .and. nufak .gt. 0) then
               if (iv.eq.ii .and. yl(neq+1).eq.1)
     .             jacelem = jacelem - nufak  !omit .and. iseqalg(iv).eq.0)
cc     .                   jacelem = jacelem - nufak*suscal(iv)/sfscal(iv)
            endif
            if (abs(jacelem*sfscal(iv)) .gt. jaccliplim) then
               if (nnz .gt. nnzmx) then
                  write(STDOUT,*)
     .             '*** jac_calc -- More storage needed for Jacobian.',
     .             ' Storage exceeded at (i,j) = (',ii,',',iv,').',
     .             ' Increase lenpfac.'
              call xerrab("")
               endif
cc               if (rdoff.ne.0.e0) jacelem=jacelem*(1.0e0+ranf()*rdoff)
               rcsc(nnz) = jacelem
               icsc(nnz) = ii
               nnz = nnz + 1
            endif

            if (istopjac.gt.0 .and. ii.eq.irstop .and. iv.eq.icstop) then
               yldot_pert(ii) = wk(ii)      ! for diagnostic only
               if (istopjac == 2) then
                 yl(iv) = yold
                 call pandf1 (xc, yc, iv, neq, t, yl, wk)
               endif
ccc               call convsr_vo (xc, yc, yl)  ! was one call to convsr
ccc               call convsr_aux (xc, yc, yl)
               call remark("***** non-zero jac_elem at irstop,icstop")
               write(*,*) 'irstop = ', irstop, ', icstop = ', icstop
               call xerrab("")
            endif
         enddo  ! end of ii loop over equations
c ... Restore dependent variable and plasma variables near its location.
         yl(iv) = yold
ccc         call convsr_vo (xc, yc, yl)  ! was one call to convsr
ccc         call convsr_aux (xc, yc, yl)
         call pandf1 (xc, yc, iv, neq, t, yl, wk)
         yldot_pert(1:neq)=wk(1:neq)
ccc 18   continue
c...  If this is the last variable before jumping to new cell, reset pandf
         if (mod(iv,numvar).eq.0 .and. isjacreset.ge.1) then
c            call DebugHelper('dump0.txt')
            call pandf1 (xc, yc, iv, neq, t, yl, wk)
         endif
         do ii=1,neq
         if (yldot_pert(ii).ne.wk(ii)) then
      write(*,'(a,i5,e20.12,e20.12)') ' *** wk modified on second call at ii=',
     . ii,yldot_pert(ii),wk(ii)
c         call DebugHelper('dump1.txt')
         call xerrab('*** Stop ***')
         if (isnan(yldot_pert(ii))) then
         call xerrab('*** Stop ***')
         endif
         endif
         enddo


c ... End loop over dependent variables and finish Jacobian storage.
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo             ! end of main loop over yl variables
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      TimeBuild=gettime(sec4)-TimeBuild
      write(*,*)'Time to build jac:',TimeBuild
      jcsc(neq+1) = nnz
       write(*,*) 'nnzcum=',nnz
c ... Convert Jacobian from compressed sparse column to compressed
c     sparse row format.
      call csrcsc (neq, 1, 1, rcsc, icsc, jcsc, jac, ja, ia)
cccc ... If desired, calculate Jacobian elements for ion continuity
cccc     equation by using INEL routines, and combine elements with
cccc     those calculated above.
ccc      if (iondenseqn .eq. "inel") then
ccc         call iondens2 (neq, nnzmx, jac, ja, ia, rcsc, icsc, jcsc)
ccc      endif
cccc ... If necessary, calculate Jacobian elements for rows corresponding
cccc     to impurity-density equations (for cells other than guard cells),
cccc     and combine elements with those calculated above.
ccc      if (isimpon .eq. 3 .or. isimpon .eq. 4) then
ccc         if (istimingon .eq. 1) tsimpjf = gettime(sec4)
ccc         call impdens2 (neq, nnzmx, jac, ja, ia, rcsc, icsc, jcsc)
ccc         if (istimingon .eq. 1) then
ccc            dtimpjf = gettime(sec4) - tsimpjf
ccc            ttimpjf = ttimpjf + dtimpjf
ccc            ttotjf = ttotjf + dtimpjf
ccc         endif
ccc      endif
c ... Accumulate cpu time spent here.
      if (istimingon .eq. 1) ttjstor = ttjstor + gettime(sec4) - tsjstor
      return
      end
c-----------------------------------------------------------------------
      subroutine jac_calc_para (neq, t, yl, yldot00, ml, mu, wk,
     .                     nnzmx, jac, ja, ia)

c ... Calculate Jacobian matrix (derivatives with respect to each
c     dependent variable of the right-hand side of each rate equation).
c     Lower and upper bandwidths are used to select for computation
c     only those Jacobian elements that may be nonzero.
c     Estimates of Jacobian elements are computed by finite differences.
c     The Jacobian is stored in compressed sparse row format.

c ... Common blocks:
            use Dim                      ! nx,ny,
                                   ! nusp[for fnorm not used here]
            use Timing                   ! istimingon,ttjstor,ttotjf,ttimpjf
            use Math_problem_size        ! neqmx,numvar
            use Grid                     ! ngrid,ig,ijac,ijactot
            use Indexes                  ! igyl,iseqalg
            use Variable_perturbation    ! del,dylconst
            use Jacobian_csc             ! rcsc,jcsc,icsc,yldot_pert
            use Ynorm                    ! suscal,sfscal
            use UEpar                    ! isphion,isnewpot,svrpkg,isbcwdt
            use Model_choice             ! iondenseqn
            use Imprad                   ! isimpon
            use Bcond                    ! isextrnpf,isextrtpf,isextrngc,
                                   ! isextrnw,isextrtw
            use Parallv                  ! nxg,nyg
            use Time_dep_nwt             ! nufak,dtreal,ylodt,dtuse
            use Selec                    ! yinc
            use Preconditioning          !lenpfac
            use OmpOptions 
c ... Functions
      implicit none

c ... Input arguments:
      integer neq      ! total number of equations (all grid points)
      real t           ! physical time
      real yl(*)       ! dependent variables
      real yldot00(neq+2) ! right-hand sides evaluated at yl
      integer ml, mu   ! lower and upper bandwidths
      integer nnzmx    ! maximum number of nonzeros in Jacobian


c ... Work-array argument:
      real wk(neq)     ! work space available to this subroutine
      integer,allocatable :: iJacConstructor(:,:)
      real,allocatable:: rJacConstructor(:,:)


c ... Output arguments:
      real jac(nnzmx)     ! nonzero Jacobian elements
      integer ja(nnzmx)   ! col indices of nonzero Jacobian elements
      integer ia(neq+1)   ! pointers to beginning of each row in jac,ja

      logical tstguardc
      real(kind=4):: gettime
cc      real(kind=4):: ranf

c ... Local variables:
      integer i,Nsize,R,thread
      integer,parameter:: Nthreadmax=256
      integer,allocatable::nnz(:),nnzcum(:)
      integer,allocatable,dimension(:,:) ::iJacCol,iJacRow
      real,allocatable,dimension(:,:) ::rJacElem
      integer ivmax(Nthreadmax),ivmin(Nthreadmax)
      real(kind=4):: sec4, tsjstor, tsimpjf, dtimpjf,time0,time1,TimeBuild
      real(kind=4):: TimeReduce
      INTEGER :: NTHREADS, TID, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
      integer:: ith,iv,neq_,nnzmxperthread
      character(len = 80) ::  filename
cJG First we get the number of threads
c$omp PARALLEL PRIVATE(TID)
      TID = OMP_GET_THREAD_NUM()
      IF (TID .EQ. 0) THEN
        NTHREADS = OMP_GET_NUM_THREADS()
      END IF
c$omp END PARALLEL

cJG Then we obtain the range of iv index per thread
       neq_=neq
       Nsize=neq_/Nthreads
       R=MOD(neq_, Nthreads)
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
                ivmax(Nthreads)=neq_
            endif
        else
            ivmin(Nthreads)=1
            ivmax(Nthreads)=neq_
        endif
*        do ith=1,Nthreads
*        totii(ith)=0
*        do iv=ivmin(ith),ivmax(ith)
*
*         ii1 = max(iv-mu, 1)
*         ii2 = min(iv+ml, neq)
*c ... Reset range if this is a potential perturbation with isnewpot=1
*ccc         if (isphion*isnewpot.eq.1 .and. mod(iv,numvar).eq.0) then
*         if (isphion*isnewpot.eq.1) then
*            ii1 = max(iv-4*numvar*nx, 1)      ! 3*nx may be excessive
*            ii2 = min(iv+4*numvar*nx, neq)    ! 3*nx may be excessive
*         endif
*c ... Reset range if extrapolation boundary conditions are used
*         if (isextrnpf+isextrtpf+isextrngc+isextrnw+isextrtw.gt.0) then
*            ii1 = max(iv-2*numvar*(nx+3), 1)      ! guess to include extrap. bc
*            ii2 = min(iv+2*numvar*(nx+3), neq)    ! guess to include extrap. bc
*         endif
*         totii(ith)=totii(ith)+ii2-ii1+1
*         enddo
*         enddo

        if (OMPDebug.gt.0) then
        write(*,*)' Ivmin,Ivmax:'
            do i=1,Nthreads
                write(*,*) ivmin(i),ivmax(i)
            enddo
        endif

      allocate(iJacRow(neq,Nthreads))
      if (Nthreads.gt.1) then
      nnzmxperthread=ceiling(real(nnzmx)/real(Nthreads-1))
      else
      nnzmxperthread=nnzmx
      endif
      write(*,*) 'Allocating size:',neq,nnzmxperthread,nnzmxperthread*8*Nthreads/1e9,'Gb'
      allocate(iJacCol(nnzmxperthread,Nthreads))
      allocate(rJacElem(nnzmxperthread,Nthreads))
      allocate(nnz(Nthreads))
      allocate(nnzcum(Nthreads))
      nnz(1:Nthreads)=-1
      nnzcum(1:Nthreads)=-1
      write(*,*) 'Allocatation done....'
      iJacCol(1:nnzmxperthread,1:Nthreads)=0
      write(*,*) 'Inialization iJacCol done....'
      rJacElem(1:nnzmxperthread,1:Nthreads)=0.0
      write(*,*) 'Inialization rJacElem done....'
      iJacRow(1:neq,1:Nthreads)=0
      write(*,*) 'Inialization iJacRow done....'
      nnzmxperthread=nnzmxperthread



ccc      save

c ... Get initial value of system cpu timer.
      if (istimingon .eq. 1) tsjstor = gettime(sec4)

c ... Pause from BASIS if a ctrl_c is typed
c     call ruthere

c ... Count Jacobian evaluations, both for total and for this case
      ijactot = ijactot + 1
      ijac(ig) = ijac(ig) + 1
      if (svrpkg.eq.'nksol') write(*,*) ' Updating Jacobian, npe =  ',
     .                                                          ijac(ig)

c ... Set up diagnostic arrays for debugging
      do iv = 1, neq
        yldot_unpt(iv) = yldot00(iv)  ! for diagnostic only
        yldot_pert(iv) = 0.
      enddo

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c ... Begin loop over dependent variables.
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cJG Running with  openmp target requires to have -fopen flag and foffload flag added to both gfortran and gcc.
c JG Let test it with a gpu (maybe summit 1 at GA?)
      write(*,*) 'OMP construction of Jacobian with Nthreads=',Nthreads
      TimeBuild=gettime(sec4)
        call OMPJacBuilder(ivmin,ivmax,neq, t, yl,yldot00, ml, mu,wk,
     .                      iJacCol,rJacElem,iJacRow,nnz,nnzmxperthread,
     .                      nthreads)

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      TimeBuild=gettime(sec4)-TimeBuild
      write(*,*)'Time to build jac:',TimeBuild
      if (DebugTime.gt.0) then
      write(*,*)'Time in var loop calc_jac:',gettime(sec4)-tsjstor
      endif
c !reduction of jacobian data into single vector
      TimeReduce=gettime(sec4)
      nnzcum(1)=nnz(1)-1

      do ith=2,Nthreads
        nnzcum(ith)=nnzcum(ith-1)+nnz(ith)-1
      enddo
      if (OMPDebug.gt.0) then
      write(*,*) 'nnzcum=',nnzcum
      endif
      if (nnzcum(Nthreads).gt.nnzmx) then
        write(*,*),'Houston we have a problem with nnz...'
      endif
c     !jcsc(iv) = nnz
*      if (OMPDebug.gt.0) then
*      write(*,*) 'size iJacRow=',size(iJacRow)
*      endif
*       write(*,*) '**** iJacRow *****'
*       do iv=1,40
*            write(*,*) iv,(iJacRow(iv,ith),ith=1,Nthreads)
*       enddo
*       write(*,*) '*******************'
*        jcsc(1:neq+1)=0
*        rcsc(1:nnzmx)=0.0
*        icsc(1:nnzmx)=0.0

        jcsc(ivmin(1):ivmax(1))= iJacRow(ivmin(1):ivmax(1),1)
      do ith=2,Nthreads
        jcsc(ivmin(ith):ivmax(ith))= iJacRow(ivmin(ith):ivmax(ith),ith)
     .                                +nnzcum(ith-1)
      enddo

c     !rcsc(nnz)= JacElem
*      if (OMPDebug.gt.0) then
*      write(*,*) 'size rJacElem=',size(rJacElem)
*      write(*,*) 'size rcsc=',size(rcsc)
*      write(*,*) 'nnz(1)=',nnz(1)
*      endif
      rcsc(1:nnz(1)-1)= rJacElem(1:nnz(1)-1,1)
      icsc(1:nnz(1)-1)= iJacCol(1:nnz(1)-1,1)
      do ith=2,Nthreads
       rcsc(nnzcum(ith-1)+1:nnzcum(ith))=
     .                            rJacElem(1:nnz(ith)-1,ith)
       icsc(nnzcum(ith-1)+1:nnzcum(ith))=
     .                            iJacCol(1:nnz(ith)-1,ith)
      enddo

*      write (filename, fmt='(a,i2.2)') 'datatest_',Nthreads
*      open (unit = 777, file = trim(filename))
*      do iv=1,nnzmx
*      if (iv.le.neq+1) then
*      write(777,*) rcsc(iv),icsc(iv),jcsc(iv)
*      else
*      write(777,*) rcsc(iv),icsc(iv),-1
*      endif
*      enddo
*
*      close(777)




      jcsc(neq+1) = nnzcum(Nthreads)+1
*      write(*,*) '********** jcsc **************'
*      do iv=1,40
*            write(*,*) iv,jcsc(iv)
*      enddo
      TimeReduce=gettime(sec4)-TimeReduce
      write(*,*)'**** Time to reduce jac:',TimeReduce
*      write(*,*) '********** rcsc,icsc **************'
*      do iv=1,200
*            write(*,*) iv,rcsc(iv),icsc(iv)
*      enddo
      deallocate(iJacRow)
      deallocate(iJacCol)
      deallocate(rJacElem)
      deallocate(nnz)
      deallocate(nnzcum)

c ... Convert Jacobian from compressed sparse column to compressed
c     sparse row format.
      time1=gettime(sec4)

      write(*,*) 'size rcsc=',size(rcsc),size(jac),size(icsc),size(ja),size(jcsc),size(ia)
*      write(*,*) 'neq=',neq
*      write(*,*)  'nnz(thread)=',nnz(Nthreads)
*      write(*,*) 'jcsc(neq+1)',jcsc(neq+1)
*      write(*,*) '********** rcsc,icsc,jcsc'
**      do iv=1,40
**    write(*,*) rcsc(iv),icsc(iv),jcsc(iv),jac(iv),ja(iv),ia(iv)
**      enddo
      call csrcsc (neq, 1, 1, rcsc, icsc, jcsc, jac, ja, ia)
      if (DebugTime.gt.0) then
      write(*,*)'Time to convert jac:',gettime(sec4)-time1
      endif
*      write(*,*) '********** Jac,ja,ia **************'
*      do iv=1,40
*            write(*,*) iv,jac(iv),ja(iv),ia(iv)
*      enddo

cccc ... If desired, calculate Jacobian elements for ion continuity
cccc     equation by using INEL routines, and combine elements with
cccc     those calculated above.
ccc      if (iondenseqn .eq. "inel") then
ccc         call iondens2 (neq, nnzmx, jac, ja, ia, rcsc, icsc, jcsc)
ccc      endif

cccc ... If necessary, calculate Jacobian elements for rows corresponding
cccc     to impurity-density equations (for cells other than guard cells),
cccc     and combine elements with those calculated above.
ccc      if (isimpon .eq. 3 .or. isimpon .eq. 4) then
ccc         if (istimingon .eq. 1) tsimpjf = gettime(sec4)
ccc         call impdens2 (neq, nnzmx, jac, ja, ia, rcsc, icsc, jcsc)
ccc         if (istimingon .eq. 1) then
ccc            dtimpjf = gettime(sec4) - tsimpjf
ccc            ttimpjf = ttimpjf + dtimpjf
ccc            ttotjf = ttotjf + dtimpjf
ccc         endif
ccc      endif

c ... Accumulate cpu time spent here.
      if (istimingon .eq. 1) ttjstor = ttjstor + gettime(sec4) - tsjstor
      if (DebugTime.gt.0) then
      write(*,*)'Time in calc_jac:',gettime(sec4)-tsjstor
      endif
      return
      end subroutine

c-----------------------------------------------------------------------
      subroutine OMPJacBuilder(ivmin,ivmax,neq, t, yl,yldot00, ml,
     .  mu,wk,iJacCol,rJacElem,iJacRow,nnz,nnzmx,
     .                      nthreads)
              use OmpOptions 
              use OmpCopybbb 
c        Use(OmpCopygrd)
              use OmpCopycom 
c        Use(OmpCopyflx)
              use OmpCopyapi 

              use Compla 
              use Dim 
              use Selec 
              use Math_problem_size 
              use MCN_dim 
              use Indices_domain_dcl 
              use Reduced_ion_interface 
        implicit none
          integer,intent(in):: nnzmx,nthreads
          integer,intent(inout)::nnz(Nthreads)
          integer,intent(in):: ivmin(Nthreads),ivmax(Nthreads)      ! index equation
          integer,intent(in):: neq      ! total number of equations (all grid points)
          real,intent(in):: t           ! physical time
          real ::yl(*)       ! dependent variables
          real,intent(inout) :: wk(neq)
          real ::wkcopy(neq)
          real::ylcopy(neq+1)
       integer,intent(inout)::iJacCol(nnzmx,nthreads)
       integer,intent(inout):: iJacRow(neq,nthreads)
          real,intent(inout):: rJacElem(nnzmx,nthreads)
       integer ::iJacColCopy(nnzmx)
       integer :: iJacRowCopy(neq)
          real :: rJacElemCopy(nnzmx)
          real,intent(in) ::yldot00(neq+2) ! right-hand sides evaluated at yl
          integer,intent(in):: ml, mu   ! lower and upper bandwidths

c ... Work-array argument:


          integer:: ith,omp_get_thread_num,tid,nnzlocal,ithcopy




      if (OMPDebug.gt.0) then
      write(*,*) '**** Copying data....'
      endif


      call OmpCopyPointerbbb
      call OmpCopyScalarbbb
      call OmpCopyPointercom
      call OmpCopyScalarcom
c      call OmpCopyPointergrd
c      call OmpCopyScalargrd
      call OmpCopyPointerapi
      call OmpCopyScalarapi
c      call OmpCopyPointerflx
c      call OmpCopyScalarflx

      iJacColCopy(1:nnzmx)=0
      rJacElemCopy(1:nnzmx)=0.0
      iJacRowCopy(1:neq)=0
      ylcopy(1:neq+1)=yl(1:neq+1) ! a very barbarian use of yl(neq+1) is implemented as a switch in pandf... Error-prone!
      wkcopy(1:neq)=wk(1:neq)

      if (OMPDebug.gt.0) then
      write(*,*) '**** Starting parallel loop....'
      endif
      tid=-1
       nnzlocal=-10000
      write(*,*) 'nnzmx=',nnzmx
c$omp parallel do default(shared)
c$omp& firstprivate(ithcopy,ivmin,ivmax,tid,nnzlocal,ylcopy,wkcopy,ml,mu,yldot00,t,neq)
c$omp& firstprivate(nnzmx,nthreads)
c$omp& firstprivate(iJacRowCopy)
c$omp& shared(nnz,iJacCol,rJacElem,iJacRow)
c$omp& firstprivate(iJacColCopy)
c$omp& firstprivate(rJacElemCopy)
        do ith=1,Nthreads !ith from 1 to Nthread, tid from 0 to Nthread-1
        tid=omp_get_thread_num()
        ithcopy=ith
        write(*,*) '!OMP: Thread id:',tid,' <-> ith:',ithcopy
      call LocalJacBuilder(ivmin(ithcopy),ivmax(ithcopy),neq, t, ylcopy,yldot00,
     . ml,mu,wkcopy,iJacColcopy,rJacElemcopy,iJacRowcopy,ithcopy,nnzlocal,
     . nnzmx,nthreads)
      write(*,*) '!,',tid,' nzlocal:',nnzlocal
c$omp critical
      iJacCol(1:nnzlocal,ithcopy)=iJacColCopy(1:nnzlocal)
      rJacElem(1:nnzlocal,ithcopy)=rJacElemCopy(1:nnzlocal)
      iJacRow(1:neq,ithcopy)=iJacRowCopy(1:neq)
      nnz(ithcopy)=nnzlocal
c$omp end critical

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      enddo             ! end of main loop over yl variables
      enddo
c$omp END PARALLEL DO
      if (OMPDebug.gt.0) then
      write(*,*) '**** End of parallel loop....'
      endif



      end subroutine OMPJacBuilder

c ----------------------------------------------------------------------
      subroutine LocalJacBuilder(ivmin,ivmax,neq, t, yl,yldot00, ml, mu, wk,
     .                      iJacCol,rJacElem,iJacRow,ith,nnz,nnzmx,
     .                      nthreads)


c ... Calculate Jacobian matrix (derivatives with respect to each
c     dependent variable of the right-hand side of each rate equation).
c     Lower and upper bandwidths are used to select for computation
c     only those Jacobian elements that may be nonzero.
c     Estimates of Jacobian elements are computed by finite differences.
c     The Jacobian is stored in compressed sparse row format.

c ... Common blocks:
            use Dim                      ! nx,ny
            use Timing                   ! istimingon,ttjstor,ttotjf,ttimpjf
            use Math_problem_size        ! neqmx,numvar
            use Grid                     ! ngrid,ig,ijac,ijactot
            use Indexes                  ! igyl,iseqalg
            use Variable_perturbation    ! del,dylconst
            use Jacobian_clipping        ! jaccliplim,istopjac,irstop,icstop
            use Ynorm                    ! suscal,sfscal
            use UEpar                    ! isphion,isnewpot,svrpkg,isbcwdt
            use Model_choice             ! iondenseqn
            use Imprad                   ! isimpon
            use Bcond                    ! isextrnpf,isextrtpf,isextrngc,
                                   ! isextrnw,isextrtw
            use Parallv                  ! nxg,nyg
            use Time_dep_nwt             ! nufak,dtreal,ylodt,dtuse
            use Selec                    ! yinc
            use Preconditioning          !lenpfac
            use Jacobian_csc 
            use OmpOptions 
            use Compla 
c ... Functions:
      implicit none

      integer,intent(in):: ith,nnzmx,nthreads
      integer,intent(inout)::nnz
      integer,intent(in):: ivmin,ivmax      ! index equation
      integer,intent(in):: neq      ! total number of equations (all grid points)
      real,intent(in):: t           ! physical time
      real ::yl(*)       ! dependent variables
      integer,intent(inout)::iJacCol(nnzmx),iJacRow(neq)
      real,intent(inout):: rJacElem(nnzmx)
      real,intent(in) ::yldot00(neq+2) ! right-hand sides evaluated at yl
      integer,intent(in):: ml, mu   ! lower and upper bandwidths

c ... Work-array argument:
      real,intent(inout) :: wk(neq)     ! work space available to this subroutine
      logical ::tstguardc
      real(kind=4):: gettime
cc      real(kind=4):: ranf

c ... Local variables:
      integer :: ii, iv, jv,ii1, ii2, xc, yc, ix, iy

      real ::yold, dyl, jacelem
      real(kind=4):: sec4, tsjstor, tsimpjf, dtimpjf,time0
      integer:: ivcount
      integer:: omp_get_thread_num,tid
c      common /ilocals/ ii, ii1, ii2, xc, yc, ix, iy
c      common /rlocals/ yold,dyl,jacelem,time0
c       time0=gettime(sec4) !here



ccc ... Only perturb variables that are being solved for (for Daspk option)
ccc      if (iseqon(iv) .eq. 0) goto 18

c ... Set beginning and ending indices of right-hand sides that might be
c     perturbed.
       tid=omp_get_thread_num()
*       write(*,*) '!OMP:',tid, ': b2) te(1,1)',te(1,1)
*       write(*,*) '!OMP:',tid, ': b2) te(nx,ny)',te(nx,ny)
*       write(*,*) '!OMP:',tid, ': b2) j8',j8
         nnz=1
         ivcount=0
        do iv=ivmin,ivmax
           ivcount=ivcount+1
         ii1 = max(iv-mu, 1)
         ii2 = min(iv+ml, neq)
c ... Reset range if this is a potential perturbation with isnewpot=1
ccc         if (isphion*isnewpot.eq.1 .and. mod(iv,numvar).eq.0) then
         if (isphion*isnewpot.eq.1) then
            ii1 = max(iv-4*numvar*nx, 1)      ! 3*nx may be excessive
            ii2 = min(iv+4*numvar*nx, neq)    ! 3*nx may be excessive
         endif
c ... Reset range if extrapolation boundary conditions are used
         if (isextrnpf+isextrtpf+isextrngc+isextrnw+isextrtw.gt.0) then
            ii1 = max(iv-2*numvar*(nx+3), 1)      ! guess to include extrap. bc
            ii2 = min(iv+2*numvar*(nx+3), neq)    ! guess to include extrap. bc
         endif

c ... Initialize all of those right-hand sides to their unperturbed
c     values.

         do ii = ii1, ii2   ! a) below wk is reset, but only over limited range
            wk(ii) = yldot00(ii)
         enddo
*       if (ivcount.lt.3 .and. OMPDebug.gt.1) then
*       write(*,'(a,i3,a,i7,i7)')'!',ith,' : Got wk, looking for xc,yc..'
*       endif
c ... Set spatial-location indices for this dependent variable.
         xc = igyl(iv,1)
         yc = igyl(iv,2)
*      if (ivcount.lt.3 .and. OMPDebug.gt.1) then

*      endif
c ... Save value of dependent variable, then perturb it.
c     The perturbation to the variable is proportional to parameter
c     del and to a measure of the size of the variable.  That measure
c     increases with the absolute value of the variable if it exceeds
c     the typical size given by dylconst/suscal but can never be less
c     than that typical size.
         yold = yl(iv)
         dyl = delperturb * (abs(yold) + dylconst / suscal(iv))
         yl(iv) = yold + dyl
c         write(*,*) 'yl',yl(iv), 'from thread = ', TID
c ... Calculate right-hand sides near location of perturbed variable.
ccc         call convsr_vo (xc, yc, yl)  ! test new convsr placement
ccc         call convsr_aux (xc, yc, yl) ! test new convsr placement

         call pandf1 (xc, yc, iv, neq, t, yl, wk)


ccc         yl(iv) = yold - dyl    ! for 2nd order Jac
ccc         call pandf1 (xc, yc, iv, neq, t, yl, yldot0) ! for 2nd order Jac
c ... Calculate possibly nonzero Jacobian elements for this variable,
c     and store nonzero elements in compressed sparse column format.
c         write(*,*) 'store non-zero elem ', 'from thread = ', TID
         iJacRow(iv) = nnz      ! sets index for first Jac. elem. of var. iv
c         write(*,'(a,i3,a,i5,i5,i5,i5,e16.6,i5)') '!',ith,':xc,yc,ii1,ii2=',xc,yc,ii1,ii2,sfscal(iv),nnz
c!!!!!!!!!!!!!!! loop ii

         do ii = ii1, ii2
ccc            if (iondenseqn .eq. "inel") then   ! skip ion cont. eqn
ccc               ix = igyl(ii,1)               ! (but not for guard cells)
ccc               iy = igyl(ii,2)
ccc               if (isnion(1) .eq. 1 .and. mod(ii,numvar) .eq. 1 .and.
ccc     .             .not. tstguardc (ix, iy)) go to 10  !end of iv loop removed
ccc            endif
c dtuse,iseqalg,jacelem,idxphi,dtphi,sfscal,jaccliplim,nnzmx

            jacelem = (wk(ii) - yldot00(ii)) / dyl
ccc            jacelem = (wk/(ii) - yldot0(ii)) / (2*dyl)  ! for 2nd order Jac
c ...  Add diagonal 1/dt for nksol
      if (((svrpkg.eq."nksol") .or. (svrpkg.eq."petsc")) .and. iv.eq.ii)
     .then
              if (iseqalg(iv)*(1-isbcwdt).eq.0) then
                jacelem = jacelem - 1/dtuse(iv)
              endif
              ix = igyl(iv,1)
              iy = igyl(iv,2)
              if (idxphi(ix,iy)==iv .and. dtphi<1e10) then !selects phi eqn
                jacelem = jacelem - 1/dtphi
              endif
c              write(*,*) 'iv,jacelem after dt:',jacelem
       endif

c ...  Add a pseudo timestep to the diagonal !! if eqn is not algebraic
       if (svrpkg .ne. "cvode" .and. nufak .gt. 0) then
               if (iv.eq.ii .and. yl(neq+1).eq.1)
     .          jacelem = jacelem - nufak  !omit .and. iseqalg(iv).eq.0)
cc     .                   jacelem = jacelem - nufak*suscal(iv)/sfscal(iv)
               endif
c        write(*,'(a,i3,i5,i5,e16.6)') '!',ith,iv,ii,jacelem
        if (abs(jacelem*sfscal(iv)) .gt. jaccliplim) then

               if (nnz .gt. nnzmx) then
                  write(*,*) nnz,nnzmx
                  write(STDOUT,*)
     .             'ith',ith,'*** jac_calc -- More storage needed for Jacobian.',
     .             ' Storage exceeded at (i,j) = (',ii,',',iv,').',
     .             ' Increase lenpfac.'
	              call xerrab("")
               endif
cc               if (rdoff.ne.0.e0) jacelem=jacelem*(1.0e0+ranf()*rdoff)
cJG We need to modify the way the jacobian elements are counted for since this is a 'linear' implementation.
cJG We cannot use the do linear construct for nnz here because there is an if condition on nnz
cWe create an array to store value of jacelem and ii per value of iv then we reduce it after the target implementation
cJG That way, we avoid having to deal with synchronization and atomic construct which are time consuming
cJG(under the assumoption that we have way enough memory (>10Gb), why bother for a single array?)
cJG               rcsc(nnz) = jacelem
cJG               icsc(nnz) = ii
             iJacCol(nnz)=ii
            rJacElem(nnz)=jacelem

        if (OMPDebug.gt.0) then
c      write(*,'(a,i3,a,i7,i7,i7, E16.6,E16.6,E16.6)') '!', ith,' : ',iv,nnz,
c     .ii,jacelem,wk(ii),yldot00(ii)
        endif


            nnz = nnz + 1

         endif

             if (istopjac.gt.0 .and. ii.eq.irstop .and. iv.eq.icstop) then
ccc               yldot_pert(ii) = wk(ii)      ! for diagnostic only
                   if (istopjac == 2) then
                     yl(iv) = yold
                     call pandf1 (xc, yc, iv, neq, t, yl, wk)
                   endif
ccc               call convsr_vo (xc, yc, yl)  ! was one call to convsr
ccc               call convsr_aux (xc, yc, yl)
                   call remark("***** non-zero jac_elem at irstop,icstop")
                  write(*,*) 'irstop = ', irstop, ', icstop = ', icstop
                   call xerrab("")
             endif
         enddo  ! end of ii loop over equations


c ... Restore dependent variable and plasma variables near its location.
         yl(iv) = yold
ccc         call convsr_vo (xc, yc, yl)  ! was one call to convsr
ccc         call convsr_aux (xc, yc, yl)
         call pandf1 (xc, yc, iv, neq, t, yl, wk)
         yldot_pert(1:neq)=wk(1:neq)
ccc 18   continue
c...  If this is the last variable before jumping to new cell, reset pandf
ccc 18   continue
cJG1 why do we need to reset?????
c...  If this is the last variable before jumping to new cell, reset pandf
         if (mod(iv,numvar).eq.0 .and. isjacreset.ge.1) then
c         call DebugHelper('dump00.txt')
           call pandf1 (xc, yc, iv, neq, t, yl, wk)
         endif
ccc$omp critical
         do ii=ii1, ii2
         if (yldot_pert(ii).ne.wk(ii)) then

      write(*,'(a,i3,a,i5,e20.12,e20.12)') ' ith:',ith,' *** wk modified on second call at ii=',
     . ii,yldot_pert(ii),wk(ii)
c         call DebugHelper('dump11.txt')
         write(*,*) 'dumped'
         stop

         endif

         enddo
cccc$omp end critical
cJG      yldot_pert(iv)=gettime(sec4)-time0
      enddo
c ... End loop over dependent variables and finish Jacobian storage.
      end subroutine LocalJacBuilder


c-----------------------------------------------------------------------
      subroutine jac_lu_decomp (neq, jac, ja, ia, wp, iwp)

c ... Compute LU decomposition of Jacobian and return it in any one
c     of the storage formats.

c ... Common blocks:
            use Timing                   ! ttmatfac
            use Decomp                   ! lbw,ubw
            use Grid                     ! ngrid,ig,ijac
            use Condition_number         ! rconds
            use Preconditioning          ! premeth,lenplumx
            use Ilutv                     ! lfililut,tolilut
            use Nonzero_diagonals        ! lfilinel,ndiagmx,ndiag,ndiagm,
                       ! adiag(neq,ndiagmx),siginel(neq),fmuinel(neq)
                       ! iwkd1(2*neq-1),iwkd2(ndiagmx),rwkd(ndiagmx)
            use Temporary_work_arrays    ! rwk1,rwk2,iwk1,iwk2,iwk3

c ... Function:
      implicit none

c ... Input arguments:
      integer neq        ! total number of equations (all grid points)
      real jac(*)        ! nonzero Jacobian elements
      integer ja(*)      ! col indices of nonzero Jacobian elements
      integer ia(neq+1)  ! pointers to beginning of each row in jac,ja

c ... Output arguments:
      real wp(*)         ! matrix elements of LU
      integer iwp(*)     ! sizes and array indices for elements of LU

      real(kind=4):: gettime

c ... Local variables:
      integer lowd, ierr, i, idum(1)
      real rcond, dum(1)
      real(kind=4):: sec4
      real tsmatfac

c ... Convert compressed sparse row to banded format and use exact
c     factorization routine sgbco from Linpack/SLATEC.
      write(*,*) '******** jac_lu_decomp ********'
c      write(*,*) 'size wp', shape(wp)
      if (premeth .eq. 'banded') then
         lowd = 2 * lbw + ubw + 1
         write(*,*) 'lowd',lowd,lbw,ubw
         call csrbnd (neq, jac, ja, ia, 0, wp, lowd, lowd,
     .                lbw, ubw, ierr)
c             csrbnd (n,    a, ja,  ia, job,abd,nabd,lowd,ml,mu,ierr)
         if (ierr .ne. 0) then
            write(STDOUT,*)
     .         '*** jac_lu_decomp -- csrbnd returned ierr =', ierr
            call xerrab("")
         endif
         tsmatfac = gettime(sec4)
         call sgbco (wp, lowd, neq, lbw, ubw, iwp(4), rcond, rwk1)
         iwp(1) = lowd
         iwp(2) = lbw
         iwp(3) = ubw

c ... Save condition number.
         i = ijac(ig)
         if (i .le. 300) rconds(i,ig) = rcond
         go to 99
      endif

c ... If sparse Jacobian matrix is in compressed sparse row storage
c     format, ...
      if (premeth .eq. 'ilut') then

c ... Reorder Jacobian rows and columns, if desired.
         call jac_reorder (neq, jac, ja, ia, wp, iwp(neq+2), iwp)

c ... Use incomplete factorization routine ilut from SparsKit.
         tsmatfac = gettime(sec4)
         call ilut (neq,jac,ja,ia,lfililut,tolilut,wp,iwp(neq+1),
     .              iwp,lenplumx,rwk1,rwk2,iwk1,
     .              iwk2,iwk3,ierr)
         if (ierr .ne. 0) then
            write(STDOUT,*) ' Error return from ilut:  ierr = ',ierr
            write(STDOUT,9000)
 9000       format(
     ./'    ierr >  0   --> Zero pivot encountered at step number ierr.'
     ./'    ierr = -1   --> Error. input matrix may be wrong.'
     ./'                     (The elimination process has generated a'
     ./'                     row in L or U with length > n.)'
     ./'    ierr = -2   --> Matrix L overflows.'
     ./'    ierr = -3   --> Matrix U overflows.'
     ./'    ierr = -4   --> Illegal value for lfililut.'
     ./'    ierr = -5   --> Zero row encountered.'
     ./'    '
     ./'    For ierr = -2 or -3, increase the value of lenplufac or'
     ./'    decrease the value of lfililut if lenplufac cannot be'
     ./'    increased.'
     .)
            call xerrab("")
         endif

c ... Use incomplete factorization routine precond5 from INEL.
c     SparsKit routines are used in preliminary steps to convert to
c     diagonal format.
      elseif (premeth .eq. 'inel') then

c ... Get the number of nonzero diagonals and the maximum in LU.
         call infdia (neq, ja, ia, iwkd1, ndiag)
         if (ndiag .gt. ndiagmx) then
            call remark('More storage for diagonals of the Jacobian')
            call remark('is needed.  Increase value of ndiagmx.')
            call xerrab("")
         endif
         ndiagm = min(lfilinel+ndiag, ndiagmx)
         iwp(1) = ndiag
         iwp(2) = ndiagm

c ... Convert to diagonal format.
         call csrdia (neq, ndiag, 10, jac, ja, ia, neq, adiag,
     .                iwp(3), dum, idum, idum, iwkd1)

c ... Reorder rows to be in increasing column order.
         call cdiagsrt (neq, adiag, ndiag, iwp(3), iwkd1, iwkd2,
     .                  rwkd)

c ... Finally, calculate approximate LU decomposition.
         tsmatfac = gettime(sec4)
         call precond5 (neq, ndiag, ndiagm, adiag, wp, rwk2, rwk1,
     .                  iwk3, iwk2, siginel, fmuinel, iwp(3))
      endif

c ... Accumulate cpu time spent here.
 99   ttmatfac = ttmatfac + (gettime(sec4) - tsmatfac)
      return
      end
c-----------------------------------------------------------------------
      subroutine jac_reorder (neq, jac, ja, ia, awk, jwk, iwk)

c ... If desired, reorder the Jacobian matrix.

c ... Common blocks:
            use Preconditioning  ! premeth
            use Timing           ! ttjreorder
            use Jacreorder       ! perm,qperm,levels,nlev,mask,maskval,ireorder

c ... Function:
      implicit none

c ... Input arguments:
      integer neq        ! total number of equations (all grid points)
      real jac(*)        ! nonzero Jacobian elements
      integer ja(*)      ! column indices of nonzero Jacobian elements
      integer ia(neq+1)  ! indices of 1st nonzero element in each row

c ... Work-array arguments:
      real awk(*)
      integer jwk(*)
      integer iwk(neq+1)

      real(kind=4):: gettime

c ... Local variables:
      real(kind=4):: sec4
      real tsjreorder
      integer i, nfirst

c ... Get initial value of system cpu timer.
      tsjreorder = gettime(sec4)

      if ((ireorder .eq. 1) .and. (premeth .eq. 'ilut')) then

c ... Copy jac, ja, and ia to awk, jwk, and iwk.
         call atob (neq, jac, ja, ia, awk, jwk, iwk)

c ... Perform a Cuthill-McKee reordering of the Jacobian.
         nfirst = 1
         perm(1) = 0
         do i = 1, neq
            mask(i) = 1
         enddo
         maskval = 1
         qperm(1) = 1
         call bfs (neq,jwk,iwk,nfirst,perm,mask,maskval,qperm,levels,
     .             nlev)

c ... Reverse the permutation to obtain the reverse Cuthill-McKee
c     reordering.
         call reversp (neq,qperm)

c ... Calculate the inverse of qperm and put it in perm.
         do i = 1, neq
            perm(qperm(i)) = i
         enddo

c ... Permute rows and columns of Jacobian using perm.
         call dperm (neq,awk,jwk,iwk,jac,ja,ia,perm,perm,1)

c ... End of If block
      endif

c ... Accumulate cpu time spent here.
      ttjreorder = ttjreorder + (gettime(sec4) - tsjreorder)
      return
      end
c-----------------------------------------------------------------------
      subroutine jac_norm_rows (neq, jac, ja, ia)

c ... If desired, normalize each row using one of three types of norm.

c ... Common blocks:
            use Timing    ! ttjrnorm
            use Math_problem_size    ! neqmx for arrays in Jacaux not used here 
            use Jacaux    ! isrnorm,normtype,fnormnw

c ... Function:
      implicit none

c ... Input arguments:
      integer neq        ! total number of equations (all grid points)
      real jac(*)        ! nonzero Jacobian elements
      integer ja(*)      ! column indices of nonzero Jacobian elements
      integer ia(neq+1)  ! indices of 1st nonzero element in each row

      real(kind=4):: gettime

c ... Local variables:
      real(kind=4):: sec4
      real tsjrnorm

c ... Get initial value of system cpu timer.
      tsjrnorm = gettime(sec4)

      if (isrnorm .eq. 1) call roscal (neq, 0, normtype, jac, ja, ia,
     .                               fnormnw, jac, ja, ia)

c ... Accumulate cpu time spent here.
      ttjrnorm = ttjrnorm + (gettime(sec4) - tsjrnorm)
      return
      end
c-----------------------------------------------------------------------
      subroutine jac_sub_cj (neq, jac, ja, ia, cj)

c ... Loop through all Jacobian elements, skipping off-diagonal
c     elements and those that correspond to algebraic equations, and
c     subtract cj.

      implicit none

c ... Input arguments:
      integer neq        ! total number of equations (all grid points)
      integer ja(*)      ! column indices of nonzero Jacobian elements
      integer ia(neq+1)  ! indices of 1st nonzero element in each row
      real cj            ! scalar to be subtracted

c ... In-out argument:
      real jac(*)    ! nonzero Jacobian elements

c ... Function:
      integer res_algeb

c ... Local variables:
      integer i   ! index to row (i.e., equation)
      integer k   ! index to nonzero element and its column index

      do i = 1, neq
         do k = ia(i), ia(i+1)-1
            if (ja(k) .ne. i) go to 10
            if (res_algeb (i) .eq. 1) go to 10
            jac(k) = jac(k) - cj
 10         continue
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      integer function res_algeb (i)
c ... Common blocks:
            use Dim        ! nx,ny,nxpt
            use Xpoint_indices       ! ixlb,ixrb
            use Math_problem_size    ! neqmx for arrays in Indexes not used here 
            use Indexes    ! igyl,idxphi
            use UEpar      ! isphionxy

c ... Local variables:
      implicit none
c ... Return 1 if the equation index corresponds to a potential
c     equation or a boundary point, otherwise return 0.

c ... Input argument:
      integer,intent(in):: i   ! equation index

      integer ix, iy, jx

      ix = igyl(i,1)   ! get spatial location for this equation
      iy = igyl(i,2)
      if (isphionxy(ix,iy).eq. 1) then
         if (i .eq. idxphi(ix,iy)) then
            res_algeb = 1
            return
         endif
      endif
      if (ix .eq. 0 .or. ix .eq. nx+1 .or. iy .eq. 0 .or. iy .eq. ny+1)
     .   then
         res_algeb = 1
      else
         res_algeb = 0
      endif
ccc   MER NOTE: Generalize the above for interior guard cells with 'dnull'
ccc             or 'dnbot' or islimon==1;
ccc             or should we use the array iseqalg(i) instead?
ccc      do jx = 1, nxpt
ccc         if (ix==ixlb(jx) .or. ix==ixrb(jx)+1 .or. iy==0 .or. iy==ny+1)
ccc     .      then
ccc            res_algeb = 1
ccc         else
ccc            res_algeb = 0
ccc         endif
ccc      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine jacd1 (t, yl, yldot, pd, cj, rpar, ipar)

c ... Calculate matrix A in banded format for daspk solver.
c     Matrix A has the form
c        A(i,j) = dG(i)/dYL(j) + CJ*dG(i)/dYLDOT(j)
c     where G is the residual described in the comments for resid.
c     Therefore, for boundary points and for the potential equation,
c        A(i,j) = df(i)/dYL(j) = Jacobian
c     Otherwise,
c        A(i,j) = df(i)/dYL(j) - CJ*(Kronecker delta ij)

c ... Common blocks:
            use Math_problem_size    ! neqmx
            use Jacaux      ! yldot0,yldot1
c     Temporary_work_arrays cannot be used here because neq is not
c     an argument.
            use Jacobian    ! nnzmx,jac,jacj,jaci

c ... Local variables:
      implicit none

c ... Input arguments:
      real t          ! physical time
      real yl(*)      ! most recent iterate of solution vector
      real yldot(*)   ! most recent iterate of time-derivative of yl
                      ! (not used)
      real rpar(*)    ! real parameter communication (not used)
      integer ipar(3) ! integer parameter communication
      real cj         ! scalar provided by daspk

c ... Output argument:
      real pd(*)      ! matrix A in (Linpack) banded format

      integer neq, lbw, ubw
      integer lowd, ierr

c ... Get total number of equations (all grid points) and
c     lower and upper bandwidths.
      neq = ipar(1)
      lbw = ipar(2)
      ubw = ipar(3)

c ... Calculate right-hand sides at unperturbed values of variables for
c ... Jacobian calculation.
      call ffun (neq, t, yl, yldot0)

c ... Calculate Jacobian of right-hand sides of UEDGE equations.
      call jac_calc (neq, t, yl, yldot0, lbw, ubw, yldot1,
     .               nnzmx, jac, jacj, jaci)

c ... Subtract cj from appropriate Jacobian elements.
      call jac_sub_cj (neq, jac, jacj, jaci, cj)

c ... Convert Jacobian from compressed sparse row to (Linpack) banded
c     format.
      lowd = 2 * lbw + ubw + 1
      call csrbnd (neq, jac, jacj, jaci, 0, pd, lowd, lowd,
     .             lbw, ubw, ierr)
      if (ierr .ne. 0) then
         write(STDOUT,*) '*** jacd1 -- ierr =', ierr
         call xerrab("")
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine jacd2 (resid, ires, neq, t, yl, yldot, rewt, savr, wk,
     .                  h, cj, wp, iwp, ierr, rpar, ipar)

c ... Calculate preconditioning matrix that approximates A (see
c     subroutine jacd1) for daspk solver.

c ... Common blocks:
            use Math_problem_size    ! neqmx
            use Jacaux      ! yldot0,jscalcol
            use Jacobian    ! nnzmx,jac,jacj,jaci

c ... Local variables:
      implicit none

c ... Input arguments:
      external resid   ! function that evaluates residuals (not used)
      integer ires     ! error flag from resid (not used)
      integer neq      ! total number of equations (all grid points)
      real t           ! physical time
      real yl(*)       ! most recent iterate of solution vector
      real yldot(neq)  ! most recent iterate of time-derivative of yl
                       ! (not used)
      real rewt(neq)   ! reciprocal error weights
      real savr(neq)   ! residual values G(t,yl,yldot) (not used)
      real h           ! step size (not used)
      real cj          ! scalar provided by daspk
      real rpar(*)     ! real parameter communication (not used)
      integer ipar(3)  ! integer parameter communication

c ... Work-array argument:
      real wk(neq)     ! work space available to this subroutine

c ... Output arguments:
      real wp(*)       ! matrix elements of LU
      integer iwp(*)   ! array indices for elements of LU
      integer ierr     ! error flag (0 means success, else failure)

      integer lbw, ubw, ii
      real sqrtn

c ... Unpack parameters.
      lbw = ipar(2)
      ubw = ipar(3)

c ... Calculate right-hand sides at unperturbed values of variables.
c     (These values could be obtained with more coding but less
c     computation by taking savr and adding yldot for differential
c     equations.)
      call ffun (neq, t, yl, yldot0)

c ... Calculate Jacobian.
      call jac_calc (neq, t, yl, yldot0, lbw, ubw, wk,
     .               nnzmx, jac, jacj, jaci)

c ... Subtract cj from appropriate Jacobian elements.
      call jac_sub_cj (neq, jac, jacj, jaci, cj)

C ... Multiply Jacobian columns by inverse of scaling vector REWT.
C     In PSOLILU, the WGHT array equals REWT/SQRT(NEQ), so we must
C     be consistent here.  Copied from coding by Peter Brown (3/10/97)
      if (jscalcol .eq. 1) then
         sqrtn = sqrt(real(neq))
         do 10 ii = 1, neq
            wk(ii) = sqrtn / rewt(ii)
 10      continue
         call amudia (neq, 0, jac, jacj, jaci, wk, jac, jacj, jaci)
      endif

c ... Normalize Jacobian rows, if desired.
      call jac_norm_rows (neq, jac, jacj, jaci)

c ... Compute LU decomposition.
      call jac_lu_decomp (neq, jac, jacj, jaci, wp, iwp)

      ierr = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine jacnw (neq, yl, f0, dt, wk, wp, iwp)

c ... Calculate LU decomposition of the Jacobian matrix
c     for use as a preconditioner by the newton solver.

c ... Common blocks:
            use Decomp      ! lbw,ubw
            use Jacobian    ! nnzmx,jac,jacj,jaci

c ... Local variables:
      implicit none

c ... Input arguments:
      integer neq      ! total number of equations (all grid points)
      real yl(*)     ! most recent iterate of solution vector
      real f0(neq)     ! function values f(yl)
      real dt          ! false timestep to improve condition number

c ... Work-array argument:
      real wk(neq)     ! work space available to this subroutine

c ... Output arguments:
      real wp(*)       ! matrix elements of LU
      integer iwp(*)   ! array indices for elements of LU

      real tp
      write(*,*) '***** jacnw *********'
c ... Flag these calls to RHS (pandf) as for the Jacobian
      yl(neq+1) = 1.

c ... Call pandf to set terms with yl(neq+1) flag on for Jacobian
      tp = 0.
      call ffun (neq, tp, yl, f0)

c ... Calculate Jacobian matrix.
      call jac_calc (neq, tp, yl, f0, lbw, ubw, wk,
     .               nnzmx, jac, jacj, jaci)

      yl(neq+1) = -1.        ! Turn-off Jacobian flag for pandf

c ... Normalize Jacobian rows, if desired.
      call jac_norm_rows (neq, jac, jacj, jaci)

c ... Compute LU decomposition.
      call jac_lu_decomp (neq, jac, jacj, jaci, wp, iwp)

      return
      end
c-----------------------------------------------------------------------
      subroutine jacvd (f, neq, tp, yl, ylsv, rewt, fty, wk, hrl1,
     .                  wp, iwp, ierr, rpar, ipar)

c ... Calculate (identity matrix) - hrl1 * (Jacobian matrix), and return
c     exact or approximate LU decomposition for the vodpk solver.

c ... Common block:
            use Jacobian    ! nnzmx,jac,jacj,jaci

c ... Local variables:
      implicit none

c ... Input arguments:
      external f       ! function that evaluates f(tp,yl)
      integer neq      ! total number of equations (all grid points)
      real tp          ! physical time
      real yl(*)     ! most recent iterate of solution vector
      real ylsv(neq)   ! a copy of yl (not used)
      real rewt(neq)   ! reciprocol error weights (not used)
      real fty(neq)    ! function values f(tv,yl)
      real hrl1        ! scalar provided by vodpk
      real rpar(*)     ! real parameter communication (not used)
      integer ipar(2)  ! integer parameter communication

c ... Work-array argument:
      real wk(neq)     ! work space available to this subroutine

c ... Output arguments:
      real wp(*)       ! matrix elements of LU
      integer iwp(*)   ! array indices for elements of LU
      integer ierr     ! error flag (0 means success, else failure)

      integer lbw, ubw, nnz

c ... Unpack parameters.
      lbw = ipar(1)
      ubw = ipar(2)

c ... Calculate Jacobian.
      call jac_calc (neq, tp, yl, fty, lbw, ubw, wk,
     .               nnzmx, jac, jacj, jaci)

c ... Jacobian elements could be saved here for reuse if change in hrl1
c     is more significant than changes in yl.

c ... Calculate (identity matrix) - hrl1 * (Jacobian matrix).
      nnz = jaci(neq+1) - 1
      call sscal (nnz, -hrl1, jac, 1)
      call aplsca (neq, jac, jacj, jaci, 1., iwp)

c ... Normalize Jacobian rows, if desired.
      call jac_norm_rows (neq, jac, jacj, jaci)

c ... Compute LU decomposition.
      call jac_lu_decomp (neq, jac, jacj, jaci, wp, iwp)

      ierr = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine FPRECO (neq, tp, yl, f0, jok, jcur, hrl1, rewt, h,
     .                   uround, nfe, v1, v2, v3, ierr)

c ... Calculate (identity matrix) - hrl1 * (Jacobian matrix), and return
c     exact or approximate LU decomposition for the pvode solver.

c ... Common block:
            use Jacobian    ! nnzmx,jac,jacj,jaci
            use Decomp      ! lbw,ubw
            use Jac_work_arrays  ! iwwp, wwp  use instead of iwp and wp for cvode cases 

c ... Local variables:
      implicit none

c ... Input arguments:
      integer neq      ! total number of equations (all grid points)
      real tp          ! physical time
      real yl(neq)     ! most recent iterate of solution vector
      real f0(neq)     ! function values f(yl)
      integer jok      ! Jacobian ok used by pvode
      real hrl1        ! scalar provided by pvode
      real rewt(neq)   ! reciprocol error weights (not used)
      real h           ! step-size from pvode, don't change
      real uround      ! from pvode, don't change
      integer nfe      ! input/output for number of RHS evaluations

c ... Work-array argument:
      real v1(neq)     ! work space available to this subroutine
      real v2(neq)     ! work space available to this subroutine
      real v3(neq)     ! work space available to this subroutine

c ... Output arguments:
      integer jcur     ! =1 if FPRECO has updated Jacobian
      integer ierr     ! error flag (0 means success, else failure)

      integer nnz

c ... Calculate Jacobian.

      call jac_calc (neq, tp, yl, f0, lbw, ubw, v1, nnzmx,
     .               jac, jacj, jaci)

c ... Jacobian elements could be saved here for reuse if change in hrl1
c     is more significant than changes in yl.

c ... Calculate (identity matrix) - hrl1 * (Jacobian matrix).
      nnz = jaci(neq+1) - 1
      call sscal (nnz, -hrl1, jac, 1)
      call aplsca (neq, jac, jacj, jaci, 1., iwwp)

c ... Normalize Jacobian rows, if desired.
      call jac_norm_rows (neq, jac, jacj, jaci)

c ... Compute LU decomposition.
      call jac_lu_decomp (neq, jac, jacj, jaci, wwp, iwwp)

cxqx Alan suggestion
      jcur = 1

      ierr = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine jacvnk (neq, yl, f0, v, z, wp, iwp)

c ... This subroutine is present merely to satisfy the requirements of

c     the option mdif=1 is not working in nksol and will be removed in
c     the future.  The Jacobian for nksol is calculated through psetnk.

      implicit none

c ... Input arguments:
      integer neq      ! total number of equations (all grid points)
      real yl(*)       ! most recent iterate of solution vector
      real f0(neq)     ! function values f(yl)
      real v(neq)      ! arbitrary vector
      real wp(*)       ! matrix elements of LU
      integer iwp(*)   ! array indices for elements of LU

c ... Output argument:
      real z(neq)      ! (most recently calculated Jacobian) * v

      return
      end
c-----------------------------------------------------------------------
      subroutine psetnk (neq, yl, f0, su, sf, wk, f, wp, iwp, ierr)

c ... Calculate LU decomposition of the Jacobian matrix
c     for use as a preconditioner by the nksol solver.

c ... Common blocks:
            use Decomp          ! lbw,ubw
            use Jacobian        ! nnzmx,jac,jacj,jaci
            use Math_problem_size    ! neqmx
            use Dim             ! nx,ny
            use Time_dep_nwt    ! nufak,ydt_max,ydt_max0,alfnuf,expnuf,nufak0
                          ! inufaknk,dtoptx,dtoptv
            use Indexes         ! idxn,idxu,idxte,idxti,idxng,idxphi
            use UEpar           ! isnion,isupon,isteon,istion,isngon,isphion,isnionxy,
                          ! isuponxy,isteonxy,istionxy,isngonxy,isphionxy
            use Jac_work_arrays  ! iwwp,wwp,liwp,lwp  ! diagnostic arrays in this sub

c ... Local variables:
      implicit none

c ... Input arguments:
      integer neq      ! total number of equations (all grid points)
      real yl(*)       ! most recent iterate of solution vector
      real f0(neq)     ! function values f(yl)
      real su(neq)     ! scale factors for yl
      real sf(neq)     ! scale factors for function values f(yl)
      external f       ! function that evaluates residuals f(yl)

c ... Work-array argument:
      real wk(neq)     ! work space available to this subroutine

c ... Output arguments:
      real wp(*)       ! matrix elements of LU
      integer iwp(*)   ! array indices for elements of LU
      integer ierr     ! error flag (0 means success, else failure)

      real tp
      integer i
      write(*,*) '**************** psetnk ********************'
c ... Calculate maximum of f0*sf to control yl(neq+2) = nufak
      ydt_max = 1.e-100
      do i = 1, neq    ! need to avoid neq+1 and neq+2
         if (abs(f0(i)*sf(i)) .gt. ydt_max)
     .                           ydt_max = abs(f0(i)*sf(i))
      enddo
      if (ydt_max0 .eq. 0) ydt_max0 = ydt_max
      nufak = min(nufak*alfnuf*(ydt_max/ydt_max0)**expnuf, nufak0)
      if (inufaknk .eq. 1) then   ! deter. if nufak is used in Krylov step
         yl(neq+2) = nufak
      else
         yl(neq+2) = 0.
      endif
      if (expnuf.ne.0.) write(*,*) ' nufak = ', nufak
      ydt_max0 = ydt_max

c ... Flag these calls to RHS (pandf) as for the Jacobian
      yl(neq+1) = 1.

c ... Call pandf to set terms with yl(neq+1) flag on for Jacobian
      call rhsnk (neq, yl, f0)

c ... Calculate Jacobian matrix.
      tp = 0.
      if (ParallelJac.gt.0) then
      call jac_calc_para (neq, tp, yl, f0, lbw, ubw, wk,
     .               nnzmx, jac, jacj, jaci)
      else
      call jac_calc (neq, tp, yl, f0, lbw, ubw, wk,
     .               nnzmx, jac, jacj, jaci)
      endif
      yl(neq+1) = -1.             ! Turn-off Jacobian flag for pandf
      call rhsnk (neq, yl, f0)    ! Reset f0 with nufak off

c ... Multiply Jacobian columns by inverse of scaling vector su.
      do i = 1, neq
         wk(i) = 1. / su(i)
      enddo
      call amudia (neq, 0, jac, jacj, jaci, wk, jac, jacj, jaci)

c ... Multiply Jacobian rows by scaling vector sf.  Scaling the
c     columns and rows allows nksol to work with quantities of
c     order unity regardless of the size of f0 and yl.
      call diamua (neq, 0, jac, jacj, jaci, sf, jac, jacj, jaci)

c ... Normalize Jacobian rows, if desired.
      call jac_norm_rows (neq, jac, jacj, jaci)

c ... Initialize iwwp and wwp
      do i = 1, liwp
        iwwp(i) = 0
      enddo
      do i = 1, lwp
        wwp(i) = 0.
      enddo

c ... Compute LU decomposition.
      call jac_lu_decomp (neq, jac, jacj, jaci, wp, iwp)

c ... Copy LU decomp into common arrays for diagostic
      do i = 1, liwp
        iwwp(i) = iwp(i)
      enddo
      do i = 1, lwp
        wwp(i) = wp(i)
      enddo

      ierr = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine psolnk (neq, yl, f0, su, sf, f, jvsol, wk, wp, iwp, bl,
     .                    ierr)

c ... Interface between linear-system solver nksol and subroutine
c     psolbody.

      implicit none

c ... Input arguments:
      integer neq      ! total number of equations (all grid points)
      real yl(neq)     ! most recent iterate of solution vector
      real f0(neq)     ! function values f(yl)
      real su(neq)     ! scale factors for yl
      real sf(neq)     ! scale factors for function values f(yl)
      external f       ! function that evaluates residuals f(yl)
      external jvsol   ! function that evaluates jac(u)*v
      real wp(*)       ! matrix elements of LU
      integer iwp(*)   ! dimensions and array indices for elements of LU

c ... Work-array argument:
      real wk(neq)     ! work space available to this subroutine

c ... In-out argument:
      real bl(neq)     ! on input, c of P*x=c; on output, x

c ... Output argument:
      integer ierr     ! error flag (0 OK, >0 try again, <0 failure)

c ... Local variable:
      logical usingsu

      usingsu = .true.
      call psolbody (neq, usingsu, su, wk, wp, iwp, bl, ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine psold (neq, t, yl, yldot, f0, wk, cj, wght, wp, iwp,
     .                  bl, eplin, ierr, rpar, ipar)

c ... Interface between linear-system solver daspk and subroutine
c     psolbody.

      implicit none

c ... Input arguments:
      integer neq      ! total number of equations (all grid points)
      real t           ! physical time (not used here)
      real yl(neq)     ! most recent iterate of solution vector
      real yldot(neq)  ! most recent iterate of time-derivative of yl
      real f0(neq)     ! G(t,yl,yldot) (not used)
      real cj          ! scalar provided by daspk (not used)
      real wght(neq)   ! error weights
      real wp(*)       ! matrix elements of LU
      integer iwp(*)   ! dimensions and array indices for elements of LU
      real eplin       ! bound on solution error (not used)
      real rpar(*)     ! real parameter communication (not used)
      integer ipar(*)  ! integer parameter communication (not used)

c ... Work-array argument:
      real wk(neq)     ! work space available to this subroutine

c ... In-out argument:
      real bl(neq)     ! on input, c of P*x=c; on output, x

c ... Output argument:
      integer ierr     ! error flag (0 OK, >0 try again, <0 failure)

c ... Local variables:
      logical usingsu
!!      real su(1)       ! array that is not used when usingsu is false

      usingsu = .false.! this may get reset in psolbody if jscalcol=1
!!      su(1) = 0.       ! set unused array so flint will be happy
      call psolbody (neq, usingsu, wght, wk, wp, iwp, bl, ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine psolnw (neq, yl, wk, wp, iwp, bl, ierr)

c ... Interface between linear-system solver newton and subroutine
c     psolbody.

      implicit none

c ... Input arguments:
      integer neq      ! total number of equations (all grid points)
      real yl(neq)     ! most recent iterate of solution vector
      real wp(*)       ! matrix elements of LU
      integer iwp(*)   ! dimensions and array indices for elements of LU

c ... Work-array argument:
      real wk(neq)     ! work space available to this subroutine

c ... In-out argument:
      real bl(neq)     ! on input, c of P*x=c; on output, x

c ... Output argument:
      integer ierr     ! error flag (0 OK, >0 try again, <0 failure)

c ... Local variables:
      logical usingsu
      real su(1)       ! array that is not used when usingsu is false

      usingsu = .false.
      su(1) = 0.       ! set unused array so flint will be happy
      call psolbody (neq, usingsu, su, wk, wp, iwp, bl, ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine psolvd (neq, t, yl, f0, wk, hrl1, wp, iwp, bl, lr,
     .                   ierr, rpar, ipar)

c ... Interface between linear-system solver vodpk and subroutine
c     psolbody.

      implicit none

c ... Input arguments:
      integer neq      ! total number of equations (all grid points)
      real t           ! physical time (not used here)
      real yl(neq)     ! most recent iterate of solution vector
      real f0(neq)     ! function values f(yl)
      real hrl1        ! scalar (not used for lr = 1)
      real wp(*)       ! matrix elements of LU
      integer iwp(*)   ! dimensions and array indices for elements of LU
      integer lr       ! type of calc. to be done (only 1 set up here)
      real rpar(*)     ! real parameter communication (not used)
      integer ipar(*)  ! integer parameter communication (not used)

c ... Work-array argument:
      real wk(neq)     ! work space available to this subroutine

c ... In-out argument:
      real bl(neq)     ! on input, c of P*x=c; on output, x

c ... Output argument:
      integer ierr     ! error flag (0 OK, >0 try again, <0 failure)

c ... Local variables:
      logical usingsu
      real su(1)       ! array that is not used when usingsu is false

      usingsu = .false.
      su(1) = 0.       ! set unused array so flint will be happy
      call psolbody (neq, usingsu, su, wk, wp, iwp, bl, ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine FPSOL (neq, t, yl, f0, wk, hrl1, rewt, delta, nfe, bl,
     .                  lr, zl, ierr)

c ... Interface between linear-system solver pvode and subroutine
c     psolbody.

c ... Common block:
            use Jac_work_arrays  ! iwwp, wwp

c ... Local variables:
      implicit none

c ... Input arguments:
      integer neq      ! total number of equations (all grid points)
      real t           ! physical time (not used here)
      real yl(neq)     ! most recent iterate of solution vector
      real f0(neq)     ! function values f(yl)
      real hrl1        ! scalar (not used for lr = 1)
      real rewt(*)     ! reciprocal of error vector
      real delta       ! input for iterative method (not used)
      integer nfe      ! input/output for number of RHS evaluations
      integer lr       ! type of calc. to be done (only 1 set up here)
      real bl(neq)     ! RHS of P*zl=bl; P is preconditioning matrix

c ... Work-array argument:
      real wk(neq)     ! work space available to this subroutine

c ... Output argument:
      real zl(neq)     ! solution of P*zl=bl
      integer ierr     ! error flag (0 OK, >0 try again, <0 failure)

      logical usingsu
      real su(1)       ! array that is not used when usingsu is false

      usingsu = .false.
      su(1) = 0.       ! set unused array so flint will be happy
cJG      zl(1:neq)=bl(1:neq)
      call scopy (neq, bl, 1, zl, 1)
      call psolbody (neq, usingsu, su, wk, wwp, iwwp, zl, ierr)

      return
      end
c  **  End of subroutine FPSOL *********
c-----------------------------------------------------------------------
      subroutine psolbody (neq, usingsu, su, wk, wp, iwp, bl, ierr)

c ... Solve the linear system P*x=c, using elements of P loaded into
c     array wp.  Loading was done by the following subroutines:
c
c       caller   loading by
c       ------   ----------
c       psold      jacd2
c       psolnk     psetnk
c       psolnw     jacnw
c       psolvd     jacvd
c       fpsol      fpreco

c ... Common blocks:
            use Timing    ! ttmatsol
            use Math_problem_size    ! neqmx for arrays in Jacaux not used here 
            use Jacaux             ! isrnorm,fnormnw,jscalcol
            use Preconditioning    ! premeth
            use Jacreorder         ! perm,qperm,ireorder
            use Dim                ! nisp,ngsp
            use UEpar              ! svrpkg  used to reset usingsu for daspk 

c ... Function:
      implicit none

c ... Input arguments:
      integer neq      ! total number of equations (all grid points)
      logical usingsu  ! .true. if su is used (svrpkg = "nksol" only)
      real su(neq)     ! scale factors for yl
      real wp(*)       ! matrix elements of LU
      integer iwp(*)   ! dimensions and array indices for elements of LU

c ... Work-array argument:
      real wk(neq)     ! work space available to this subroutine

c ... In-out argument:
      real bl(neq)     ! on input, c of P*x=c; on output, x

c ... Output arguments:
      integer ierr     ! error flag (0 OK, >0 try again, <0 failure)

      real(kind=4):: gettime

c ... Local variables:
      integer i
      integer lowd, lbw, ubw
      integer ndiag, ndiagm
      real(kind=4):: sec4
      real tsmatsol

c ... Get initial value of system cpu timer.
      tsmatsol = gettime(sec4)

c ... Scale c by multiplying by row-normalization factors, if used,
c     and by column scaling vector su (for svrpkg='nksol' only).
      if (isrnorm .eq. 1) then
         do i = 1, neq
            bl(i) = bl(i) * fnormnw(i)
         enddo
      endif
      if (usingsu) then    ! omit for daspk column scaling - correct?
         do i = 1, neq
            bl(i) = bl(i) * su(i)
         enddo
      endif

c ... Solve P*x=c for a preconditioner stored as a banded matrix.
      if (premeth .eq. 'banded') then
         lowd = iwp(1)
         lbw = iwp(2)
         ubw = iwp(3)
         call sgbsl (wp, lowd, neq, lbw, ubw, iwp(4), bl, 0)
cJG         wk(1:neq)=bl(1:neq)
         call scopy (neq, bl, 1, wk, 1)

c ... Solve P*x=c for a preconditioner stored as a sparse matrix in
c     compressed sparse row format.
c     If rows and columns of P were reordered (permuted), permute c,
c     then use inverse permutation on x.
      elseif (premeth .eq. 'ilut') then
	 if (ireorder .eq. 1) call dvperm (neq, bl, perm)
         call lusol0 (neq, bl, wk, wp, iwp(neq+1), iwp)
	 if (ireorder .eq. 1) call dvperm (neq, wk, qperm)

c ... Solve P*x=c for a preconditioner stored as a sparse matrix in
c     diagonal storage format.
      else
         ndiag = iwp(1)
         ndiagm = iwp(2)
         call minvmul (neq, ndiag, ndiagm, wp, iwp(3), wk, bl)
      endif

c ... Divide solution x by column-scaling factors (for svrpkg='nksol').
      if (usingsu) then
         do i = 1, neq
            bl(i) = wk(i) / su(i)
         enddo
      elseif (svrpkg.eq."daspk" .and. jscalcol.eq.1) then
         do i = 1, neq
            bl(i) = wk(i) / su(i)
         enddo
      elseif (premeth .ne. 'banded') then
cJG        bl(1:neq)=wk(1:neq)
      call scopy (neq, wk, 1, bl, 1)

      endif

c ... Accumulate cpu time spent here.
      ierr = 0
      ttmatsol = ttmatsol + (gettime(sec4) - tsmatsol)
      return
      end
c-----------------------------------------------------------------------
      subroutine sfsetnk (neq, yl, su, sf)

c ... Calculate the scaling vector of function values for the nksol
c     routine.  The Jacobian at yl is calculated, and sf(i) is set
c     to the reciprocal of the norm of row i.

c ... Common blocks:
            use Decomp               ! lbw,ubw
            use Jacobian             ! nnzmx,jac,jacj,jaci
            use Math_problem_size    ! neqmx
            use Jacaux               ! yldot0,normtype
            use Dim                  ! nx,ny
            use Time_dep_nwt         ! ydt_max0,nufak0
            use Share                ! cutlo

c ... Local variables:
      implicit none

c ... Input arguments:
      integer neq      ! total number of equations (all grid points)
      real yl(*)       ! dependent variables
      real su(neq)     ! scale factors for yl

c ... Output arguments:
      real sf(neq)     ! scaling vector (after use as work array in
                       ! jac_calc and here)

      real tp
      integer i
      write(*,*) '**************** sfsetnk ********************'
c ... Flag these calls to RHS (pandf) as Jacobian calculations
      yl(neq+1) = 1.  ! with ccc, dont include nufak in scaling Jacobian

c ... Calculate right-hand sides at unperturbed values of variables.
      call rhsnk (neq, yl, yldot0)

c ... Calculate Jacobian matrix.
      tp = 0.
      if (ParallelJac.gt.0) then
      call jac_calc_para (neq, tp, yl, yldot0, lbw, ubw, sf,
     .               nnzmx, jac, jacj, jaci)
      else
      call jac_calc (neq, tp, yl, yldot0, lbw, ubw, sf,
     .               nnzmx, jac, jacj, jaci)
      endif


      yl(neq+1) = -1.      ! Turn-off Jacobian flag for pandf
c ... Compute inverse of column-scaling vector, and perform column
c     scaling.
      do i = 1, neq
         sf(i) = 1. / su(i)
      enddo
      call amudia (neq, 0, jac, jacj, jaci, sf, jac, jacj, jaci)

c ... Calculate one of three types of norm for each row.
c ... Also find initial maximum of yldot*sf = ydt_max0 for scaling nufak
      call rnrms (neq, normtype, jac, jacj, jaci, sf)
      nufak0 = nufak                    ! record initial nufak value
      ydt_max0 = cutlo
      do i = 1, neq
c      write(*,*) 'i,sf(i)',i,sf(i),jac(i),jacj(i),jaci(i)
         if (abs(sf(i)) .lt. 1e20*cutlo) then
            write(*,*) '*** Error: Jacobian row = 0 for eqn iv =', i
            call xerrab("")
         else
            sf(i) = 1./sf(i)
         endif
         if (abs(yldot0(i)*sf(i)) .gt. ydt_max0)
     .                            ydt_max0 = abs(yldot0(i)*sf(i))
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_dt(neq, yl, f0)

c ... Calculates the time step to be used with svrpkg="nksol" based
c ... on dtreal and various cases of (yl/yldot)

c ... Common blocks:
            use Math_problem_size    ! neqmx
            use Dim             ! nx,ny
            use Time_dep_nwt    ! nufak,ydt_max,ydt_max0,alfnuf,expnuf,nufak0
                          ! inufaknk,dtoptx,dtoptv
            use Indexes         ! idxn,idxu,idxte,idxti,idxng,idxphi,iseqalg
            use UEpar           ! isbcwdt,isnionxy,isuponxy,isteonxy,istionxy,
                          ! isngonxy,isphionxy
            use Share           ! geometry,nxc,isnonog,cutlo
            use Selec           ! ixm1,ixp1
            use Indices_domain_dcl     ! ixmnbcl,ixmxbcl,iymnbcl,iymxbcl

c ... Local variables:
      implicit none

c ... Input arguments:
      integer neq      ! total number of equations (all grid points)
      real yl(*)       ! most recent iterate of solution vector
      real f0(neq)     ! function values f(yl)

      integer iv,iv1,iv2,iv3,iv4,ifld,igsp,ix,iy,iym1,iyp1,ixm1u,ixp1u
      real up_5ca


      call rhsnk (neq, yl, f0)    ! Reset f0 with nufak off

c ... special new section for adjustable timesteps
      do iy = 1-iymnbcl, ny+iymxbcl
         iym1 = max(1-iymnbcl,iy-1)
	 iyp1 = min(ny+iymxbcl,iy+1)
         do ix = 1-ixmnbcl, nx+ixmxbcl
            do ifld = 1, nisp
               if(isnionxy(ix,iy,ifld) .eq. 1) then
                  iv = idxn(ix,iy,ifld)
                  dtoptv(iv) = deldt*abs(ylodt(iv)/(f0(iv)+cutlo))
                  dtoptx(ix,iy) = min(dtoptv(iv), dtoptx(ix,iy))
                  if (model_dt .eq. 0) then
                    dtuse(iv) = dtreal
                  elseif (model_dt .eq. 1) then
                    dtuse(iv) = dtreal*dtoptv(iv)/(dtreal+dtoptv(iv))
                  elseif (model_dt .eq. 2) then
                    dtuse(iv) = dtoptv(iv)
                  elseif (model_dt .eq. 3) then
                    dtuse(iv) = sqrt(dtreal*dtoptv(iv))
                  endif
               endif
            enddo

            if(ix.ne.nx+2*isbcwdt) then
                           ! nx test - for algebr. eq. unless isbcwdt=1
               do ifld = 1, nusp
                  if(isuponxy(ix,iy,ifld).eq.1) then
                    ixm1u = max(1-ixmnbcl, ixm1(ix,iy))
                    ixp1u = min(nx+ixmxbcl, ixp1(ix,iy))
                    iv = idxu(ix,iy,ifld)
                    iv1 = idxu(ixm1u,iy,ifld)
                    iv2 = idxu(ixp1u,iy,ifld)
                    iv3 = idxu(ix,iyp1,ifld)
                    iv4 = idxu(ix,iym1,ifld)
                    up_5ca = ( abs(ylodt(iv)) + abs(ylodt(iv1)) +
     .                       abs(ylodt(iv2))+ abs(ylodt(iv3)) +
     .                       abs(ylodt(iv4)) )/5
                    if (abs(f0(iv)).gt.cutlo) dtoptv(iv) =
     .                                     deldt*abs(up_5ca/(f0(iv)))
                    if (model_dt .eq. 0) then
                      dtuse(iv) = dtreal
                    elseif (model_dt .eq. 1) then
                      dtuse(iv) = dtreal*dtoptv(iv)/(dtreal+dtoptv(iv))
                    elseif (model_dt .eq. 2) then
                      dtuse(iv) = dtoptv(iv)
                    elseif (model_dt .eq. 3) then
                      dtuse(iv) = sqrt(dtreal*dtoptv(iv))
                    endif
                  endif
               enddo
            endif

            if(isteonxy(ix,iy).eq.1) then
               iv =  idxte(ix,iy)
               dtoptv(iv) = deldt*abs(ylodt(iv)/(f0(iv)+cutlo))
               dtoptx(ix,iy) = min(dtoptv(iv), dtoptx(ix,iy))
               if (model_dt .eq. 0) then
                 dtuse(iv) = dtreal
               elseif (model_dt .eq. 1) then
                 dtuse(iv) = dtreal*dtoptv(iv)/(dtreal+dtoptv(iv))
               elseif (model_dt .eq. 2) then
                 dtuse(iv) = dtoptv(iv)
               elseif (model_dt .eq. 3) then
                 dtuse(iv) = sqrt(dtreal*dtoptv(iv))
               endif
            endif

            if(istionxy(ix,iy).eq.1) then
               iv = idxti(ix,iy)
               dtoptv(iv) = deldt*abs(ylodt(iv)/(f0(iv)+cutlo))
               dtoptx(ix,iy) = min(dtoptv(iv), dtoptx(ix,iy))
               if (model_dt .eq. 0) then
                 dtuse(iv) = dtreal
               elseif (model_dt .eq. 1) then
                 dtuse(iv) = dtreal*dtoptv(iv)/(dtreal+dtoptv(iv))
               elseif (model_dt .eq. 2) then
                 dtuse(iv) = dtoptv(iv)
               elseif (model_dt .eq. 3) then
                 dtuse(iv) = sqrt(dtreal*dtoptv(iv))
               endif
            endif

            do igsp = 1, ngsp
               if(isngonxy(ix,iy,igsp).eq.1) then
                  iv = idxg(ix,iy,igsp)
                  dtoptv(iv) = deldt*abs(ylodt(iv)/(f0(iv)+cutlo))
                  dtoptx(ix,iy) = min(dtoptv(iv), dtoptx(ix,iy))
                  if (model_dt .eq. 0) then
                    dtuse(iv) = dtreal
                  elseif (model_dt .eq. 1) then
                    dtuse(iv) = dtreal*dtoptv(iv)/(dtreal+dtoptv(iv))
                  elseif (model_dt .eq. 2) then
                    dtuse(iv) = dtoptv(iv)
                  elseif (model_dt .eq. 3) then
                    dtuse(iv) = sqrt(dtreal*dtoptv(iv))
                  endif
               endif
            enddo

            do igsp = 1, ngsp
               if(istgonxy(ix,iy,igsp).eq.1) then
                  iv = idxtg(ix,iy,igsp)
                  dtoptv(iv) = deldt*abs(ylodt(iv)/(f0(iv)+cutlo))
                  dtoptx(ix,iy) = min(dtoptv(iv), dtoptx(ix,iy))
                  if (model_dt .eq. 0) then
                    dtuse(iv) = dtreal
                  elseif (model_dt .eq. 1) then
                    dtuse(iv) = dtreal*dtoptv(iv)/(dtreal+dtoptv(iv))
                  elseif (model_dt .eq. 2) then
                    dtuse(iv) = dtoptv(iv)
                  elseif (model_dt .eq. 3) then
                    dtuse(iv) = sqrt(dtreal*dtoptv(iv))
                  endif
               endif
            enddo

            if (isphionxy(ix,iy).eq.1 .and. isbcwdt.eq.1) then
               iv = idxphi(ix,iy)
               dtoptv(iv) = dtoptv(idxte(ix,iy))  ! same as for Te
               if (model_dt .eq. 0) then
                    dtuse(iv) = dtreal
                  elseif (model_dt .eq. 1) then
                    dtuse(iv) = dtreal*dtoptv(iv)/(dtreal+dtoptv(iv))
                  elseif (model_dt .eq. 2) then
                    dtuse(iv) = dtoptv(iv)
                  elseif (model_dt .eq. 3) then
                    dtuse(iv) = sqrt(dtreal*dtoptv(iv))
               endif
            endif
         enddo
      enddo
ccc
c ... If model_dt < 4, then jump over this to 23; otherwise use this
c ... section to define time-step based on cell minimum-step, dtoptx
ccc
      if (model_dt .lt. 4) goto 23
      do iy = 1-iymnbcl, ny+iymxbcl
         do ix = 1-ixmnbcl, nx+ixmxbcl
            do ifld = 1, nisp
               if(isnionxy(ix,iy,ifld) .eq. 1) then
                  iv = idxn(ix,iy,ifld)
                  if (model_dt .eq. 4) then
                    dtuse(iv) = dtreal*dtoptx(ix,iy)/
     .                          (dtreal+dtoptx(ix,iy))
                  elseif (model_dt .eq. 5) then
                    dtuse(iv) = dtoptx(ix,iy)
                  elseif (model_dt .eq. 6) then
                    dtuse(iv) = sqrt(dtreal*dtoptx(ix,iy))
                  endif
               endif
            enddo
            if(ix.ne.nx+2*isbcwdt) then
                           ! nx test - for algebr. eq. unless isbcwdt=1
               do ifld = 1, nusp
                  if(isuponxy(ix,iy,ifld).eq.1) then

                  iv = idxu(ix,iy,ifld)
                  if (model_dt .eq. 4) then
                    dtuse(iv) = dtreal*dtoptx(ix,iy)/
     .                          (dtreal+dtoptx(ix,iy))
                  elseif (model_dt .eq. 5) then
                    dtuse(iv) = dtoptx(ix,iy)
                  elseif (model_dt .eq. 6) then
                    dtuse(iv) = sqrt(dtreal*dtoptx(ix,iy))
                  endif

                  endif
               enddo
            endif
            if(isteonxy(ix,iy).eq.1) then
               iv =  idxte(ix,iy)
               if (model_dt .eq. 4) then
                 dtuse(iv) = dtreal*dtoptx(ix,iy)/
     .                       (dtreal+dtoptx(ix,iy))
               elseif (model_dt .eq. 5) then
                 dtuse(iv) = dtoptx(ix,iy)
               elseif (model_dt .eq. 6) then
                 dtuse(iv) = sqrt(dtreal*dtoptx(ix,iy))
               endif
            endif
            if(istionxy(ix,iy).eq.1) then
            endif
            do igsp = 1, ngsp
               if(isngonxy(ix,iy,igsp).eq.1) then
                  iv = idxg(ix,iy,igsp)
               endif
               if (model_dt .eq. 4) then
                 dtuse(iv) = dtreal*dtoptx(ix,iy)/
     .                       (dtreal+dtoptx(ix,iy))
               elseif (model_dt .eq. 5) then
                 dtuse(iv) = dtoptx(ix,iy)
               elseif (model_dt .eq. 6) then
                 dtuse(iv) = sqrt(dtreal*dtoptx(ix,iy))
               endif
            enddo
            do igsp = 1, ngsp
               if(istgonxy(ix,iy,igsp).eq.1) then
                  iv = idxtg(ix,iy,igsp)
               endif
               if (model_dt .eq. 4) then
                 dtuse(iv) = dtreal*dtoptx(ix,iy)/
     .                       (dtreal+dtoptx(ix,iy))
               elseif (model_dt .eq. 5) then
                 dtuse(iv) = dtoptx(ix,iy)
               elseif (model_dt .eq. 6) then
                 dtuse(iv) = sqrt(dtreal*dtoptx(ix,iy))
               endif
            enddo
            if (isphionxy(ix,iy).eq.1 .and. isbcwdt.eq.1) then
               iv = idxphi(ix,iy)
               if (model_dt .eq. 4) then
                 dtuse(iv) = dtreal*dtoptx(ix,iy)/
     .                       (dtreal+dtoptx(ix,iy))
               elseif (model_dt .eq. 5) then
                 dtuse(iv) = dtoptx(ix,iy)
               elseif (model_dt .eq. 6) then
                 dtuse(iv) = sqrt(dtreal*dtoptx(ix,iy))
               endif
            endif
         enddo
      enddo
 23   continue

c ... Be sure that algebraic equations (esp. limiter) have large dt
      if (isbcwdt .eq. 0) then
        do iv = 1, neq
          if (iseqalg(iv).eq.1) dtuse(iv) = 1.e20
        enddo
      endif

      return
      end
c-----------------------------------------------------------------------

      subroutine jacmap

c ... Common blocks:
            use Math_problem_size    ! neqmx for arrays in Lsode not used here 
            use Lsode            ! neq
            use Jacobian         ! jac,jacj,jaci
            use Jacobian_full    ! jacfull
            use Jacreorder      ! ireorder

c ... Local variables:
      implicit none
c_mpi      include 'mpif.h'

c ... Output a map of the Jacobian matrix.

c ... All information is passed through common for convenience in
c     calling this routine from the UEDGE> prompt.

      integer ierr
      integer us
c_mpi      integer my_pe
      character*24 filename

c ... Allocate full Jacobian for jacmap; warning of size
      call remark("*** CAUTION: allocating large jacfull(neq,neq)***")
      call gallot("Jacobian_full",0)
      write (STDOUT,*) '*** Full Jacobian size is neq**2 = ', neq*neq

c ... Issue a warning if reordering is on which may rearrange Jacobian
      if (ireorder .eq. 1) then
         write (STDOUT,*) '***ireorder=1, Jacobian may be rearranged***'
      endif

cccc ... Attempt to allocate space for Jacobian in full storage format;
cccc     tell user and abort if space is insufficient.
ccc      if (allot('jacfull', neq*neq) .ne. 0) then
ccc         write (STDOUT,*)
ccc     .      '*** jacmap cannot get space for full Jacobian storage.'
ccc         call xerrab("")
ccc      endif

c ... Convert Jacobian matrix to full storage format.
      call csrdns (neq, neq, jac, jacj, jaci, jacfull, neq, ierr)
      if (ierr .ne. 0) then
         write (STDOUT,*)
     .      '*** jacmap got error return ierr =', ierr,
     .      ' from csrdns.'
         call xerrab("")
      endif

c ... Open a file, and output the map.
      call freeus (us)
      filename = 'Jacobian_map.dat'
c_mpi      if(MY_PE().eq.0) then
c_mpi        us = 59
c_mpi        filename = 'Jacobian_map.dat0'
c_mpi      else
c_mpi        us = 69
c_mpi        filename = 'Jacobian_map.dat1'
c_mpi      endif
      open(unit=us, file=filename, status='unknown')
      call jmap (neq, jacfull, us)

c ... Close file, and report file name.
      close(us)
      write (STDOUT,*) ' Jacobian map in data file:  ', filename

      return
      end
c-----------------------------------------------------------------------
      subroutine jacout

c ... Common blocks:
            use Dim         ! nusp for array fnorm in Ynorm not used here 
            use Math_problem_size    ! neqmx
            use Lsode       ! neq,yldot
            use Ynorm       ! sfscal
            use Jacobian    ! jac,jacj,jaci
            use UEpar       ! svrpkg

c ... Local variables:
      implicit none

c ... Output Jacobian matrix and right-hand side in Boeing-Harwell
c     format.

c ... All information is passed through common for convenience in
c     calling this routine from the UEDGE> prompt.

      integer i, us, ifmt
      character*24 filename
      character*72 title

c ... For the nksol solver, scale the right-hand side.
      if ((svrpkg .eq. 'nksol') .or. (svrpkg.eq.'petsc')) then
         do i = 1, neq
            yldot(i) = yldot(i) * sfscal(i)
         enddo
      endif

c ... Open a file, and output data.
      call freeus (us)
      filename = 'Uedge_Test_Matrix.dat'
      open(unit=us, file=filename, status='unknown')
      title = ' UEDGE Test Matrix '
      ifmt = 15
      call prtmt (neq,neq,jac,jacj,jaci,yldot,
     .   'NN',title,'SPARSKIT','RUA',ifmt,3,us)

c ... Close file, and report file name.
      close(us)
      write (STDOUT,*) ' Jacobian matrix in data file:  ', filename

      return
      end
******* end of subroutine jacout *******
c-----------------------------------------------------------------------
      subroutine engbal(pwrin)

c ... Common blocks:
            use Dim               !nx,ny,nzspt,nzsp,ngsp,nfsp,nxpt
            use Xpoint_indices    !ixlb,ixrb
            use Comgeo            !gx,gy,vol,sx,sy
            use Selec             !ixm1,ixp1
            use Comflo            !feex,feix,feey,feiy,fqx,fqy,fnix,fniy
            use Compla            !up,mi,ng,v2,vy
            use UEpar             !ediss,eion,ebind,iigsp,ishosor,fsprd,ishymol
            use Phyvar            !ev,qe
            use Share             !cutlo
            use Rhsides           !psor,psorg,psorxr,erliz,erlrc,psor_tmpov
            use Coefeq            !cfvisx,cfvisy,cnsor,cngmom
            use Bcond             !ckinfl
            use Parallv           ! nxg,nyg
            use Conduc            !visx,visy,eeli,nuiz,nucx
            use Imprad            !prad,pradz
            use Postproc          !complete group
            use Indices_domain_dcl   ! ixmnbcl,ixmxbcl
            use Noggeo            ! angfx

c ... Local variables:
      implicit none

c ... Calculates various components of the 2-D energy flow and the
c ... ionization and radiation for use in the postprocessing file
c ... balancee to determine energy balance; these 2-D loops become
c ... expensive when done from the parser.

c ... Input arguments:
      real pwrin            !total input power for normalization

      integer ix,iy,ix1,ix2,ix3,ix4,iimp,id,igsp,jz,jx,ixt,ixt1,ixr,ixr1
      real ave,a1,a2,t1,t2,thetaix,thetaix2,eta_dup2dy

c ... Implicit function:
      ave(a1,a2) = a1 * a2 / (a1 + a2 + cutlo)
cccac$omp

      if (ishosor .eq. 1) then  ! averge power terms is ishosor=1
         call volavenv(nx, ny, 1, ny, 1, nx, ixp1(0,0), ixm1(0,0),
     .                               fsprd, psor_tmpov(0,0), prad)
         do igsp = nhgsp+1, ngsp
           jz = igsp - nhgsp
           do iimp = 0, nzsp(jz)
             call volavenv(nx, ny, 1, ny, 1, nx, ixp1(0,0), ixm1(0,0),
     .                     fsprd, psor_tmpov(0,0), pradz(0,0,iimp,jz))
           enddo
        enddo
      endif

! Set arrays to check energy conserv; add ion parallel drift and visc heat

      do 10 ix=0,nx
         do 9 iy=0,ny
            ix1 = ixm1(ix,iy)
            ix2 = ixp1(ix,iy)
            ix3 = ixm1(ix,iy-1)
            ix4 = ixm1(ix,iy+1)
            fety(ix,iy) = 0.
            fetx(ix,iy) = 0.
            do id = 1, nusp
	       thetaix =  0.5*(angfx(ix1,iy) + angfx(ix,iy))
               thetaix2 = 0.5*(angfx(ix,iy) + angfx(ix2,iy))
               eta_dup2dy = 0.25*( visy(ix,iy+1,id)*
     .                       (upi(ix1,iy+1,id)+upi(ix,iy+1,id))**2 -
     .                             visy(ix,iy  ,id)*
     .                       (upi(ix1,iy  ,id)+upi(ix,iy  ,id))**2 )
               fety(ix,iy) = fety(ix,iy) + (mi(id)/32)*( upi(ix1,iy,id)+
     .                         upi(ix,iy,id)+upi(ix1,iy+1,id)+
     .                         upi(ix,iy+1,id) )**2*fniy(ix,iy,id) -
     .                         cfvisy*0.5*sy(ix,iy)*gyf(ix,iy)*eta_dup2dy
               fetx(ix,iy) = fetx(ix,iy) + 0.5*mi(id)*upi(ix,iy,id)**2*
     .                          fnix(ix,iy,id) - cfvisx*0.25*sx(ix,iy)*(
     .                         visx(ix ,iy,id)*gx(ix ,iy)*cos(thetaix)*
     .                          ( upi(ix,iy,id)**2 - upi(ix1,iy,id)**2 ) +
     .                        visx(ix2,iy,id)*gx(ix2,iy)*cos(thetaix2)*
     .                          ( upi(ix2,iy,id)**2 - upi(ix,iy,id)**2 ) )
               fetx(ix,iy) = fetx(ix,iy) - upi(ix,iy,id)*fmixy(ix,iy,id)
            enddo
            fety(ix,iy) = fety(ix,iy) + feey(ix,iy) + feiy(ix,iy)
            fetx(ix,iy) = fetx(ix,iy) + feex(ix,iy) + feix(ix,iy)
  9      continue
 10   continue

! Now correct the boundary x-fluxes if non-unity ckinfl
      if (abs(ckinfl-1.) > 1.e-10) then
       do jx = 1, nxpt
         ixt  = ixlb(jx)
         ixt1 = ixt + 1
         ixr  = ixrb(jx)
         ixr1 = ixr - 1
         do 15 iy = 0, ny
            fetx(ixt,iy) = 0.
            fetx(ixr,iy) = 0.
            do id = 1, nfsp
               fetx(ixt,iy) = fetx(ixt,iy) +
     .                        0.5*mi(id)*up(ixt,iy,id)**2*fnix(ixt,iy,id) -
     .                   ckinfl*0.5*sx(ixt,iy)*visx(ixt1,iy,id)*gx(ixt1,iy)*
     .                          ( up(ixt1,iy,id)**2 - up(ixt,iy,id)**2 )
               fetx(ixr,iy) = fetx(ixr,iy) +
     .                        0.5*mi(id)*up(ixr,iy,id)**2*fnix(ixr,iy,id) -
     .                   ckinfl*0.5*sx(ixr,iy)*visx(ixr,iy,id)*gx(ixr,iy)*
     .                          ( up(ixr,iy,id)**2 - up(ixr1,iy,id)**2 )
            enddo
            fetx(ixt,iy) = fetx(ixt,iy) + feex(ixt,iy) + feix(ixt,iy)
            fetx(ixr,iy) = fetx(ixr,iy) + feex(ixr,iy) + feix(ixr,iy)
 15      continue
       enddo  ! end do-loop over nxpt mesh regions
      endif   ! test on ckinfl-1

      pvmomcx = 0.e0
      ptjdote = 0.e0

      do 31 ix=1,nx
         do 30 iy=1,ny
            ix1 = ixm1(ix,iy)
            ix2 = ixp1(ix,iy)
! Here peirad includes sum of electron and ion energy losses; note that binding
! energy is included in eeli term, but it is carried by the ions.
! Note also that eion and ediss generally balance in the next line
! because should have ediss=2*eion - transfer from electron to ion energy
            peirad(ix,iy) = cnsor*( erliz(ix,iy) + erlrc(ix,iy) +
     .                              ebind*ev*psor(ix,iy,1) -
     .                              ebind*ev*psorrg(ix,iy,1) +
     .                              ediss*ev*(0.5*psordis(ix,iy)) -
     .                        ceisor*eion*ev*(psordis(ix,iy)) )

! other energy diagnostics are given below

cc            jdote(ix,iy) = -   ! this energy is included in resee, not lost
cc     .                  0.5 * fqx(ix ,iy)*(phi(ix2,iy  )+phi(ix ,iy)) +
cc     .                  0.5 * fqx(ix1,iy)*(phi(ix ,iy  )+phi(ix1,iy)) -
cc     .                  0.5 * fqy(ix ,iy)*(phi(ix ,iy+1)+phi(ix ,iy)) +
cc     .                  0.5 * fqy(ix,iy-1)*(phi(ix,iy)+phi(ix,iy-1))
cc            ptjdote = ptjdote + jdote(ix,iy)
            ptjdote = ptjdote + wjdote(ix,iy)

            if (isupgon(1) .eq. 0) then
               pmomv(ix,iy)=cngmom(1)*up(ix,iy,1)*sx(ix,iy)*rrv(ix,iy)*
     .                           ( ng(ix2,iy,1)*tg(ix2,iy,1)-
     .                             ng(ix ,iy,1)*tg(ix ,iy,1) ) +
     .             cmwall(1)*0.125*mi(1)*(up(ix,iy,1)+up(ix1,iy,1))**2*
     .                ng(ix,iy,1)*nucx(ix,iy,1)*vol(ix,iy)
            elseif (isupgon(1) .eq. 1) then    ! inertial neutrals
               pmomv(ix,iy) = 0.  ! coupled back to therm eng for inertial neut
            endif
            pvmomcx = pvmomcx + pmomv(ix,iy)

            engerr(ix,iy) = ( fetx(ix1,iy  )-fetx(ix,iy)+
     .                        fety(ix ,iy-1)-fety(ix,iy)-
     .                        peirad(ix,iy)-png2ni(ix,iy) ) /
     .                         abs(pwrin)
            if (isimpon.ne.0) then  ! prad allocated only if isimpon.ne.0
               engerr(ix,iy) = engerr(ix,iy) - prad(ix,iy)*vol(ix,iy)/
     .                                                     abs(pwrin)
            endif
 30      continue
 31   continue

! ionization and background sources

      iion_tot = 0.e0
      irecomb_tot = 0.e0
      icxgas_tot = 0.e0
      pradrc = 0.e0
      pradiz = 0.e0
      pradht = 0.e0
      prdiss = 0.e0
      pibirth = 0.e0
      pbinde = 0.e0
      pbindrc = 0.e0
      pradzbind = 0.e0
      do igsp = 1, ngsp
         iion(igsp) = 0.e0
         irecomb(igsp) = 0.e0
         icxgas(igsp) = 0.e0
      enddo
      do igsp = 1, max(1, ngsp-nhgsp)
         pradimpt(igsp) = 0.
         do iimp = 0, nzsp(igsp)
            pradimp(iimp,igsp) = 0.
         enddo
      enddo
      pradfft = 0.

      do iy = 1, ny
         do ix = 1, nx
            do igsp = 1, ngsp
              if (ishymol.eq.0 .or. igsp.ne.2) then
               iion(igsp) = iion(igsp) - cnsor*qe*psorg(ix,iy,igsp)
               irecomb(igsp) = irecomb(igsp) -cnsor*qe*psorrg(ix,iy,igsp)
               icxgas(igsp) = icxgas(igsp) - qe*psorcxg(ix,iy,igsp)
              endif
            enddo
         enddo
      enddo

      do igsp = 1, ngsp
         iion_tot = iion_tot + iion(igsp)
         irecomb_tot = irecomb_tot + irecomb(igsp)
         icxgas_tot = icxgas_tot + icxgas(igsp)
      enddo

      do iy = 1, ny
         do ix = 1, nx
            ix1 = ixm1(ix,iy)
            pradrc = pradrc + cnsor*erlrc(ix,iy)
            pradiz = pradiz + (eeli(ix,iy)-ebind*ev) * psor(ix,iy,1)
            pbinde = pbinde + ebind*ev * psor(ix,iy,1)
            pbindrc = pbindrc + ebind*ev*psorrg(ix,iy,1)
            prdiss = prdiss + ediss * ev * (0.5*psordis(ix,iy))
            pibirth = pibirth + ceisor* eion*ev * (psordis(ix,iy)) -
     .                ccoldsor*ng(ix,iy,1)*(1.5*ti(ix,iy)-eion*ev)*
     .                                          nucx(ix,iy,1)*vol(ix,iy)
         enddo
      enddo
      pradht = pradiz + pradrc

      if (isimpon .eq. 2 .or. isimpon .eq. 7) then !fixed fraction model
         do iy = 1, ny
            do ix = 1, nx
               pradfft = pradfft + pradcff(ix,iy)*vol(ix,iy)
            enddo
         enddo
      endif

      if (isimpon .gt. 2) then  !separate impurity species
         do 39 iy = 1, ny
            do 38 ix = 1, nx
               do igsp = nhgsp+1, ngsp
                  jz = igsp - nhgsp
                  do iimp = 0, nzsp(jz)
                     pradimp(iimp,jz) = pradimp(iimp,jz) +
     .                             pradz(ix,iy,iimp,jz)*vol(ix,iy)
                  enddo
               enddo
            pradzbind = pradzbind + (pwrze(ix,iy)-prad(ix,iy))*
     .                                vol(ix,iy) ! only total pradzbind calc
 38         continue
 39      continue
         do igsp = nhgsp+1, ngsp
            jz = igsp - nhgsp
            do iimp = 0, nzsp(jz)
               pradimpt(jz) = pradimpt(jz) + pradimp(iimp,jz)
            enddo
         enddo
      endif
cccac$omp end parallel
      return
      end

******* end of subroutine engbal *******

c----------------------------------------------------------------------c
      subroutine pradpltwl

            use Dim             ! nx,ny
            use Share           ! nxomit
            use Imprad          ! isimpon, prad
            use Comgeo          ! sx, sy, vol
            use Noggeo          ! angfx
            use RZ_grid_info    ! rm, zm
            use Phyvar          ! pi, ev
            use Xpoint_indices  ! ixlb, ixrb
            use Conduc          ! eeli
            use Rhsides         ! psor
            use UEpar           ! ebind
            use Postproc 	  !pwr_pltz,pwr_plth,pwr_wallz,pwr_wallh
                          !pwr_pfwallz, pwr_pfwallh

c ... Local variables
      Implicit none

c ... Calc radiation power to divertor and outer wall surfaces
c ... Use as a diagnostic to call from the BASIS parser

      real prdu(0:nx+1,0:ny+1)
      real theta_ray1, theta_ray2, dthgy, dthgx, sxo, frth
      integer ixv, iyv, nj, ix, iy, ip, iodd

! Initialize arrays
      call sfill ((ny+2)*2*nxpt, 0., pwr_pltz, 1)
      call sfill ((ny+2)*2*nxpt, 0., pwr_plth, 1)
      call sfill ((nx+2), 0., pwr_wallz, 1)
      call sfill ((nx+2), 0., pwr_wallh, 1)
      call sfill ((nx+2)*nxpt, 0., pwr_pfwallz, 1)
      call sfill ((nx+2)*nxpt, 0., pwr_pfwallh, 1)

      if (isimpon > 0) then   ! use prdu since prad might not be dimensioned
         call s2copy (nx+2,ny+2, prad,1,nx+2, prdu,1,nx+2) !prad --> prdu
      else
         call s2fill (nx+2, ny+2, 0., prdu, 1, nx+2)
      endif

      nj = nxomit

c ... First do the divertor plate surfaces
      do ip = 1, 2*nxpt
        iodd = (ip+1)/2 - ip/2  ! =1 if ip odd, =0 if ip even
	if (iodd==1) then
          ixv = ixlb(ip/2+1)    ! viewing ix position
	else
	  ixv = ixrb(ip/2) + 1
	endif
        do iyv = 1, ny		! loop over viewing ix
          do iy = 1, ny	        ! loop over source iy
            do ix = 1, nx	! loop over source ix
              theta_ray1 = atan2( zm(ixv+nj,iyv,1)-zm(ix+nj,iy,0),
     .                            rm(ixv+nj,iyv,1)-rm(ix+nj,iy,0) )
              theta_ray2 = atan2( zm(ixv+nj,iyv,3)-zm(ix+nj,iy,0),
     .                            rm(ixv+nj,iyv,3)-rm(ix+nj,iy,0) )
              dthgy = abs(theta_ray1-theta_ray2)
              frth = min(dthgy, 2*pi-dthgy)/(2*pi)  ! frac.; need angle < pi
              sxo = sx(ixv,iyv)/(cos(angfx(ixv,iyv)))
              pwr_pltz(iyv,ip) = pwr_pltz(iyv,ip) +
     .                                 prdu(ix,iy)*vol(ix,iy)*frth/sxo
              pwr_plth(iyv,ip) = pwr_plth(iyv,ip) + (
     .                           (eeli(ix,iy)-ebind*ev)*psor(ix,iy,1)
     .                                        + erlrc(ix,iy))*frth/sxo
            enddo
          enddo
        enddo
c ... Set corner values
      pwr_pltz(0,ip)    = pwr_pltz(1,ip)
      pwr_pltz(ny+1,ip) = pwr_pltz(ny,ip)
      pwr_plth(0,ip)    = pwr_plth(1,ip)
      pwr_plth(ny+1,ip) = pwr_plth(ny,ip)

      enddo             ! end of ip loop over divertor plates


c ... Now do the "outer" wall surface, i.e., iy=ny+1
      iyv = ny+1                ! viewing iy position
      do ixv = 1, nx		! loop over viewing ix
        do iy = 1, ny	        ! loop over source iy
          do ix = 1, nx	        ! loop over source ix
            theta_ray1 = atan2( zm(ixv+nj,iyv,1)-zm(ix+nj,iy,0),
     .                          rm(ixv+nj,iyv,1)-rm(ix+nj,iy,0) )
            theta_ray2 = atan2( zm(ixv+nj,iyv,2)-zm(ix+nj,iy,0),
     .                          rm(ixv+nj,iyv,2)-rm(ix+nj,iy,0) )
            dthgx = abs(theta_ray1-theta_ray2)
            frth = min(dthgx, 2*pi-dthgx)/(2*pi)  ! frac; need angle < pi
            pwr_wallz(ixv) = pwr_wallz(ixv) +
     .                           prdu(ix,iy)*vol(ix,iy)*frth/sy(ixv,iyv)
            pwr_wallh(ixv) = pwr_wallh(ixv) + (
     .                       (eeli(ix,iy)-ebind*ev)*psor(ix,iy,1) +
     .                          erlrc(ix,iy) )*frth/sy(ixv,iyv)
          enddo
        enddo
      enddo
      pwr_wallz(0) = pwr_wallz(1)	! Because prad(0,) is not calculated
      pwr_wallz(nx+1) = pwr_wallz(nx)
      pwr_wallh(0) = pwr_wallh(1)
      pwr_wallh(nx+1) = pwr_wallh(nx)

c ... Finally do the private-flux wall surfaces, i.e., iy=0
      iyv = 0                        ! viewing iy position
      do ip = 1, nxpt                ! loop over number of x-points
        do ixv = 1, nx   	     ! loop over viewing ix
          do iy = 1, ny	             ! loop over source iy
            do ix = 1, nx	     ! loop over source ix
              theta_ray1 = atan2( zm(ixv+nj,iyv,1)-zm(ix+nj,iy,0),
     .                          rm(ixv+nj,iyv,1)-rm(ix+nj,iy,0) )
              theta_ray2 = atan2( zm(ixv+nj,iyv,2)-zm(ix+nj,iy,0),
     .                          rm(ixv+nj,iyv,2)-rm(ix+nj,iy,0) )
              dthgx = abs(theta_ray1-theta_ray2)
              frth = min(dthgx, 2*pi-dthgx)/(2*pi)  ! frac; need angle < pi
              pwr_pfwallz(ixv,ip) = pwr_pfwallz(ixv,ip) +
     .                           prdu(ix,iy)*vol(ix,iy)*frth/sy(ixv,iyv)
              pwr_pfwallh(ixv,ip) = pwr_pfwallh(ixv,ip) + (
     .                       (eeli(ix,iy)-ebind*ev)*psor(ix,iy,1) +
     .                          erlrc(ix,iy) )*frth/sy(ixv,iyv)
            enddo
          enddo
          if(ixv>ixpt1(ip) .and. ixv <ixpt2(ip)+1) then  ! 0 in core region
            pwr_pfwallh(ixv,ip) = 0.
            pwr_pfwallh(ixv,ip) = 0.
            pwr_pfwallz(ixv,ip) = 0.
            pwr_pfwallz(ixv,ip) = 0.
          endif
        enddo
        pwr_pfwallz(0,ip) = pwr_pfwallz(1,ip)  ! prad(0,) is not calculated
        pwr_pfwallz(nx+1,ip) = pwr_pfwallz(nx,ip)
        pwr_pfwallh(0,ip) = pwr_pfwallh(1,ip)
        pwr_pfwallh(nx+1,ip) = pwr_pfwallh(nx,ip)
      enddo

      return
      end
******* end of subroutine pradpltwl *******
c-------------------------------------------------------------------------c

      subroutine plateflux

            use Dim             ! nx,ny
            use Comgeo          ! sx,vol,gx,gxf
            use Noggeo          ! angfx
            use Phyvar          ! pi, ev
            use Postproc        ! pwr_plth,Pwr_pltz,sdel,rb;sdil,rb; sdbindl,rb,
                          ! sdtl,rb;gdil,rd
            use Xpoint_indices  ! ixlb, ixrb
            use Compla          ! mi,ni,up,te,ti
            use Comflo          ! fnix,feex,feix
            use Conduc          ! visx,hcxn
            use UEpar           ! ebind
            use Bcond           ! ckinfl
            use Parallv         ! nxg,nyg
            use Poten           ! phi0l,phi0r

      Implicit none

c ... Calc particle & heat fluxes to divertor plates
c ... Use as diagnostic called from BASIS/Python parser (alt. to balancee)

c  Local variables
      integer iu,jx,id,ixi,ixo,ix,iy
      real tempvo,tempvi
      real pdivilb(1:nxpt),pdivirb(1:nxpt),pdivelb(1:nxpt),
     .     pdiverb(1:nxpt)
      real pbindlb(1:nxpt),pbindrb(1:nxpt),pdivnlb(1:nxpt),
     .     pdivnrb(1:nxpt)
      real sdilbd(0:ny+1,1:nfsp,1:nxpt),sdirbd(0:ny+1,1:nfsp,1:nxpt)
      real sdnlb(0:ny+1,1:nxpt),sdnrb(0:ny+1,1:nxpt)
      real sxi(0:ny+1,1:nxpt),sxo(0:ny+1,1:nxpt)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First do the particle flux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do jx=1,nxpt
        ixi=ixlb(jx)	! ixi=0
        ixo=ixrb(jx)	! ixo=nx
        do iy=1,ny+1
          sxo(iy,jx) = sx(ixo,iy)/(cos(angfx(ixo,iy)))
          sxi(iy,jx) = sx(ixi,iy)/(cos(angfx(ixi,iy)))
          do id = 1, nfsp
	    gdilb(iy,id,jx) = -fnix(ixi,iy,id)/sxi(iy,jx)
	    gdirb(iy,id,jx) =  fnix(ixo,iy,id)/sxo(iy,jx)
            engilb(iy,id,jx) = (2.*ti(ixi,iy)/ev + zi(id)*
     .                          (phi(ixi,iy)-phi0l(iy,jx)) )
            engirb(iy,id,jx) = 2.*ti(ixo+1,iy)/ev + zi(id)*
     .                          (phi(ixo+1,iy)-phi0r(iy,jx))
          enddo
        enddo
      enddo

! Fix corner boundary values
      do jx = 1, nxpt
        do id = 1, nfsp
          gdilb(0,id,jx) =     gdilb(1,id,jx)
          gdilb(ny+1,id,jx) =  gdilb(ny,id,jx)
          gdirb(0,id,jx) =     gdirb(1,id,jx)
          gdirb(ny+1,id,jx) =  gdirb(ny,id,jx)
          engilb(0,id,jx) =    engilb(1,id,jx)
          engilb(ny+1,id,jx) = engilb(ny,id,jx)
          engirb(0,id,jx) =    engirb(1,id,jx)
          engirb(ny+1,id,jx) = engirb(ny,id,jx)
        enddo
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now do the heat flux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

! Get radiation power in pwr_plth and pwr_plt
      call pradpltwl

! Fill radiation heat flux arrays and initialize tot powers
      do jx = 1, nxpt
        pdivirb(jx) = 0.
        pdiverb(jx) = 0.
        pdivilb(jx) = 0.
        pdivelb(jx) = 0.
        pbindrb(jx) = 0.
        pbindlb(jx) = 0.
        do iy = 0, ny+1
          iu = 2*(jx/2)+1  ! iu/iu+1 gives "odd/even" plates (in/out)
          sdrlb(iy,jx) = pwr_plth(iy,iu)+pwr_pltz(iy,iu)
          sdrrb(iy,jx) = pwr_plth(iy,iu+1)+pwr_pltz(iy,iu+1)
        enddo
      enddo

! here the sds are ion and electron poloidal power fluxes in W/m**2
      do jx=1,nxpt
        ixi=ixlb(jx)	! ixi=0
        ixo=ixrb(jx)	! ixo=nx
        do iy=1,ny+1
          sdirb(iy,jx) = 0.
          sdilb(iy,jx) = 0.
          do id = 1, nfsp
           if (zi(id) .gt. 0) then
	    sdirbd(iy,id,jx) = ( 0.5*mi(id)*upi(ixo,iy,id)**2*
     .                           fnix(ixo,iy,id) )/sxo(iy,jx)
	    sdilbd(iy,id,jx) = (-0.5*mi(id)*upi(ixi,iy,id)**2*
     .                           fnix(ixi,iy,id) )/sxi(iy,jx)
           else    ! note: upi=0 for neutral species; use up instead
	    sdirbd(iy,id,jx) = ( 0.5*mi(id)*up(ixo,iy,id)**2*
     .                           fnix(ixo,iy,id) )/sxo(iy,jx)
	    sdilbd(iy,id,jx) = (-0.5*mi(id)*up(ixi,iy,id)**2*
     .                           fnix(ixi,iy,id) )/sxi(iy,jx)
           endif
           sdirb(iy,jx) = sdirb(iy,jx) + sdirbd(iy,id,jx)
           sdilb(iy,jx) = sdilb(iy,jx) + sdilbd(iy,id,jx)
         enddo
         do id = 1, nusp      ! note: up neutral species in nonzero
           tempvo =  - ckinfl*0.5*sx(ixo,iy)*visx(ixo,iy,id)*gx(ixo,iy)*
     .               ( up(ixo,iy,id)**2 - up(ixo-1,iy,id)**2 ) /sxo(iy,jx)
	   tempvi =  + ckinfl*0.5*sx(ixi,iy)*visx(ixi+1,iy,id)*gx(ixi+1,iy)*
     .                ( up(ixi+1,iy,id)**2 - up(ixi,iy,id)**2 ) / sxi(iy,jx)
	   sdirbd(iy,id,jx) = sdirbd(iy,id,jx) + tempvo
	   sdirb(iy,jx) = sdirb(iy,jx) + tempvo
	   sdilbd(iy,id,jx) = sdilbd(iy,id,jx) + tempvi
	   sdilb(iy,jx) = sdilb(iy,jx) + tempvi
         enddo
         sdirb(iy,jx) = sdirb(iy,jx) + feix(ixo,iy)/sxo(iy,jx)
         sbindrb(iy,jx) =  fnix(ixo,iy,1) * ebind*ev / sxo(iy,jx)
         sderb(iy,jx) =  ( feex(ixo,iy)+fqx(ixo,iy)*
     .                     (phi(ixo+1,iy)-phi0r(iy,jx)) )/sxo(iy,jx)
         sdtrb(iy,jx) = sderb(iy,jx) + sdirb(iy,jx) + sbindrb(iy,jx)
         sdilb(iy,jx) = sdilb(iy,jx) - feix(ixi,iy)/sxi(iy,jx)
         sbindlb(iy,jx) = -fnix(ixi,iy,1) * ebind*ev / sxi(iy,jx)
         sdelb(iy,jx) = -( feex(ixi,iy)+fqx(ixi,iy)*
     .                     (phi(ixi  ,iy)-phi0l(iy,jx)) )/sxi(iy,jx)
         sdtlb(iy,jx) = sdelb(iy,jx) + sdilb(iy,jx) + sbindlb(iy,jx)
         pdivirb(jx) = pdivirb(jx) + sdirb(iy,jx)*sxo(iy,jx)
         pdiverb(jx) = pdiverb(jx) + sderb(iy,jx)*sxo(iy,jx)
         pdivilb(jx) = pdivilb(jx) + sdilb(iy,jx)*sxi(iy,jx)
         pdivelb(jx) = pdivelb(jx) + sdelb(iy,jx)*sxi(iy,jx)
         pbindrb(jx) = pbindrb(jx) + sbindrb(iy,jx)*sxo(iy,jx)
         pbindlb(jx) = pbindlb(jx) + sbindlb(iy,jx)*sxi(iy,jx)
         if (isupgon(1).eq.1) then    ! Approx. neutral energy flux
	   sdnrb(iy,jx)=sdirbd(iy,jx,2) + ( sx(ixo,iy)*hcxn(ixo,iy)*
     .                    gxf(ixo,iy)*(ti(ixo,iy)-ti(ixo+1,iy)) +
     .   		  2.5*fnix(ixo,iy,2)*ti(ixo+1,iy) ) / sxo(iy,jx)
	   sdnlb(iy,jx)=sdilbd(iy,jx,2) - ( sx(ixi,iy)*hcxn(ixi,iy)*
     .                    gxf(ixi,iy)*(ti(ixi,iy)-ti(ixi+1,iy)) +
     .                    2.5*fnix(ixi,iy,2)*ti(ixi  ,iy) ) / sxi(iy,jx)
	   pdivnrb(jx) = pdivnrb(jx) + sdnrb(iy,jx)*sxo(iy,jx)
	   pdivnlb(jx) = pdivnlb(jx) + sdnlb(iy,jx)*sxi(iy,jx)
        endif
        enddo  ! end do-loop for iy=1,ny+1
      enddo  ! end do-loop for jx=1,nxpt

! Fix corner boundary values
      do jx = 1, nxpt
         sdtlb(0,jx) = sdtlb(1,jx)
         sdtlb(ny+1,jx) = sdtlb(ny,jx)
         sdelb(0,jx) = sdelb(1,jx)
         sdelb(ny+1,jx) = sdelb(ny,jx)
         sdilb(0,jx) = sdilb(1,jx)
         sdilb(ny+1,jx) = sdilb(ny,jx)
         sdrlb(ny+1,jx) = sdrlb(ny,jx)
         sdtrb(0,jx) = sdtrb(1,jx)
         sdtrb(ny+1,jx) = sdtrb(ny,jx)
         sderb(0,jx) = sderb(1,jx)
         sderb(ny+1,jx) = sderb(ny,jx)
         sdirb(0,jx) = sdirb(1,jx)
         sdirb(ny+1,jx) = sdirb(ny,jx)
         sdrrb(ny+1,jx) = sdrrb(ny,jx)
      enddo

      return
      end
******* end of subroutine plateflux *******
c-------------------------------------------------------------------------c

      subroutine wallflux

            use Dim             ! nx,ny,nxpt
            use Comgeo          ! sy
            use Phyvar          ! pi, ev
            use Postproc        ! pwr_wallh, Pwr_wallz
                          ! swallr, swalli, swalle, swbind, swallt
                          ! spfwallr
            use Compla          ! mi,ni,up,te,ti
            use Comflo          ! fniy,feey,feiy
            use UEpar           ! ebind

      Implicit none

c ... Calc particle & heat fluxes to outer wall surfaces; only PF rad flux?
c ... Use as diagnostic called from BASIS/Python parser (alt. to balancee)


c  Local variables
      integer jx,id,ixi,ixo,ip,ix,iy
      real pwallr,pwalli,pwalle,pwbind
      real swallid(0:nx+1,nfsp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First do the wall particle flux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do ix = 1, nx
        do id = 1, nfsp
	  gwalli(ix,id) = fniy(ix,ny,id)/sy(ix,ny)
          engwalli(ix,id) = 2.*ti(ix,ny+1)/ev + zi(id)*phi(ix,ny+1)
        enddo
      enddo

! Fix corner boundary values
      do id = 1, nfsp
        gwalli(0,id) = gwalli(1,id)
        gwalli(nx+1,id) = gwalli(nx,id)
        engwalli(0,id) = engwalli(1,id)
        engwalli(nx+1,id) = engwalli(nx,id)
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now do the wall heat flux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Initialize total powers (local use only)
      pwallr = 0.
      pwalli = 0.
      pwalle = 0.
      pwbind = 0.

! Get radiation power in pwr_plth and pwr_plt
      call pradpltwl

! Fill radiation heat flux arrays
      do ix = 0, nx+1
	swallr(ix) = pwr_wallh(ix) + pwr_wallz(ix)
        do ip = 1, nxpt
          spfwallr(ix,ip) = pwr_pfwallz(ix,ip)+pwr_pfwallh(ix,ip)
        enddo
      enddo

! swalls are ion and electron radial power fluxes in W/m**2
      do ix = 1, nx
        swalli(ix) = 0.
        do id = 1, nfsp
          if (zi(id) .gt. 0) then
            swallid(ix,id) = ( 0.5*mi(id)*upi(ix,ny,id)**2*
     .                        fniy(ix,ny,id) )/sy(ix,ny)
          else    ! note: upi=0 for neutral species; use up instead
	    swallid(ix,id) = ( 0.5*mi(id)*up(ix,ny,id)**2*
     .                         fniy(ix,ny,id) )/sy(ix,ny)
          endif
          swalli(ix) = swalli(ix) + swallid(ix,id)
        enddo
        swalli(ix) = swalli(ix) + feiy(ix,ny)/sy(ix,ny)
        swbind(ix) = fniy(ix,ny,1) * ebind*ev / sy(ix,ny)
        swalle(ix) = (feey(ix,ny)+fqy(ix,ny)*phi(ix,ny+1) )/sy(ix,ny)
        swallt(ix) = swalle(ix) + swalli(ix) + swbind(ix) + swallr(ix)
        pwallr = pwallr + swallr(ix)*sy(ix,ny)
        pwalle = pwalle + swalle(ix)*sy(ix,ny)
        pwalli = pwalli + swalli(ix)*sy(ix,ny)
        pwbind = pwbind + swbind(ix)*sy(ix,ny)
      enddo  ! end do-loop for ix=1,nx

! Fix corner boundary values
      swallr(0) = swallr(1)
      swalli(0) = swalli(1)
      swalle(0) = swalle(1)
      swbind(0) = swbind(1)
      swallt(0) = swallt(1)
      swallr(nx+1) = swallr(nx)
      swalli(nx+1) = swalli(nx)
      swalle(nx+1) = swalle(nx)
      swbind(nx+1) = swbind(nx)
      swallt(nx+1) = swallt(nx)

      return
      end
******* end of subroutine wallflux *******
c-------------------------------------------------------------------------c

      subroutine bbb2wdf
      use Dim 		! nx,ny
      use Share 		! nycore,nysol,nxcore,nxleg,igrid,geometry
      use RZ_grid_info 	! br,bz,bpol,bphi,b
      use Bfield              ! rbfbt
      use Bcond 		! fngysi,fngyso
      use Parallv             ! nxg,nyg
      use Selec 		! ixm1
      use Comflo 		! fnix
      use Compla 		! ni,uu,up,v2,vy,te,ti,ne
      use Linkbbb 		!

c     Local variables --
      implicit none
      integer ix,iy,nunit
      real uuc,upc,vyc,v2c

c     Compute data for output to wdf package

      nxbbb=nx
      nybbb=ny
      nycorebbb=nycore(igrid)
      nysolbbb=nysol(igrid)
      nxleg1bbb=nxleg(igrid,1)
      nxcore1bbb=nxcore(igrid,1)
      nxleg2bbb=nxleg(igrid,2)
      nxcore2bbb=nxcore(igrid,2)
      geometrybbb=geometry

      do ix=0,nx+1
         do iy=0,ny+1
            nibbb(ix,iy)=ni(ix,iy,1)
            tibbb(ix,iy)=ti(ix,iy)
            nebbb(ix,iy)=ne(ix,iy)
            tebbb(ix,iy)=te(ix,iy)
            fnixbbb(ix,iy)=fnix(ix,iy,1)
         enddo
         fngysibbb(ix)=fngysi(ix,1)
         fngysobbb(ix)=fngyso(ix,1)
      enddo

c     Compute ion velocity at cell centers :
      do ix=1,nx
	 do iy=1,ny
	    uuc=0.5*(uu(ix,iy,1)+uu(ixm1(ix,iy),iy,1))
	    vyc=0.5*(vy(ix,iy,1)+vy(ix,iy-1,1))
	    vflowxbbb(ix,iy) =   uuc*br(ix,iy,0)/bpol(ix,iy,0)
     &                       - vyc*bz(ix,iy,0)/bpol(ix,iy,0)
	    vflowzbbb(ix,iy) =   uuc*bz(ix,iy,0)/bpol(ix,iy,0)
     &                       + vyc*br(ix,iy,0)/bpol(ix,iy,0)

	    v2c=0.5*(v2(ix,iy,1)+v2(ixm1(ix,iy),iy,1))
	    upc=0.5*(up(ix,iy,1)+up(ixm1(ix,iy),iy,1))
	    vflowybbb(ix,iy) = - v2c*bpol(ix,iy,0)/b(ix,iy,0)
     &                       + upc*rbfbt(ix,iy)
	 enddo
      enddo

c ********* Write the file to be read by the wdf package *************

      call freeus (nunit)
      open (nunit, file='bbb-wdf', form='unformatted', status='unknown')
      write(nunit)nxbbb,nybbb,nycorebbb,nysolbbb,nxleg1bbb,nxcore1bbb,
     &            nxleg2bbb,nxcore2bbb
      write(nunit)nibbb,tibbb,nebbb,tebbb,vflowxbbb,vflowybbb,vflowzbbb,
     &            fnixbbb,fngysibbb,fngysobbb
      write(nunit)geometrybbb
      close (nunit)

      return
      end

c----------------------------------------------------------------------c

      subroutine scale_mcn
      use Dim 	! nx,ny,nisp,nusp
      use UEpar       ! isupgon,iigsp
      use Comflo 	! fnix
      use MCN_dim 
      use MCN_sources 

      implicit none
      integer ix,iy,ifld,istra

      external remark

c     This subroutine scales plasma source terms obtained from the
c     Monte-Carlo-Neutrals code.  These sources are assumed to scale
c     with the neutral source currents at the divertor plates.
c     The flag that controls the scaling is ismcnvar :
c         ismcnvar=0  -->  sources are fixed
c         ismcnvar=1  -->  sources scale with current
c
c     In the terminology of the EIRENE code, a stratum is a surface
c     or volume element where neutral particles originate.  Here
c     we consider only two possible strata: the inboard and outboard
c     divertor plates.  We neglect strata associated with gas puffing
c     or recombination.
c
c     Plasma source terms from EIRENE are :
c         mcnsor_ni
c         mcnsor_up
c         mcnsor_te
c         mcnsor_ti
c     where mcnsor_ni(ix,iy,ifld,istra) is the particle source
c     for ion fluid 'ifld' at cell (ix,iy) due to neutral source
c     stratum 'istra'.
c
c     mcncurr(istra) is the incident ion current that
c     normalizes the plasma sources due to neutrals
c     that originate at stratum 'istra'.
c
c     The scaled plasma source terms for UEDGE are :
c         uesor_ni * cmneutsor_ni
c         uesor_up * cmneutsor_mi
c         uesor_te * cmneutsor_ee
c         uesor_ti * cmneutsor_ei
c     where uesor_ni is obtained by summing the scaled mcnsor_ni
c     over all strata, and similarly for up,te,ti.

ccc=================================================================
ccc     NOTE: For the cmod-box problem there is only ONE strata !!!
ccc     (in general, we need to devise a test to identify strata)
ccc=================================================================

c     Compute scale factors for each strata :
      do istra=1,nstra
         strascal(istra)=1.	! default is no scaling
      enddo
      if (ismcnvar==1) then	! sources scale with current
         do istra=1,nstra
            uecurr(istra)=0.
            if (istra==1) then	! for east target plate only
               do iy=1,ny
                  uecurr(istra)=uecurr(istra)+fnix(nx,iy,1)
               enddo
            else
               call remark("***")
               call remark("***    subroutine scale_mcn    ***")
               call remark("***  not defined for nstra > 1  ***")
               call remark("***")
            endif
            if(mcncurr(istra) > 0) then
               strascal(istra)=uecurr(istra)/mcncurr(istra)
            endif
         enddo
      endif

c     Scale source terms and sum over strata :
      do iy=0,ny+1
         do ix=0,nx+1
            do ifld=1,nisp
               uesor_ni(ix,iy,ifld)=0.
            enddo
            do ifld=1,nusp
               uesor_up(ix,iy,ifld)=0.
            enddo
            uesor_te(ix,iy)=0.
            uesor_ti(ix,iy)=0.
         enddo
      enddo
      do istra=1,nstra
         do iy=0,ny+1
            do ix=0,nx+1
               do ifld=1,nisp
                  uesor_ni(ix,iy,ifld)=uesor_ni(ix,iy,ifld)
     &                    +mcnsor_ni(ix,iy,ifld,istra)*strascal(istra)
               enddo
               do ifld=1,nusp
                  uesor_up(ix,iy,ifld)=uesor_up(ix,iy,ifld)
     &                    +mcnsor_up(ix,iy,ifld,istra)*strascal(istra)
               enddo
               uesor_te(ix,iy)=uesor_te(ix,iy)
     &                    +mcnsor_te(ix,iy,istra)*strascal(istra)
               uesor_ti(ix,iy)=uesor_ti(ix,iy)
     &                    +mcnsor_ti(ix,iy,istra)*strascal(istra)
            enddo
         enddo
      enddo

      return
      end

!----------------------------------------------------------------------!

      subroutine write30 (fname, runid)
      use Dim             ! nx,ny
      use UEint           ! mhdgeo
      use Share           ! geometry,nxomit
      use Xpoint_indices  ! ixpt1,ixpt2,iysptrx
      use RZ_grid_info    ! rm,zm

      implicit none
      character*(*) fname, runid

      integer nunit, ix, iy, i, iv(4)
      integer nncut,nniso
      integer nxcut1(2),nxcut2(2),nycut1(2),nycut2(2)
      integer nxiso1(2),nxiso2(2),nyiso1(2),nyiso2(2)

      external remark

      data nncut/0/,nniso/0/
      data nxcut1/2*0/,nxcut2/2*0/,nycut1/2*0/,nycut2/2*0/
      data nxiso1/2*0/,nxiso2/2*0/,nyiso1/2*0/,nyiso2/2*0/

! This subroutine writes the geometry data file that is read by the
! EIRENE code, using the B2 code variable definitions.

cccMER NOTE: generalize the following for 2 or more x-points
! Convert UEDGE variables to B2/EIRENE data -
      if (mhdgeo .eq. 1) then
         nncut=2
         nxcut1(1)=ixpt1(1)
         nxcut2(1)=ixpt2(1)+1
         nycut1(1)=0
         nycut2(1)=iysptrx1(1)
         nxcut1(2)=ixpt2(1)
         nxcut2(2)=ixpt1(1)+1
         nycut1(2)=0
         nycut2(2)=iysptrx1(1)
      endif
! UEDGE cell vertex iv(i) corresponds to B2 cell vertex i :
      iv(1)=2
      iv(2)=4
      iv(3)=3
      iv(4)=1

! Write the data -
      call freeus (nunit)
      open (nunit, file=fname, form='formatted', status='unknown')
      write (nunit,*) runid
      write (nunit,*) ' '
      write (nunit,*) nx,ny,nncut
      write (nunit,*)
     &              (nxcut1(i),nxcut2(i),nycut1(i),nycut2(i),i=1,nncut)
      if (nncut .gt. 2) then
         write (nunit,*) nniso
         write (nunit,*)
     &              (nxiso1(i),nxiso2(i),nyiso1(i),nyiso2(i),i=1,nniso)
      endif
      write (nunit,*) ' '
      do ix=1,nx
         do iy=1,ny
            if (mhdgeo .eq. 1) then
               write (nunit,'(4e15.7)') (rm(ix,iy,iv(i)),i=1,4)
               write (nunit,'(4e15.7)') (zm(ix,iy,iv(i)),i=1,4)
            else
               write (nunit,'(4e15.7)') (zm(ix,iy,iv(i)),i=1,4)
               write (nunit,'(4e15.7)') (rm(ix,iy,iv(i)),i=1,4)
            endif
         enddo
      enddo

      close (nunit)
      call remark(" *** geometry file written for EIRENE ***")

      return
      end

!----------------------------------------------------------------------!


      subroutine write31 (fname, runid)
      use Dim             ! nx,ny,nxm,nisp
      use UEpar           ! isupgon
      use Compla          ! ni,uu,up,vy,te,ti,pr,zi
      use Comflo          ! fnix,fniy,feix,feiy,feex,feey
      use RZ_grid_info    ! b
      use Comgeo          ! vol,rr

      implicit none
      character*(*) fname, runid

      integer nunit,nxd,nyd,ifld
      external remark

! This subroutine writes the plasma data file that is read by the
! EIRENE code, using the B2 code variable definitions.  This file only
! includes data for charged species; inertial neutrals are excluded.

      nxd = nx
      nyd = ny

      call freeus (nunit)
      open (nunit, file=fname, form='formatted', status='unknown')
      do ifld = 1, nisp
         if (zi(ifld) .gt. 0)
     .      call gfsub3 (nunit,nx,ny,nxd,nyd,1,ni(0,0,ifld))
      enddo
      do ifld = 1, nisp
         if (zi(ifld) .gt. 0)
     .      call gfsub3 (nunit,nx,ny,nxd,nyd,1,uu(0,0,ifld))
      enddo
      do ifld = 1, nisp
         if (zi(ifld) .gt. 0)
     .      call gfsub3 (nunit,nx,ny,nxd,nyd,1,vy(0,0,ifld))
      enddo
      call gfsub3 (nunit,nx,ny,nxd,nyd,1,te(0,0))
      call gfsub3 (nunit,nx,ny,nxd,nyd,1,ti(0,0))
      call gfsub3 (nunit,nx,ny,nxd,nyd,1,pr(0,0))
      do ifld = 1, nisp
         if (zi(ifld) .gt. 0)
     .      call gfsub3 (nunit,nx,ny,nxd,nyd,1,upi(0,0,ifld))
      enddo
      call gfsub3 (nunit,nx,ny,nxd,nyd,1,rr(0,0))
      do ifld = 1, nisp
         if (zi(ifld) .gt. 0)
     .      call gfsub3 (nunit,nx,ny,nxd,nyd,1,fnix(0,0,ifld))
      enddo
      do ifld = 1, nisp
         if (zi(ifld) .gt. 0)
     .      call gfsub3 (nunit,nx,ny,nxd,nyd,1,fniy(0,0,ifld))
      enddo
      call gfsub3 (nunit,nx,ny,nxd,nyd,1,feix(0,0))
      call gfsub3 (nunit,nx,ny,nxd,nyd,1,feiy(0,0))
      call gfsub3 (nunit,nx,ny,nxd,nyd,1,feex(0,0))
      call gfsub3 (nunit,nx,ny,nxd,nyd,1,feey(0,0))
      call gfsub3 (nunit,nx,ny,nxd,nyd,1,vol(0,0))

      nxd = nxm		! note dimension nxm rather than nx for b(,,)
      call gfsub3 (nunit,nx,ny,nxd,nyd,1,b(0,0,0))

      close (nunit)
      call remark(" *** background plasma file written for EIRENE ***")

      return
      end

!----------------------------------------------------------------------!

      subroutine gfsub3 (nunit,nx,ny,ndimx,ndimy,ndims,dummy)
      implicit none
      integer nunit,nx,ny,ndimx,ndimy,ndims
      real dummy(0:ndimx+1,0:ndimy+1,1:ndims)

      integer nd1,lim,is,iy,ix,iii
      real eps
      data eps/1.0e-90/

! This subroutine writes an array 'dummy' into an output file.
! (Copied from the B2 code by MER 95/12/12; modified 97/03/11)
! The additive constant, eps, ensures that extremely small non-zero
! values of array elements will not corrupt the e-format in the output file.

      nd1 = nx + 2
      lim = (nd1/5)*5 - 4
      do is=1,ndims
         do iy=0,ny+1
            do ix=1,lim,5
               write(nunit,910) (dummy(ix-2+iii,iy,is)+eps,iii=1,5)
            enddo
            if ( (lim+4) .lt. nd1 ) ! write partial line
     .         write(nunit,910) (dummy(ix-1,iy,is)+eps,ix=lim+5,nd1)
         enddo
      enddo

 910  format(5(e16.8))

      return
      end

!----------------------------------------------------------------------!

      subroutine write_eirene
      implicit none

! This subroutine writes the geometry (fort.30) and plasma (fort.31)
! files that supply information to the EIRENE Monte Carlo neutrals code.


      call write30("fort.30","UEDGE geometry data")
      call write31("fort.31","UEDGE plasma data")

      return
      end

!----------------------------------------------------------------------!

      subroutine read32(fname)
      use Dim 		! nx,ny
      use MCN_dim 
      use MCN_sources 
      implicit none
      character*(*) fname

      integer nunit,ix,iy,ifl,istra

c     Read data from EIRENE code output data file 'fort.32':

c     NOTE: Dimensions nxf and nyf should be read from
c            EIRENE code output data file 'fort.44'
c     Here, we assume nxf=nx and nyf=ny elements in ix and iy.
c     Also, nstra and nfl need to be set properly.

      call freeus (nunit)
      open (nunit, file=fname, form='formatted', status='unknown')

      do istra=1,nstra

c     Normalization constants:
         read (nunit,*) wsor(istra), esor(istra)

c     Normalized plasma source terms:
         do ifl=1,nfl
            read (nunit,*) ((sni(ix,iy,ifl,istra), ix=1,nx), iy=1,ny)
            read (nunit,*) ((smo(ix,iy,ifl,istra), ix=1,nx), iy=1,ny)
         enddo
         read (nunit,*) ((see(ix,iy,istra), ix=1,nx), iy=1,ny)
         read (nunit,*) ((sei(ix,iy,istra), ix=1,nx), iy=1,ny)

      enddo

      close (nunit)
      call remark(" *** plasma sources read from file fort.32 ***")

      return
      end

c----------------------------------------------------------------------c

      subroutine writemcnfile(fname, runidtag)

      use Dim              ! nx,ny,nxpt,nisp
      use Xpoint_indices   ! ixlb,ixpt1,ixmdp,ixpt2,ixrb,iysptrx1,iysptrx2
      use RZ_grid_info     ! rm,zm,br,bz,bphi
      use Ext_neutrals  		 ! ext_verbose
      implicit none
      character*(*) fname, runidtag

c     This subroutine writes the UEDGE mesh and plasma data file that
c     is read by the DEGAS2 code.  Plasma velocities are converted to
c     Cartesian or cylindrical components at density cell centers.

      integer nunit,ix,iy,n,jx

      call freeus (nunit)
      open (nunit, file=fname, form='formatted', status='unknown')

c write geometry data
      write (nunit,*) runidtag
      write (nunit,*) nx,ny,nxpt
      do jx=1,nxpt
         write (nunit,*) iysptrx1(jx),iysptrx2(jx)
         write (nunit,*) ixlb(jx),ixpt1(jx),ixmdp(jx),ixpt2(jx),ixrb(jx)
      enddo
      write (nunit,*) (((rm(ix,iy,n),ix=1,nx),iy=1,ny),n=0,4)
      write (nunit,*) (((zm(ix,iy,n),ix=1,nx),iy=1,ny),n=0,4)
      write (nunit,*) (((br(ix,iy,n),ix=1,nx),iy=1,ny),n=0,4)
      write (nunit,*) (((bz(ix,iy,n),ix=1,nx),iy=1,ny),n=0,4)
      write (nunit,*) (((bphi(ix,iy,n),ix=1,nx),iy=1,ny),n=0,4)

c write plasma data
      call gchange("MCN_bkgd",0)
      call writemcnbkgd(nunit)

      close (nunit)
      if (ext_verbose)
     . call remark(" *** data file "//fname//" written for DEGAS2 ***")

      return
      end

c----------------------------------------------------------------------c

      subroutine writemcnbkgd(nunit)
      use Dim              ! nx,ny,nxpt,nisp
      use Xpoint_indices   ! ixlb,ixrb
      use RZ_grid_info     ! rm,zm,br,bz,bphi
      use Bfield           ! rbfbt
      use Comgeo           ! rr
      use Compla           ! zi,ni,upi,uu,vy,v2,te,ti
      use Comflo           ! fnix
      use MCN_bkgd         ! various velocities
      use UEint            ! mhdgeo
      use Selec            ! ixm1,ixp1

      implicit none
      integer nunit

c local variables --
      integer ix,iy,ifld,jx,ixt,ixt1
      real sinb1,cosb1


c compute UEDGE flow velocity at cell centers:
      do ifld=1,nisp
         do ix=1,nx
            do iy=1,ny
               ! 3 orthogonal components of the velocity are:
               v2c(ix,iy,ifld)=(v2(ixm1(ix,iy),iy,ifld)+v2(ix,iy,ifld))/2.
               vyc(ix,iy,ifld)=(vy(ix,iy-1,ifld)+vy(ix,iy,ifld))/2.
               upc(ix,iy,ifld)=(upi(ixm1(ix,iy),iy,ifld)+upi(ix,iy,ifld))/2.
               ! from these we construct the poloidal and toroidal components:
               uuc(ix,iy,ifld)=upc(ix,iy,ifld)*rr(ix,iy)
     .                         +v2c(ix,iy,ifld)*rbfbt(ix,iy)
               utc(ix,iy,ifld)=upc(ix,iy,ifld)*rbfbt(ix,iy)
     .                         -v2c(ix,iy,ifld)*rr(ix,iy)
            enddo
         enddo
      enddo
      if (mhdgeo .eq. 1) then		! convert to cylindrical coordinates:
         do ifld=1,nisp
            do ix=1,nx
               do iy=1,ny
	          cosb1=br(ix,iy,0)/bpol(ix,iy,0)
	          sinb1=bz(ix,iy,0)/bpol(ix,iy,0)
	          vr(ix,iy,ifld)=uuc(ix,iy,ifld)*cosb1-vyc(ix,iy,ifld)*sinb1
	          vz(ix,iy,ifld)=uuc(ix,iy,ifld)*sinb1+vyc(ix,iy,ifld)*cosb1
	          vphi(ix,iy,ifld)=utc(ix,iy,ifld)
               enddo
            enddo
         enddo
      elseif (mhdgeo .eq. -1) then	! cartesian slab coordinates:
         do ifld=1,nisp
            do ix=1,nx
               do iy=1,ny
                  vr(ix,iy,ifld)=vyc(ix,iy,ifld)
	          vz(ix,iy,ifld)=uuc(ix,iy,ifld)
	          vphi(ix,iy,ifld)=utc(ix,iy,ifld)
               enddo
            enddo
         enddo
      endif

c write plasma data (charged species only)
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) ((ni(ix,iy,ifld),ix=1,nx),iy=1,ny)
         endif
      enddo
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) ((vr(ix,iy,ifld),ix=1,nx),iy=1,ny)
         endif
      enddo
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) ((vz(ix,iy,ifld),ix=1,nx),iy=1,ny)
         endif
      enddo
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) ((vphi(ix,iy,ifld),ix=1,nx),iy=1,ny)
         endif
      enddo
      write (nunit,*) ((ti(ix,iy),ix=1,nx),iy=1,ny)
      write (nunit,*) ((te(ix,iy),ix=1,nx),iy=1,ny)

c loop over nxpt mesh regions
      do jx=1,nxpt
c compute UEDGE flow velocity at left target plate:
      do ifld=1,nisp
         do iy=1,ny
            ixt=ixlb(jx)       ! analog of ix=0
            ixt1=ixp1(ixt,iy)  ! analog of ix=1
            ! 3 orthogonal components of the velocity are:
            v2tg1(iy,ifld)=v2(ixt,iy,ifld)
            vytg1(iy,ifld)=(vy(ixt1,iy-1,ifld)+vy(ixt1,iy,ifld))/2. ! <- NOTE
            uptg1(iy,ifld)=upi(ixt,iy,ifld)
            ! from these we construct the poloidal and toroidal components:
            uutg1(iy,ifld)=uptg1(iy,ifld)*rr(ixt,iy)
     .                      +v2tg1(iy,ifld)*rbfbt(ixt,iy)
            uttg1(iy,ifld)=uptg1(iy,ifld)*rbfbt(ixt,iy)
     .                      -v2tg1(iy,ifld)*rr(ixt,iy)
         enddo
      enddo
      if (mhdgeo .eq. 1) then		! convert to cylindrical coordinates:
         do ifld=1,nisp
            do iy=1,ny
               cosb1=br(ixt,iy,0)/bpol(ixt,iy,0)
               sinb1=bz(ixt,iy,0)/bpol(ixt,iy,0)
               vrtg1(iy,ifld)=uutg1(iy,ifld)*cosb1-vytg1(iy,ifld)*sinb1
               vztg1(iy,ifld)=uutg1(iy,ifld)*sinb1+vytg1(iy,ifld)*cosb1
               vphitg1(iy,ifld)=uttg1(iy,ifld)
            enddo
         enddo
      elseif (mhdgeo .eq. -1) then	! cartesian slab coordinates:
         do ifld=1,nisp
            do iy=1,ny
               vrtg1(iy,ifld)=vytg1(iy,ifld)
               vztg1(iy,ifld)=uutg1(iy,ifld)
               vphitg1(iy,ifld)=uttg1(iy,ifld)
            enddo
         enddo
      endif

c write inner target data
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) (ni(ixt,iy,ifld),iy=1,ny)
         endif
      enddo
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) (vrtg1(iy,ifld),iy=1,ny)
         endif
      enddo
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) (vztg1(iy,ifld),iy=1,ny)
         endif
      enddo
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) (vphitg1(iy,ifld),iy=1,ny)
         endif
      enddo
      write (nunit,*) (ti(ixt,iy),iy=1,ny)
      write (nunit,*) (te(ixt,iy),iy=1,ny)
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) (-fnix(ixt,iy,ifld),iy=1,ny)
         endif
      enddo

c compute UEDGE flow velocity at right target plate:
      do ifld=1,nisp
         do iy=1,ny
            ixt=ixrb(jx)+1     ! analog of nx+1
            ixt1=ixm1(ixt,iy)  ! analog of nx
            ! 3 orthogonal components of the velocity are:
            v2tg2(iy,ifld)=v2(ixt1,iy,ifld)
            vytg2(iy,ifld)=(vy(ixt1,iy-1,ifld)+vy(ixt1,iy,ifld))/2. ! <- NOTE
            uptg2(iy,ifld)=upi(ixt1,iy,ifld)
            ! from these we construct the poloidal and toroidal components:
            uutg2(iy,ifld)=uptg2(iy,ifld)*rr(ixt,iy)
     .                      +v2tg2(iy,ifld)*rbfbt(ixt,iy)
            uttg2(iy,ifld)=uptg2(iy,ifld)*rbfbt(ixt,iy)
     .                      -v2tg2(iy,ifld)*rr(ixt,iy)
         enddo
      enddo
      if (mhdgeo .eq. 1) then		! convert to cylindrical coordinates:
         do ifld=1,nisp
            do iy=1,ny
               cosb1=br(ixt,iy,0)/bpol(ixt,iy,0)
               sinb1=bz(ixt,iy,0)/bpol(ixt,iy,0)
               vrtg2(iy,ifld)=uutg2(iy,ifld)*cosb1-vytg2(iy,ifld)*sinb1
               vztg2(iy,ifld)=uutg2(iy,ifld)*sinb1+vytg2(iy,ifld)*cosb1
               vphitg2(iy,ifld)=uttg2(iy,ifld)
            enddo
         enddo
      elseif (mhdgeo .eq. -1) then	! cartesian slab coordinates:
         do ifld=1,nisp
            do iy=1,ny
               vrtg2(iy,ifld)=vytg2(iy,ifld)
               vztg2(iy,ifld)=uutg2(iy,ifld)
               vphitg2(iy,ifld)=uttg2(iy,ifld)
            enddo
         enddo
      endif

c write outer target data
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) (ni(ixt,iy,ifld),iy=1,ny)
         endif
      enddo
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) (vrtg2(iy,ifld),iy=1,ny)
         endif
      enddo
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) (vztg2(iy,ifld),iy=1,ny)
         endif
      enddo
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) (vphitg2(iy,ifld),iy=1,ny)
         endif
      enddo
      write (nunit,*) (ti(ixt,iy),iy=1,ny)
      write (nunit,*) (te(ixt,iy),iy=1,ny)
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) (fnix(ixt1,iy,ifld),iy=1,ny)
         endif
      enddo

      enddo   ! end do-loop over nxpt mesh regions

      return
      end

c----------------------------------------------------------------------c

      subroutine readmcnsor(fname)
      use Dim 		! nx,ny
      use MCN_dim 
      use MCN_sources 
      use Ext_neutrals  		 ! ext_verbose
      implicit none
      character*(*) fname

      integer nunit,ix,iy,ifl,istra

c     Read data from DEGAS2 code output data file 'sources.out':

      call freeus (nunit)
      open (nunit, file=fname, form='formatted', status='unknown')

      do istra=1,nstra

c     Integrated ion particle source:
         read (nunit,*) wsor(istra)

c     Plasma source terms:
         do ifl=1,nfl
            read (nunit,*) ((sni(ix,iy,ifl,istra), ix=1,nx), iy=1,ny)
            read (nunit,*) ((smor(ix,iy,ifl,istra), ix=1,nx), iy=1,ny)
            read (nunit,*) ((smophi(ix,iy,ifl,istra), ix=1,nx), iy=1,ny)
            read (nunit,*) ((smoz(ix,iy,ifl,istra), ix=1,nx), iy=1,ny)
         enddo
         read (nunit,*) ((see(ix,iy,istra), ix=1,nx), iy=1,ny)
         read (nunit,*) ((sei(ix,iy,istra), ix=1,nx), iy=1,ny)

      enddo

      close (nunit)
      if (ext_verbose)
     . call remark(" *** plasma sources read from DEGAS2 file "
     .                                                  //fname//" ***")

      return
      end

!----------------------------------------------------------------------!

      subroutine readmcntest(fname)
      use Dim 
      use MCN_dim 

      implicit none
      character*(*) fname

      integer nunit
c     Read data from DEGAS2 code output file 'testdata.out':

      call freeus (nunit)
      open (nunit, file=fname, form='formatted', status='unknown')

c     Dimensions:
      read (nunit,*) nxf, nyf, nmcsp

      if (nmcsp .gt. nmcmx) then
         call remark('***')
         call remark('*** READMCNTEST: nmcsp > nmcmx')
         call remark('                 re-compile with larger nmcmx')
         call remark('***')
         call xerrab("")
      endif

      call gchange("MCN_test",0)

      call readmcntesta(nunit)


      close(nunit)
      call remark(" *** neutral diagnostics read from DEGAS2 file "
     .                                                  //fname//" ***")

      return
      end

!----------------------------------------------------------------------!

      subroutine readmcntesta(nunit)
      use Dim 
      use MCN_dim 
      use MCN_test 

      implicit none
      integer nunit

      integer ix,iy,id

c Species names:
      read (nunit,*) (labelmc(id), id=1,nmcsp)

c Densities and 'temperatures':
      do id=1,nmcsp
         read (nunit,*) ((nmc(ix,iy,id),ix=1,nxf),iy=1,nyf)
         read (nunit,*) ((tmc(ix,iy,id),ix=1,nxf),iy=1,nyf)
      enddo

c Particle fluxes:
      do id=1,nmcsp
         read (nunit,*) ((fnmcx(ix,iy,id),ix=1,nxf),iy=1,nyf)
         read (nunit,*) ((fnmcy(ix,iy,id),ix=1,nxf),iy=1,nyf)
         read (nunit,*) ((fnmcz(ix,iy,id),ix=1,nxf),iy=1,nyf)
      enddo

      return
      end

!----------------------------------------------------------------------!

      subroutine read44(fname)
      use Dim 
      use MCN_dim 
      use MCN_sources 

      implicit none
      character*(*) fname

      integer nunit

c     Read data from EIRENE code output data file 'fort.44':

      call freeus (nunit)
      open (nunit, file=fname, form='formatted', status='unknown')

c     Dimensions:
      read (nunit,*) nxf, nyf
      read (nunit,*) natmi, nmoli, nioni

c     NOTE:
c     Character arrays (labela,labelm,labeli) in MCN_sources have fixed
c     dimension; Dynamic allocation not allowed on C90 (MER 97/03/11).
      if ( (natmi .gt. nmcmx) .or. (nmoli .gt. nmcmx)
     .                         .or. (nioni .gt. nmcmx) ) then
         call remark('***')
         call remark('*** READ44: natmi or nmoli or nioni > nmcmx')
         call remark('            re-compile with larger nmcmx')
         call remark('***')
         call xerrab("")
      endif

      call gchange("MCN_sources",0)

      call read44a(nunit)

      close (nunit)
      call remark(' *** neutral diagnostics read from file fort.44 ***')

      return
      end

!----------------------------------------------------------------------!

      subroutine read44a(nunit)
      use Dim 
      use MCN_dim 
      use MCN_sources 

      implicit none
      integer nunit

      integer ix,iy,id

c     Species names:
      read (nunit,*) (labela(id), id=1,natmi)
      read (nunit,*) (labelm(id), id=1,nmoli)
      read (nunit,*) (labeli(id), id=1,nioni)

c     Densities and 'temperatures' of atoms, molecules, test ions:
      read (nunit,*) (((naf(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,natmi)
      read (nunit,*) (((taf(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,natmi)
      read (nunit,*) (((nmf(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,nmoli)
      read (nunit,*) (((tmf(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,nmoli)
      read (nunit,*) (((ntf(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,nioni)
      read (nunit,*) (((ttf(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,nioni)

c     Radial particle fluxes:
      read (nunit,*) (((fnay(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,natmi)
      read (nunit,*) (((fnmy(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,nmoli)

c     'parallel' particle flux:
      read (nunit,*) (((fnax(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,natmi)
      read (nunit,*) (((fnmx(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,nmoli)

c     Radial energy flux:
      read (nunit,*) (((feay(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,natmi)
      read (nunit,*) (((femy(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,nmoli)

c     'parallel' energy flux:
      read (nunit,*) (((feax(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,natmi)
      read (nunit,*) (((femx(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,nmoli)

c     H-alpha emissivity from atoms and molecules:
      read (nunit,*) ((hatm(ix,iy),ix=1,nxf),iy=1,nyf)
      read (nunit,*) ((hmol(ix,iy),ix=1,nxf),iy=1,nyf)

      return
      end

!----------------------------------------------------------------------!

      subroutine volave(n1, n2, j2, j5, i2, i5, ixp, ixm, fsprd,
     .                                             vol, ps_tmp, ps)
c ... This subroutine does a volume integral, or average, of cell source-like
c ... quantities by including interpolation based on adjacent-cell quantities

c *** Input:
c        n1 is nx for one of the two array dimensions
c        n2 is ny for the second array dimension
c        j2 and j5 are lower and upper limits of iy loop
c        i2 and i5 are lower and upper limits of ix loop
c        ixp is the 2-D index array ixp1
c        ixm is the 2-D index array ixm1
c        fsprd is the fraction spread to each of 4 adjacent cells
c        vol is the 2-D real array of cell volume
c        pstmp is an (0:n1+1,0:n2+1) array used as a local work array
c        ps is an (0:n1+1,0:n2+1) array containing the input source-like array

c *** Output:
c        ps is an (n1,n2) array containing the result of volume-averaging ps

      implicit none !39

c *** Input and output variables
      integer n1, n2, j2, j5, i2, i5
      integer ixm(0:n1+1,0:n2+1), ixp(0:n1+1,0:n2+1)
      real vol(0:n1+1,0:n2+1), ps_tmp(0:n1+1,0:n2+1), ps(0:n1+1,0:n2+1)
      real fsprd

c *** Local variables
      integer iy, ix, ix1, ix2, iy1, iy2
      integer omp_get_thread_num
      real fs0, signps
       if (omp_get_thread_num().ne.0) call xerrab('non-master omp thread not supposed to run here')
      fs0 = 1. - 4*fsprd    ! fraction to central cell

            do iy = j2, j5
             do ix = i2, i5
               ix1 = ixm(ix,iy)
               ix2 = ixp(ix,iy)
               iy1 = max(0,iy-1)
               iy2 = min(n2+1,iy+1)
               if (abs( ps(ix,iy)*ps(ix1,iy)*ps(ix2,iy)*
     .                  ps(ix,iy1)*ps(ix,iy2) ) .gt. 1.e-200) then
                 signps = 1.
cc                 if (ps(ix,iy) .lt. 0.) signps = -1.
cc                 ps_tmp(ix,iy) =
cc     .                      fs0*  log(abs(ps(ix ,iy )/vol(ix ,iy)))
cc     .                    + fsprd*log(abs(ps(ix1,iy )/vol(ix1,iy)))
cc     .                    + fsprd*log(abs(ps(ix2,iy )/vol(ix2,iy)))
cc     .                    + fsprd*log(abs(ps(ix ,iy1)/vol(ix,iy1)))
cc     .                    + fsprd*log(abs(ps(ix ,iy2)/vol(ix,iy2)))
                  ps_tmp(ix,iy) = fs0*ps(ix,iy) + fsprd*(ps(ix1,iy)+
     .                                ps(ix2,iy)+ps(ix,iy1)+ps(ix,iy2))
               endif
             enddo
           enddo
           do iy = j2, j5
             do ix = i2, i5
               ix1 = ixm(ix,iy)
               ix2 = ixp(ix,iy)
               iy1 = max(0,iy-1)
               iy2 = min(n2+1,iy+1)
               ps(ix,iy) = ps_tmp(ix,iy)
cc               if (abs( ps(ix,iy)*ps(ix1,iy)*ps(ix2,iy)*
cc     .                  ps(ix,iy1)*ps(ix,iy2) ) .gt. 1.e-200) then
cc               signps = 1.
cc                 if (ps(ix,iy) .lt. 0.) signps = -1.
cc                 ps(ix,iy) = sign( exp(ps_tmp(ix,iy)) * vol(ix,iy),
cc     .                             signps )
cc               endif
             enddo
           enddo
ccdcac$ompend parallel
      return
      end
!----------------------------------------------------------------------!

      subroutine volavenv(n1, n2, j2, j5, i2, i5, ixp, ixm, fsprd,
     .                    ps_tmp, ps)
c ... This subroutine does a volume integral, or average, of cell source-like
c ... quantities by including interpolation based on adjacent-cell quantities
c ... like volave, except here per m**3 quantities are averaged directly,
c ... so volumes are not involved (note absence of vol(ix,iy))

c *** Input:
c        n1 is nx for one of the two array dimensions
c        n2 is ny for the second  array dimension
c        j2 and j5 are lower and upper limits of iy loop
c        i2 and i5 are lower and upper limits of ix loop
c        ixp is the 2-D index array ixp1
c        ixm is the 2-D index array ixm1
c        fsprd is the fraction spread to each of 4 adjacent cells
c        pstmp is an (0:n1+1,0:n2+1) array used as a local work array
c        ps is an (0:n1+1,0:n2+1) array containing the input source-like array
c        pstmp is an (0:n1+1,0:n2+1) array used as a local work array
c        ps is an (0:n1+1,0:n2+1) array containing the input source-like array

c *** Output:
c        ps is an (n1,n2) array containing the result of volume-averaging ps

      implicit none !38

c *** Input and output variables
      integer n1, n2, j2, j5, i2, i5
      integer ixm(0:n1+1,0:n2+1), ixp(0:n1+1,0:n2+1)
      real ps_tmp(0:n1+1,0:n2+1), ps(0:n1+1,0:n2+1)
      real fsprd

c *** Local variables
      integer iy, ix, ix1, ix2, iy1, iy2
      real fs0, signps

ccdcac$ompparallel default(firstprivate)
      fs0 = 1. - 4*fsprd    ! fraction to central cell

            do iy = j2, j5
             do ix = i2, i5
               ix1 = ixm(ix,iy)
               ix2 = ixp(ix,iy)
               iy1 = max(0,iy-1)
               iy2 = min(n2+1,iy+1)
               if (abs( ps(ix,iy)*ps(ix1,iy)*ps(ix2,iy)*
     .                  ps(ix,iy1)*ps(ix,iy2) ) .gt. 1.e-200) then
                 signps = 1.
                 if (ps(ix,iy) .lt. 0.) signps = -1.
                 ps_tmp(ix,iy) =   fs0*  log(abs(ps(ix ,iy)))
     .                           + fsprd*log(abs(ps(ix1,iy)))
     .                           + fsprd*log(abs(ps(ix2,iy)))
     .                           + fsprd*log(abs(ps(ix ,iy1)))
     .                           + fsprd*log(abs(ps(ix ,iy2)))
               endif
             enddo
           enddo
           do iy = j2, j5
             do ix = i2, i5
               ix1 = ixm(ix,iy)
               ix2 = ixp(ix,iy)
               iy1 = max(0,iy-1)
               iy2 = min(n2+1,iy+1)
               if (abs( ps(ix,iy)*ps(ix1,iy)*ps(ix2,iy)*
     .                  ps(ix,iy1)*ps(ix,iy2) ) .gt. 1.e-200) then
                 signps = 1.
                 if (ps(ix,iy) .lt. 0.) signps = -1.
                 ps(ix,iy) = sign( exp(ps_tmp(ix,iy)), signps)
               endif
             enddo
           enddo
ccdcac$ompend parallel
      return
      end
c ======================================================================
c
      subroutine coneq

c ... This subroutine calculates the fluxes needed for the ion continuity
c ... equations

            use Dim       ! nx,ny,nhsp,nisp,ngsp
            use Share     ! nxpt,geometry,nxc,cutlo
            use Xpoint_indices       ! ixpt1,ixpt2,iysptrx
            use UEpar     ! methn,nlimix,nlimiy,nlimiy
            use Coefeq    ! cnfx,cnfy
            use Selec     ! i1,i2,i3,i4,i5,i6,i7,i8,j1,j2,j3,j4,j5,j6,j7,j8
                    ! xlinc,xrinc,yinc,ixm1,ixp1
            use Comgeo    ! vol, gx, gy, sx ,xy
            use Compla    ! mi, zi, ni, uu, up, v2, v2ce, vygtan, mg
            use Comflo    ! flnix,flniy
            use Indices_domain_dcl    ! ixmxbcl

      implicit none !37

*  -- local variables
      integer methnx,methny,ifld
      integer ix,iy,igsp,iv,iv1,iv2,iv3,ix2,ix3,ix4,ix5,ix6
      real t,t0,t1,t2,a
c ------------------
      do 104 ifld = 1, nfsp
*  ---------------------------------------------------------------------
*     compute flux, residual
*     The residual is: res := snic + sniv * ni - outflow(ni).
*  ---------------------------------------------------------------------

*  -- compute flnix, flox, conx --

         methnx = mod(methn, 10)
         methny = methn/10
         do 81 iy = j4, j8
            do 80 ix = i1, i5
              if (zi(ifld).eq.0. .and. ineudif <= 2) then
                 flnix(ix,iy,ifld) = fngx(ix,iy,1)
              else
               ix2 = ixp1(ix,iy)

               if (methnx .eq. 2) then   ! central differencing

                  t2 = ( lni(ix, iy,ifld) + lni(ix2,iy,ifld) ) / 2

               elseif (methnx .eq. 3) then   ! upwind differencing

                  if( uu(ix,iy,ifld) .ge. 0.) then
                     t2 = lni(ix,iy,ifld)
                  else
                     t2 = lni(ix2,iy,ifld)
                  endif

               else   ! interp. ave or harmonic ave depending on wind*grad

                  t0 = ( lni(ix, iy,ifld)*gx(ix, iy) +
     .                   lni(ix2,iy,ifld)*gx(ix2,iy) ) /
     .                                      ( gx(ix,iy)+gx(ix2,iy) )
                  t1 = ( gx(ix,iy)+gx(ix2,iy) ) * lni(ix,iy,ifld) *
     .                   lni(ix2,iy,ifld) / ( cutlo + lni(ix,iy,ifld)*
     .                   gx(ix2,iy) + lni(ix2,iy,ifld)*gx(ix,iy) )
                  if( uu(ix,iy,ifld)*(lni(ix,iy,ifld)-lni(ix2,iy,ifld))
     .                                                     .ge. 0.) then
                     t2 = t0
                  else
                     t2 = t1
                  endif

               endif

               flnix(ix,iy,ifld) = cnfx*uu(ix,iy,ifld) * sx(ix,iy) * t2
cc               if (uu(ix,iy,ifld)*(lni(ix,iy,ifld)-lni(ix2,iy,ifld))
cc     .                                                   .lt. 0.) then
                  flnix(ix,iy,ifld) = flnix(ix,iy,ifld)/( 1 +
     .                             (nlimix(ifld)/lni(ix2,iy,ifld))**2 +
     .                             (nlimix(ifld)/lni(ix ,iy,ifld))**2 )
               endif
cfw              if (zi(ifld).eq.0.) then
cfw                 flnix(ix,iy,ifld) = flnix(ix,iy,ifld) - fngxy(ix,iy,1)
cfw              end if
cc              endif
   80      continue
           if ((isudsym==1.or.geometry.eq.'dnXtarget')
     .                               .and. nxc > 1) flnix(nxc,iy,ifld)=0.
           if (islimon.ne.0 .and. iy.ge.iy_lims) flnix(ix_lim,iy,ifld)=0.
           if (nxpt==2.and.ixmxbcl==1) flnix(ixrb(1)+1,iy,ifld)=0.
   81    continue


*  -- compute flniy, floy, cony --

         do 83 iy = j1, j5
            do 82 ix = i4, i8
               if (zi(ifld).eq.0.) then
                  flniy(ix,iy,ifld) = fngy(ix,iy,1)
               else
                  if (methny .eq. 2) then   ! central differencing

                     t2 = ( niy0(ix,iy,ifld) + niy1(ix,iy,ifld) ) / 2

                  elseif (methny .eq. 3) then    ! upwind differencing

                     if( vy(ix,iy,ifld) .ge. 0.) then
                        t2 = niy0(ix,iy,ifld)
                     else
                        t2 = niy1(ix,iy,ifld)
                     endif

                  else    ! interp. ave or harmonic ave depending on wind*grad

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

                  flniy(ix,iy,ifld) = cnfy*vy(ix,iy,ifld)*sy(ix,iy)*t2
                  if (vy(ix,iy,ifld)*(lni(ix,iy,ifld)-lni(ix,iy+1,ifld))
     .                                                      .lt. 0.) then
                     flniy(ix,iy,ifld) = flniy(ix,iy,ifld)/( 1 +
     .                               (nlimiy(ifld)/lni(ix,iy+1,ifld))**2 +
     .                               (nlimiy(ifld)/lni(ix,iy  ,ifld))**2 )
!             To prevent bad behavior for the nonorthognal mesh at y-bndry
c                  if (iy .eq. 0) flniy(ix,iy,ifld) = cnfy*vy(ix,iy,ifld)*
c     .                                     sy(ix,iy)*lni(ix,1,ifld)
c                  if (iy .eq. ny) flniy(ix,iy,ifld) = cnfy*vy(ix,iy,ifld)*
c     .                                     sy(ix,iy)*lni(ix,ny,ifld)
                  endif
c...  Note: nonorthogonality comes in through calc. of vy
               endif
 82         continue
 83      continue

         do ix = i4, i8
            flniy(ix,ny+1,ifld) = 0.0e0
         enddo

 104  continue

      return
      end
c ***** End of subroutine coneq **********
c ======================================================================
c
      subroutine upvisneo

c ... This subroutine calculates the total up ion viscosity term with
c ... neoclassical effects

            use Dim       ! nx,ny,nhsp,nisp,ngsp,nusp
            use Xpoint_indices       ! ixlb,ixpt1,ixpt2,ixrb,iysptrx1,iysptrx2
                               ! iysptrx
            use Comgeo    ! vol, gx, gy, sx ,xy
            use RZ_grid_info   ! bpol,b12,b32,bsqr
            use Bfield    !rbfbt
            use Share     ! nxpt,geometry,nxc,cutlo
            use Selec     ! i1,i2,i3,i4,i5,i6,i7,i8,xlinc,xrinc,ixs1,
                    ! j1,j1p,j2,j2p,j3,j4,j5,j6,j5p,j6p,j7,j8,
            use Compla    ! ni,up,te,ti,ng,phi,v2cd,
            use Comflo    ! qipar
            use Conduc    ! visxneo,nuii,alfneo,visvol_v,visvol_q

      implicit none !35

*  -- local variables
      integer ifld
      real tempp,tempm,diffp,diffm
      integer ix,iy,igsp,iv,iv1,iv2,iv3,ix1,ix2,ix3,ix4,ix5,ix6
      real tv,t0,t1,t2,a
      do ifld = 1, nfsp
       if(zi(ifld) > 0.) then
        do iy = j2, j5
          do ix = i2, i5
            ix1 = ixp1(ix,iy)
            ix2 = ixm1(ix,iy)
c     First do the viscosity driven by particle flux
            tempp = b12(ix1,iy)*( up(ix1,iy,ifld) +
     .                (v2cd(ix1,iy,ifld)+v2ce(ix1,iy,ifld))*
     .                                       rbfbt(ix1,iy)/rrv(ix1,iy) )
            tempm = b12(ix,iy)*( up(ix,iy,ifld) +
     .                (v2cd(ix,iy,ifld)+v2ce(ix,iy,ifld))*
     .                                       rbfbt(ix,iy)/rrv(ix,iy) )
            if(ix < nx) then
              diffp = (visxneo(ix1,iy,ifld)/rr(ix1,iy))*(tempp-tempm)*
     .                                          gx(ix1,iy)/bsqr(ix1,iy)
            else
              diffp = (visxneo(ix1,iy,ifld)/rr(ix1,iy))*(tempp-tempm)*
     .                                       2.*gx(nx,iy)/bsqr(ix1,iy)
            endif
            tempp = tempm
            tempm = b12(ix2,iy)*( up(ix2,iy,ifld) +
     .               (v2cd(ix2,iy,ifld)+v2ce(ix2,iy,ifld))*
     .                                       rbfbt(ix2,iy)/rrv(ix2,iy) )
            diffm = (visxneo(ix,iy,ifld)/rr(ix,iy))*(tempp-tempm)*
     .                                             gx(ix,iy)/bsqr(ix,iy)
            visvol_v(ix,iy,ifld)=(4./3.)*rrv(ix,iy)*
     .                                 b32(ix,iy)*(diffp-diffm)*
     .                                   gxf(ix,iy)*volv(ix,iy)
c     Now do the viscosity driven by heat flux
            tempp = b12(ix1,iy)*( qipar(ix1,iy,ifld) +
     .                   q2cd(ix1,iy,ifld)*rbfbt(ix1,iy)/rrv(ix1,iy) )
            tempm = b12(ix,iy)*( qipar(ix,iy,ifld) +
     .                   q2cd(ix,iy,ifld)*rbfbt(ix,iy)/rrv(ix,iy) )
            if(ix < nx) then
              diffp = alfneo(ix1,iy,ifld)*rr(ix1,iy)*(tempp-tempm)*
     .                  gx(ix1,iy) / (nuii(ix1,iy,ifld)*bsqr(ix1,iy))
            else
              diffp = alfneo(ix1,iy,ifld)*rr(ix1,iy)*(tempp-tempm)*
     .                 2.*gx(nx,iy) / (nuii(ix1,iy,ifld)*bsqr(ix1,iy))
            endif
            tempp = tempm
            tempm = b12(ix2,iy)*( qipar(ix2,iy,ifld) +
     .                   q2cd(ix2,iy,ifld)*rbfbt(ix2,iy)/rrv(ix2,iy) )
            diffm = alfneo(ix,iy,ifld)*rr(ix,iy)*(tempp-tempm)*
     .                   gx(ix,iy) / (nuii(ix,iy,ifld)*bsqr(ix,iy))
            visvol_q(ix,iy,ifld) = rrv(ix,iy)*b32(ix,iy)*
     .                              (diffp-diffm)*gxf(ix,iy)*volv(ix,iy)
          enddo
        enddo
       endif  ! if-test on zi(ifld)
      enddo
      return
      end
c ***** End of subroutine upvisneo **********
c ======================================================================
c
      subroutine jvisneo

c ... This subroutine calculates the plasma current from
c ... neoclassical viscosity effects

            use Dim       ! nx,ny,nhsp,nisp,ngsp,nusp
            use Xpoint_indices       ! ixlb,ixpt1,ixpt2,ixrb,iysptrx1,iysptrx2
                               ! iysptrx
            use Comgeo    ! vol, gx, gy, sx ,xy
            use RZ_grid_info   ! bpol,b12,b32,bsqr,b12ctr
            use Bfield    !rbfbt
            use Share     ! nxpt,geometry,nxc,cutlo
            use Selec     ! i1,i2,i3,i4,i5,i6,i7,i8,xlinc,xrinc,ixs1,
                    ! j1,j1p,j2,j2p,j3,j4,j5,j6,j5p,j6p,j7,j8,
            use Compla    ! ni,up,te,ti,ng,phi,v2cd,
            use Comflo    ! qipar
            use Conduc    ! visxneo,nuii,alfneo,visvol_v,visvol_q

      implicit none !34

*  -- local variables
      integer ifld
      real tempp,tempm,temp0,diffp,diffm,gradx_vpiv,gradx_vpiq
      integer ix,iy,igsp,iv,iv1,iv2,iv3,ix1,ix2,ix3,ix4,ix5,ix6
      real tv,t0,t1,t2,a
      do ifld = 1, 1  ! limit to main species for now
       if(zi(ifld) > 0.) then !place-holder for when ifld > 1 used
        do iy = j2, j5
          do ix = i2, i5
            ix1 = ixp1(ix,iy)
            ix2 = ixm1(ix,iy)
c     First do the current driven by neoclassical particle flux
            tempp = up(ix1,iy,ifld) +
     .                        (v2cd(ix1,iy,ifld)+v2ce(ix1,iy,ifld))*
     .                                     rbfbt(ix1,iy)/rrv(ix1,iy)
            temp0 = up(ix,iy,ifld) +
     .                        (v2cd(ix2,iy,ifld)+v2ce(ix2,iy,ifld))*
     .                                     rbfbt(ix2,iy)/rrv(ix2,iy)
            tempm = up(ix2,iy,ifld) +
     .                        (v2cd(ix,iy,ifld)+v2ce(ix,iy,ifld))*
     .                                     rbfbt(ix,iy)/rrv(ix,iy)
            gradx_vpiv = (visxneo(ix,iy,ifld)*b12(ix,iy)*rrv(ix,iy)/3.)
     .                    *( (tempp+temp0)*b12ctr(ix1,iy)-
     .                       (temp0+tempm)*b12ctr(ix,iy) )*0.5*gxf(ix,iy)
            fq2pneo(ix,iy) = -gradx_vpiv*dbm2dy(ix1,iy)
            fqypneo(ix,iy) = 0.5*(rbfbt(ix,iy)+rbfbt(ix1,iy))*
     .                            gradx_vpiv*dbm2dx(ix1,iy)

c     Now do the current driven by neoclassical heat flux
            tempp = qipar(ix1,iy,ifld) +
     .                   q2cd(ix1,iy,ifld)*rbfbt(ix1,iy)/rrv(ix1,iy)
            temp0 = qipar(ix,iy,ifld) +
     .                   q2cd(ix,iy,ifld)*rbfbt(ix,iy)/rrv(ix,iy)
            tempm = qipar(ix2,iy,ifld) +
     .                   q2cd(ix2,iy,ifld)*rbfbt(ix2,iy)/rrv(ix2,iy)
            gradx_vpiq =( alfneo(ix,iy,ifld)*0.24*
     .                     b12(ix,iy)*rrv(ix,iy)/nuii(ix,iy,ifld) )*
     .                     ( (tempp+temp0)*b12ctr(ix1,iy)-
     .                       (temp0+tempm)*b12ctr(ix,iy) )*0.5*gxf(ix,iy)
            fq2qneo(ix,iy) = -gradx_vpiq*dbm2dy(ix,iy)
            fqyqneo(ix,iy) = 0.5*(rbfbt(ix,iy)+rbfbt(ix1,iy))*
     .                            gradx_vpiq*dbm2dx(ix1,iy)

          enddo
        enddo
       endif  ! if-test on zi(ifld)
      enddo

      return
      end
c ***** End of subroutine jvisneo **********
c
c-----------------------------------------------------------------------
c--  Used by ANL for PETSc development -------------------
c-------------------------------------------------------------
      subroutine jacwrite(n, jac, jacj, jaci)

c  This function serves to output the jacobian for viewing purposes
c  The output is the file jacwrite.txt

      integer n,j,k
      real jac(*)
      integer jacj(*), jaci(n+1)

      open(UNIT=88,FILE="jacwrite.txt",STATUS='REPLACE')
 77   format(/)

      write(88,*)"This is the jacobian after some scaling"
      do j=1,n
        do k=jaci(j),jaci(j+1)-1
          write(88,*)j,"  ",jacj(k),"  ",jac(k)
        end do
      end do

      close(88)
      write(*,*)"Jacobian written successfully to jacwrite.txt"
      end
c ***** End of subroutine jacwrite **********

      SUBROUTINE scopy(N,SX,INCX,SY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
        INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
        REAL SX(*),SY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
        INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
        INTRINSIC mod
*     ..
       IF (n.LE.0) RETURN
       IF (incx.EQ.1 .AND. incy.EQ.1) THEN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
              m = mod(n,7)
              IF (m.NE.0) THEN
                 DO i = 1,m
                    sy(i) = sx(i)
                 END DO
                 IF (n.LT.7) RETURN
              END IF
              mp1 = m + 1
              DO i = mp1,n,7
                 sy(i) = sx(i)
                 sy(i+1) = sx(i+1)
                 sy(i+2) = sx(i+2)
                 sy(i+3) = sx(i+3)
                 sy(i+4) = sx(i+4)
                 sy(i+5) = sx(i+5)
                 sy(i+6) = sx(i+6)
              END DO
           ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
          ix = 1
          iy = 1
          IF (incx.LT.0) ix = (-n+1)*incx + 1
          IF (incy.LT.0) iy = (-n+1)*incy + 1
          DO i = 1,n
             sy(iy) = sx(ix)
             ix = ix + incx
             iy = iy + incy
          END DO
       END IF
       RETURN
       END


       SUBROUTINE saxpy(N,SA,SX,INCX,SY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
       REAL SA
       INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
       REAL SX(*),SY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
       INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
       INTRINSIC mod
*     ..
       IF (n.LE.0) RETURN
       IF (sa.EQ.0.0) RETURN
       IF (incx.EQ.1 .AND. incy.EQ.1) THEN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
          m = mod(n,4)
          IF (m.NE.0) THEN
             DO i = 1,m
                sy(i) = sy(i) + sa*sx(i)
             END DO
          END IF
          IF (n.LT.4) RETURN
          mp1 = m + 1
          DO i = mp1,n,4
             sy(i) = sy(i) + sa*sx(i)
             sy(i+1) = sy(i+1) + sa*sx(i+1)
             sy(i+2) = sy(i+2) + sa*sx(i+2)
             sy(i+3) = sy(i+3) + sa*sx(i+3)
          END DO
       ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
          ix = 1
          iy = 1
          IF (incx.LT.0) ix = (-n+1)*incx + 1
          IF (incy.LT.0) iy = (-n+1)*incy + 1
          DO i = 1,n
           sy(iy) = sy(iy) + sa*sx(ix)
           ix = ix + incx
           iy = iy + incy
          END DO
       END IF
       RETURN
       END

             subroutine sgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info
      real a(lda,1)
c
c     sgefa factors a real matrix by gaussian elimination.
c
c     sgefa is usually called by sgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for sgeco) = (1 + 9/n)*(time for sgefa) .
c
c     on entry
c
c        a       real(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that sgesl or sgedi will divide by zero
c                     if called.  use  rcond  in sgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas saxpy,sscal,isamax
c
c     internal variables
c
      real t
      integer isamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = isamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0e0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0e0/a(k,k)
            call sscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call saxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0e0) info = n
      return
      end

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

      subroutine DebugHelper(FileName)
            use Bcond 
            use Cfric 
            use Coefeq 
            use Comflo 
            use Comgeo 
            use Compla 
            use Comtra 
            use Conduc 
            use Dim 
            use Gradients 
            use Imprad 
            use Locflux 
            use MCN_sources 
            use PNC_params 
            use Rhsides 
            use Save_terms 
            use Selec 
            use Time_dep_nwt 
            use Timing 
            use UEpar 
            use Wkspace 
      implicit none
      integer,parameter:: iunit=777
      character(len = *) ::  filename

c$ omp  parallel default(firstprivate)
      open (unit = iunit, file = trim(filename))
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
c$ omp end  parallel
      end subroutine DebugHelper

