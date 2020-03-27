Package bbb 
============
.. py:attribute:: neqmx
   
   number of math. eqns to be solved/integrated
   
   :Default: 
   :Dimension: None
   :Group: Math_problem_size
   :Type: integer
   :Unit: 
.. py:attribute:: numvar
   
   number of physical variables per cell
   
   :Default: 
   :Dimension: None
   :Group: Math_problem_size
   :Type: integer
   :Unit: 
.. py:attribute:: numvarbwpad
   
   add to numvar for bandwidth calc;safety param
   
   :Default: 
   :Dimension: None
   :Group: Math_problem_size
   :Type: integer
   :Unit: 
.. py:attribute:: csfacti
   
   Bohm speed = sqrt((te+csfacti*ti)/mi)
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: cslim
   
   frac of cs used for limiter Bohm sheath b.c.
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: dcslim
   
   reduce sonic flow at limiter by the factor
   cslim*[1-exp(-(iy-iy_lims+1)/dcslim)]
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: lnlam
   
   Coulomb log;shouldn't be constant
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: methe
   
   elec. eng. eqn: 22-cd, 33-uw, 44-hyb, 55-p-law
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: methu
   
   ion mom. eqn: 22-cd, 33-uw, 44-hyb, 55-p-law
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: methn
   
   ion cont. eqn: 22-harmonic average, 33-uw
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: methi
   
   ion eng. eqn: 22-cd, 33-uw, 44-hyb, 55-p-law
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: methg
   
   neut. gas eqn: 22-cd, 33-uw, 44-hyb, 55-p-law
   66 nonorth. log intrp, 77 nonorth. 1/ng intrp
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: methp
   
   potential eqn: 22-cd, 33-uw, 44-hyb, 55-p-law
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isgxvon
   
   =0 uses gx in fmix; =1 for harmonic ave of gxf
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: ishavisy
   
   =1 uses harmonic ave for conxi up
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: cfaccony
   
   scales conxi for up
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: isvylog
   
   =0 for vy~(1/n)dn/dy; =1 for vy~d(log(n))/dy
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: isintlog
   
   nonog logrithm interp for remaining terms
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: concap
   
   =1 fixes Te and Ti to afix for thermal cond.
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: convis
   
   =1 fixes Te to afix for ion viscosity
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: cniatol
   
   multiplier for atol for ni
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: cngatol
   
   multiplier for atol for ng
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: cupatol
   
   multiplier for atol for up
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: cteatol
   
   multiplier for atol for te
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: ctiatol
   
   multiplier for atol for ti
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: cphiatol
   
   multiplier for atol for phi
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: tolbf
   
   multiplier for atol&rtol for the boundary eqns
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: tadj
   
   reduces time step by 1/tadj if iopts=1
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: icnuiz
   
   =1 constant ioniz. freq., cnuiz; =2 freezes
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: icnucx
   
   =0, var nucx;=1 const. nucx=cnucx;
   =2, use sigcx, so nucx~(Tg)**.5
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: cnuiz
   
   constant ioniz. freq. for icnuiz=1
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 1/s
.. py:attribute:: cnucx
   
   constant charge exhange freq. for icnucx=1
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 1/s
.. py:attribute:: isrecmon
   
   flag to turn-on recombination (yes=1); use
   cfrecom to turn-off recomb after isrecmon was on
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: cfrecom
   
   scale factor multiplying recombination freq.
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: igas
   
   =1 invokes local rate eqn. for ng
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: ingb
   
   background gas source=nuiz*ngbackg*
   (.9+.1*(ngbackg/ng)**ingb)
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: inflbg
   
   expon to force flalfg large near ng~ngback
   ex:flalfgx,y*(1.+(cflgb*ngbackg/ng)**inflbg)
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: cflbg
   
   scaling fac for flalfgx,y using inflbg
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: inzb
   
   background impurity source=nuiz*nzbackg*
   (.9+.1*(nzbackg/nzi)**ingb)
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: tebg
   
   backgrd elec eng sor to limit te~tebg
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: eV
.. py:attribute:: tibg
   
   backgrd ion eng sor to limit te~tebg
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: eV
.. py:attribute:: iteb
   
   exponent of (tebg*ev/te)**iteb for bkg sor
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: temin
   
   min value of te allow; if less, reset to
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: eV
.. py:attribute:: temin2
   
   soft floor with te=sqrt[te**2+(temin2*ev)**2]
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: eV
.. py:attribute:: pwrbkg_c
   
   const background factor in pwrebkg express
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: W/m**3
.. py:attribute:: pwribkg_c
   
   const background factor in pwribkg express
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: W/m**3
.. py:attribute:: cfwjdotelim
   
   factor scaling reduction of wjdote if te<tebg
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: nlimgx
   
   factor to prevent ion density pump out in x
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: nlimgy
   
   factor to prevent ion density pump out in y
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: is1D_gbx
   
   =1 turns on 1-D gas-box model
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: xgbx
   
   poloidal location of 1-D gas box
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: m
.. py:attribute:: ixgb
   
   poloidal index of xgbx 1-D gas box (calc)
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: agdc
   
   exp. decay factor ng from gas-box edge
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: pcolwid
   
   width of plasma column for 1-D gas-box model
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: m
.. py:attribute:: eion
   
   energy that ionized ion is born with
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: eV
.. py:attribute:: ediss
   
   elec eng lost by mol. dissoc; should = 2*eion
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: eV
.. py:attribute:: ebind
   
   binding energy carried by hydrogen ion
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: eV
.. py:attribute:: tfcx
   
   NO LONGER USED; instead see tgas
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: eV
.. py:attribute:: tfcy
   
   NO LONGER USED; instead see tgas
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: eV
.. py:attribute:: afix
   
   Te,i for fixed cond.(concap), visc.(convis)
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: e
.. py:attribute:: coef
   
   factor for ion viscosity: was 1.92 ???
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: ce
   
   factor for electron thermal conductivity
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: ci
   
   factor for ion thermal conductivity
   The zeff dependence of ce has been explicitly added in zcoef,
   thus ce should always be left as 3.16 even if zeff is not 1,
   provided zeff is less than or equal to 4.
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: ncrhs
   
   
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: istep
   
   
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: iter
   
   
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: dp1
   
   
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: qfl
   
   
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: csh
   
   
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: qsh
   
   
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: mfl
   
   
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: msh
   
   
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: ro
   
   
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: cs
   
   
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: fxe
   
   
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: ctaue
   
   
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: fxi
   
   
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: ctaui
   
   
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: zcoef
   
   factor (calc) give zeff dependence of elec thermal c.
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: coef1
   
   factor (calc) for energy equipartion rate
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: cnurn
   
   scales nurlx rate for ion continuity eqn.
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: cnuru
   
   scales nurlx rate for ion mom. eqn.
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: cnure
   
   scales nurlx rate for elec. eng. eqn.
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: cnuri
   
   scales nurlx rate for ion eng. eqn.
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: cnurg
   
   scales nurlx rate for gas eqn.
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: cnurp
   
   scales nurlx rate for potential eqn.
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: nurlxn
   
   
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: nurlxu
   
   
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: nurlxe
   
   
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: nurlxi
   
   
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: nurlxg
   
   
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: nurlxp
   
   
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: rnewpot
   
   mixture of fqy=(1-rnewpot)*fqy_old+rnewpot*fqy_new
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: r0slab
   
   effect. major radius for isnewpot j_r calc in slab
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: m
.. py:attribute:: ishymol
   
   =1 turns on hydr. mol; requires nhgsp=2
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: te_s_dis
   
   Te shift of ioniz curve to approx dissociation curve
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: isfqpave
   
   =0 for lin interp for fqp terms; =1 for simple ave.
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isgpye
   
   change -vy*dP/dy eng. terms; =1 for old B2; =2 for Knoll
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: iigsp
   
   Ion index for neutrals when isupgon=1 (zi(iigsp)=0)
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: fnnuiz
   
   fraction of new nuiz used for Jacobian
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: itrap_negni
   
   flag to trap negative ni condition
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: itrap_negt
   
   flag to trap negative Te,i condition
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: itrap_negng
   
   flag to trap negative ng condition
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isybdryog
   
   =1 sets fx0, fmx stencil to orthog values at iy=0 & ny
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isybdrywd
   
   =1 vy diffusion-only for iy=0 & ny if matwalli,o=1
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isxmpog
   
   =1 sets fy0, fmy stencil to orthog values at ix=nxc-1
   and ix=nxc+1 for geometry='dnbot'
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: iexclnxc1
   
   if=0; include nxc+1 for fee,iytotc if geometry=dnbot;
   if=1; exclude nxc+1 for fee,iytotc
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: ineudif
   
   =1 gas sub. neudif uses ng, tg for gas vel & fngx->fnix
   =2 gas sub. neudifgp uses pg for gas vel & fngx->fnix
   =3 gas sub. neudifl use log_ng, tg for gas vel
   otherwise, old case has ug=ui (strong cx coupling)
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: thetar
   
   rotate (R,Z) coordinates by angle theta (degrees)
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: isbcwdt
   
   include dtreal in B.C. if isbcwdt=1
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: ishosor
   
   if=1, integrate hydr. sources over cell; full RHS only
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: iseesorave
   
   cell ave factor; 0 ctr only; 1 5 pt ave elec eng sors
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: ispsorave
   
   cell ave factor; 0 ctr only; 1 5 pt ave of psorg,psor,etc.
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: fsprd
   
   fraction of eng. sor. spread to each of 4 neighbors
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: issyvxpt0
   
   if=1, set syv=0 around x-point; ambig. rad. mom. flux
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isrrvave
   
   if=0, rrv from vertex B's; if=1, rrv=0.5*(rr_1+rr_2);
   if=2, average of cases 0 and 1
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: rr_fac
   
   scale factor to multiple rr and rrv
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: rrmin
   
   min rr used in calc of u_tor & fqy for potential calc.
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: isdtsfscal
   
   if=1, dt is included in sfscal Jac scaling factor
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: frfqpn
   
   frac. of new fqp at ix=0,nx using grad at ix=1,nx-1
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: cffqpsat
   
   factor by which fqp can exceed fqpsatlb,rb (sat. cur)
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: isplflxl
   
   =0, flalfe,i not active at ix=0 & nx;=1 active all ix
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isplflxlv
   
   =0, flalfv not active at ix=0 & nx;=1 active all ix
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isplflxlgx
   
   =0, flalfgx not active at ix=0 & nx;=1 active all ix
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isplflxlgxy
   
   =0, flalfgxy not active at ix=0 & nx;=1 active all ix
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: iswflxlgy
   
   =0, flalfgy not active at iy=0 & ny;=1 active all iy
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isplflxlvgx
   
   =0, flalfvgx not active at ix=0 & nx;=1 active all ix
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isplflxlvgxy
   
   =0, flalfvgxy not active at ix=0 & nx;=1 active all ix
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: iswflxlvgy
   
   =0, flalfvgy not active at iy=0 & ny;=1 active all iy
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isplflxltgx
   
   =0, flalfvgx not active at ix=0 & nx;=1 active all ix
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isplflxltgxy
   
   =0, flalfvgxy not active at ix=0 & nx;=1 active all ix
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: iswflxltgy
   
   =0, flalfvgy not active at iy=0 & ny;=1 active all iy
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: flalfipl
   
   ion therm flux lim factor on plates when isplflxl=0
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: flalfepl
   
   elec therm flux lim factor on plates when isplflxl=0
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: isfeexpl0
   
   if=1, feex cannot be out of inner/outer plates
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isfeixpl0
   
   if=1, feix cannot be out of inner/outer plates
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isofric
   
   If =1, use old (B2) interspecies up drag expression
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: del_te_ro
   
   te width in eV of tanh which turns off pwrze below 1 eV
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: iskaboom
   
   =1 turns on ijmgetmr 'k' or 'kaboom' stopping option
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isnglf
   
   =1 gives ng=nglfix at ix=0
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: nglfix
   
   value of ng at ix=0 if isnglf=1
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: isngrf
   
   =1 gives ng=nglfix at ix=nx+1
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: ngrfix
   
   value of ng at ix=nx+1 if isnglf=1
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: isup1up2
   
   =1 sets up2=rup21*up1
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: rup21
   
   rup21=up2/up1 if isup1up2
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: isteon
   
   user:turns on (=1) electron energy eqn.
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: istion
   
   user:turns on (=1) ion enegy eqn.
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isphion
   
   user:turns on (=1) potential eqn.
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isphiofft
   
   user:=1 leaves old cur. on & ex=-d(phis)/dx; must be used
   with isphion=0
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isnewpot
   
   user:turns on (=1) new potential; J_r from tor. mom. bal.
   =-2 sets phi constant on core boundary with
   total core current = icoreelec
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isugfm1side
   
   =0, use pol ave gas vels in par up eqn
   =1, use 1-sided vals for domain decomp
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isnupdot1sd
   
   =0, use 2-pt ndot for (n*up)_dot;
   =1, use 1-sided n_dot for (n*up)_dot
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isphicore0
   
   =1 sets phi=0 in core if isphion=1
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: is_z0_imp_const
   
   =0 use hydr Keilhacker;=1 z0_imp_const
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: z0_imp_const
   
   z0 in therm force if is_z0_imp_const=1
   
   :Default: 
   :Dimension: None
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: cnfx
   
   X-flux coef for conv. in n-eq.
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cnfy
   
   Y-flux coef for conv. in n-eq.
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cnsor
   
   Coef for particle src. in n-eq.
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfneut
   
   Coef for fluid neutrals contrib's to resid's
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfnidh
   
   Coef for neutral-ion drift heating
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfupcx
   
   Coef for nucx*(up_ion - up_gas) momentum coupling
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfticx
   
   Coef for nucx*(up_ion-up_gas)**2 heating in Ti Eq
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfupimpg
   
   Coef for impur up Cx/elast drag on up=0 imp gas
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cftiimpg
   
   Coef for Ti cooling CX/elast loss to cold imp gas
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cmneut
   
   Coef for Monte Carlo neutrals contrib's to resid's
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: chradi
   
   Coef for hyd. ioniz. rad. loss in elec. eng. eq.
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: chradr
   
   Coef for hyd. recomb. rad. loss in elec. eng. eq.
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: chioniz
   
   Coef for hydrogen ionization in elec. eng. eq.
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: ifxnsgi
   
   =1 sets ne for <sig*v>_i to cne_sgvi
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: integer
   :Unit: 
.. py:attribute:: cne_sgvi
   
   ne for <sig*v>_i if ifxnsgi=1
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 1/m**3
.. py:attribute:: ctsor
   
   Coef for eng. src. in Ti eq. 0.5*mi*up**2*psor
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: ceisor
   
   scale fac for ion energy source term (nu_i & eion)
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: ccoldsor
   
   scale fac for ion eng loss from cold cx
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: sxgsol
   
   stretches x-coord. for gas in sol & core regions
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: sxgpr
   
   stretches x-coord. for gas in private flux region
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: xstscal
   
   scale-length with stretch-coord decays from plates
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: m
.. py:attribute:: cngsor
   
   Coef for part. src. in ng-eq.
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: lgvmax
   
   max gas scale length for calc. viscous D_g
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cmfx
   
   X-flux coef for conv. in up-eq.
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cmfy
   
   Y-flux coef for conv. in up-eq.
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cpgx
   
   Coef for Grad(p) in up-eq.
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfvisx
   
   Coef. for x-visc. in ti-eq.
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfvisy
   
   Coef. for y-visc. in ti-eq.
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfanomvisxg
   
   Coef. for neut x-visc ~travis(2)
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfanomvisyg
   
   Coef. for neut y-visc ~travis(2)
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfvisxn
   
   Coef. for neutral x-visc. in up(,,iispg) eqn
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: isvisxn_old
   
   =1 uses sigcx,rrfac=1; =0 uses kelhihg, rrfac=rr**2
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: integer
   :Unit: 
.. py:attribute:: cfvxnrr
   
   =1 gives rr**2 in visx gas; =0 gives old 1 factor
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfvisyn
   
   Coef. for neutral y-visc. in up(,,iispg) eqn
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: vboost
   
   previously scaled eqp; no longer in use
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cvgp
   
   Coef for v.Grad(p) terms.
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfbgt
   
   Coef for the B x Grad(T) terms.
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfjhf
   
   Coef for convective cur (fqp) heat flow
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: jhswitch
   
   Coef for the Joule-heating terms
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: integer
   :Unit: 
.. py:attribute:: cf2ef
   
   Coef for ExB drift in 2-direction
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfyef
   
   Coef for ExB drift in y-direction
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cftef
   
   Coef for ExB drift in toroidal direction
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cf2bf
   
   Coef for Grad B drift in 2-direction
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfybf
   
   Coef for Grad B drift in y-direction
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfcbti
   
   Coef for adding fnixcb & fniycb to Ti eqn.
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfcurv
   
   Coef for curvature part of Grad_B drift
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfgradb
   
   Coef for p_perp part of Grad_B drift
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfq2bf
   
   Coef for Grad_B current in 2-direction
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfqybf
   
   Coef for Grad_B current in y-direction
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfqyn
   
   Coef for cx coll. rad current in y-direction
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfqym
   
   Coef for spatial inertial rad current in y-dir.
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfqydt
   
   Coef for time-dep inertial rad current in y-dir.
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cf2dd
   
   Coef for diamagnetic drift in 2-direction
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfydd
   
   Coef for diamagnetic drift in y-direction
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cftdd
   
   Coef for diamagnetic drift in toroidal direction
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfrd
   
   Coef for resistive cross-field drift
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfvisxneov
   
   Coef for v-driven parallel viscosity
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfvisxneoq
   
   Coef for q-driven parallel viscosity
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfvycr
   
   Coef for thermal force class. vel. vycr
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfvycf
   
   Coef for visc. force class. vel. vycf
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfvyavis
   
   Coef for vy from anom perp viscosity
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfjve
   
   Coef for J-contribution to ve.
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfjp2
   
   Coef for B x gradP terms in div(J) eqn
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfjpy
   
   Coef for B x gradP terms in div(J) eqn
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: isnfmiy
   
   diff fmiy for symmetry for vel. cells touching x-pt.
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: integer
   :Unit: 
.. py:attribute:: cfnfmiy
   
   Coef for new fmiy for vel. cells touching x-pt.
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cnimp
   
   Coef for impurity radiation loss
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: fac2sp
   
   factor to test 2-species model; for
   equal densities, set fac2sp=2
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cftnm
   
   Coef for neutral cx in toroidal mom. eq for fqy
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfupjr
   
   coef to include u_par in Jr calc.
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfcximp1
   
   coef multi. kcxrz for imp(+1)+D(0)->imp(0)+D(+1)
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfcximp2
   
   coef mult. kcxrz;imp(+p)+D(0)->imp(p-1)+D(+1),p>1
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfnetap
   
   coef mult. netap*fqp term in frice express.
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: fcdif
   
   coef mult all constant anomal diff coef
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfmsor
   
   coef mult msor and msorxr in up eqn.
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfloyi
   
   coef mult ion radial convective energy flow
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfloye
   
   coef mult elec radial convective energy flow
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfcvte
   
   coef mult elec poloidal convect(~5/2) energy flow
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfcvti
   
   coef mult ion & neut pol convect(~5/2) energy flow
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfcvtg
   
   coef mult gas pol convect(~5/2) energy flow
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfloxiplt
   
   coef mult neutral convect engy from plates
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfloygwall
   
   coef mult neutral convect engy from walls
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: exjbdry
   
   exponent pwr to limit fqp < fqpsat at plates
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfgpijr
   
   scalar factor for grad_Pi term in fqya
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: iszeffcon
   
   if =1, zeff=zeffcon
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: integer
   :Unit: 
.. py:attribute:: zeffcon
   
   value of zeff if iszeffcon=1
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: alftng
   
   neutral thermal force coeff; careful of sign
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfqya
   
   Coef for anomalous current in y-direction (new model)
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfqyao
   
   Coef for anomalous current in y-direction (old model)
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfqyae
   
   Coef for anomalous electron current in y-direction
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfqyai
   
   Coef for anomalous ion current in y-direction
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cftgcond
   
   Coef for gas thermal cond (usually molecules)
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cftgeqp
   
   Coef for gas thermal equipartion (usually molecules)
   
   :Default: 
   :Dimension: None
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: ibctepl
   
   Switch for ix=0 energy flux bc's
   =0, fixed te (see tepltl)
   =1, standard sheath transmission b.c.
   =2, zero poloidal gradients for te
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: ibctipl
   
   Same as ibctepl, with te --> ti
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: ibctepr
   
   Switch for ix=nx+1 energy flux bc's
   =0, fixed te (see tepltr)
   =1, standard sheath transmission b.c.
   =2, zero poloidal gradients for te
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: ibctipr
   
   Same as ibctepr, with te --> ti
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isphilbc
   
   Switch for ix=0 b.c. on phi
   =0, phi = phi0l + kappal * te
   =1, phi = phi0l
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isphirbc
   
   Switch for ix=nx+1 b.c. on phi
   =0, phi = phi0r + kappar * te
   =1, phi = phi0r
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: iphibcc
   
   core BC at iy=1 when isnewpot=1;iy=0 
   =1, d^2(ey)/dy^2=0
   =2, te=constant & ey(ixmp,0)=eycore
   =3, phi=constant & ey(ixmp,0)=eycore
   >3, dphi(ix,1)=dphi_iy1,isutcore ctrls ix=ixmp
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: iphibcwi
   
   =0, d(ey)/dy=0
   =1, phi(ix,0) = phintewi*te(ix,0)/ev
   =3, d(phi)/dy/phi = 1/lyphi(1)
   =4, phi(ix,0)=phiwi(ix) in PF region
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: iphibcwo
   
   =0, d(ey)/dy=0
   =1, phi(ix,ny+1) = phintewi*te(ix,ny+1)/ev
   =3, d(phi)/dy/phi = 1/lyphi(2)
   =4, phi(ix,ny+1)=phiwo(ix)
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: phintewi
   
   phi/te on inner wall if iphibcwi=1
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: phintewo
   
   phi/te on outer wall if iphibcwo=1
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: icoreelec
   
   electrical current from core
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: A
.. py:attribute:: eycore
   
   rad E-field core BC is iphibcc=3
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: V/m
.. py:attribute:: cfniybbo
   
   factor to includ. vycb in fniy,feiy at iy=0 only
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: cfniydbo
   
   factor to includ. vycp in fniy,feiy at iy=0 only
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: cfeeybbo
   
   factor to includ. vycb in feey at iy=0 only
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: cfeeydbo
   
   factor to includ. vycp in feey at iy=0 only
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: cfeexdbo
   
   factor includ v2cde & BxgradTe in BC at ix=0,nx
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: cfeixdbo
   
   factor includ v2cdi & BxgradTi in BC at ix=0,nx
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: cfqybbo
   
   factor to includ. fqyb in core current B.C. only
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: cfqydbo
   
   factor to includ. fqyd in core current B.C. only
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: nfqya0core
   
   num iy core cells beyond iy=0 where force fqya=0
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: nfqya0pf
   
   num. iy pf cells beyond iy=0 where force fqya=0
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: nfqya0ow
   
   num iy outer wall cell below iy=ny+1 with fqya=0
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: ixfixnc
   
   ix where ni=ncore if isnicore=2 begins
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: incixc
   
   ix range for ni=ncore from ixfixnc if isnicore=2
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: tcoree
   
   core elecron temp if iflcore=0
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: eV
.. py:attribute:: tcorei
   
   core ion temp if iflcore=0
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: eV
.. py:attribute:: tedge
   
   edge ion and elec. temp; used for te,iwalli,o
   arrays if last element zero (as in interpolation)
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: eV
.. py:attribute:: tepltl
   
   left plate Te B.C. if ibctepl=0
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: eV
.. py:attribute:: tipltl
   
   left plate Ti B.C. if ibctipl=0
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: eV
.. py:attribute:: tepltr
   
   right plate Te B.C. if ibcteplr=0
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: eV
.. py:attribute:: tipltr
   
   right plate Ti B.C. if ibctiplr=0
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: eV
.. py:attribute:: tbmin
   
   min. wall & pf temp for extrap. b.c.(isextrt..)
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: eV
.. py:attribute:: nbmin
   
   min. wall & pf den for extrap. b.c.(isextrn..)
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: m**-3
.. py:attribute:: ngbmin
   
   min. core gas den for extrap. b.c.(isextrngc)
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: m**-3
.. py:attribute:: istewc
   
   switch for outer-wall BC on Te
   =0, set zero energy flux
   =1, set fixed temp to tedge or tewallo
   =2, use extrapolation BC
   =3, set Te scale length to lyte
   =4, set feey = bceew*fniy*te
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: istiwc
   
   switch for outer-wall BC on Ti, see istewc detail
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: istepfc
   
   switch for priv.-flux BC on Te
   =0, set zero energy flux
   =1, set fixed temp to tedge or tewalli
   =2, use extrapolation BC
   =3, set Te scale length to lyte
   =4, set feey = bceew*fniy*te
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: istipfc
   
   switch for priv.-flux BC on Ti, see istewc detail
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isulytex
   
   if=0, lytex filled with lyte
   if=1, user values of lytex used
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isulytix
   
   if=0, lytix filled with lyti
   if=1, user values of lytex used
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isulyphix
   
   if=0, lyphix filled with lyphix
   if=1, user values of lyphix used
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isextrnp
   
   =1 sets extrap. b.c. at div. plate bound'y for ni
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isextrnpf
   
   =1 sets extrap. b.c. at p.f. bound'y for ni
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isextrtpf
   
   =1 sets extrap. b.c. at p.f. bound'y for Te & Ti
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isextrngc
   
   =1 sets extrap. b.c. on core bdry for ng
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isextrnw
   
   =1 sets extrap. b.c. at outer wall for ni
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isextrtw
   
   =1 sets extrap. b.c. at outer wall for Te & Ti
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: iflcore
   
   =0, core Te,i=tcoree,i; =1 core power=pcoree,i;
   =-1, core d(Te,i)/dy=0
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: pcoree
   
   electron power from core if iflcore=1
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: W
.. py:attribute:: pcorei
   
   ion power from core if iflcore=1
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: W
.. py:attribute:: ifluxni
   
   flag for setting iy=0,ny+1 dens flux to 0 (=1,yes)
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: ckinfl
   
   includes kinetic viscosity in energy bound. cond.
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: isbohmms
   
   =0 for single-species Bohm; =1 for multispecies B
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isulynix
   
   if=0, lynix filled with lyni
   if=1, user values of lynix used
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isulyupx
   
   if=0, lyupx filled with lyup
   if=1, user values of lynup used
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: nwsor
   
   number of sources on wall; must be < 10
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: sinphi
   
   sine of angle between side wall and flux surf.
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: teb
   
   left plate electron temp for isfixlb=1
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: eV
.. py:attribute:: tib
   
   left plate ion temp for isfixlb=1
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: eV
.. py:attribute:: yylb0
   
   radial shift in LHB profiles for isfixlb=1
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: ywnii
   
   inner Gaussian radial width of nib
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: ywnio
   
   outer Gaussian radial width of nib
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: ywupi
   
   inner Gaussian radial width of upb
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: ywupo
   
   outer Gaussian radial width of upb
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: ywtei
   
   inner Gaussian radial width of teb
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: ywteo
   
   outer Gaussian radial width of teb
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: ywtii
   
   inner Gaussian radial width of tib
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: ywtio
   
   outer Gaussian radial width of tib
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: engbsr
   
   energy factor for backscattered neutrals to Ti
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: epsbs
   
   small fac added (substracted) from Rbs (Rfc)
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: isguardc
   
   using guard cells? (=1 yes, =0 no)
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: rlimiter
   
   position of limiter at ix=0 for isfixlb=2
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: islimsor
   
   =1 extends sources into limiter region
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isutcore
   
   Used for ix=ixcore phi BC ONLY IF iphibcc > 3
   =0, tor mom=lzcore on core;
   =1, d<uz>/dy=0;
   >1, d^2(Ey)/dy^2=0 at outer midplane
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: islbcn
   
   b.c. for ni at limiter guard cells;
   =0,1 set ni in 2 cells
   =2 set ni in 1 cell, fnix at interface
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: islbcu
   
   b.c. for up at limiter guard cells;
   =0,1 set up in 3 cells
   =2 set up in 2 cells, fmix at interface
   =3,4,6 set fmix at interface
   =5 set fmix-fmixy at interface
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: islbce
   
   b.c. for te at limiter guard cells;
   =0,1 set te in 2 cells
   =2 set te in 1 cell, feex at interface
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: islbci
   
   b.c. for ti at limiter guard cells;
   =0,1 set ti in 2 cells
   =2 set ti in 1 cell, feix at interface
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: islbcg
   
   b.c. for ng at limiter guard cells;
   =0,1 set ng in 2 cells
   =2 set ng in 1 cell, fngx at interface
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: islbcp
   
   b.c. for phi at limiter guard cells;
   =0,1 set phi in 2 cells
   =2 set phi in 1 cell, fqx at interface
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: matt
   
   output flag from syld96 for sputt. target mat.
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: matp
   
   output flag from syld96 for sputt. plasma
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: cion
   
   input to syld96; atom num. of sputt. target
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: cizb
   
   input to syld96; max charge state of plasma
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: crmb
   
   input to syld96; mass of plasma ions
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: AMU
.. py:attribute:: eincid
   
   incident energy of ion or neut. for chem sputt
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: eV
.. py:attribute:: t_wall
   
   temp. of side wall; now use tvwallo,i
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: K
.. py:attribute:: t_plat
   
   temp. of divertor plate; now use tvplatlb,rb
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: K
.. py:attribute:: flux_in
   
   incident ion or neutral flux for chem sputt
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: ychem
   
   chem sputt. yield output from sputchem
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: fchemywi
   
   deprecated var; use fchemygwi; no harm if=1
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: fchemywo
   
   deprecated var; use fchemygwo; no harm if=1
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: isexunif
   
   =1 forces ex ~ uniform at div. plates
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: xcnearlb
   
   =TRUE if Jac'n 'box' overlaps a left boundary
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: xcnearrb
   
   =TRUE if Jac'n 'box' overlaps a right boundary
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: openbox
   
   =TRUE if Jac'n 'box' is wide open
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: kappa0
   
   modified sheath drop (allows j>jsat) for kappa > kappa0
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: kappamx
   
   maximum kappa value
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: cfueb
   
   scale factor for ueb in plate b.c.'s
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: ikapmod
   
   =1 for new kappa model; =0 for qpfac model
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: cfvytanbc
   
   factor for adding vytan to plate B.C.
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: cgpl
   
   scale fac atom eng plate loss; experim.
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: cgpld
   
   scale fac disso eng loss; experim.
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: cgengpl
   
   new scale fac atom eng plate loss; old cgpl
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: cgengw
   
   new scale fac atom eng wall loss
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: cgmompl
   
   scale fac atom par mom plate loss
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: vgmomp
   
   vel used in exp factor of atom mom loss
   
   :Default: 
   :Dimension: None
   :Group: Bcond
   :Type: double
   :Unit: m/s
.. py:attribute:: recycm
   
   momentum recycling/Rp for inertial gas;
   if recycm betwn -9.9 & -10.1 d(up)/dx=0
   if recycm < -10.1, therm mom flux used
   
   :Default: 
   :Dimension: None
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: recyce
   
   energy recycling/Rp for inertial gas
   
   :Default: 
   :Dimension: None
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: recycl
   
   recycling coef. at a limiter (ix_lim)
   
   :Default: 
   :Dimension: None
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: recycml
   
   momentum recycling/Rp for gas at limtr
   
   :Default: 
   :Dimension: None
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: albedo_by_user
   
   if=1, user fills albedoo,i & albdlb,rb
   
   :Default: 
   :Dimension: None
   :Group: Rccoef
   :Type: integer
   :Unit: 
.. py:attribute:: isrefluxclip
   
   =1 prohib outward gas for inward ion
   
   :Default: 
   :Dimension: None
   :Group: Rccoef
   :Type: integer
   :Unit: 
.. py:attribute:: gamsec
   
   secondary elec emiss coeff on plates
   
   :Default: 
   :Dimension: None
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: sputtr
   
   sputtering coef. at plates
   
   :Default: 
   :Dimension: None
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: ipsputt_s
   
   start dens-index phys sputt species
   
   :Default: 
   :Dimension: None
   :Group: Rccoef
   :Type: integer
   :Unit: 
.. py:attribute:: ipsputt_e
   
   end dens-index of phys sputt species
   
   :Default: 
   :Dimension: None
   :Group: Rccoef
   :Type: integer
   :Unit: 
.. py:attribute:: npltsor
   
   number sources on plates; must be <= 10
   
   :Default: 
   :Dimension: None
   :Group: Rccoef
   :Type: integer
   :Unit: 
.. py:attribute:: isextpltmod
   
   =1 use ext gas plate fluxes fngxextlb,rb
   and feixextlb,rb
   
   :Default: 
   :Dimension: None
   :Group: Rccoef
   :Type: integer
   :Unit: 
.. py:attribute:: isextwallmod
   
   =1 use ext gas wall fluxes fngyexti,o
   and feiyexti,o
   
   :Default: 
   :Dimension: None
   :Group: Rccoef
   :Type: integer
   :Unit: 
.. py:attribute:: isoutwall
   
   =1 call outwallflux to export wall fluxes
   
   :Default: 
   :Dimension: None
   :Group: Rccoef
   :Type: integer
   :Unit: 
.. py:attribute:: ifixsrc
   
   =1 turns on fixed Gaussian source
   
   :Default: 
   :Dimension: None
   :Group: Fixsrc
   :Type: integer
   :Unit: 
.. py:attribute:: ifixpsor
   
   =1 freezes part. source to initial val.
   
   :Default: 
   :Dimension: None
   :Group: Fixsrc
   :Type: integer
   :Unit: 
.. py:attribute:: xxsrc
   
   
   
   :Default: 
   :Dimension: None
   :Group: Fixsrc
   :Type: double
   :Unit: 
.. py:attribute:: yysrc
   
   
   
   :Default: 
   :Dimension: None
   :Group: Fixsrc
   :Type: double
   :Unit: 
.. py:attribute:: c1n
   
   
   
   :Default: 
   :Dimension: None
   :Group: Fixsrc
   :Type: double
   :Unit: 
.. py:attribute:: c1e
   
   
   
   :Default: 
   :Dimension: None
   :Group: Fixsrc
   :Type: double
   :Unit: 
.. py:attribute:: c1i
   
   
   
   :Default: 
   :Dimension: None
   :Group: Fixsrc
   :Type: double
   :Unit: 
.. py:attribute:: a1n
   
   
   
   :Default: 
   :Dimension: None
   :Group: Fixsrc
   :Type: double
   :Unit: 
.. py:attribute:: a1e
   
   
   
   :Default: 
   :Dimension: None
   :Group: Fixsrc
   :Type: double
   :Unit: 
.. py:attribute:: a1i
   
   
   
   :Default: 
   :Dimension: None
   :Group: Fixsrc
   :Type: double
   :Unit: 
.. py:attribute:: b1n
   
   
   
   :Default: 
   :Dimension: None
   :Group: Fixsrc
   :Type: double
   :Unit: 
.. py:attribute:: b1e
   
   
   
   :Default: 
   :Dimension: None
   :Group: Fixsrc
   :Type: double
   :Unit: 
.. py:attribute:: b1i
   
   
   
   :Default: 
   :Dimension: None
   :Group: Fixsrc
   :Type: double
   :Unit: 
.. py:attribute:: i1
   
   
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: i2
   
   
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: i2p
   
   used for 4th-order diffusion in x
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: i3
   
   
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: i4
   
   
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: i5
   
   
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: i5m
   
   same as i5, except restricted to ix<nx
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: i6
   
   
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: i7
   
   
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: i8
   
   
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: j1
   
   
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: j2
   
   
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: j3
   
   
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: j4
   
   
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: j5
   
   
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: j5m
   
   same as j5, except restricted to iy<ny
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: j6
   
   
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: j7
   
   
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: j8
   
   
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: j1p
   
   y-index lower range for potential eqn
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: j2p
   
   y-index lower range for potential eqn
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: j5p
   
   y-index upper range for potential eqn
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: j6p
   
   y-index upper range for potential eqn
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: ixs1
   
   
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: ixf6
   
   
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: iys1
   
   
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: iyf6
   
   
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: xlinc
   
   
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: xrinc
   
   
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: yinc
   
   
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: isjaccorall
   
   if=1 uses all ix cells for iy=0 Jac
   
   :Default: 
   :Dimension: None
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: ix
   
   
   
   :Default: 
   :Dimension: None
   :Group: Aux
   :Type: integer
   :Unit: 
.. py:attribute:: iy
   
   
   
   :Default: 
   :Dimension: None
   :Group: Aux
   :Type: integer
   :Unit: 
.. py:attribute:: igsp
   
   
   
   :Default: 
   :Dimension: None
   :Group: Aux
   :Type: integer
   :Unit: 
.. py:attribute:: iv
   
   
   
   :Default: 
   :Dimension: None
   :Group: Aux
   :Type: integer
   :Unit: 
.. py:attribute:: iv1
   
   
   
   :Default: 
   :Dimension: None
   :Group: Aux
   :Type: integer
   :Unit: 
.. py:attribute:: iv2
   
   
   
   :Default: 
   :Dimension: None
   :Group: Aux
   :Type: integer
   :Unit: 
.. py:attribute:: iv3
   
   
   
   :Default: 
   :Dimension: None
   :Group: Aux
   :Type: integer
   :Unit: 
.. py:attribute:: ix1
   
   
   
   :Default: 
   :Dimension: None
   :Group: Aux
   :Type: integer
   :Unit: 
.. py:attribute:: ix2
   
   
   
   :Default: 
   :Dimension: None
   :Group: Aux
   :Type: integer
   :Unit: 
.. py:attribute:: ix3
   
   
   
   :Default: 
   :Dimension: None
   :Group: Aux
   :Type: integer
   :Unit: 
.. py:attribute:: ix4
   
   
   
   :Default: 
   :Dimension: None
   :Group: Aux
   :Type: integer
   :Unit: 
.. py:attribute:: ix5
   
   
   
   :Default: 
   :Dimension: None
   :Group: Aux
   :Type: integer
   :Unit: 
.. py:attribute:: ix6
   
   
   
   :Default: 
   :Dimension: None
   :Group: Aux
   :Type: integer
   :Unit: 
.. py:attribute:: ixmp
   
   poloidal index of outer midplane; for yyc,f
   
   :Default: 
   :Dimension: None
   :Group: Aux
   :Type: integer
   :Unit: 
.. py:attribute:: tv
   
   
   
   :Default: 
   :Dimension: None
   :Group: Aux
   :Type: double
   :Unit: 
.. py:attribute:: t0
   
   
   
   :Default: 
   :Dimension: None
   :Group: Aux
   :Type: double
   :Unit: 
.. py:attribute:: t1
   
   
   
   :Default: 
   :Dimension: None
   :Group: Aux
   :Type: double
   :Unit: 
.. py:attribute:: t2
   
   
   
   :Default: 
   :Dimension: None
   :Group: Aux
   :Type: double
   :Unit: 
.. py:attribute:: a
   
   
   
   :Default: 
   :Dimension: None
   :Group: Aux
   :Type: double
   :Unit: 
.. py:attribute:: errmsgflag
   
   =0 turns off error messages, =1 turns them on
   
   :Default: 
   :Dimension: None
   :Group: Err_msg_out
   :Type: integer
   :Unit: 
.. py:attribute:: errunit
   
   output unit for error messages
   (nksol ignores errunit, sending everything to 6)
   
   :Default: 
   :Dimension: None
   :Group: Err_msg_out
   :Type: integer
   :Unit: 
.. py:attribute:: inopt
   
   resets iopts for solvers (vodpk, daspk)
   
   :Default: 
   :Dimension: None
   :Group: Opt_input
   :Type: integer
   :Unit: 
.. py:attribute:: runtim
   
   time of first output; total time=runtim*trange
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: double
   :Unit: 
.. py:attribute:: trange
   
   factor multiplying runtim to give total sim. time
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: double
   :Unit: 
.. py:attribute:: neq
   
   total number of equations over whole domain
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: integer
   :Unit: 
.. py:attribute:: jacflg
   
   flag for computing Jacobian in vodpk
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: integer
   :Unit: 
.. py:attribute:: jpre
   
   flag for using the preconditioning step in vodpk
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: integer
   :Unit: 
.. py:attribute:: itol
   
   
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: integer
   :Unit: 
.. py:attribute:: itask
   
   
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: integer
   :Unit: 
.. py:attribute:: istate
   
   
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: integer
   :Unit: 
.. py:attribute:: iopts
   
   internally set to inopt, an input variable
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: integer
   :Unit: 
.. py:attribute:: mf
   
   vodpk flag: mf=21, full user J; mf=22,full lsode J
   mf=24, banded user J; mf=25, banded lsode J
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: integer
   :Unit: 
.. py:attribute:: idid
   
   
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: integer
   :Unit: 
.. py:attribute:: ires
   
   
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: integer
   :Unit: 
.. py:attribute:: ts
   
   start time for ODE solvers
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: double
   :Unit: 
.. py:attribute:: tout
   
   output times for ODE solvers
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: double
   :Unit: 
.. py:attribute:: dtmax
   
   maximum allowed dt for daspk if info(7)=1
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: double
   :Unit: 
.. py:attribute:: dtinit
   
   starting dt for daspk if info(8)=1
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: double
   :Unit: 
.. py:attribute:: maxpoly
   
   maximum polynomial power used in daspk timestepping
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: integer
   :Unit: 
.. py:attribute:: eplidpk
   
   optional input for daspk when info(13) = 1
   tolerance for linear Krylov iteration.
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: double
   :Unit: 
.. py:attribute:: epnldpk
   
   optional input for daspk when info(13) = 1
   tolerance for Newton iteration convergence.
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: double
   :Unit: 
.. py:attribute:: srtolpk
   
   del=srtolpk*rtol for num. diff. (daspk,vodpk)
   Now set internally as srtolpk=del/rtolv
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: double
   :Unit: 
.. py:attribute:: efacn
   
   scaling factor for Newton error test in vodpk
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: double
   :Unit: 
.. py:attribute:: ftol
   
   stop tolerance of su*f for nksol (=epsmch**(1/3))
   ( maxnorm(f) .le. ftol to stop. )
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: double
   :Unit: 
.. py:attribute:: stptol
   
   stop tolerance of yl(k)-yl(k-1) in nksol
   ( maxnorm(yl(k)-yl(k-1)) .le. stptol to stop. )
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: double
   :Unit: 
.. py:attribute:: epscon1
   
   linear solve tolerance in nksol, epsfac =
   epscon1*min(epscon2,frnm)
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: double
   :Unit: 
.. py:attribute:: epscon2
   
   linear solve tolerance in nksol, epsfac =
   epscon1*min(epscon2,frnm)
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: double
   :Unit: 
.. py:attribute:: iterm
   
   output flag for nksol
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: integer
   :Unit: 
.. py:attribute:: mdif
   
   nksol flag for user-supplied j*v product (0=internal)
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: integer
   :Unit: 
.. py:attribute:: ipflag
   
   nksol flag to precondition (1=yes)
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: integer
   :Unit: 
.. py:attribute:: mfnksol
   
   nksol method flag; =1 means dogleg strategy,
   =2 means linesearch with Arnoldi method,
   =3 means linesearch with GMRES method.
   negative mfnksol ignores global constaints
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: integer
   :Unit: 
.. py:attribute:: iprint
   
   nksol optional statistics flag.
   =0 means no optional statistics are printed.
   =1 means iteration count, norm of F(u) and
   no. of F evaluations are printed.
   =2 means irpint=1 statistics are printed, and
   statistics regarding the convergence of the
   Krylov iteration, dogleg strategy. See
   nksol documentation for more details.
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: integer
   :Unit: 
.. py:attribute:: itermx
   
   maximum number of nonlinear iterations for nksol.
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: integer
   :Unit: 
.. py:attribute:: stepmx
   
   maximum length of a Newton step for nksol.
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: double
   :Unit: 
.. py:attribute:: del2nksol
   
   if nonzero, size of del**2 for diff. quot. Jac
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: double
   :Unit: 
.. py:attribute:: taunksol
   
   initial size of trust region for dogleg strategy
   (mfnksol = 1) in nksol.
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: double
   :Unit: 
.. py:attribute:: incpset
   
   maximum number of nonlinear iterations before
   the preconditioner is reevaluated within nksol.
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: integer
   :Unit: 
.. py:attribute:: ismmaxuc
   
   =1 for calc. mmaxu internally from nx and ny
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: integer
   :Unit: 
.. py:attribute:: mmaxu
   
   maximum Krylov subspace dimension.
   currently, only used in nksol			 # If ismmaxuc=1, calc. internally; ismmaxuc=0 use input
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: integer
   :Unit: 
.. py:attribute:: icntnunk
   
   nksol continuation call flag.
   =1 tells nksol not to call the preconditioner routine
   pset on the current call. In this case, nksol
   assumes that the preconditioner was evaulated
   on an earlier call, and is to be used for as
   many steps as it is successful on this call.
   =0 tells nksol that this is not a continuation call.
   The preconditioner routine pset is called to
   evaluate and factor the Jacobian matrix.
   
   :Default: 
   :Dimension: None
   :Group: Lsode
   :Type: integer
   :Unit: 
.. py:attribute:: nlocal
   
   number of equations on given processor
   
   :Default: 
   :Dimension: None
   :Group: Parallv
   :Type: integer
   :Unit: 
.. py:attribute:: neqg
   
   total number of equations over all processors
   
   :Default: 
   :Dimension: None
   :Group: Parallv
   :Type: integer
   :Unit: 
.. py:attribute:: nxg
   
   number of global poloidal mesh points = nxg+2
   
   :Default: 
   :Dimension: None
   :Group: Parallv
   :Type: integer
   :Unit: 
.. py:attribute:: nyg
   
   number of global radial mesh points = nyg+2g
   
   :Default: 
   :Dimension: None
   :Group: Parallv
   :Type: integer
   :Unit: 
.. py:attribute:: meth
   
   input for fpvmalloc; spec. method (lmm)
   
   :Default: 
   :Dimension: None
   :Group: Parallv
   :Type: integer
   :Unit: 
.. py:attribute:: itmeth
   
   input for fpvmalloc; spec. interation method (iter)
   
   :Default: 
   :Dimension: None
   :Group: Parallv
   :Type: integer
   :Unit: 
.. py:attribute:: iatol
   
   input for fpvmalloc; spec. error array type
   
   :Default: 
   :Dimension: None
   :Group: Parallv
   :Type: integer
   :Unit: 
.. py:attribute:: igs
   
   input for fcvspgrm2; Gram-Schmidt process
   
   :Default: 
   :Dimension: None
   :Group: Parallv
   :Type: integer
   :Unit: 
.. py:attribute:: maxkd
   
   maximum Krylov dimension for kinsol
   
   :Default: 
   :Dimension: None
   :Group: Parallv
   :Type: integer
   :Unit: 
.. py:attribute:: maxlrst
   
   for kinsol
   
   :Default: 
   :Dimension: None
   :Group: Parallv
   :Type: integer
   :Unit: 
.. py:attribute:: msbpre
   
   preconditioner flag for kinsol
   
   :Default: 
   :Dimension: None
   :Group: Parallv
   :Type: integer
   :Unit: 
.. py:attribute:: globalstrat
   
   global strategy flag for kinsol
   
   :Default: 
   :Dimension: None
   :Group: Parallv
   :Type: integer
   :Unit: 
.. py:attribute:: rtol_pv
   
   relative tol. for parallel pvode
   
   :Default: 
   :Dimension: None
   :Group: Parallv
   :Type: double
   :Unit: 
.. py:attribute:: atol_pv
   
   relative tol. for parallel pvode
   
   :Default: 
   :Dimension: None
   :Group: Parallv
   :Type: double
   :Unit: 
.. py:attribute:: delt_pv
   
   linear converg. error-test param. for pvode
   
   :Default: 
   :Dimension: None
   :Group: Parallv
   :Type: double
   :Unit: 
.. py:attribute:: icflag
   
   flag to use constraint that ni, etc. not < 0
   =1 turns on for nksol(with rlx) and vodpk(no rlx)
   =2 adds rlx constraint to vodpk
   
   :Default: 
   :Dimension: None
   :Group: Constraints
   :Type: integer
   :Unit: 
.. py:attribute:: rlx
   
   fractional change allowed per iteration
   
   :Default: 
   :Dimension: None
   :Group: Constraints
   :Type: double
   :Unit: 
.. py:attribute:: rlxv
   
   fractional change in up allowed for svrpkg=newton
   
   :Default: 
   :Dimension: None
   :Group: Constraints
   :Type: double
   :Unit: 
.. py:attribute:: adjf1
   
   if mfnksol=3 glob strat, frnm_new/adjf1>=fnrm_old
   
   :Default: 
   :Dimension: None
   :Group: Constraints
   :Type: double
   :Unit: 
.. py:attribute:: iscolnorm
   
   =0 for no implicit scaling (suscal=1)
   =1 for scaling by normalization constants
   =2 for scaling by max(abs(yl),floors)
   =3 combination of global scaling with nnorm,
   etc, followed by local scaling by each yl
   
   :Default: 
   :Dimension: None
   :Group: Ynorm
   :Type: integer
   :Unit: 
.. py:attribute:: var_scale_floor
   
   factor from norm_c to floor_c except for up
   factor multiplied by normalization constants to get floors for scaling
   
   :Default: 
   :Dimension: None
   :Group: Ynorm
   :Type: double
   :Unit: 
.. py:attribute:: vsf_up
   
   var_scale_floor factor for up eqns
   
   :Default: 
   :Dimension: None
   :Group: Ynorm
   :Type: double
   :Unit: 
.. py:attribute:: vsf_phi
   
   var_scale_floor factor for phi eqns
   
   :Default: 
   :Dimension: None
   :Group: Ynorm
   :Type: double
   :Unit: 
.. py:attribute:: temp0
   
   normalization temperature
   
   :Default: 
   :Dimension: None
   :Group: Ynorm
   :Type: double
   :Unit: eV
.. py:attribute:: isflxvar
   
   sets variables for ODE, Jacobian
   =1 for yl=n,nv,nT; =0 for yl=n,v,T
   =2 for yl=n,v,nT
   
   :Default: 
   :Dimension: None
   :Group: Ynorm
   :Type: integer
   :Unit: 
.. py:attribute:: isrscalf
   
   rescales ODE rhs if isflxvar.ne.1
   
   :Default: 
   :Dimension: None
   :Group: Ynorm
   :Type: integer
   :Unit: 
.. py:attribute:: dx0
   
   norm. grid spacing factor for phi eqn
   
   :Default: 
   :Dimension: None
   :Group: Ynorm
   :Type: double
   :Unit: m
.. py:attribute:: nnorm
   
   normalization density(calc)
   
   :Default: 
   :Dimension: None
   :Group: Ynorm
   :Type: double
   :Unit: m**-3
.. py:attribute:: ennorm
   
   normalization energy density(calc)
   
   :Default: 
   :Dimension: None
   :Group: Ynorm
   :Type: double
   :Unit: J/m**3
.. py:attribute:: sigbar0
   
   normalization parallel cond. (calc)
   
   :Default: 
   :Dimension: None
   :Group: Ynorm
   :Type: double
   :Unit: Mho/m
.. py:attribute:: vpnorm
   
   normalization ion paral.velocity(calc)
   
   :Default: 
   :Dimension: None
   :Group: Ynorm
   :Type: double
   :Unit: m/s
.. py:attribute:: rdoff
   
   ranf-induced roundoff error compared to unity
   
   :Default: 
   :Dimension: None
   :Group: Ynorm
   :Type: double
   :Unit: 
.. py:attribute:: isyloext
   
   =1 allows d(yl)/dt using ext. yloext
   
   :Default: 
   :Dimension: None
   :Group: Ynorm
   :Type: integer
   :Unit: 
.. py:attribute:: pi
   
   Pi
   
   :Default: 
   :Dimension: None
   :Group: Phyvar
   :Type: double
   :Unit: 
.. py:attribute:: me
   
   Electron mass
   
   :Default: 
   :Dimension: None
   :Group: Phyvar
   :Type: double
   :Unit: 
.. py:attribute:: mp
   
   Proton mass
   
   :Default: 
   :Dimension: None
   :Group: Phyvar
   :Type: double
   :Unit: 
.. py:attribute:: ev
   
   1 electron volt
   
   :Default: 
   :Dimension: None
   :Group: Phyvar
   :Type: double
   :Unit: 
.. py:attribute:: qe
   
   Elementary charge
   
   :Default: 
   :Dimension: None
   :Group: Phyvar
   :Type: double
   :Unit: 
.. py:attribute:: mu0
   
   Vac. magnetic perm.
   
   :Default: 
   :Dimension: None
   :Group: Phyvar
   :Type: double
   :Unit: 
.. py:attribute:: eps0
   
   Vac. dielectric perm.
   
   :Default: 
   :Dimension: None
   :Group: Phyvar
   :Type: double
   :Unit: 
.. py:attribute:: rt8opi
   
   sqrt(8/pi)
   
   :Default: 
   :Dimension: None
   :Group: Phyvar
   :Type: double
   :Unit: 
.. py:attribute:: cdifnit
   
   =1 for all turb., =0 all fixed, D coef
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: none
.. py:attribute:: inbtdif
   
   if isbohmcalc=3, D,chi ~1/Bt**inbtdif
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: inbpdif
   
   if isbohmcalc=3, D,chi ~1/Bp**inbpdif
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: integer
   :Unit: 
.. py:attribute:: ixbpmin
   
   isbohmcalc=3, min bpol(ixpt2-ixbpmin,
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: integer
   :Unit: 
.. py:attribute:: isbohmcalc
   
   if=1, calc Bohm diff if facb... > 0
   if=2, harmonic ave of Bohm, difni, etc.
   if=3, D=difniv*(B0/B)**inbdif, etc
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: integer
   :Unit: 
.. py:attribute:: facbni
   
   factor for Bohm density y-diff. coeff.
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: facbup
   
   factor for Bohm parll v y-diff. coeff.
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: facbni2
   
   factor for Bohm density 2-diff. coeff.
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: facbee
   
   factor for Bohm Te diff. coeff.
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: facbei
   
   factor for Bohm Ti diff. coeff.
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: difcng
   
   constant gas diff. coeff if isgasdc=1
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: m**2/s
.. py:attribute:: isgasdc
   
   switch to turn on constant gas dif coef
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: integer
   :Unit: 
.. py:attribute:: flalfe
   
   || heat flux limit factor for elec.
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalfi
   
   || heat flux limit factor for ions
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: lxtemax
   
   max pol. scale len of elec heat-flux lim
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: m
.. py:attribute:: lxtimax
   
   max pol. scale len of ion heat-flux lim
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: m
.. py:attribute:: lxtgmax
   
   max pol. scale len of gas heat-flux lim
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: m
.. py:attribute:: flalftf
   
   elec. thermal force flux-lim factor
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flgam
   
   exponent for ion flux-limit expression
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flgamv
   
   exponent for vel flux-limit expression
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flgamg
   
   exponent for gas dens flux-limit
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flgamvg
   
   exponent for gas visc flux-limit
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flgamtg
   
   exponent for gas temp flux-limit
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: fricflf
   
   flux-limiting factor for inputs to
   multispecies friction (and upi) calc
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: isflxlde
   
   =1,elec flux limit diff;=0, conv/diff
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: integer
   :Unit: 
.. py:attribute:: isflxldi
   
   =1,ion flux limit diff;=0, conv/diff
   =2, diff on individ hxcij
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: integer
   :Unit: 
.. py:attribute:: kxe
   
   pol elec heat conduc factor; 1.35->Balescu
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: alfkxi
   
   reduces ion thermal conduc, K_||, if
   |ti(ix+1)-ti(ix)|<alfkxi*ti(ix)
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: alfkxe
   
   reduces elec thermal conduc, K_||, if
   |te(ix+1)-te(ix)|<alfkxe*te(ix)
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: rkxecore
   
   pol elec heat diff. reduc fac in core
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: inkxc
   
   expon on yyf/yyf(0) fac for core kxe
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: integer
   :Unit: 
.. py:attribute:: kye
   
   radial electron heat diffusivity
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kye4order
   
   4th order Te radial diff. coef.
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kyet
   
   turb. radial elec. heat diff. multiplier
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: none
.. py:attribute:: ckyet
   
   =1 for all turb., =0 all fixed, chi_e
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: none
.. py:attribute:: kxi
   
   poloidal ion heat diff. multi. fac
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: kxicore
   
   poloidal ion heat diff. factor in core
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: kxn
   
   poloidal cx-neutral heat diff. factor
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: kyi
   
   radial ion heat diffusivity
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kyi4order
   
   4th order Ti radial diff. coef.
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kyit
   
   turb. radial ion heat diff. multiplier
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: none
.. py:attribute:: ckyit
   
   =1 for all turb., =0 all fixed, chi_i
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: none
.. py:attribute:: kyn
   
   radial cx-neutral heat diff. factor
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: feqp
   
   (Te-Ti) equipartition multiplier
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: alfeqp
   
   reduces equipart. term if te~ti
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalfgnx
   
   flux-limit on total fngx;for safety
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalfgny
   
   flux-limit on total fngy;for safety
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: gcfacgx
   
   mult tot conv gas x-flux at ix=0 & nx
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: gcfacgy
   
   mult tot conv gas y-flux at iy=0 & ny
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: gcfacgtx
   
   mult grad Ti conv gas x-flux ix=0 & nx
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: gcfacgty
   
   mult grad Ti conv gas y-flux iy=0 & ny
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: isdifxg_aug
   
   =1 enhances D_xgas with flx-lim factor
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: integer
   :Unit: 
.. py:attribute:: isdifyg_aug
   
   =1 enhances D_ygas with flx-lim factor
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: integer
   :Unit: 
.. py:attribute:: flalfv
   
   parallel velocity flux limit factor
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: isupdrag
   
   =1 adds nonunif B-field drag on v_||
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: integer
   :Unit: 
.. py:attribute:: con_leng
   
   connect length used for coll trans fac
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: m
.. py:attribute:: frac_pt
   
   invers aspect ratio for colliless drag
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalftgx
   
   poloidal atom temp. diff. flux limit
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalftgy
   
   radial atom temp. diff. flux limit
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalftmx
   
   poloidal mol temp. diff. flux limit
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalftmy
   
   radial mol temp. diff. flux limit
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalfvgx
   
   poloidal gas parall viscosity flux lim
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalfvgy
   
   radial gas parall viscosity flux limit
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalfvgxy
   
   FL for neutral fmixy
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalftxy
   
   FL for feixy (addition to hcy FL)
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: rnn2cx
   
   ratio of neut-neut coll. to cx coll.
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: rscat2cx
   
   fraction of cx counted as visx scatt
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: sigcx
   
   cx cross-sect if icnucx=2
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: m**2
.. py:attribute:: kelhihg
   
   elastic coll. coeff:hyd_ion+hyd_atom
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: m**3/s
.. py:attribute:: kelhghg
   
   elastic coll. coeff:hyd_atm+hyd_atom
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: m**3/s
.. py:attribute:: kelhmhg
   
   elastic coll. coeff:hyd_mol+hyd_atom
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: m**3/s
.. py:attribute:: kelhmhm
   
   elastic coll. coeff:hyd_mol+hyd_mol
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: m**3/s
.. py:attribute:: cfmassfac
   
   scales elas scat factor 16mi/(3mg+mi)
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: sigvi_floor
   
   minimum of ioniz. rates allowed(1e-18)
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: m**3/s
.. py:attribute:: fupe_cur
   
   =1 fixes cur err to upe for isimpon=6
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: cfnus_e
   
   factor mult nu_star_e for elec coll_fe
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: cfnus_i
   
   factor mult nu_star_i for ion coll_fi
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: tibsep
   
   Ion temp on sep for banana width in lconi
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: eV
.. py:attribute:: tebsep
   
   Elec temp on sep for banana width in lcone
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: eV
.. py:attribute:: cfelecbwd
   
   Factor for elec banana width in lcone; makes
   elec banana width not to small for mesh
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: fluxfacy
   
   multiples y-fluxes & v(dP/dy) for 1D sims
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: isdifbetap
   
   =1 turns on betap-dependent & difniv diffusion
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: integer
   :Unit: 
.. py:attribute:: iexpbp
   
   exponent for diff ~ betap**iexpbp
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: dfacbp
   
   diff. coeff dens *betap**iexpbp;dif_use,dif_use
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: m**2/s
.. py:attribute:: trfacbp
   
   diff. coeff up *betap**iexpbp;tray_use,trax_use
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kefacbp
   
   diff. coeff Te *betap**iexpbp;kye_use,kxe_use
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kifacbp
   
   diff. coeff Ti *betap**iexpbp;kyi_use,kxi_use
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: m**2/s
.. py:attribute:: cfmolcool
   
   scale factor for molec cooling if ishymol=1
   
   :Default: 
   :Dimension: None
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: isdifuseinterp
   
   =1 use dif_int in dif_use, etc, only
   
   :Default: 
   :Dimension: None
   :Group: Interprettrans
   :Type: integer
   :Unit: 
.. py:attribute:: isadjsolprof
   
   adjs SOL nis,tes,tis to smooth connect
   
   :Default: 
   :Dimension: None
   :Group: Interprettrans
   :Type: integer
   :Unit: 
.. py:attribute:: denrdrop
   
   rel drop in nis rad prof for adjsolprof
   
   :Default: 
   :Dimension: None
   :Group: Interprettrans
   :Type: double
   :Unit: 
.. py:attribute:: terdrop
   
   rel drop in tes rad prof for adjsolprof
   
   :Default: 
   :Dimension: None
   :Group: Interprettrans
   :Type: double
   :Unit: 
.. py:attribute:: tirdrop
   
   rel drop in tis rad prof for adjsolprof
   
   :Default: 
   :Dimension: None
   :Group: Interprettrans
   :Type: double
   :Unit: 
.. py:attribute:: maxchgdiff
   
   max change allowed in dif_int per cell
   
   :Default: 
   :Dimension: None
   :Group: Interprettrans
   :Type: double
   :Unit: 
.. py:attribute:: dif_use_max
   
   max used values of diff. coeff
   
   :Default: 
   :Dimension: None
   :Group: Interprettrans
   :Type: double
   :Unit: m**2/s
.. py:attribute:: dif_use_min
   
   max used values of diff. coeff
   
   :Default: 
   :Dimension: None
   :Group: Interprettrans
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kye_use_max
   
   max used values of kye coeff
   
   :Default: 
   :Dimension: None
   :Group: Interprettrans
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kye_use_min
   
   max used values of kye coeff
   
   :Default: 
   :Dimension: None
   :Group: Interprettrans
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kyi_use_max
   
   max used values of kyeicoeff
   
   :Default: 
   :Dimension: None
   :Group: Interprettrans
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kyi_use_min
   
   max used values of kye coeff
   
   :Default: 
   :Dimension: None
   :Group: Interprettrans
   :Type: double
   :Unit: m**2/s
.. py:attribute:: difni_sol
   
   added to SOL dif_use if interp mode
   
   :Default: 
   :Dimension: None
   :Group: Interprettrans
   :Type: double
   :Unit: m**2/s
.. py:attribute:: difni_pf
   
   added to PF dif_use if interp mode
   
   :Default: 
   :Dimension: None
   :Group: Interprettrans
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kye_sol
   
   added to SOL kye_use if interp mode
   
   :Default: 
   :Dimension: None
   :Group: Interprettrans
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kye_pf
   
   added to PF kye_use if interp mode
   
   :Default: 
   :Dimension: None
   :Group: Interprettrans
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kyi_sol
   
   added to SOL kyi_use if interp mode
   
   :Default: 
   :Dimension: None
   :Group: Interprettrans
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kyi_pf
   
   added to PF kye_use if interp mode
   #taudndt real [s] /1.e10/ # global density rise-time
   
   :Default: 
   :Dimension: None
   :Group: Interprettrans
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kappabar
   
   field-line avg'd curvature
   
   :Default: 
   :Dimension: None
   :Group: Turbulence
   :Type: double
   :Unit: 1/m
.. py:attribute:: lambdan
   
   dens(divertor) / dens(midplane)
   
   :Default: 
   :Dimension: None
   :Group: Turbulence
   :Type: double
   :Unit: none
.. py:attribute:: lambdat
   
   temp(midplane) / temp(divertor)
   
   :Default: 
   :Dimension: None
   :Group: Turbulence
   :Type: double
   :Unit: none
.. py:attribute:: gammasi
   
   secondary emis. coef. from ion bombard.
   
   :Default: 
   :Dimension: None
   :Group: Turbulence
   :Type: double
   :Unit: none
.. py:attribute:: lambdap
   
   e (dPhi0 / dr0) / (dTed / dr0)
   
   :Default: 
   :Dimension: None
   :Group: Turbulence
   :Type: double
   :Unit: none
.. py:attribute:: suppress_lmode
   
   =1 to suppress L-mode turbulence in SOL
   
   :Default: 
   :Dimension: None
   :Group: Turbulence
   :Type: integer
   :Unit: 
.. py:attribute:: maxmag_lmode
   
   max magn. of step in bracketing kymax
   
   :Default: 
   :Dimension: None
   :Group: Turbulence
   :Type: double
   :Unit: 
.. py:attribute:: nky
   
   number of ky values in search for kymax
   
   :Default: 
   :Dimension: None
   :Group: Turbulence
   :Type: integer
   :Unit: 
.. py:attribute:: kybeg
   
   lower limit of acceptable kymax
   
   :Default: 
   :Dimension: None
   :Group: Turbulence
   :Type: double
   :Unit: none
.. py:attribute:: kyend
   
   upper limit of acceptable kymax
   
   :Default: 
   :Dimension: None
   :Group: Turbulence
   :Type: double
   :Unit: none
.. py:attribute:: kya
   
   one initial point in search for kymax
   
   :Default: 
   :Dimension: None
   :Group: Turbulence
   :Type: double
   :Unit: none
.. py:attribute:: kyb
   
   other initial point in search for kymax
   
   :Default: 
   :Dimension: None
   :Group: Turbulence
   :Type: double
   :Unit: none
.. py:attribute:: iprint_lmode
   
   =1 for diagnostic output, =2 for more output
   
   :Default: 
   :Dimension: None
   :Group: Turbulence
   :Type: integer
   :Unit: 
.. py:attribute:: tol_lmode
   
   abs & rel tolerance in search for kymax
   
   :Default: 
   :Dimension: None
   :Group: Turbulence
   :Type: double
   :Unit: none
.. py:attribute:: isturbnloc
   
   =1 to turn on nonlocal dependence of D,chi
   
   :Default: 
   :Dimension: None
   :Group: Turbulence
   :Type: integer
   :Unit: 
.. py:attribute:: isturbcons
   
   =1 to make turbulent D,chi const within SOL
   =2 to apply radial digital filter to D,chi
   
   :Default: 
   :Dimension: None
   :Group: Turbulence
   :Type: integer
   :Unit: 
.. py:attribute:: diffusrange
   
   radial range of digital filter
   
   :Default: 
   :Dimension: None
   :Group: Turbulence
   :Type: double
   :Unit: m
.. py:attribute:: diffuslimit
   
   no. of surfaces in half-range of digital filter
   
   :Default: 
   :Dimension: None
   :Group: Turbulence
   :Type: integer
   :Unit: 
.. py:attribute:: islmodebeta
   
   =1 to turn on finite-beta correction
   
   :Default: 
   :Dimension: None
   :Group: Turbulence
   :Type: integer
   :Unit: 
.. py:attribute:: gradvconst
   
   factor involving rad. grad. of v(parallel)
   
   :Default: 
   :Dimension: None
   :Group: Turbulence
   :Type: double
   :Unit: none
.. py:attribute:: epsilon
   
   rhos / lte
   
   :Default: 
   :Dimension: None
   :Group: Turbulence_comm
   :Type: double
   :Unit: none
.. py:attribute:: turbdelta
   
   param. depending on temp. & length ratios
   
   :Default: 
   :Dimension: None
   :Group: Turbulence_comm
   :Type: double
   :Unit: none
.. py:attribute:: ssqthsqavg
   
   field-line avg of s**2 * theta**2
   
   :Default: 
   :Dimension: None
   :Group: Turbulence_comm
   :Type: double
   :Unit: none
.. py:attribute:: kxconst
   
   constant appearing in calc. of kx
   
   :Default: 
   :Dimension: None
   :Group: Turbulence_comm
   :Type: double
   :Unit: none
.. py:attribute:: cubrtnu
   
   cube root of parameter nu
   
   :Default: 
   :Dimension: None
   :Group: Turbulence_comm
   :Type: double
   :Unit: none
.. py:attribute:: bcoef0
   
   
   
   :Default: 
   :Dimension: None
   :Group: Turbulence_comm
   :Type: double complex
   :Unit: none
.. py:attribute:: ccoef1
   
   
   
   :Default: 
   :Dimension: None
   :Group: Turbulence_comm
   :Type: double complex
   :Unit: none
.. py:attribute:: ccoef2
   
   
   
   :Default: 
   :Dimension: None
   :Group: Turbulence_comm
   :Type: double
   :Unit: none
.. py:attribute:: ccoef3
   
   
   
   :Default: 
   :Dimension: None
   :Group: Turbulence_comm
   :Type: double
   :Unit: none
.. py:attribute:: nsteps
   
   number of logarithmically spaced output times
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: integer
   :Unit: 
.. py:attribute:: n_stor
   
   number of storage pts for solution
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: integer
   :Unit: 
.. py:attribute:: dt_init_rundt
   
   init dtreal for first call to rundt
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: double
   :Unit: 
.. py:attribute:: istep_nk
   
   array index for time-dep. nksol
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: integer
   :Unit: 
.. py:attribute:: nsteps_nk
   
   number of nksol time-steps
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: integer
   :Unit: 
.. py:attribute:: rdtphidtr
   
   ratio dtphi/dtreal
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: double
   :Unit: 
.. py:attribute:: ismfnkauto
   
   if =1, mfnksol=3 for dtreal<dtmfnk3, otherwise=-3
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: integer
   :Unit: 
.. py:attribute:: dtmfnk3
   
   dtreal for mfnksol sign change if ismfnkauto=1
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: double
   :Unit: 
.. py:attribute:: mult_dt
   
   factor expanding dtreal after ii2max steps
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: double
   :Unit: 
.. py:attribute:: ii1max
   
   number of changes to dtreal
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: integer
   :Unit: 
.. py:attribute:: ii2max
   
   number of timesteps at current dtreal
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: integer
   :Unit: 
.. py:attribute:: itermxrdc
   
   value of itermx used by rdcontdt
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: integer
   :Unit: 
.. py:attribute:: itermx_dt
   
   sets itermx for rundt; max iters per dt try
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: integer
   :Unit: 
.. py:attribute:: ftol_dt
   
   fnrm tolerance for the time-dependent steps
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: double
   :Unit: 
.. py:attribute:: ftol_min
   
   value of fnrm where time advance will stop
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: double
   :Unit: 
.. py:attribute:: dt_tot
   
   total time accumulated for run (output, not input)
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: double
   :Unit: 
.. py:attribute:: t_stop
   
   value of dt_tot (sec) where calculation will stop
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: double
   :Unit: 
.. py:attribute:: dt_max
   
   max time step for dtreal
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: double
   :Unit: 
.. py:attribute:: dt_kill
   
   min allowed time step; rdcontdt stops if reached
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: double
   :Unit: 
.. py:attribute:: deldt_min
   
   min relative change allowed for model_dt > 0
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: double
   :Unit: 
.. py:attribute:: initjac
   
   if=1, calc initial Jac upon reading rdcontdt
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: integer
   :Unit: 
.. py:attribute:: irev
   
   flag allows reduced dt advance after cutback
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: integer
   :Unit: 
.. py:attribute:: numrevjmax
   
   num dt reducts before Jac recalculated
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: integer
   :Unit: 
.. py:attribute:: numfwdjmax
   
   num dt increases before Jac recalculated
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: integer
   :Unit: 
.. py:attribute:: numrev
   
   count dt reducts in rdcontdt
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: integer
   :Unit: 
.. py:attribute:: numfwd
   
   count dt increases in rdcontdt
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: integer
   :Unit: 
.. py:attribute:: numrfcum
   
   number of cumulative reducts/increase in dt
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: integer
   :Unit: 
.. py:attribute:: tstor_s
   
   beginning time for storing solution
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: double
   :Unit: 
.. py:attribute:: tstor_e
   
   ending time for storing solution
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: double
   :Unit: 
.. py:attribute:: ipt
   
   index of variable printed at each out timestep
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: integer
   :Unit: 
.. py:attribute:: iprtrundt
   
   =1, then print rundt diag; .ne.1, no printing
   
   :Default: 
   :Dimension: None
   :Group: Timary
   :Type: integer
   :Unit: 
.. py:attribute:: rtauxfac
   
   fac*rtaux, Ly-a optic depth to plate
   =1 standard; <=0 skips rtau calc.
   
   :Default: 
   :Dimension: None
   :Group: Compla
   :Type: double
   :Unit: 
.. py:attribute:: rtauyfac
   
   fac*rtauy, Ly-a optic depth to wall
   
   :Default: 
   :Dimension: None
   :Group: Compla
   :Type: double
   :Unit: 
.. py:attribute:: rt_scal
   
   factor to scale rtaux,y & thus rtau
   
   :Default: 
   :Dimension: None
   :Group: Compla
   :Type: double
   :Unit: 
.. py:attribute:: fracvgpgp
   
   frac of vgp in vgradp eng terms
   
   :Default: 
   :Dimension: None
   :Group: Compla
   :Type: double
   :Unit: 
.. py:attribute:: ptjdote
   
   sum of J.E heating
   
   :Default: 
   :Dimension: None
   :Group: Postproc
   :Type: double
   :Unit: W
.. py:attribute:: ptigas
   
   sum of png2ni heating
   
   :Default: 
   :Dimension: None
   :Group: Postproc
   :Type: double
   :Unit: W/m^3
.. py:attribute:: pvmomcx
   
   sum of pmomv heating
   
   :Default: 
   :Dimension: None
   :Group: Postproc
   :Type: double
   :Unit: W/m^3
.. py:attribute:: iion_tot
   
   total ioniz. current over all isotopes
   
   :Default: 
   :Dimension: None
   :Group: Postproc
   :Type: double
   :Unit: A
.. py:attribute:: irecomb_tot
   
   total recomb. current over all isotopes
   
   :Default: 
   :Dimension: None
   :Group: Postproc
   :Type: double
   :Unit: A
.. py:attribute:: icxgas_tot
   
   total cx current over all gas isotopes
   
   :Default: 
   :Dimension: None
   :Group: Postproc
   :Type: double
   :Unit: A
.. py:attribute:: pradht
   
   total H photon rad. loss (- binding eng.)
   
   :Default: 
   :Dimension: None
   :Group: Postproc
   :Type: double
   :Unit: W
.. py:attribute:: pradiz
   
   ionization radiation energy loss
   
   :Default: 
   :Dimension: None
   :Group: Postproc
   :Type: double
   :Unit: W
.. py:attribute:: pradrc
   
   recombination radiation energy loss
   
   :Default: 
   :Dimension: None
   :Group: Postproc
   :Type: double
   :Unit: W
.. py:attribute:: pradfft
   
   net radiation loss via fixed-fraction impurity
   
   :Default: 
   :Dimension: None
   :Group: Postproc
   :Type: double
   :Unit: W
.. py:attribute:: pradzbind
   
   elec loss at imp ioniz carried as bind eng
   
   :Default: 
   :Dimension: None
   :Group: Postproc
   :Type: double
   :Unit: W
.. py:attribute:: pbinde
   
   pwr stored in ion binding pot. eng. from ioniz.
   
   :Default: 
   :Dimension: None
   :Group: Postproc
   :Type: double
   :Unit: W
.. py:attribute:: pbindrc
   
   binding eng. pwr released to elec. from recomb.
   
   :Default: 
   :Dimension: None
   :Group: Postproc
   :Type: double
   :Unit: W
.. py:attribute:: prdiss
   
   net photon pwr lost in dissoc.
   
   :Default: 
   :Dimension: None
   :Group: Postproc
   :Type: double
   :Unit: W
.. py:attribute:: pibirth
   
   net ion energy gain from dissoc.
   
   :Default: 
   :Dimension: None
   :Group: Postproc
   :Type: double
   :Unit: W
.. py:attribute:: jcvsor
   
   total core-region current for voljcsor
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: double
   :Unit: A
.. py:attribute:: ix_sjcsor
   
   if nonzero, beginning ix for voljcsor
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: integer
   :Unit: 
.. py:attribute:: ix_ejcsor
   
   if nonzero, ending ix for voljcsor
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: integer
   :Unit: 
.. py:attribute:: iy_sjcsor
   
   if nonzero, beginning iy for voljcsor
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: integer
   :Unit: 
.. py:attribute:: iy_ejcsor
   
   if nonzero, ending iy for voljcsor
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: integer
   :Unit: 
.. py:attribute:: pvole
   
   total power into electrons
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: double
   :Unit: W
.. py:attribute:: pvoli
   
   total power into ions
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: double
   :Unit: W
.. py:attribute:: z0pe
   
   axial or x loc. of elec. power profile
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: double
   :Unit: m
.. py:attribute:: z0pi
   
   axial or x loc. of ion power profile
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: double
   :Unit: m
.. py:attribute:: r0pe
   
   radial or y loc. of elec. power profile
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: double
   :Unit: m
.. py:attribute:: r0pi
   
   radial or y loc. of ion power profile
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: double
   :Unit: m
.. py:attribute:: zwpe
   
   axial or x Gauss. 1/2 width of e-power
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: double
   :Unit: m
.. py:attribute:: zwpi
   
   axial or y Gaussian 1/2 width ion power
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: double
   :Unit: m
.. py:attribute:: rwpe
   
   rad. or x Gaussian 1/2 width e-power
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: double
   :Unit: m
.. py:attribute:: rwpi
   
   rad. or y Gaussian 1/2 width ion power
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: double
   :Unit: m
.. py:attribute:: z0ni
   
   axial or x loc. of ion particle profile
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: double
   :Unit: m
.. py:attribute:: r0ni
   
   rad. or y loc. of ion particle profile
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: double
   :Unit: m
.. py:attribute:: zwni
   
   axial or y Gaussian 1/2 width ion prtcl
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: double
   :Unit: m
.. py:attribute:: rwni
   
   rad. or y Gaussian 1/2 width ion prtcl
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: double
   :Unit: m
.. py:attribute:: z0up
   
   axial or x loc. of ion mom. profile
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: double
   :Unit: m
.. py:attribute:: r0up
   
   rad. or y loc. of ion mom. profile
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: double
   :Unit: m
.. py:attribute:: zwup
   
   axial or y Gaussian 1/2 width ion mom.
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: double
   :Unit: m
.. py:attribute:: rwup
   
   rad. or y Gaussian 1/2 width ion mom.
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: double
   :Unit: m
.. py:attribute:: thetarot
   
   rotation angle for R,Z with effec. R,Z
   R_e= R0+(R-R0)cos(th)+(Z-Z0)sin(th),
   Z_e= Z0-(R-R0)sin(th)+(Z-Z0)cos(th),
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: double
   :Unit: rad
.. py:attribute:: rcutmin
   
   source zero if R<rcutmin
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: double
   :Unit: m
.. py:attribute:: zcutmin
   
   source zero if Z<zcutmin
   
   :Default: 
   :Dimension: None
   :Group: Volsrc
   :Type: double
   :Unit: m
.. py:attribute:: b0
   
   scale factor for magnetic field
   
   :Default: 
   :Dimension: None
   :Group: Bfield
   :Type: double
   :Unit: 
.. py:attribute:: isrozhfac
   
   =0 sets bfacx,yrozh=1; =1 computes ave
   
   :Default: 
   :Dimension: None
   :Group: Bfield
   :Type: integer
   :Unit: 
.. py:attribute:: cfparcur
   
   scale fac fqp=cfparcur*parcurrent if	
   isimpon=5 (fmombal from Hirshman)
   
   :Default: 
   :Dimension: None
   :Group: Comflo
   :Type: double
   :Unit: 
.. py:attribute:: lacor
   
   location in rwork where error vector begins
   
   :Default: 
   :Dimension: None
   :Group: Stat
   :Type: integer
   :Unit: 
.. py:attribute:: lewt
   
   location in rwork where ewt**-1 begins
   
   :Default: 
   :Dimension: None
   :Group: Stat
   :Type: integer
   :Unit: 
.. py:attribute:: iskaplex
   
   =1 if kappal is set externally (from parser)
   
   :Default: 
   :Dimension: None
   :Group: Poten
   :Type: integer
   :Unit: 
.. py:attribute:: iskaprex
   
   =1 if kappar is set externally (from parser)
   
   :Default: 
   :Dimension: None
   :Group: Poten
   :Type: integer
   :Unit: 
.. py:attribute:: bcee
   
   electron sheath energy trans. factor(newbc=0)
   
   :Default: 
   :Dimension: None
   :Group: Poten
   :Type: double
   :Unit: 
.. py:attribute:: bcei
   
   ion sheath energy trans. factor(newbc=0)
   
   :Default: 
   :Dimension: None
   :Group: Poten
   :Type: double
   :Unit: 
.. py:attribute:: bceew
   
   elec wall energy trans factor
   
   :Default: 
   :Dimension: None
   :Group: Poten
   :Type: double
   :Unit: 
.. py:attribute:: bceiw
   
   ion wall energy trans factor
   
   :Default: 
   :Dimension: None
   :Group: Poten
   :Type: double
   :Unit: 
.. py:attribute:: bcen
   
   neut energy trans. factor on plates 
   For combined neutral+ion energy equation
   
   :Default: 
   :Dimension: None
   :Group: Poten
   :Type: double
   :Unit: 
.. py:attribute:: bcenw
   
   neut eng trans fac on walls
   
   :Default: 
   :Dimension: None
   :Group: Poten
   :Type: double
   :Unit: 
.. py:attribute:: isfdiax
   
   switch to turn on diamagnetic drift for sheath
   potential calculation
   
   :Default: 
   :Dimension: None
   :Group: Poten
   :Type: double
   :Unit: 
.. py:attribute:: cthe
   
   electron thermal force coeff.
   
   :Default: 
   :Dimension: None
   :Group: Poten
   :Type: double
   :Unit: 
.. py:attribute:: cthi
   
   ?ion thermal force coeff.?; enters if zeff.ne.1
   
   :Default: 
   :Dimension: None
   :Group: Poten
   :Type: double
   :Unit: 
.. py:attribute:: sigma1
   
   parallel conductivity coeff.
   sigma1=1/(5.19e-5*Z*ln_lambda) & Z=1, ln_lambda=12.9
   
   :Default: 
   :Dimension: None
   :Group: Poten
   :Type: double
   :Unit: 1/(eV**1.5_Ohm_m)
.. py:attribute:: cfsigm
   
   scale factor for parallel cond. sigma1
   
   :Default: 
   :Dimension: None
   :Group: Poten
   :Type: double
   :Unit: 
.. py:attribute:: rsigpl
   
   ad hoc radial electrical conductivity - global
   
   :Default: 
   :Dimension: None
   :Group: Poten
   :Type: double
   :Unit: 
.. py:attribute:: rsigplcore
   
   ad hoc radial electrical conduct - core only
   ratio of perp to parallel conductivity
   
   :Default: 
   :Dimension: None
   :Group: Poten
   :Type: double
   :Unit: 
.. py:attribute:: cfkincor
   
   factor for kincorlb,rb denom. factor
   Variables for the grid-sequencing.
   yet to be defined?
   
   :Default: 
   :Dimension: None
   :Group: Poten
   :Type: double
   :Unit: 
.. py:attribute:: einduc
   
   inductive tor. E-field - input
   
   :Default: 
   :Dimension: None
   :Group: Gradients
   :Type: double
   :Unit: V/m
.. py:attribute:: cfgti
   
   scale factor for ion thermal force
   
   :Default: 
   :Dimension: None
   :Group: Cfric
   :Type: double
   :Unit: 
.. py:attribute:: cfgte
   
   scale factor for elec. thermal force
   
   :Default: 
   :Dimension: None
   :Group: Cfric
   :Type: double
   :Unit: 
.. py:attribute:: cftaud
   
   scale factor for ion-ion drag time
   
   :Default: 
   :Dimension: None
   :Group: Cfric
   :Type: double
   :Unit: 
.. py:attribute:: ngrid
   
   
   
   :Default: 
   :Dimension: None
   :Group: Grid
   :Type: integer
   :Unit: 
.. py:attribute:: ig
   
   counter for mesh-seq number
   
   :Default: 
   :Dimension: None
   :Group: Grid
   :Type: integer
   :Unit: 
.. py:attribute:: imeth
   
   imeth=inewton(igrid)
   
   :Default: 
   :Dimension: None
   :Group: Grid
   :Type: integer
   :Unit: 
.. py:attribute:: nurlx
   
   rate coeff. to relax to boundary cond.
   
   :Default: 
   :Dimension: None
   :Group: Grid
   :Type: double
   :Unit: 1/s
.. py:attribute:: ijactot
   
   tot Jac calcs, used as check when icntnunk=1
   
   :Default: 
   :Dimension: None
   :Group: Grid
   :Type: integer
   :Unit: 
.. py:attribute:: fnuizx
   
   fraction of nuiz in nuix (see nuix)
   
   :Default: 
   :Dimension: None
   :Group: Conduc
   :Type: double
   :Unit: 
.. py:attribute:: fnucxx
   
   fraction of nucx in nuix (see nuix)
   
   :Default: 
   :Dimension: None
   :Group: Conduc
   :Type: double
   :Unit: 
.. py:attribute:: cfvlh
   
   scal fac for hyd rate in nuvl
   
   :Default: 
   :Dimension: None
   :Group: Conduc
   :Type: double
   :Unit: 
.. py:attribute:: l_parloss
   
   parall length for nuvl loss rate
   
   :Default: 
   :Dimension: None
   :Group: Conduc
   :Type: double
   :Unit: m
.. py:attribute:: tdiflim
   
   lim on hcxe/ne; reduces hcxe if >0
   
   :Default: 
   :Dimension: None
   :Group: Conduc
   :Type: double
   :Unit: s
.. py:attribute:: lmfplim
   
   hcxe,i -> hcxe,i/(1+lmfp/lmfelim)
   
   :Default: 
   :Dimension: None
   :Group: Conduc
   :Type: double
   :Unit: m
.. py:attribute:: cfeta1
   
   scale factor for eta1
   
   :Default: 
   :Dimension: None
   :Group: Conduc
   :Type: double
   :Unit: 
.. py:attribute:: cfrtaue
   
   scale factor for cfrtaue
   
   :Default: 
   :Dimension: None
   :Group: Conduc
   :Type: double
   :Unit: 
.. py:attribute:: cfcl_e
   
   scale fac for dclass_e
   
   :Default: 
   :Dimension: None
   :Group: Conduc
   :Type: double
   :Unit: 
.. py:attribute:: cfcl_i
   
   scale fac for dclass_i
   
   :Default: 
   :Dimension: None
   :Group: Conduc
   :Type: double
   :Unit: 
.. py:attribute:: omgci_taui
   
   ion gy_freq*coll_rate for cl_model
   
   :Default: 
   :Dimension: None
   :Group: Conduc
   :Type: double
   :Unit: 
.. py:attribute:: omgce_taue
   
   elec gy_freq*coll_rate for cl_model
   
   :Default: 
   :Dimension: None
   :Group: Conduc
   :Type: double
   :Unit: 
.. py:attribute:: nuneo
   
   neoclass pol. damping rate; for fqyn
   
   :Default: 
   :Dimension: None
   :Group: Conduc
   :Type: double
   :Unit: 
.. py:attribute:: nstra
   
   number of 'strata' or 'source groups' in Monte-Carlo-Neutrals model;
   i.e., a surface or volume element where neutrals originate;
   for multi-species neutrals, each is a separate source group.
   
   :Default: 
   :Dimension: None
   :Group: MCN_dim
   :Type: integer
   :Unit: 
.. py:attribute:: nfl
   
   number of plasma fluids recognized by Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: None
   :Group: MCN_dim
   :Type: integer
   :Unit: 
.. py:attribute:: natmi
   
   number of atomic neutral species in EIRENE code
   
   :Default: 
   :Dimension: None
   :Group: MCN_dim
   :Type: integer
   :Unit: 
.. py:attribute:: nmoli
   
   number of molecular neutral species in EIRENE code
   
   :Default: 
   :Dimension: None
   :Group: MCN_dim
   :Type: integer
   :Unit: 
.. py:attribute:: nioni
   
   number of molecular ion species in EIRENE code
   
   :Default: 
   :Dimension: None
   :Group: MCN_dim
   :Type: integer
   :Unit: 
.. py:attribute:: nxf
   
   ix dimension from EIRENE file fort.44 or DEGAS2 file testdata.out
   
   :Default: 
   :Dimension: None
   :Group: MCN_dim
   :Type: integer
   :Unit: 
.. py:attribute:: nyf
   
   iy dimension from EIRENE file fort.44 or DEGAS2 file testdata.out
   
   :Default: 
   :Dimension: None
   :Group: MCN_dim
   :Type: integer
   :Unit: 
.. py:attribute:: nmcsp
   
   number of Monte Carlo species
   
   :Default: 
   :Dimension: None
   :Group: MCN_dim
   :Type: integer
   :Unit: 
.. py:attribute:: ismcnon
   
   flag for turning on plasma source terms from Monte-Carlo-Neutrals
   ismcnon=0 --> MCN plasma source terms are OFF (default)
   ismcnon=1 --> MCN-only is used for both Jac'n and RHS in pandf
   ismcnon=2 --> MCN-only is used for RHS, fluid-only is used for Jac'n
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: integer
   :Unit: 
.. py:attribute:: ismcnvar
   
   flag for scaling plasma source terms from Monte-Carlo-Neutrals
   ismcnvar=0 --> MCN plasma source terms are constant (default)
   ismcnvar=1 --> MCN plasma source terms scale with plate currents
   Special case of neutral atoms emitted with finite energy from walls
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: integer
   :Unit: 
.. py:attribute:: eedisspl
   
   energy loss for prompt dissociation at left plate
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: eV
.. py:attribute:: eedisspr
   
   energy loss for prompt dissociation at right plate
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: eV
.. py:attribute:: eidisspl
   
   energy gain for prompt dissociation at left plate
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: eV
.. py:attribute:: eidisspr
   
   energy gain for prompt dissociation at right plate
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: eV
.. py:attribute:: cmntgpl
   
   coeff. for neutral energy at left plate: cmntipl*ti
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: cmntgpr
   
   coeff. for neutral energy at right plate: cmntipr*ti
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: edisswo
   
   energy for prompt dissociation loss at outer wall
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: eV
.. py:attribute:: edisswi
   
   energy for prompt dissociation loss at private flux wall
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: eV
.. py:attribute:: cmntgwo
   
   coeff. for neutral energy at outer plate: cmntiwo*ti
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: cmntgwi
   
   coeff. for neutral energy at private flux plate: cmntiwi*ti
   cfneut 		real /1./ #Coef to turn on all fluid neutrals contrib's to resid's
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: cfneutsor_ni
   
   coeff. for fluid neutral particle source in resco
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: cfneutsor_mi
   
   coeff. for fluid neutral momentum source in resmo
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: cfneutsor_ei
   
   coeff. for fluid neutral energy source in resei
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: cfneutsor_ee
   
   coeff. for fluid neutral energy source in resee
   cmneut 		real /0./ #Coef to turn on all Monte Carlo neutral sources
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: cmneutsor_ni
   
   coeff. for MC neutral particle source in resco
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: cmneutsor_mi
   
   coeff. for MC neutral momentum source in resmo
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: cmneutsor_ei
   
   coeff. for MC neutral energy source in resei
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: cmneutsor_ee
   
   coeff. for MC neutral energy source in resee
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: cfneutdiv
   
   coeff. to turn on divergence of all fluid neutral fluxes
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: cfneutdiv_fng
   
   coeff. for div. fluid neutral particle flux in resng
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: cfneutdiv_fmg
   
   coeff. for div. fluid neutral momentum flux in resmo
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: cfneutdiv_feg
   
   coeff. for div. fluid neutral energy flux in resei
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: cmneutdiv
   
   coeff. to turn on divergence of all MC neutral fluxes
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: cmneutdiv_fng
   
   coeff. for div. fluid neutral particle flux in resng
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: cmneutdiv_fmg
   
   coeff. for div. fluid neutral momentum flux in resmo
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: cmneutdiv_feg
   
   coeff. for div. fluid neutral energy flux in resei
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: mcalpha_ng
   
   coeff. for blending kinetic and fluid ng
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: mcalpha_pg
   
   coeff. for blending kinetic and fluid pg
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: mcalpha_fng
   
   coeff. for blending kinetic and fluid fng
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: mcalpha_fmg
   
   coeff. for blending kinetic and fluid fmg
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: mcalpha_feg
   
   coeff. for blending kinetic and fluid feg
   ## Scalars ###
   
   :Default: 
   :Dimension: None
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: extneutopt
   
   specifies which external neutral program to use
   
   :Default: 
   :Dimension: None
   :Group: Ext_neutrals
   :Type: integer
   :Unit: 
.. py:attribute:: isextneuton
   
   whether to use external neutrals implicitly within exmain
   
   :Default: 
   :Dimension: None
   :Group: Ext_neutrals
   :Type: integer
   :Unit: 
.. py:attribute:: extneutmeth
   
   method for external neutrals: default=sources, 1=div. fluxes
   
   :Default: 
   :Dimension: None
   :Group: Ext_neutrals
   :Type: integer
   :Unit: 
.. py:attribute:: ext_verbose
   
   whether to print system call commands
   
   :Default: 
   :Dimension: None
   :Group: Ext_neutrals
   :Type: integer
   :Unit: 
.. py:attribute:: istimecmdon
   
   whether to time system call commands
   
   :Default: 
   :Dimension: None
   :Group: Ext_neutrals
   :Type: integer
   :Unit: 
.. py:attribute:: ismpicmdon
   
   whether to use MPI for external system call
   
   :Default: 
   :Dimension: None
   :Group: Ext_neutrals
   :Type: integer
   :Unit: 
.. py:attribute:: npext
   
   number of procs for external system call
   
   :Default: 
   :Dimension: None
   :Group: Ext_neutrals
   :Type: integer
   :Unit: 
.. py:attribute:: get_neutral_sources
   
   whether to use neutral source data
   
   :Default: 
   :Dimension: None
   :Group: Ext_neutrals
   :Type: integer
   :Unit: 
.. py:attribute:: get_neutral_moments
   
   whether to use neutral moment data
   
   :Default: 
   :Dimension: None
   :Group: Ext_neutrals
   :Type: integer
   :Unit: 
.. py:attribute:: pnc_verbose
   
   print diagnostic info
   
   :Default: 
   :Dimension: None
   :Group: PNC_params
   :Type: integer
   :Unit: 
.. py:attribute:: pnc_opt
   
   specifies choice of plasma-neutral coupling
   
   :Default: 
   :Dimension: None
   :Group: PNC_params
   :Type: integer
   :Unit: 
.. py:attribute:: pnc_step
   
   step count for plasma-neutral coupling
   
   :Default: 
   :Dimension: None
   :Group: PNC_params
   :Type: integer
   :Unit: 
.. py:attribute:: pnc_maxstep
   
   maximum number of coupled plasma+neutral steps
   
   :Default: 
   :Dimension: None
   :Group: PNC_params
   :Type: integer
   :Unit: 
.. py:attribute:: pnc_time
   
   time since beginning of coupled run
   pnc_ftol	real			/1e-4/					# ftol for PNC
   
   :Default: 
   :Dimension: None
   :Group: PNC_params
   :Type: double
   :Unit: 
.. py:attribute:: dtneut
   
   time step for neutrals
   
   :Default: 
   :Dimension: None
   :Group: PNC_params
   :Type: double
   :Unit: s
.. py:attribute:: dtplasma
   
   time step for plasma-neutral coupling
   
   :Default: 
   :Dimension: None
   :Group: PNC_params
   :Type: double
   :Unit: s
.. py:attribute:: dtold
   
   old time step
   
   :Default: 
   :Dimension: None
   :Group: PNC_params
   :Type: double
   :Unit: s
.. py:attribute:: relax_p
   
   relaxation parameter for plasma
   
   :Default: 
   :Dimension: None
   :Group: PNC_params
   :Type: double
   :Unit: 
.. py:attribute:: relax_g
   
   relaxation parameter for neutral gas
   
   :Default: 
   :Dimension: None
   :Group: PNC_params
   :Type: double
   :Unit: 
.. py:attribute:: pnc_ngs_mc
   
   replace fluid density with MC value
   
   :Default: 
   :Dimension: None
   :Group: PNC_params
   :Type: integer
   :Unit: 
.. py:attribute:: pnc_upgs_mc
   
   replace fluid parallel velocity with MC value
   
   :Default: 
   :Dimension: None
   :Group: PNC_params
   :Type: integer
   :Unit: 
.. py:attribute:: pnc_tgs_mc
   
   replace fluid temperature with MC value
   
   :Default: 
   :Dimension: None
   :Group: PNC_params
   :Type: integer
   :Unit: 
.. py:attribute:: pnc_fp
   
   file unit
   
   :Default: 
   :Dimension: None
   :Group: PNC_params
   :Type: integer
   :Unit: 
.. py:attribute:: pnc_print_norm
   
   choice of normalization: 0=absolute, 1=relative to max
   
   :Default: 
   :Dimension: None
   :Group: PNC_params
   :Type: integer
   :Unit: 
.. py:attribute:: pnc_nsave
   
   number of steps before saving data
   
   :Default: 
   :Dimension: None
   :Group: PNC_params
   :Type: integer
   :Unit: 
.. py:attribute:: pnc_dobalance
   
   run dobalance function MUST BE READ FIRST!!!
   
   :Default: 
   :Dimension: None
   :Group: PNC_params
   :Type: integer
   :Unit: 
.. py:attribute:: res_ni
   
   standard deviation of relative change in density
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: res_up
   
   standard deviation of relative change in parallel velocity
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: res_ti
   
   standard deviation of relative change in ion temperature
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: res_te
   
   standard deviation of relative change in electron temperature
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: res_phi
   
   standard deviation of relative change in electric potential
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: res_ng
   
   standard deviation of relative change in neutral density
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: res_upg
   
   standard deviation of relative change in neutral parallel velocity
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: res_tg
   
   standard deviation of relative change in neutral temperature
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: res_sng
   
   standard deviation of relative change in neutral density source
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: res_smg
   
   standard deviation of relative change in neutral parallel momentum source
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: res_seg
   
   standard deviation of relative change in neutral energy source
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: res_sni
   
   standard deviation of relative change in ion particle source
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: res_smor
   
   standard deviation of relative change in radial ion momentum source
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: res_smophi
   
   standard deviation of relative change in toroidal ion momentum source
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: res_smoz
   
   standard deviation of relative change in vertical ion momentum source
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: res_sei
   
   standard deviation of relative change in ion energy source
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: res_see
   
   standard deviation of relative change in electron energy source
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: del_ni
   
   maximum absolute change in density
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: del_up
   
   maximum absolute change in parallel velocity
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: del_ti
   
   maximum absolute change in ion temperature
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: del_te
   
   maximum absolute change in electron temperature
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: del_phi
   
   maximum absolute change in electric potential
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: del_ng
   
   maximum absolute change in neutral density
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: del_upg
   
   maximum absolute change in neutral parallel velocity
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: del_tg
   
   maximum absolute change in neutral temperature
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: del_sng
   
   maximum absolute change in neutral density source
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: del_smg
   
   maximum absolute change in neutral parallel momentum source
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: del_seg
   
   maximum absolute change in neutral energy source
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: del_sni
   
   maximum absolute change in ion particle source
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: del_smor
   
   maximum absolute change in radial ion momentum source
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: del_smophi
   
   maximum absolute change in toroidal ion momentum source
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: del_smoz
   
   maximum absolute change in vertical ion momentum source
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: del_sei
   
   maximum absolute change in ion energy source
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: del_see
   
   maximum absolute change in electron energy source
   
   :Default: 
   :Dimension: None
   :Group: PNC_data
   :Type: double
   :Unit: 
.. py:attribute:: nufak
   
   pseudo freq. on precond.-Jac diag for nksol
   
   :Default: 
   :Dimension: None
   :Group: Time_dep_nwt
   :Type: double
   :Unit: 1/s
.. py:attribute:: nufak0
   
   initial value of nufak0 saved (calc)
   
   :Default: 
   :Dimension: None
   :Group: Time_dep_nwt
   :Type: double
   :Unit: 1/s
.. py:attribute:: inufaknk
   
   flag for using nufak in Krylov step of nksol
   
   :Default: 
   :Dimension: None
   :Group: Time_dep_nwt
   :Type: integer
   :Unit: 
.. py:attribute:: dtreal
   
   real timestep (both Jac and RHS) for nksol
   Do not use large nufak and small dtreal simult
   
   :Default: 
   :Dimension: None
   :Group: Time_dep_nwt
   :Type: double
   :Unit: s
.. py:attribute:: dtdamp
   
   mix old/new as frac=1/(1+(dtdamp/dtreal)**itdamp)
   
   :Default: 
   :Dimension: None
   :Group: Time_dep_nwt
   :Type: double
   :Unit: s
.. py:attribute:: itdamp
   
   exponent for mix of old/new dt solutions
   
   :Default: 
   :Dimension: None
   :Group: Time_dep_nwt
   :Type: double
   :Unit: 
.. py:attribute:: dtreal_old
   
   previous value of dtreal
   
   :Default: 
   :Dimension: None
   :Group: Time_dep_nwt
   :Type: double
   :Unit: s
.. py:attribute:: dtphi
   
   additional dt to relax phi equation
   
   :Default: 
   :Dimension: None
   :Group: Time_dep_nwt
   :Type: double
   :Unit: s
.. py:attribute:: ydt_max
   
   maximum of yldot*sfscal
   
   :Default: 
   :Dimension: None
   :Group: Time_dep_nwt
   :Type: double
   :Unit: 
.. py:attribute:: ydt_max0
   
   old value of ydt_max
   
   :Default: 
   :Dimension: None
   :Group: Time_dep_nwt
   :Type: double
   :Unit: 
.. py:attribute:: alfnuf
   
   dtnewt->dtnewt*alfdtn*exp(ydt_max0/ydt_max)
   **expdtn
   
   :Default: 
   :Dimension: None
   :Group: Time_dep_nwt
   :Type: double
   :Unit: 
.. py:attribute:: expnuf
   
   see alfdtn
   
   :Default: 
   :Dimension: None
   :Group: Time_dep_nwt
   :Type: double
   :Unit: 
.. py:attribute:: deldt
   
   frac. of var. change per cell for var. dt
   
   :Default: 
   :Dimension: None
   :Group: Time_dep_nwt
   :Type: double
   :Unit: 
.. py:attribute:: model_dt
   
   determines dtuse for nksol time-step:
   =0, use dtreal
   =1, use dtreal*dtoptv/(dtreal+dtoptv)
   =2, use dtoptv
   =3, use sqrt(dtreal*dtoptv);
   =4, use dtreal*dtoptx/(dtreal+dtoptx)
   =5, use dtoptx
   =6, use sqrt(dtreal*dtoptx)
   
   :Default: 
   :Dimension: None
   :Group: Time_dep_nwt
   :Type: integer
   :Unit: 
.. py:attribute:: ubw
   
   
   
   :Default: 
   :Dimension: None
   :Group: Decomp
   :Type: integer
   :Unit: 
.. py:attribute:: lbw
   
   
   
   :Default: 
   :Dimension: None
   :Group: Decomp
   :Type: integer
   :Unit: 
.. py:attribute:: neqp1
   
   Dimension (=neq+1) of jaci
   
   :Default: 
   :Dimension: None
   :Group: Jacobian
   :Type: integer
   :Unit: 
.. py:attribute:: nnzmx
   
   Maximum no. of nonzeros in Jacobian matrix.
   
   :Default: 
   :Dimension: None
   :Group: Jacobian
   :Type: integer
   :Unit: 
.. py:attribute:: isjacstnlon
   
   Compute 9-pt stencil in ivl2gstnl - serial
   
   :Default: 
   :Dimension: None
   :Group: Jacobian
   :Type: integer
   :Unit: 
.. py:attribute:: nnz1mx
   
   Length of arrays in Jacobian_part
   
   :Default: 
   :Dimension: None
   :Group: Jacobian_part
   :Type: integer
   :Unit: 
.. py:attribute:: nlev
   
   Number of levels in levels array.
   See subroutine bfs for more details.
   
   :Default: 
   :Dimension: None
   :Group: Jacreorder
   :Type: integer
   :Unit: 
.. py:attribute:: maskval
   
   Scalar used with mask.
   
   :Default: 
   :Dimension: None
   :Group: Jacreorder
   :Type: integer
   :Unit: 
.. py:attribute:: ireorder
   
   Flag used to determine if a reordering
   of the Jacobian matrix is desired.
   = 1 means a reverse Cuthill-McKee
   reordering of the rows and columns
   of the Jacobian is done.
   = 0 means no reordering.
   
   :Default: 
   :Dimension: None
   :Group: Jacreorder
   :Type: integer
   :Unit: 
.. py:attribute:: lenpfac
   
   fudge factor to multiply neq by to get an
   estimate for the number of nonzeros in the
   preconditioner matrix.
   
   :Default: 
   :Dimension: None
   :Group: Preconditioning
   :Type: integer
   :Unit: 
.. py:attribute:: lenplufac
   
   fudge factor to multiply neq by to get an
   estimate for the number of nonzeros in the
   factored preconditioner matrix.
   
   :Default: 
   :Dimension: None
   :Group: Preconditioning
   :Type: integer
   :Unit: 
.. py:attribute:: lenplumx
   
   maximum number of nonzeros in the
   factored preconditioner matrix
   lenplumx = nnzmx + lenplufac*neq.
   
   :Default: 
   :Dimension: None
   :Group: Preconditioning
   :Type: integer
   :Unit: 
.. py:attribute:: tolilut
   
   threshold tolerance for ILUT.
   
   :Default: 
   :Dimension: None
   :Group: Ilutv
   :Type: double
   :Unit: 
.. py:attribute:: lfililut
   
   fill-in parameter used in ILUT. ILUT
   will allow up to lfililut additional nonzeros
   in each row of L and U.
   
   :Default: 
   :Dimension: None
   :Group: Ilutv
   :Type: integer
   :Unit: 
.. py:attribute:: ndiagmx
   
   maximum number of nonzero diagonals in the
   Jacobian matrix
   
   :Default: 
   :Dimension: None
   :Group: Nonzero_diagonals
   :Type: integer
   :Unit: 
.. py:attribute:: ndiag
   
   actual number of nonzero diagonals in the
   Jacobian matrix
   
   :Default: 
   :Dimension: None
   :Group: Nonzero_diagonals
   :Type: integer
   :Unit: 
.. py:attribute:: lfilinel
   
   fill-in parameter used in INEL preconditioner
   lfilinel= number of additional diagonals
   used in the INEL ILU preconditioner
   lfilinel+ndiag .le. ndiagmx.
   
   :Default: 
   :Dimension: None
   :Group: Nonzero_diagonals
   :Type: integer
   :Unit: 
.. py:attribute:: ndiagm
   
   number of nonzero diagonals stored in the
   INEL ILU preconditioner
   = min(lfilinel+ndiag,ndiagmx)
   
   :Default: 
   :Dimension: None
   :Group: Nonzero_diagonals
   :Type: integer
   :Unit: 
.. py:attribute:: newgeo
   
   flag to calculate new grid (1=yes)
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: mhdgeo
   
   flag for grid geometry
   mhdgeo = 2 ==> toroidal circular limiter
   mhdgeo = 1 ==> toroidal MHD equilibrium
   mhdgeo = 0 ==> cylindrical geometry
   mhdgeo = -1 ==> cartesian geometry
   mhdgeo = -2 ==> mag mirror (FRC-annulus)
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: gengrid
   
   flag to generate grid, else read from file gridue
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: isgindx
   
   =1 for interpolating grid based on indices
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: nfmax
   
   
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: restart
   
   flag for restart from previous case(yes=1)
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: initsol
   
   flag to initially solve algebraic eqns for
   DASPK (yes=1)
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: ttbeg
   
   initial Te in Joules = tinit/ev (calc)
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: double
   :Unit: 
.. py:attribute:: tinit
   
   initial electron temperature Te in eV
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: double
   :Unit: 
.. py:attribute:: tscal
   
   ratio of initial Ti & Tg to Te
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: double
   :Unit: 
.. py:attribute:: xgscal
   
   exponential scale of initial gas (m)
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: double
   :Unit: 
.. py:attribute:: isallloc
   
   =1 for local process. allocation with mpi
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: newaph
   
   =1 calls aphread for hyd. atomic data;=0 not
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: newapi
   
   =1, call readmc for new imp. data;=0, no
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: read_diffs
   
   =0,a flag to signal whether to read diffusivities
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: dif_io
   
   =0,a flag to signal whether to read/write dif_use
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: tra_io
   
   =0,a flag to signal whether to read/write tra_use
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: dutm_io
   
   =0,a flag to signal whether to read/write dutm_use
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: kye_io
   
   =0,a flag to signal whether to read/write kye_use
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: kyi_io
   
   =0,a flag to signal whether to read/write kyi_use
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: vy_io
   
   =0,a flag to signal whether to read/write kyi_use
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: vyup_io
   
   =0,a flag to signal whether to read/write kyi_use
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: vyte_io
   
   =0,a flag to signal whether to read/write kyi_use
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: vyti_io
   
   =0,a flag to signal whether to read/write kyi_use
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: fniyos_io
   
   =0,a flag to signal whether to read/write kyi_use
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: feeyosn_io
   
   =0,a flag to signal whether to read/write kyi_use
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: feiyosn_io
   
   =0,a flag to signal whether to read/write kyi_use
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: isvolsorext
   
   volsor sources if =0; or user sors if =1
   
   :Default: 
   :Dimension: None
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: isnintp
   
   switch to turn on new interpol. (=1)
   also check isgindx switch in UEint
   
   :Default: 
   :Dimension: None
   :Group: Interp
   :Type: integer
   :Unit: 
.. py:attribute:: isimesh
   
   flag for initial mesh => must copy
   save variables and not interpolate
   
   :Default: 
   :Dimension: None
   :Group: Interp
   :Type: integer
   :Unit: 
.. py:attribute:: isumesh2
   
   for parallel vers;=1, interp new mesh
   
   :Default: 
   :Dimension: None
   :Group: Interp
   :Type: integer
   :Unit: 
.. py:attribute:: nxold
   
   
   
   :Default: 
   :Dimension: None
   :Group: Interp
   :Type: integer
   :Unit: 
.. py:attribute:: nyold
   
   
   
   :Default: 
   :Dimension: None
   :Group: Interp
   :Type: integer
   :Unit: 
.. py:attribute:: nxoldg
   
   
   
   :Default: 
   :Dimension: None
   :Group: Interp
   :Type: integer
   :Unit: 
.. py:attribute:: nyoldg
   
   
   
   :Default: 
   :Dimension: None
   :Group: Interp
   :Type: integer
   :Unit: 
.. py:attribute:: iysptrxo
   
   prev. grid value for iysptrx
   
   :Default: 
   :Dimension: None
   :Group: Interp
   :Type: integer
   :Unit: 
.. py:attribute:: npes
   
   total number of processors
   
   :Default: 
   :Dimension: None
   :Group: Npes_mpi
   :Type: integer
   :Unit: 
.. py:attribute:: mype
   
   processor number of local processor (domain)
   
   :Default: 
   :Dimension: None
   :Group: Npes_mpi
   :Type: integer
   :Unit: 
.. py:attribute:: ismpion
   
   flag to indicate using MPI (if=1)
   
   :Default: 
   :Dimension: None
   :Group: Npes_mpi
   :Type: integer
   :Unit: 
.. py:attribute:: hascomm
   
   flag indicates communicator has been set (if=1)
   
   :Default: 
   :Dimension: None
   :Group: Npes_mpi
   :Type: integer
   :Unit: 
.. py:attribute:: isparmultdt
   
   =1 for multistep parallel beyond 1st step
   
   :Default: 
   :Dimension: None
   :Group: Npes_mpi
   :Type: integer
   :Unit: 
.. py:attribute:: isddcon
   
   switch to turn on domain decomposition
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: ndxt
   
   total number of x-domains
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: ndyt
   
   total number of y-domains
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: ndomain
   
   total number of domains
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: ndomain_orig
   
   tot num orig domains before par_data gather
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: nvrsend
   
   size of global real send/recv array for MPI
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: nvisend
   
   size of global integer send/recv array for MPI
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: neq_locgmx
   
   maximum of neq_locg
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: nx_loc
   
   number of ix cells for given processor
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: ny_loc
   
   number of iy cells for given processor
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: nvrsendl
   
   size of local real send/recv array for MPI
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: nvisendl
   
   size of local integer send/recv array for MPI
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: ixmnbcl
   
   B.C. type at ix=ixmin bdry;=0 intern,=1 extern
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: ixmxbcl
   
   B.C. type at ix=ixmax bdry;=0 intern,=1 extern
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: iymnbcl
   
   B.C. type at ix=iymin bdry;=0 intern,=1 extern
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: iymxbcl
   
   B.C. type at iy=iymax bdry;=0 intern,=1 extern
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: idxp1
   
   domain to the right of given domain (ix+1)
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: idxm1
   
   domain to the left of given domain (ix-1)
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: idyp1
   
   domain to the above given domain (iy+1)
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: idym1
   
   domain to the below given domain (iy-1)
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: neq_locl
   
   number of variables on local processor
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: numvarl
   
   =numvar global via MPI_BCAST for parallel
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: ispwrbcl
   
   =1 if domain has cell for core power BC
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: ixpt1l
   
   local ixpt1 before par_data gather
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: ixpt2l
   
   local ixpt2 before par_data gather
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: iysptrx1l
   
   local iysptrx1 before par_data gather
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: ixlbl
   
   local ixlb before par_data gather
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: ixrbl
   
   local ixrb before par_data gather
   
   :Default: 
   :Dimension: None
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: scrit
   
   
   
   :Default: 
   :Dimension: None
   :Group: Jacaux
   :Type: double
   :Unit: 
.. py:attribute:: normtype
   
   0,1,2 for max-norm, 1-norm, or 2-norm row scaling
   
   :Default: 
   :Dimension: None
   :Group: Jacaux
   :Type: integer
   :Unit: 
.. py:attribute:: issfon
   
   =1 calc sfscal for row scaling (norml.) by nksol
   
   :Default: 
   :Dimension: None
   :Group: Jacaux
   :Type: integer
   :Unit: 
.. py:attribute:: isrnorm
   
   =1 causes row normaliza. of Jac. (see normtype)
   
   :Default: 
   :Dimension: None
   :Group: Jacaux
   :Type: integer
   :Unit: 
.. py:attribute:: jscalcol
   
   =1 causes column scaling for daspk
   
   :Default: 
   :Dimension: None
   :Group: Jacaux
   :Type: integer
   :Unit: 
.. py:attribute:: del
   
   fractional change for finite diffs
   
   :Default: 
   :Dimension: None
   :Group: Variable_perturbation
   :Type: double
   :Unit: 
.. py:attribute:: delpy
   
   Forthon del; used to set del if > 0
   
   :Default: 
   :Dimension: None
   :Group: Variable_perturbation
   :Type: double
   :Unit: 
.. py:attribute:: dylconst
   
   factor in floor term in dyl
   
   :Default: 
   :Dimension: None
   :Group: Variable_perturbation
   :Type: double
   :Unit: 
.. py:attribute:: isjacreset
   
   if=1, pandf1 reset for last variable
   
   :Default: 
   :Dimension: None
   :Group: Variable_perturbation
   :Type: integer
   :Unit: 
.. py:attribute:: jaccliplim
   
   rel. value of elements to be retained
   
   :Default: 
   :Dimension: None
   :Group: Jacobian_clipping
   :Type: double
   :Unit: 
.. py:attribute:: istopjac
   
   flag to stop if non-zero elem at irstop,icstop
   
   :Default: 
   :Dimension: None
   :Group: Jacobian_clipping
   :Type: integer
   :Unit: 
.. py:attribute:: irstop
   
   row (or eqn) index of non-zero stopping test
   
   :Default: 
   :Dimension: None
   :Group: Jacobian_clipping
   :Type: integer
   :Unit: 
.. py:attribute:: icstop
   
   column (or var.) index of n-z stopping test
   
   :Default: 
   :Dimension: None
   :Group: Jacobian_clipping
   :Type: integer
   :Unit: 
.. py:attribute:: icsum
   
   
   
   :Default: 
   :Dimension: None
   :Group: Newtaux
   :Type: integer
   :Unit: 
.. py:attribute:: rwmin
   
   value of sumnew1 to stop Newton iter.
   
   :Default: 
   :Dimension: None
   :Group: Newtaux
   :Type: double
   :Unit: 
.. py:attribute:: saux
   
   
   
   :Default: 
   :Dimension: None
   :Group: Newtaux
   :Type: double
   :Unit: 
.. py:attribute:: saux1
   
   
   
   :Default: 
   :Dimension: None
   :Group: Newtaux
   :Type: double
   :Unit: 
.. py:attribute:: sumnew
   
   
   
   :Default: 
   :Dimension: None
   :Group: Newtaux
   :Type: double
   :Unit: 
.. py:attribute:: sumrdy
   
   
   
   :Default: 
   :Dimension: None
   :Group: Newtaux
   :Type: double
   :Unit: 
.. py:attribute:: nmaxnewt
   
   max number of Newton iterations
   
   :Default: 
   :Dimension: None
   :Group: Newtaux
   :Type: integer
   :Unit: 
.. py:attribute:: ifexmain
   
   scalar to indicate if subroutine allocate
   is called by exmain.
   =1 means allocate is called by exmain,
   =0 means it is not.
   
   :Default: 
   :Dimension: None
   :Group: Cdv
   :Type: integer
   :Unit: 
.. py:attribute:: iallcall
   
   flag to signal first call to allocate
   
   :Default: 
   :Dimension: None
   :Group: Cdv
   :Type: integer
   :Unit: 
.. py:attribute:: exmain_aborted
   
   Set to .true. in Python version on control-C abort
   
   :Default: 
   :Dimension: None
   :Group: Cdv
   :Type: integer
   :Unit: 
.. py:attribute:: isimpon
   
   switch for impurity model:
   0 for no impurities
   2 for fixed-fraction model
   3 for average-impurity-ion model(disabled)
   4 for INEL multi-charge-state model(disabled)
   5 for Hirshman's reduced-ion model
   6 for force-balance model or nusp_imp > 0;
   see also isofric for full-Z drag term
   7 for simultaneous fixed-fraction and
   multi-charge-state (isimpon=6) models
   
   :Default: 
   :Dimension: None
   :Group: Imprad
   :Type: integer
   :Unit: 
.. py:attribute:: nusp_imp
   
   fixes nusp for total num. of par. mom. eqns.
   
   :Default: 
   :Dimension: None
   :Group: Imprad
   :Type: integer
   :Unit: 
.. py:attribute:: isupimpap
   
   =1 includes imp atm phys in up eqn; =0, omits
   
   :Default: 
   :Dimension: None
   :Group: Imprad
   :Type: integer
   :Unit: 
.. py:attribute:: ismctab
   
   Determines which data is used for multi-charge-state rates.
   =1 tables originally generated by R. Campbell for D. Knoll,
   data file name is specified by inelmc=....
   corresponding rate evaluation routines are imprates and radimpmc.
   =2 tables generated by code from B. Braams,
   data file name is specified by mcfilename=...,
   corresponding rate evaluation routines are mcrates and radmc.
   
   :Default: 
   :Dimension: None
   :Group: Imprad
   :Type: integer
   :Unit: 
.. py:attribute:: misotope
   
   number of isotopes (including electrons)
   
   :Default: 
   :Dimension: None
   :Group: Reduced_ion_interface
   :Type: integer
   :Unit: 
.. py:attribute:: nchstate
   
   maximum charge state among all isotopes
   
   :Default: 
   :Dimension: None
   :Group: Reduced_ion_interface
   :Type: integer
   :Unit: 
.. py:attribute:: ib_idiv
   
   begin inner divertor
   
   :Default: 
   :Dimension: None
   :Group: Bdy_indexlims
   :Type: integer
   :Unit: 
.. py:attribute:: ie_idiv
   
   end inner divertor
   
   :Default: 
   :Dimension: None
   :Group: Bdy_indexlims
   :Type: integer
   :Unit: 
.. py:attribute:: ib_comwall
   
   begin common flux wall
   
   :Default: 
   :Dimension: None
   :Group: Bdy_indexlims
   :Type: integer
   :Unit: 
.. py:attribute:: ie_comwall
   
   end common flux wall
   
   :Default: 
   :Dimension: None
   :Group: Bdy_indexlims
   :Type: integer
   :Unit: 
.. py:attribute:: ib_odiv
   
   begin outer divertor
   
   :Default: 
   :Dimension: None
   :Group: Bdy_indexlims
   :Type: integer
   :Unit: 
.. py:attribute:: ie_odiv
   
   end outer divertor
   
   :Default: 
   :Dimension: None
   :Group: Bdy_indexlims
   :Type: integer
   :Unit: 
.. py:attribute:: ib_opfwall
   
   begin outer part of pf wall
   
   :Default: 
   :Dimension: None
   :Group: Bdy_indexlims
   :Type: integer
   :Unit: 
.. py:attribute:: ie_opfwall
   
   end outer part of pf wall
   
   :Default: 
   :Dimension: None
   :Group: Bdy_indexlims
   :Type: integer
   :Unit: 
.. py:attribute:: ib_ipfwall
   
   begin inner part of pf wall
   
   :Default: 
   :Dimension: None
   :Group: Bdy_indexlims
   :Type: integer
   :Unit: 
.. py:attribute:: ie_ipfwall
   
   end inner part of pf wall
   
   :Default: 
   :Dimension: None
   :Group: Bdy_indexlims
   :Type: integer
   :Unit: 
.. py:attribute:: liw
   
   length of iwork
   
   :Default: 
   :Dimension: None
   :Group: Solver_work_arrays
   :Type: integer
   :Unit: 
.. py:attribute:: lrw
   
   length of rwork
   
   :Default: 
   :Dimension: None
   :Group: Solver_work_arrays
   :Type: integer
   :Unit: 
.. py:attribute:: liwp
   
   length of iwwp
   
   :Default: 
   :Dimension: None
   :Group: Jac_work_arrays
   :Type: integer
   :Unit: 
.. py:attribute:: lwp
   
   length of wwp
   
   :Default: 
   :Dimension: None
   :Group: Jac_work_arrays
   :Type: integer
   :Unit: 
.. py:attribute:: xpzag
   
   Zagorski's xp
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: xszag
   
   Zagorski's xs
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: xkzag
   
   Zagorski's xk
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: ypzag
   
   Zagorski's yp
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: ykzag
   
   Zagorski's yk
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: inh
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: integer
   :Unit: 
.. py:attribute:: inz
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: integer
   :Unit: 
.. py:attribute:: ihf
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: integer
   :Unit: 
.. py:attribute:: istype
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: integer
   :Unit: 
.. py:attribute:: ibound
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: integer
   :Unit: 
.. py:attribute:: iboundz
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: integer
   :Unit: 
.. py:attribute:: zs
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: spuff
   
   PARAM common block
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: lst
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: integer
   :Unit: 
.. py:attribute:: wx
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: wy
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: ht
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: rmach
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: recyc
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: zrecyc
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: imap
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: integer
   :Unit: 
.. py:attribute:: snz
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: sn
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: sqe
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: sqi
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: fe0
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: fi0
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: flime
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: flimi
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: wkr
   
   OUT common block
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: xsi
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: zxsi
   
   
   
   :Default: 
   :Dimension: None
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: nxx
   
   number of mesh points
   
   :Default: 
   :Dimension: None
   :Group: Convdiffeqn
   :Type: integer
   :Unit: 
.. py:attribute:: alfz
   
   
   
   :Default: 
   :Dimension: None
   :Group: Convdiffeqn
   :Type: double
   :Unit: 
.. py:attribute:: vrfac
   
   scales convective velocity
   
   :Default: 
   :Dimension: None
   :Group: Convdiffeqn
   :Type: double
   :Unit: 
.. py:attribute:: sp
   
   
   
   :Default: 
   :Dimension: None
   :Group: Convdiffeqn
   :Type: double
   :Unit: 
.. py:attribute:: courant
   
   
   
   :Default: 
   :Dimension: None
   :Group: Convdiffeqn
   :Type: double
   :Unit: 
.. py:attribute:: tend
   
   final time
   
   :Default: 
   :Dimension: None
   :Group: Convdiffeqn
   :Type: double
   :Unit: 
.. py:attribute:: ndtmax
   
   max number of timesteps allowed
   
   :Default: 
   :Dimension: None
   :Group: Convdiffeqn
   :Type: integer
   :Unit: 
.. py:attribute:: ntim
   
   number of output times
   
   :Default: 
   :Dimension: None
   :Group: Convdiffeqn
   :Type: integer
   :Unit: 
.. py:attribute:: ito
   
   
   
   :Default: 
   :Dimension: None
   :Group: Convdiffeqn
   :Type: integer
   :Unit: 
.. py:attribute:: csfaclb
   
   frac of cs used for Bohm sheath b.c.
   
   :Default: 
   :Dimension: (31,2)
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: csfacrb
   
   frac of cs used for Bohm sheath b.c.
   
   :Default: 
   :Dimension: (31,2)
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: ngbackg
   
   background gas density
   
   :Default: 
   :Dimension: (6)
   :Group: UEpar
   :Type: double
   :Unit: 1/m**3
.. py:attribute:: facngbackg2ngs
   
   fraction of ngbackg add to initial ngs
   
   :Default: 
   :Dimension: (6)
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: nzbackg
   
   background impurity density
   
   :Default: 
   :Dimension: (31)
   :Group: UEpar
   :Type: double
   :Unit: 1/m**3
.. py:attribute:: facnzbackg2nis
   
   fraction of nzbackg add to initial nis
   
   :Default: 
   :Dimension: (31)
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: upclng
   
   max ion vel at beginning of iteration
   
   :Default: 
   :Dimension: (31)
   :Group: UEpar
   :Type: double
   :Unit: m/s
.. py:attribute:: facupclng2ups
   
   fraction of upclng subtract from initial ups
   
   :Default: 
   :Dimension: (31)
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: nlimix
   
   factor to prevent ion density pump out in x
   
   :Default: 
   :Dimension: (31)
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: nlimiy
   
   factor to prevent ion density pump out in y
   
   :Default: 
   :Dimension: (31)
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: label
   
   code name and run time, date, and machine
   
   :Default: 
   :Dimension: 
   :Group: UEpar
   :Type: character(72)
   :Unit: 
.. py:attribute:: svrpkg
   
   use solver pkg daspk,vodpk,nksol,newton
   reset to newton if inewton=1
   
   :Default: 
   :Dimension: 
   :Group: UEpar
   :Type: character(8)
   :Unit: 
.. py:attribute:: petscoptfile
   
   specify options file for petsc code
   
   :Default: 
   :Dimension: 
   :Group: UEpar
   :Type: character(80)
   :Unit: 
.. py:attribute:: petscopts
   
   specify options for petsc code
   
   :Default: 
   :Dimension: 
   :Group: UEpar
   :Type: character(2048)
   :Unit: 
.. py:attribute:: isnion
   
   user:turns on (=1) ion continuity eqn.
   
   :Default: 
   :Dimension: (31)
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isupon
   
   user:turns on (=1) parallel vel. eqn.
   
   :Default: 
   :Dimension: (31)
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isupgon
   
   user:=1 for par neutral vel. eqn.; index igsp
   
   :Default: 
   :Dimension: (6)
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isngon
   
   user:turns on (=1) neutral eqn.; index igsp
   
   :Default: 
   :Dimension: (6)
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: istgon
   
   user:turns on (=1) gas enegy eqn.
   
   :Default: 
   :Dimension: (6)
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isnionxy
   
   calc:=1 for ni eqn on; =0 for eqn off
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nisp)
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isuponxy
   
   calc:=1 for up eqn on; =0 for eqn off
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nisp)
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isteonxy
   
   calc:=1 for te eqn on; =0 for eqn off
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: istionxy
   
   calc:=1 for ti eqn on; =0 for eqn off
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isngonxy
   
   calc:=1 for ng eqn on; =0 for eqn off
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,ngsp)
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: istgonxy
   
   calc:=1 for tg eqn on; =0 for eqn off
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,ngsp)
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isphionxy
   
   calc:=1 for phi eqn on; =0 for eqn off
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isnioffxy
   
   user:=1, ni eqn off; =0 for eqn on
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nisp)
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isupoffxy
   
   user:=1, up eqn off; =0 for eqn on
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nisp)
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isteoffxy
   
   user:=1, te eqn off; =0 for eqn on
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: istioffxy
   
   user:=1, ti eqn off; =0 for eqn on
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isngoffxy
   
   user:=1, ng eqn off; =0 for eqn on
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,ngsp)
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: istgoffxy
   
   user:=1, tg eqn off; =0 for eqn on
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,ngsp)
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: isphioffxy
   
   user:=1, phi eqn off;=0 for eqn on
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: UEpar
   :Type: integer
   :Unit: 
.. py:attribute:: fdtnixy
   
   user:=1 for ni eqn off; =0 for eqn on
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nisp)
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: fdtupxy
   
   user:=1 for up eqn off; =0 for eqn on
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nisp)
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: fdttexy
   
   user:=1 for te eqn off; =0 for eqn on
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: fdttixy
   
   user:=1 for ti eqn off; =0 for eqn on
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: fdtngxy
   
   user:=1 for ng eqn off; =0 for eqn on
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,ngsp)
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: fdttgxy
   
   user:=1 for tg eqn off; =0 for eqn on
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,ngsp)
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: fdtphixy
   
   user:=1 for phi eqn off; =0 for eqn on
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: UEpar
   :Type: double
   :Unit: 
.. py:attribute:: iondenseqn
   
   ion continuity equation
   
   :Default: 
   :Dimension: 
   :Group: Model_choice
   :Type: character(8)
   :Unit: 
.. py:attribute:: cnflux
   
   coef for particle flux in n-eq. (resco)
   
   :Default: 
   :Dimension: (6)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cngfx
   
   scale fac for flux from grad_x T_g in gas eqn.
   
   :Default: 
   :Dimension: (6)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cngfy
   
   scale fac for flux from grad_y T_g in gas eqn.
   
   :Default: 
   :Dimension: (6)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cngmom
   
   mom. cx-loss coeff for diffusve-neut hydr only
   
   :Default: 
   :Dimension: (31)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cmwall
   
   mom. wall-loss coeff for diff-neut hydr only
   
   :Default: 
   :Dimension: (31)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cngtgx
   
   X-flux coef for gas comp. of Ti eqn.
   
   :Default: 
   :Dimension: (6)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cngtgy
   
   Y-flux coef for gas comp. of Ti eqn.
   
   :Default: 
   :Dimension: (6)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: rld2dxg
   
   ratio of gas decay-length to dx via artificial diff.
   
   :Default: 
   :Dimension: (6)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: rld2dyg
   
   ratio of gas decay-length to dy via artificial diff.
   
   :Default: 
   :Dimension: (6)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cngflox
   
   fac for x-flux from convection in ng-eqn.
   
   :Default: 
   :Dimension: (6)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cngfloy
   
   fac for y-flux from convection in ng-eqn.
   
   :Default: 
   :Dimension: (6)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cngniflox
   
   fac for rel ion-neut x-vel in ng-eqn
   
   :Default: 
   :Dimension: (31,6)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cngnifloy
   
   fac for rel ion-neut y-vel in ng-eqn
   
   :Default: 
   :Dimension: (31,6)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cdifg
   
   scale factor for gas diffusion coeff.
   
   :Default: 
   :Dimension: (6)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: lgmax
   
   max gas scale length for calc particle D_g
   
   :Default: 
   :Dimension: (6)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: lgtmax
   
   max gas scale length for calc. thermal D_g
   
   :Default: 
   :Dimension: (6)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: rtg2ti
   
   ratio of gas temp to ion temp
   
   :Default: 
   :Dimension: (6)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: tgas
   
   value of tg if istgcon=1
   
   :Default: 
   :Dimension: (6)
   :Group: Coefeq
   :Type: double
   :Unit: eV
.. py:attribute:: cfvisxy
   
   Coef. mult fmixy(ifld)
   
   :Default: 
   :Dimension: (1:31)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfvcsx
   
   Coefs for x-visc. in ti-eq. with ismcnon>0
   
   :Default: 
   :Dimension: (1:31)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfvcsy
   
   Coefs for y-visc. in ti-eq. with ismcnon>0
   
   :Default: 
   :Dimension: (1:31)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfvgpx
   
   Coefs for x components of v.grad(p) in ti-eq
   
   :Default: 
   :Dimension: (1:31)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfvgpy
   
   Coefs for y components of v.grad(p) in ti-eq
   
   :Default: 
   :Dimension: (1:31)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cpiup
   
   mult. press. grad term in up eqn
   
   :Default: 
   :Dimension: (31)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cftgdiss
   
   coef mult tg*nu_diss eng loss
   
   :Default: 
   :Dimension: (6)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfhcxgc
   
   Coef constant pol heat conduct (chixg_use)
   
   :Default: 
   :Dimension: (6)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: cfhcygc
   
   Coef constant rad heat conduct (chiyg_use)
   
   :Default: 
   :Dimension: (6)
   :Group: Coefeq
   :Type: double
   :Unit: 
.. py:attribute:: phiwi
   
   /(nx+2)*0./
   PF wall phi profile if iphicwi=4;user set
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Bcond
   :Type: double
   :Unit: eV
.. py:attribute:: phiwo
   
   /(nx+2)*0./
   outer wall phi profile if iphicwo=4;user set
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Bcond
   :Type: double
   :Unit: eV
.. py:attribute:: ncore
   
   core ion dens if isnicore=1
   
   :Default: 
   :Dimension: (31)
   :Group: Bcond
   :Type: double
   :Unit: m**-3
.. py:attribute:: upcore
   
   core ion parall vel. if isupcore=0
   
   :Default: 
   :Dimension: (31)
   :Group: Bcond
   :Type: double
   :Unit: m/s
.. py:attribute:: ngcore
   
   core gas dens if isngcore=1
   
   :Default: 
   :Dimension: (6)
   :Group: Bcond
   :Type: double
   :Unit: m**-3
.. py:attribute:: tgcore
   
   core gas temp if istgcore=1
   
   :Default: 
   :Dimension: (6)
   :Group: Bcond
   :Type: double
   :Unit: eV
.. py:attribute:: isnicore
   
   switch for ion-density core B.C.
   =1, set uniform, fixed density, ncore
   =0, set flux to curcore/sy locally in ix
   =2, set flux & ni over range
   =3, set icur=curcore-recycc*fngy, const ni
   =4, use impur. source terms (impur only)
   =5, set d(ni)/dy=-ni/lynicore at midp &
   ni constant poloidally
   
   :Default: 
   :Dimension: (31)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isupcore
   
   =0 sets up=upcore on core bdry
   =1 sets d(up)/dy=0 on the core bdry
   =2 sets d^2(up)/dy^2 = 0
   =3 sets fmiy = 0
   =4 sets tor. ang mom flux = lzflux & n*up/R=const
   =5 sets ave tor vel = utorave & n*up/R=const
   
   :Default: 
   :Dimension: (31)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isngcore
   
   switch for neutral-density core B.C.
   =0, set loc flux= -(1-albedoc)*ng*vtg/4
   =1, set uniform, fixed density, ngcore
   =2, not available
   =3, extrapolation, but limited
   =anything else, set zero deriv which was
   prev default inert hy
   anything else same as =0
   
   :Default: 
   :Dimension: (6)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: istgcore
   
   switch for neutral-density core B.C.
   =0, set tg(ixcore,0,igsp)=ti(ixcore,0)
   =1, set fixed temp tgcore(igsp)
   if > 1, set zero grad; tg(,0,)=tg(,1,)
   
   :Default: 
   :Dimension: (6)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: curcore
   
   value of current from core if isnicore=0
   
   :Default: 
   :Dimension: (1:31)
   :Group: Bcond
   :Type: double
   :Unit: A
.. py:attribute:: lzcore
   
   tor. ang. mom dens core bdry; phi eqn
   
   :Default: 
   :Dimension: (1:31)
   :Group: Bcond
   :Type: double
   :Unit: kg/ms
.. py:attribute:: lzflux
   
   tor. ang. mom dens flux core bdry; up eqn
   
   :Default: 
   :Dimension: (1:31)
   :Group: Bcond
   :Type: double
   :Unit: kg/s**2
.. py:attribute:: utorave
   
   ave tor ave vel = utorave; up eqn
   
   :Default: 
   :Dimension: (1:31)
   :Group: Bcond
   :Type: double
   :Unit: m/s
.. py:attribute:: istgwc
   
   switch for outer-wall BC on Tg(,0,igsp)
   =0, set fixed temp to tgwall
   =1, use extrapolation BC
   =2, set Tg scale length to lytg(2,
   >2, report error in input
   
   :Default: 
   :Dimension: (6)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: istgpfc
   
   switch for PF BC on Tg(,0,igsp)
   =0, set fixed temp to tgwall
   =1, use extrapolation BC
   =2, set Tg scale length to lytg(1,
   >2, report error in input
   
   :Default: 
   :Dimension: (6)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: tewalli
   
   /(nx+2)*0./
   inner wall Te for istepfc=1.; = tedge if not set
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Bcond
   :Type: double
   :Unit: eV
.. py:attribute:: tiwalli
   
   /(nx+2)*0./
   inner wall Ti for istipfc=1.; = tedge if not set
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Bcond
   :Type: double
   :Unit: eV
.. py:attribute:: tewallo
   
   /(nx+2)*0./
   outer wall Te for istewc=1.; = tedge if not set
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Bcond
   :Type: double
   :Unit: eV
.. py:attribute:: tiwallo
   
   /(nx+2)*0./
   outer wall Ti for istiwc=1.; = tedge if not set
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Bcond
   :Type: double
   :Unit: eV
.. py:attribute:: tgwall
   
   Wall gas temp BC
   
   :Default: 
   :Dimension: (6)
   :Group: Bcond
   :Type: double
   :Unit: eV
.. py:attribute:: lyte
   
   decaying rad Te grad leng;(1,2) istepfc,wc=3
   
   :Default: 
   :Dimension: (1:2)
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: lytex
   
   pol dep radial te grad length if set < 1e5
   istepfc,wc=3: 1:2=i:o, 2nd dim ix
   
   :Default: 
   :Dimension: (2,0:nx+1)
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: lyti
   
   decaying rad Ti grad leng;(1,2) istipfc,wc=3
   
   :Default: 
   :Dimension: (1:2)
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: lytg
   
   rad tg scale length: PF (1,; Outer(2,
   
   :Default: 
   :Dimension: (1:2,6)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: lytix
   
   pol dep radial ti grad length if set < 1e5
   istipfc,wc=3: 1:2=i:o, 2nd dim ix
   
   :Default: 
   :Dimension: (2,0:nx+1)
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: lyphi
   
   decaying rad phi grad leng;(1,2) iphibcwi,o=3
   
   :Default: 
   :Dimension: (1:2)
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: lyphix
   
   pol dep radial phi grad length if set < 1e5
   isphipfc,wc=3: 1:2=i:o, 2nd dim ix
   
   :Default: 
   :Dimension: (2,0:nx+1)
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: isupss
   
   =0, up=cs; =1, up>=1; =-1, dup/dx=0 at plates
   
   :Default: 
   :Dimension: (31)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isnwconi
   
   switch for private-flux wall (iy=0) density B.C.
   =0, old case; if ifluxni=0, dn/dy=0; if ifluxni=1, fniy=0
   =1, fixed density to nwalli(ix) array
   =2, extrapolation B.C.
   =3, approx grad-length lyni, but limited by nwimin
   
   :Default: 
   :Dimension: (1:31)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isnwcono
   
   switch for outer wall (iy=ny+1) density B.C.
   =0, old case; if ifluxni=0, dn/dy=0; if ifluxni=1, fniy=0
   =1, fixed density to nwallo(ix) array
   =2, extrapolation B.C.
   =3, approx grad-length lyni, but limited by nwomin
   
   :Default: 
   :Dimension: (1:31)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: nwalli
   
   inner wall dens set by isnwconi
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Bcond
   :Type: double
   :Unit: m**-3
.. py:attribute:: nwallo
   
   outer wall dens set by isnwcono
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Bcond
   :Type: double
   :Unit: m**-3
.. py:attribute:: nwimin
   
   min inner wall dens if isnwconi=3
   
   :Default: 
   :Dimension: (31)
   :Group: Bcond
   :Type: double
   :Unit: m**-3
.. py:attribute:: nwomin
   
   min outer wall dens if isnwcono=3
   
   :Default: 
   :Dimension: (31)
   :Group: Bcond
   :Type: double
   :Unit: m**-3
.. py:attribute:: ncoremin
   
   min ncore for isnicore=5
   
   :Default: 
   :Dimension: (31)
   :Group: Bcond
   :Type: double
   :Unit: m**-3
.. py:attribute:: lyni
   
   rad dens grad length -isnwconi,o=3
   
   :Default: 
   :Dimension: (2)
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: lynix
   
   pol dep radial dens grad length if set < 1e5
   isnwconi,o=3: 1:2=i:o, 2nd dim ix, 3rd spec
   
   :Default: 
   :Dimension: (2,0:nx+1,nisp)
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: lynicore
   
   ni core BC rad scale-length if
   isnicore=5
   
   :Default: 
   :Dimension: (31)
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: lyup
   
   radial up grad length if isupwi,o=3: 1:2=i:o
   
   :Default: 
   :Dimension: (2)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: lyupx
   
   pol dep radial up grad length if set < 1e5
   isupwi,o=3: indices,1:2=i:o, 2nd dim ix, 3rd spec
   
   :Default: 
   :Dimension: (2,0:nx+1,nusp)
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: igasi
   
   Gas currents from inner wall (iy=0)
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: double
   :Unit: Amp
.. py:attribute:: igaso
   
   Gas currents from outer wall (iy=ny+1)
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: double
   :Unit: Amp
.. py:attribute:: igspsori
   
   index of gas species for inner wall sources
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: igspsoro
   
   index of gas species for outer wall sources
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: issorlb
   
   flag for coord. origin of source.;=1, left plate;
   =0, right plate
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: jxsori
   
   xgasi=0. is located at left boundary
   of mesh region jxsori
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: jxsoro
   
   xgaso=0. is located at left boundary
   of mesh region jxsoro
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: xgasi
   
   location of inner wall sources; if issorlb(i)=1,0
   measured from left plate, right plate
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: xgaso
   
   location of outer wall sources; if issorlb(i)=1,0,
   measured from left plate, right plate
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: wgasi
   
   total cosine widths of inner wall gas sources
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: wgaso
   
   total cosine widths of outer wall gas sources
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: albdsi
   
   albedos at inner gas source locations
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: albdso
   
   albedos at outer gas source locations
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: double
   :Unit: m
.. py:attribute:: chemsputi
   
   chem sputt coeff, priv flux surface, flux(i)=
   sum(chemsputi(i,j)*ng(j)*vt*sy)
   
   :Default: 
   :Dimension: (10,10)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: chemsputo
   
   chem sputt coeff, outer wall - see chemsputi def.
   
   :Default: 
   :Dimension: (10,10)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: matwsi
   
   material wall at inner gas source locations
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: matwso
   
   material wall at outer gas source locations
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: issori
   
   starting ix cell index for inner source
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: iesori
   
   ending ix cell index for inner source
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: issoro
   
   starting ix cell index for outer source
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: iesoro
   
   ending ix cell index for outer source
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: iwalli
   
   current from inner source region isor for coupling
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: iwallo
   
   current from outer source region isor for coupling
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: ncpli
   
   flag for coupling between inner srce isor & ncpli
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: ncplo
   
   flag for coupling between outer srce isor & ncpli
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: cplsori
   
   coeff. giving coupling from inner isor to ncpli
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: cplsoro
   
   coeff. giving coupling from outer isor to ncpli
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: iscpli
   
   (=1) => ix pt involved in inner bndry coupling
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: iscplo
   
   (=1) => ix pt involved in outer bndry coupling
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: fwsori
   
   profile of inner wall source isor (missing igasi)
   
   :Default: 
   :Dimension: (0:nx+1,10)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: fwsoro
   
   profile of outer wall source isor (missing igasi)
   
   :Default: 
   :Dimension: (0:nx+1,10)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: fngysi
   
   gas input flux from igasi on inner wall (calc)
   
   :Default: 
   :Dimension: (0:nx+1,ngsp)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: fngyi_use
   
   user supplied gas input flux*area
   
   :Default: 
   :Dimension: (0:nx+1,ngsp)
   :Group: Bcond
   :Type: double
   :Unit: 1/m**3s
.. py:attribute:: fngysig
   
   global value of fngysi if domain decomp (parll)
   
   :Default: 
   :Dimension: (0:nxg+1,ngsp)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: fngyso
   
   gas input flux from igaso on outer wall (calc)
   
   :Default: 
   :Dimension: (0:nx+1,ngsp)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: fngyo_use
   
   user supplied gas input flux*area
   
   :Default: 
   :Dimension: (0:nx+1,ngsp)
   :Group: Bcond
   :Type: double
   :Unit: 1/m**3s
.. py:attribute:: fngysog
   
   global value of fngyso if domain-decomp (parll)
   
   :Default: 
   :Dimension: (0:nxg+1,ngsp)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: albedoo
   
   albedo outer iy=ny+1 surface for neutrals (calc)
   
   :Default: 
   :Dimension: (0:nx+1,ngsp)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: albedoog
   
   global val albedoo if domain-decomp (parll)
   
   :Default: 
   :Dimension: (0:nxg+1,ngsp)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: albedoi
   
   albedo of inner iy=0 surface for neutrals (calc)
   
   :Default: 
   :Dimension: (0:nx+1,ngsp)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: albedoig
   
   global val albedoi if domain-decomp (parll)
   
   :Default: 
   :Dimension: (0:nxg+1,ngsp)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: matwallo
   
   flag (=1) denoting outer material side wall
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: matwallog
   
   global val matwallo if domain-decomp (parll)
   
   :Default: 
   :Dimension: (0:nxg+1)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: matwalli
   
   flag (=1) denoting inner material side wall (pf)
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: matwallig
   
   global val matwalli if domain-decomp (parll)
   
   :Default: 
   :Dimension: (0:nxg+1)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isfixlb
   
   =1 fixes values left bndry;=2 for symm. pt.
   
   :Default: 
   :Dimension: (2)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isfixrb
   
   =2 for symmetry pt. at ix=nx+1
   
   :Default: 
   :Dimension: (2)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: nib
   
   left plate density for isfixlb=1
   
   :Default: 
   :Dimension: (20)
   :Group: Bcond
   :Type: double
   :Unit: m**-3
.. py:attribute:: upb
   
   left plate parallel velocity for isfixlb=1
   
   :Default: 
   :Dimension: (20)
   :Group: Bcond
   :Type: double
   :Unit: m/s
.. py:attribute:: nibprof
   
   radial profile of nib
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: upbprof
   
   radial profile of upb
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: tebprof
   
   radial profile of teb
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: tibprof
   
   radial profile of tib
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: isupwi
   
   =0 sets up=0 on inner wall
   =1 sets fmiy=0 (parallel mom-dens y-flux)
   =2 sets dup/dy=0 on inner wall
   =3 sets (1/up)dup/dy=1/lyup(1) scale length
   
   :Default: 
   :Dimension: (31)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isupwo
   
   =0 sets up=0 on outer wall
   =1 sets fmiy=0 (parallel mom-dens y-flux)
   =2 sets dup/dy=0 on outer wall
   =3 sets (1/up)dup/dy=1/lyup(2) scale length
   
   :Default: 
   :Dimension: (31)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isph_sput
   
   flag for plate sputtering;
   0=old fixed case; 1=DIVIMP/JET phys sputt fits
   =2 adds h-ion chem sputt;=3 adds h-neut c_sput
   
   :Default: 
   :Dimension: (6)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isi_sputw
   
   flag for outer wall ion-based sputter;
   =0, no ion sputtering
   =1 adds phys ion sputt; =2 adds chem ion sputt
   
   :Default: 
   :Dimension: (6)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isi_sputpf
   
   flag for priv flux ion-based sputter;
   =0, no ion sputtering
   =1 adds phys ion sputt; =2 adds chem ion sputt
   
   :Default: 
   :Dimension: (6)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: isch_sput
   
   chem sputt. opt; 0=old;
   5=Roth,G-R; 6=Haasz97; 7=Haasz97+Davis at low E
   
   :Default: 
   :Dimension: (6)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: tvwallo
   
   user outer wall temp if iswalltempc=0
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Bcond
   :Type: double
   :Unit: K
.. py:attribute:: tvwalli
   
   user inner wall temp if iswalltempc=0
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Bcond
   :Type: double
   :Unit: K
.. py:attribute:: tvplatlb
   
   user left plate temp if isplttempc=0
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Bcond
   :Type: double
   :Unit: K
.. py:attribute:: tvplatrb
   
   user left plate temp if isplttempc=0
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Bcond
   :Type: double
   :Unit: K
.. py:attribute:: yld_carbi
   
   chem sputt. yield, inner wall if isch_sput=5,6
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: yld_carbo
   
   chem sputt. yield, outer wall if isch_sput=5,6
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: fchemygwi
   
   fac mult pf wall gas chem yield if isch_sput>0
   
   :Default: 
   :Dimension: (6)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: fchemygwo
   
   fac mult outer wall gas chem yield; isch_sput>0
   
   :Default: 
   :Dimension: (6)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: fchemyiwi
   
   fac mult pf wall ion chem yield if isch_sput>0
   
   :Default: 
   :Dimension: (6)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: fchemyiwo
   
   fac mult outer wall ion chem yield; isch_sput>0
   
   :Default: 
   :Dimension: (6)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: fphysyiwi
   
   fac mult pf wall ion phys yield if isch_sput>0
   
   :Default: 
   :Dimension: (6)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: fphysyiwo
   
   fac mult outer wall ion phys yield; isch_sput>0
   
   :Default: 
   :Dimension: (6)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: fchemylb
   
   fac*inner plt gas chem yield; isch_sput>0
   
   :Default: 
   :Dimension: (6,2)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: fchemyrb
   
   fac*outer plt gas chem yield; isch_sput>0
   
   :Default: 
   :Dimension: (6,2)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: fphysylb
   
   fac*inner plt ion phys sp yield;isch_sput>0
   
   :Default: 
   :Dimension: (6,2)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: fphysyrb
   
   fac*outer plt ion phys sp yield;isch_sput>0
   
   :Default: 
   :Dimension: (6,2)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: fqpsatlb
   
   ion saturation current at left boundary
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: fqpsatrb
   
   ion saturation current at right boundary
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: fvapi
   
   scale factor for inner evap vapor source
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: avapi
   
   linear coeff. for inner evap vapor source
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: bvapi
   
   exponent coeff. for inner evap vapor source
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: fvapo
   
   scale factor for outer evap vapor source
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: avapo
   
   linear coeff. for outer evap vapor source
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: bvapo
   
   exponent coeff. for outer evap vapor source
   
   :Default: 
   :Dimension: (10)
   :Group: Bcond
   :Type: double
   :Unit: 
.. py:attribute:: tvapi
   
   inner wall temp for evap; input after alloc
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Bcond
   :Type: double
   :Unit: K
.. py:attribute:: tvapo
   
   outer wall temp for evap; input after alloc
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Bcond
   :Type: double
   :Unit: K
.. py:attribute:: totfeexl
   
   elec polod energy flux*area on 'left' plate
   
   :Default: 
   :Dimension: (0:ny+1,nxpt)
   :Group: Bcond
   :Type: double
   :Unit: W
.. py:attribute:: totfeexr
   
   elec polod energy flux*area on 'right' plate
   
   :Default: 
   :Dimension: (0:ny+1,nxpt)
   :Group: Bcond
   :Type: double
   :Unit: W
.. py:attribute:: totfeixl
   
   elec polod energy flux*area on 'left' plate
   
   :Default: 
   :Dimension: (0:ny+1,nxpt)
   :Group: Bcond
   :Type: double
   :Unit: W
.. py:attribute:: totfeixr
   
   elec polod energy flux*area on 'right' plate
   
   :Default: 
   :Dimension: (0:ny+1,nxpt)
   :Group: Bcond
   :Type: double
   :Unit: W
.. py:attribute:: istglb
   
   =0 for tg=tgwall; =1 for extrap
   
   :Default: 
   :Dimension: (6)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: istgrb
   
   =0 for tg=tgwall; =1 for extrap
   
   :Default: 
   :Dimension: (6)
   :Group: Bcond
   :Type: integer
   :Unit: 
.. py:attribute:: ue_part_fluxelb
   
   inner plt elec flux
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Outpwall
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: ue_part_fluxerb
   
   outer plt elec flux
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Outpwall
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: ue_part_fluxeyi
   
   inner (PF) wall elec part flux
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Outpwall
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: ue_part_fluxeyo
   
   outer (PF) wall elec part flux
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Outpwall
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: ue_part_fluxh2p1lb
   
   inner plt deut ion part flux
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Outpwall
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: ue_part_fluxh2p1rb
   
   outer plt deut ion part flux
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Outpwall
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: ue_part_fluxh2p1yi
   
   inner (PF) wall deut ion part flux
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Outpwall
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: ue_part_fluxh2p1yo
   
   outer (PF) wall deut ion part flux
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Outpwall
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: ue_part_fluxh2lb
   
   inner plt deut neut part flux
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Outpwall
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: ue_part_fluxh2rb
   
   outer plt deut neut part flux
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Outpwall
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: ue_part_fluxh2yi
   
   inner (PF) wall deut neut part flux
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Outpwall
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: ue_part_fluxh2yo
   
   outer (PF) wall deut neut part flux
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Outpwall
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: ue_heat_fluxelb
   
   inner plt elec heat flux
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Outpwall
   :Type: double
   :Unit: W/m**2
.. py:attribute:: ue_heat_fluxerb
   
   outer plt elec heat flux
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Outpwall
   :Type: double
   :Unit: W/m**2
.. py:attribute:: ue_heat_fluxeyi
   
   inner (PF) wall elec heat flux
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Outpwall
   :Type: double
   :Unit: W/m**2
.. py:attribute:: ue_heat_fluxeyo
   
   outer (PF) wall elec heat flux
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Outpwall
   :Type: double
   :Unit: W/m**2
.. py:attribute:: ue_heat_fluxh2p1lb
   
   inner plt deut ion heat flux
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Outpwall
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: ue_heat_fluxh2p1rb
   
   outer plt deut ion heat flux
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Outpwall
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: ue_heat_fluxh2p1yi
   
   inner (PF) wall deut ion heat flux
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Outpwall
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: ue_heat_fluxh2p1yo
   
   outer (PF) wall deut ion heat flux
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Outpwall
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: ue_heat_fluxh2lb
   
   inner plt deut neut heat flux
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Outpwall
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: ue_heat_fluxh2rb
   
   outer plt deut neut heat flux
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Outpwall
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: ue_heat_fluxh2yi
   
   inner (PF) wall deut neut heat flux
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Outpwall
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: ue_heat_fluxh2yo
   
   outer (PF) wall deut neut heat flux
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Outpwall
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: ue_mean_engelb
   
   inner plt elec mean energy
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Outpwall
   :Type: double
   :Unit: J
.. py:attribute:: ue_mean_engerb
   
   outer plt elec mean energy
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Outpwall
   :Type: double
   :Unit: J
.. py:attribute:: ue_mean_engeyi
   
   inner (PF) wall elec mean energy
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Outpwall
   :Type: double
   :Unit: J
.. py:attribute:: ue_mean_engeyo
   
   outer (PF) wall elec mean energy
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Outpwall
   :Type: double
   :Unit: J
.. py:attribute:: ue_mean_engh2p1lb
   
   inner plt deut ion mean energy
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Outpwall
   :Type: double
   :Unit: J
.. py:attribute:: ue_mean_engh2p1rb
   
   outer plt deut ion mean energy
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Outpwall
   :Type: double
   :Unit: J
.. py:attribute:: ue_mean_engh2p1yi
   
   inner (PF) wall deut ion mean energy
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Outpwall
   :Type: double
   :Unit: J
.. py:attribute:: ue_mean_engh2p1yo
   
   outer (PF) wall deut ion mean energy
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Outpwall
   :Type: double
   :Unit: J
.. py:attribute:: ue_mean_engh2lb
   
   inner plt deut neut mean energy
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Outpwall
   :Type: double
   :Unit: J
.. py:attribute:: ue_mean_engh2rb
   
   outer plt deut neut mean energy
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Outpwall
   :Type: double
   :Unit: J
.. py:attribute:: ue_mean_engh2yi
   
   inner (PF) wall deut neut mean energy
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Outpwall
   :Type: double
   :Unit: J
.. py:attribute:: ue_mean_engh2yo
   
   outer (PF) wall deut neut mean energy
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Outpwall
   :Type: double
   :Unit: J
.. py:attribute:: ue_pot_engh2p1lb
   
   inner plt deut ion pot energy
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Outpwall
   :Type: double
   :Unit: J
.. py:attribute:: ue_pot_engh2p1rb
   
   outer plt deut ion pot energy
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Outpwall
   :Type: double
   :Unit: J
.. py:attribute:: ue_pot_engh2p1yi
   
   inner (PF) wall deut ion pot energy
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Outpwall
   :Type: double
   :Unit: J
.. py:attribute:: ue_pot_engh2p1yo
   
   outer (PF) wall deut ion pot energy
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Outpwall
   :Type: double
   :Unit: J
.. py:attribute:: recylb
   
   tot inner plate recycling coeff. (calc)
   if recylb > 0, recycling coeff
   if in range [-1,0], acts as albedo
   if in range (-2,-1), gives ng=nglfix
   if recylb <= -2, gives ng(1)=ng(0)
   
   :Default: 
   :Dimension: (0:ny+1,6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: recyrb
   
   tot outer plate recycling coeff. (calc)
   if recyrb > 0, recycling coeff
   if in range [-1,0], acts as albedo
   if in range (-2,-1), gives ng=ngrfix
   if recyrb <= -2, gives ng(nx+1)=ng(nx)
   
   :Default: 
   :Dimension: (0:ny+1,6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: recylb_use
   
   inner plate recycling coeff. user input
   
   :Default: 
   :Dimension: (0:ny+1,6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: recyrb_use
   
   outer plate recycling coeff. user input
   
   :Default: 
   :Dimension: (0:ny+1,6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: recycp
   
   recycling coef at plates if ndatlb,rb=0
   
   :Default: 
   :Dimension: (6)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: recycflb
   
   extra factor for recycling at ix=0
   
   :Default: 
   :Dimension: (6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: recycfrb
   
   extra factor for recycling at ix=nx+1
   
   :Default: 
   :Dimension: (6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: recycmlb_use
   
   inner plt mom-recycl coeff user input
   
   :Default: 
   :Dimension: (0:ny+1,6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: recycmrb_use
   
   outer plt mom-recycl coeff user input
   
   :Default: 
   :Dimension: (0:ny+1,6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: recycmlb
   
   total inner plt mom recycling coeff
   
   :Default: 
   :Dimension: (0:ny+1,6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: recycmrb
   
   total outer plt mom recycling coeff
   
   :Default: 
   :Dimension: (0:ny+1,6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: recycc
   
   core recycling coeff. if isnicore=3
   
   :Default: 
   :Dimension: (6)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: albedoc
   
   core neut albedo for isngcore=0
   
   :Default: 
   :Dimension: (6)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: albedolb
   
   albedo at inner plate if ndatlb=0
   
   :Default: 
   :Dimension: (6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: albedorb
   
   albedo at outer plate if ndatrb=0
   
   :Default: 
   :Dimension: (6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: ndatlb
   
   number of recycp data pts on inner plt
   
   :Default: 
   :Dimension: (6,2)
   :Group: Rccoef
   :Type: integer
   :Unit: 
.. py:attribute:: ndatrb
   
   number of recycp data pts on outer plt
   
   :Default: 
   :Dimension: (6,2)
   :Group: Rccoef
   :Type: integer
   :Unit: 
.. py:attribute:: ydatlb
   
   inner data pt location from sep.
   
   :Default: 
   :Dimension: (6,50,2)
   :Group: Rccoef
   :Type: double
   :Unit: m
.. py:attribute:: ydatrb
   
   outer data pt location from sep.
   
   :Default: 
   :Dimension: (6,50,2)
   :Group: Rccoef
   :Type: double
   :Unit: m
.. py:attribute:: rdatlb
   
   inner recycp data for each ydatlb
   
   :Default: 
   :Dimension: (6,50,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: rdatrb
   
   outer recycp data for each ydatrb
   
   :Default: 
   :Dimension: (6,50,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: alblb
   
   inner plate albedo; used if <1 (calc)
   
   :Default: 
   :Dimension: (0:ny+1,6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: albrb
   
   outer plate albedo; used if <1 (calc)
   
   :Default: 
   :Dimension: (0:ny+1,6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: fngxslb
   
   inner plt liq vapor gas sour. if sputtlb>0
   
   :Default: 
   :Dimension: (0:ny+1,6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 1/s
.. py:attribute:: fngxsrb
   
   outer plt liq vapor gas sour. if sputtlb>0
   
   :Default: 
   :Dimension: (0:ny+1,6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 1/s
.. py:attribute:: fngxlb_use
   
   user external left plate source
   
   :Default: 
   :Dimension: (0:ny+1,6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 1/s
.. py:attribute:: fngxrb_use
   
   user external left plate source
   
   :Default: 
   :Dimension: (0:ny+1,6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 1/s
.. py:attribute:: adatlb
   
   inner albdedo data for each ydati
   
   :Default: 
   :Dimension: (6,50,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: adatrb
   
   outer albdedo data for each ydati
   
   :Default: 
   :Dimension: (6,50,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: recycw
   
   recycling coef. at side walls
   
   :Default: 
   :Dimension: (6)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: recypf_use
   
   priv flux recycling coef; user input
   
   :Default: 
   :Dimension: (0:nx+1,6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: recywall_use
   
   outer wall recycling coef; user input
   
   :Default: 
   :Dimension: (0:nx+1,6)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: recycwit
   
   tot recyc coeff on PF wall
   
   :Default: 
   :Dimension: (0:nx+1,6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: recycwot
   
   tot recyc coeff on outer wall
   
   :Default: 
   :Dimension: (0:nx+1,6)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: sputtlb
   
   set sputt coef. inner plate (iy,igsp)
   
   :Default: 
   :Dimension: (0:ny+1,6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: sputtrb
   
   set sputt coef. outer plate (iy,igsp)
   
   :Default: 
   :Dimension: (0:ny+1,6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: sputflxlb
   
   calc sput flux inner plate (iy,igsp)
   
   :Default: 
   :Dimension: (0:ny+1,6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: sputflxrb
   
   calc sput flux outer plate (iy,igsp)
   
   :Default: 
   :Dimension: (0:ny+1,6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: sputflxw
   
   calc sput flux outer wall (ix,igsp)
   
   :Default: 
   :Dimension: (0:nx+1,6)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: sputflxpf
   
   calc sput flux PF wall (ix,igsp)
   
   :Default: 
   :Dimension: (0:nx+1,6)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: ngplatlb
   
   ng on inner plate if sputti < -9.9
   
   :Default: 
   :Dimension: (6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: ngplatrb
   
   ng on outer plate if sputto < -9.9
   
   :Default: 
   :Dimension: (6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: igaslb
   
   Gas cur from left-hand plate(s) (ix=0)
   
   :Default: 
   :Dimension: (10,2)
   :Group: Rccoef
   :Type: double
   :Unit: Amp
.. py:attribute:: igasrb
   
   Gas cur from right-hand plate(s) (ix=nx)
   
   :Default: 
   :Dimension: (10,2)
   :Group: Rccoef
   :Type: double
   :Unit: Amp
.. py:attribute:: igspsorlb
   
   gas species index, left-hand plate sources
   
   :Default: 
   :Dimension: (10,2)
   :Group: Rccoef
   :Type: integer
   :Unit: 
.. py:attribute:: igspsorrb
   
   gas species index, right-hand plate sources
   
   :Default: 
   :Dimension: (10,2)
   :Group: Rccoef
   :Type: integer
   :Unit: 
.. py:attribute:: ygaslb
   
   loc of left-plate sources wrt strike pt
   
   :Default: 
   :Dimension: (10,2)
   :Group: Rccoef
   :Type: double
   :Unit: m
.. py:attribute:: ygasrb
   
   loc of right-plate sources wrt strike pt
   
   :Default: 
   :Dimension: (10,2)
   :Group: Rccoef
   :Type: double
   :Unit: m
.. py:attribute:: wgaslb
   
   total cos width of left-plate gas sources
   
   :Default: 
   :Dimension: (10,2)
   :Group: Rccoef
   :Type: double
   :Unit: m
.. py:attribute:: wgasrb
   
   total cos width of right-plate gas sources
   
   :Default: 
   :Dimension: (10,2)
   :Group: Rccoef
   :Type: double
   :Unit: m
.. py:attribute:: fvaplb
   
   scale factor left-plate evap vapor source
   
   :Default: 
   :Dimension: (6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: avaplb
   
   lin coeff left-plate evapor sor
   
   :Default: 
   :Dimension: (6,2)
   :Group: Rccoef
   :Type: double
   :Unit: k**.5/(m**2s)
.. py:attribute:: bvaplb
   
   expon. coeff. left-plate evap vapor source
   
   :Default: 
   :Dimension: (6,2)
   :Group: Rccoef
   :Type: double
   :Unit: K
.. py:attribute:: fvaprb
   
   scale factor right-plate evap vapor source
   
   :Default: 
   :Dimension: (6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 
.. py:attribute:: avaprb
   
   lin coeff right-plate evapor sor
   
   :Default: 
   :Dimension: (6,2)
   :Group: Rccoef
   :Type: double
   :Unit: k**.5/(m**2s)
.. py:attribute:: bvaprb
   
   expon coeff. right-plate evap vapor source
   
   :Default: 
   :Dimension: (6,2)
   :Group: Rccoef
   :Type: double
   :Unit: K
.. py:attribute:: tvaplb
   
   left-plate temp for evap; input after alloc
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Rccoef
   :Type: double
   :Unit: K
.. py:attribute:: tvaprb
   
   right-plate temp for evap; input after alloc
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Rccoef
   :Type: double
   :Unit: K
.. py:attribute:: fngxextlb
   
   inner plt external particle flux*A
   
   :Default: 
   :Dimension: (0:ny+1,6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 1/s
.. py:attribute:: fngxextrb
   
   outer plt external particle flux*A
   
   :Default: 
   :Dimension: (0:ny+1,6,2)
   :Group: Rccoef
   :Type: double
   :Unit: 1/s
.. py:attribute:: fngyexti
   
   inner wall external particle flux*A
   
   :Default: 
   :Dimension: (0:nx+1,6)
   :Group: Rccoef
   :Type: double
   :Unit: 1/s
.. py:attribute:: fngyexto
   
   outer wall external particle flux*A
   
   :Default: 
   :Dimension: (0:nx+1,6)
   :Group: Rccoef
   :Type: double
   :Unit: 1/s
.. py:attribute:: feixextlb
   
   inner plt external energy flux*A
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Rccoef
   :Type: double
   :Unit: J/s
.. py:attribute:: feixextrb
   
   outer plt external energy flux*A
   
   :Default: 
   :Dimension: (0:ny+1,2)
   :Group: Rccoef
   :Type: double
   :Unit: J/s
.. py:attribute:: feiyexti
   
   inner wall external energy flux*A
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Rccoef
   :Type: double
   :Unit: J/s
.. py:attribute:: feiyexto
   
   outer wall external energy flux*A
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Rccoef
   :Type: double
   :Unit: J/s
.. py:attribute:: ixm1
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: ixp1
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: iym1a
   
   for mdsplus use only
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: iyp1a
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Selec
   :Type: integer
   :Unit: 
.. py:attribute:: stretcx
   
   array for stretching gas x-coord.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Selec
   :Type: double
   :Unit: 
.. py:attribute:: iworkin
   
   optional integer array input for vodpk
   
   :Default: 
   :Dimension: (1:25)
   :Group: Opt_input
   :Type: integer
   :Unit: 
.. py:attribute:: rworkin
   
   optional real array input for vodpk
   
   :Default: 
   :Dimension: (1:25)
   :Group: Opt_input
   :Type: double
   :Unit: 
.. py:attribute:: rtolv
   
   relative tol. vector used in convert.m
   
   :Default: 
   :Dimension: (30)
   :Group: Lsode
   :Type: double
   :Unit: 
.. py:attribute:: ipar
   
   integer parameter communication
   
   :Default: 
   :Dimension: (3)
   :Group: Lsode
   :Type: integer
   :Unit: 
.. py:attribute:: rpar
   
   real parameter communication
   
   :Default: 
   :Dimension: (1)
   :Group: Lsode
   :Type: double
   :Unit: 
.. py:attribute:: rtol
   
   relative tol. used by solvers; rtol=rtolv(igrid)*tolbf
   
   :Default: 
   :Dimension: (neqmx)
   :Group: Lsode
   :Type: double
   :Unit: 
.. py:attribute:: atol
   
   absolute tolerance for lsode-like routines
   
   :Default: 
   :Dimension: (neqmx)
   :Group: Lsode
   :Type: double
   :Unit: 
.. py:attribute:: yl
   
   primary variables for ODE's
   
   :Default: 
   :Dimension: (neqmx)
   :Group: Lsode
   :Type: double
   :Unit: 
.. py:attribute:: yldot
   
   time-derivatives of yl's; RHS of ODE's
   
   :Default: 
   :Dimension: (neqmx)
   :Group: Lsode
   :Type: double
   :Unit: 
.. py:attribute:: delta
   
   
   
   :Default: 
   :Dimension: (neqmx)
   :Group: Lsode
   :Type: double
   :Unit: 
.. py:attribute:: iextra
   
   
   
   :Default: 
   :Dimension: (5)
   :Group: Lsode
   :Type: integer
   :Unit: 
.. py:attribute:: info
   
   set for daspk
   info(1) =0, indicates start of new problem to initializd code
   info(2) =0, rtol,atol scalar; = 1, rtol,atol vectors
   info(3) =0, output only at tout; =1, output at intermed. t
   info(4) =0, no restrict on t; =1, stop at t=tstop=rwork(1)
   info(5) =0, daspk gen Jac, only if info(12)=0; =1, ext Jac
   info(6) =0, full mat sol, only if info(12)=0; =1, band sol
   iwork(1)=lbw, iwork(2)=ubw
   info(7) =0, code sets dt_max; =1, rwork(2)=dt_max
   info(8) =0, code sets dt_init; =1, rwork(3)=dt_init
   info(9) =0, maxord=5; =1 iwork(3)=maxord (=<5)
   info(10)=0, no constraints; =1, constr initial cond.;
   =2, constr Y>0; =3, options 1 & 2 added
   info(11)=0, consist init cond; =1, calc init cond(alg + der)
   =2, calc init cond (use Y', calc Y)
   info(12)=0, direct mat sol; =1, Krylov method
   info(13)=0, default Krylov param; =1, iwork(24)=maxl,
   iwork(25)=kmp, iwork(26)=nrmax; rwork(10)=eplidpk
   default Krylov and direct pram: rwork(16)=epnldpk
   info(14)=0, proc after init; =1, stop after init, then
   reset info(1)=0 to avoid recalc init
   (or info(14)=0 & info(11)=0 to avoid recalc init?)
   info(15)=0, no Jac routine; =1 Jac provided
   info(16)=0, err chk on all Y; =1, no err on alg Y
   info(17)=0, init calc control default; =1,input init control
   info(18)=0. no init print; =1, min print; =2, full print
   info(19) not specified
   info(20) not specified
   
   :Default: 
   :Dimension: (20)
   :Group: Lsode
   :Type: integer
   :Unit: 
.. py:attribute:: iopt
   
   opt. input/output array
   
   :Default: 
   :Dimension: (40)
   :Group: Parallv
   :Type: integer
   :Unit: 
.. py:attribute:: ropt
   
   opt. input/output array
   
   :Default: 
   :Dimension: (40)
   :Group: Parallv
   :Type: double
   :Unit: 
.. py:attribute:: icnstr
   
   nksol constraint array; =1 means must be > 0
   
   :Default: 
   :Dimension: (neqmx)
   :Group: Constraints
   :Type: integer
   :Unit: 
.. py:attribute:: constr
   
   kinsol constraint array; =1 means must be > 0
   
   :Default: 
   :Dimension: (neqmx)
   :Group: Constraints
   :Type: double
   :Unit: 
.. py:attribute:: ylprevc
   
   yl vector from previous call from vodpk
   
   :Default: 
   :Dimension: (neqmx)
   :Group: Constraints
   :Type: double
   :Unit: 
.. py:attribute:: ylchng
   
   change in yl from prev. step, yl-ylprevc
   
   :Default: 
   :Dimension: (neqmx)
   :Group: Constraints
   :Type: double
   :Unit: 
.. py:attribute:: norm_cons
   
   normalization constants (calculated)
   
   :Default: 
   :Dimension: (numvar)
   :Group: Ynorm
   :Type: double
   :Unit: 
.. py:attribute:: floor_cons
   
   floor constants (calculated)
   
   :Default: 
   :Dimension: (numvar)
   :Group: Ynorm
   :Type: double
   :Unit: 
.. py:attribute:: n0
   
   normalization ion density
   
   :Default: 
   :Dimension: (31)
   :Group: Ynorm
   :Type: double
   :Unit: m**-3
.. py:attribute:: n0g
   
   normalization gas density
   
   :Default: 
   :Dimension: (6)
   :Group: Ynorm
   :Type: double
   :Unit: m**-3
.. py:attribute:: fnorm
   
   normalization momentum flux(calc)
   
   :Default: 
   :Dimension: (1:nusp)
   :Group: Ynorm
   :Type: double
   :Unit: kg/m**2 s
.. py:attribute:: suscal
   
   scale factors for yl's in nksol routine (nominally=1)
   
   :Default: 
   :Dimension: (neqmx)
   :Group: Ynorm
   :Type: double
   :Unit: 
.. py:attribute:: sfscal
   
   scale factors for f's in nksol routine (nominally=1)
   
   :Default: 
   :Dimension: (neqmx)
   :Group: Ynorm
   :Type: double
   :Unit: 
.. py:attribute:: yloext
   
   last var for d(yl)/dt set externally
   
   :Default: 
   :Dimension: (neqmx)
   :Group: Ynorm
   :Type: double
   :Unit: 
.. py:attribute:: parvis
   
   factor times parallel visc.
   
   :Default: 
   :Dimension: (1:31)
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: travis
   
   value of perp. visc.
   
   :Default: 
   :Dimension: (1:31)
   :Group: Comtra
   :Type: double
   :Unit: m**2/s
.. py:attribute:: difni
   
   value of density radial diff. coef.
   
   :Default: 
   :Dimension: (1:31)
   :Group: Comtra
   :Type: double
   :Unit: m**2/s
.. py:attribute:: dif4order
   
   4th ord ion density radial diff. coef.
   
   :Default: 
   :Dimension: (1:31)
   :Group: Comtra
   :Type: double
   :Unit: m**4/s
.. py:attribute:: difgy4order
   
   4th ord gas density radial diff. coef.
   
   :Default: 
   :Dimension: (1:6)
   :Group: Comtra
   :Type: double
   :Unit: m**4/s
.. py:attribute:: difgx4order
   
   4th ord gas density poloid diff. coef.
   
   :Default: 
   :Dimension: (1:6)
   :Group: Comtra
   :Type: double
   :Unit: m**4/s
.. py:attribute:: difax
   
   poloid. diff coeff.,scaled with dn/dx
   
   :Default: 
   :Dimension: (1:31)
   :Group: Comtra
   :Type: double
   :Unit: m*2/s
.. py:attribute:: difnit
   
   turb. radial diff. multiplier
   
   :Default: 
   :Dimension: (1:31)
   :Group: Comtra
   :Type: double
   :Unit: none
.. py:attribute:: difpr
   
   value of pressure radial diff. coef.
   
   :Default: 
   :Dimension: (1:31)
   :Group: Comtra
   :Type: double
   :Unit: m**2/s
.. py:attribute:: difni2
   
   value of density e_ll x e_r diff. coef.
   
   :Default: 
   :Dimension: (1:31)
   :Group: Comtra
   :Type: double
   :Unit: m**2/s
.. py:attribute:: difpr2
   
   value of pressure e_ll x e_r diff. coef
   
   :Default: 
   :Dimension: (1:31)
   :Group: Comtra
   :Type: double
   :Unit: m**2/s
.. py:attribute:: difutm
   
   value of toroidal mom. diff. coef
   
   :Default: 
   :Dimension: (1:31)
   :Group: Comtra
   :Type: double
   :Unit: m**2/s
.. py:attribute:: difniv
   
   dens diff. if isbohmcalc=3, varys w/B
   
   :Default: 
   :Dimension: (0:ny+1,1:nisp)
   :Group: Comtra
   :Type: double
   :Unit: m**2/s
.. py:attribute:: difprv
   
   press diff if isbohmcalc=3, varys w/B
   
   :Default: 
   :Dimension: (0:ny+1,1:nisp)
   :Group: Comtra
   :Type: double
   :Unit: m**2/s
.. py:attribute:: difniv2
   
   dens2 diff. if isbohmcalc=3, varys w/B
   
   :Default: 
   :Dimension: (0:ny+1,1:nisp)
   :Group: Comtra
   :Type: double
   :Unit: m**2/s
.. py:attribute:: travisv
   
   viscosity. if isbohmcalc=3, varys w/B
   
   :Default: 
   :Dimension: (0:ny+1,1:nisp)
   :Group: Comtra
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kyev
   
   elec eng chi if isbohmcalc=3, vary(0:ny+1)s w/B
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Comtra
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kyiv
   
   ion eng chi if isbohmcalc=3, varys w/B
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Comtra
   :Type: double
   :Unit: m**2/s
.. py:attribute:: difutmv
   
   dens diff. if isbohmcalc=3, varys w/B
   
   :Default: 
   :Dimension: (0:ny+1,1:nisp)
   :Group: Comtra
   :Type: double
   :Unit: m**2/s
.. py:attribute:: vconyv
   
   convec radial vel if isbohmcalc=3
   
   :Default: 
   :Dimension: (0:ny+1,1:nisp)
   :Group: Comtra
   :Type: double
   :Unit: m/s
.. py:attribute:: vcony
   
   value of constant radial velocity
   
   :Default: 
   :Dimension: (1:31)
   :Group: Comtra
   :Type: double
   :Unit: m/s
.. py:attribute:: flalfgx
   
   poloidal gas diff. flux limit
   
   :Default: 
   :Dimension: (10)
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalfgy
   
   radial gas diff. flux limit
   
   :Default: 
   :Dimension: (10)
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalfgxy
   
   nonorthog pol-face gas flux limit
   
   :Default: 
   :Dimension: (10)
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: sigcxms
   
   cx x-sect for (ifld,igsp) coll
   
   :Default: 
   :Dimension: (nisp,ngsp)
   :Group: Comtra
   :Type: double
   :Unit: m**2
.. py:attribute:: rcxighg
   
   ratio of charge-exchange rate for
   ng_imp+ni_hydrn -> ng_hydrn+ni_Z=1_imp
   
   :Default: 
   :Dimension: (6)
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: kelighi
   
   elastic coll. coeff:imp_gas+hyd_ion
   
   :Default: 
   :Dimension: (6)
   :Group: Comtra
   :Type: double
   :Unit: m**3/s
.. py:attribute:: kelighg
   
   elastic coll. coeff:imp_gas+hyd_gas
   
   :Default: 
   :Dimension: (6)
   :Group: Comtra
   :Type: double
   :Unit: m**3/s
.. py:attribute:: keligii
   
   elastic coll. coeff:imp_gas+imp_ion
   
   :Default: 
   :Dimension: (6)
   :Group: Comtra
   :Type: double
   :Unit: m**3/s
.. py:attribute:: keligig
   
   elastic coll. coeff:imp_gas+imp_gas
   
   :Default: 
   :Dimension: (6)
   :Group: Comtra
   :Type: double
   :Unit: m**3/s
.. py:attribute:: diffusivity
   
   anomalous (turbulent) diffusivity (calculated during rhs eval)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: diffusivwrk
   
   anomalous (turbulent) diffusivity (mixed w/difni using cdifnit)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: diffusivloc
   
   anomalous (turbulent) diffusivity (local values for isturbcons=2)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: coll_fe
   
   nu_star_e/(1+nu_star_e) for elec CF drifts
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: coll_fi
   
   nu_star_i/(1+nu_star_i) for ion CF drifts
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalfea
   
   calc:elec thermal flux-limit array (see flalfe)
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalfia
   
   calc:ion thermal flux-limit array (see flalfi)
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalfva
   
   calc:ion visc flux-limit array (see flalfv)
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalfgxa
   
   calc:neut pol flux-limit array (see flalfgx)
   
   :Default: 
   :Dimension: (0:nx+1,10)
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalfgxya
   
   calc:neut xy flux-limit array (see flalfgxy)
   
   :Default: 
   :Dimension: (0:nx+1,10)
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalfgya
   
   calc:neut rad flux-limit array (see flalfgy)
   
   :Default: 
   :Dimension: (0:ny+1,10)
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalfvgxa
   
   calc:neut part pol flux-limit array (see flalfgx)
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalfvgxya
   
   calc:neut part xy flux-limit array (see flalfgxy)
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalfvgya
   
   calc:neut part rad flux-limit array (see flalfgy)
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalftgxa
   
   calc:neut part pol flux-limit array (see flalfgx)
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalftgxya
   
   calc:neut part xy flux-limit array (see flalfgxy)
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: flalftgya
   
   calc:neut part rad flux-limit array (see flalfgy)
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Comtra
   :Type: double
   :Unit: 
.. py:attribute:: del_sp
   
   dN/dt for ion source term
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: 1/m**3s
.. py:attribute:: del_witot
   
   Total dWi/dt; not used
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: W/m**3
.. py:attribute:: del_wetot
   
   Total dWe/dt; not used
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: W/m**3
.. py:attribute:: del_dndt
   
   dN/dt deduced from data (input)
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: 1/m**3s
.. py:attribute:: del_deedt
   
   dNTe/dt deduced from data (input)
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: W/m**3s
.. py:attribute:: del_deidt
   
   dNTi/dt deduced from data (input)
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: W/m**3s
.. py:attribute:: del_wicd
   
   dWi/dt for ions rad diff heat flux
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: W/m**3
.. py:attribute:: del_wicv
   
   dWi/dt for ions rad conv heat flux
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: W/m**3
.. py:attribute:: del_wecd
   
   dWe/dt for elec rad diff heat flux
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: W/m**3
.. py:attribute:: del_wecv
   
   dWe/dt for elec rad conv heat flux
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: W/m**3
.. py:attribute:: del_wicdd
   
   diag: dWi/dt ion rad diff heat flux
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: W/m**3
.. py:attribute:: del_wecdd
   
   diag: dWe/dt elec rad diff heat flux
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: W/m**3
.. py:attribute:: del_cei
   
   dWe/dt=-dWi/dt elec-ion coll exchange
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: W/m**3
.. py:attribute:: del_wivdp
   
   dWi/dt p-v type pressure-work terms
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: W/m**3
.. py:attribute:: del_wevdp
   
   dWe/dt p-v type pressure-work terms
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: W/m**3
.. py:attribute:: dif_int
   
   interp particle diff, D
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kyi_int
   
   interp ion radial conduc., Chi_i
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kye_int
   
   interp elec radial conduc., Chi_e
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: m**2/s
.. py:attribute:: vyn_int
   
   interp ion radial particle drift vel
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: m/s
.. py:attribute:: vyei_int
   
   interp ion radial energy drift vel
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: m/s
.. py:attribute:: vyee_int
   
   interp elec radial energy drift vel
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: m/s
.. py:attribute:: gamp
   
   radial ion particle flux
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: gamei
   
   ion radial energy flux
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: W/m**2s
.. py:attribute:: gamee
   
   elec radial energy flux
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: W/m**2s
.. py:attribute:: pfmpg
   
   radial neutral particle flux
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: facgam
   
   geom factor for flux-surf averaging
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: 
.. py:attribute:: floyd
   
   equiv radial ion particle current
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Interprettrans
   :Type: double
   :Unit: 1/s
.. py:attribute:: diffuswgts
   
   weights for radial digital filter
   
   :Default: 
   :Dimension: (-9:9)
   :Group: Turbulence
   :Type: double
   :Unit: none
.. py:attribute:: chinorml
   
   norm. anom. diffusivity (L-mode turbulence)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Turbulence_diagnostics
   :Type: double
   :Unit: 
.. py:attribute:: chinormh
   
   norm. anom. diffusivity (H-mode turbulence)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Turbulence_diagnostics
   :Type: double
   :Unit: 
.. py:attribute:: nist1
   
   density for ODE output times
   
   :Default: 
   :Dimension: (nsteps,0:nx+1,0:ny+1,nisp)
   :Group: Timary
   :Type: double
   :Unit: 1/m**3
.. py:attribute:: upst1
   
   parallel vel for ODE output times
   
   :Default: 
   :Dimension: (nsteps,0:nx+1,0:ny+1,nisp)
   :Group: Timary
   :Type: double
   :Unit: m/s
.. py:attribute:: test1
   
   Te at for ODE output times
   
   :Default: 
   :Dimension: (nsteps,0:nx+1,0:ny+1)
   :Group: Timary
   :Type: double
   :Unit: eV
.. py:attribute:: tist1
   
   Ti at for ODE output times
   
   :Default: 
   :Dimension: (nsteps,0:nx+1,0:ny+1)
   :Group: Timary
   :Type: double
   :Unit: eV
.. py:attribute:: ngst1
   
   ng at for ODE output times
   
   :Default: 
   :Dimension: (nsteps,0:nx+1,0:ny+1,ngsp)
   :Group: Timary
   :Type: double
   :Unit: 1/m**3
.. py:attribute:: phist1
   
   phi at various times
   
   :Default: 
   :Dimension: (nsteps,0:nx+1,0:ny+1)
   :Group: Timary
   :Type: double
   :Unit: V
.. py:attribute:: toutlsod
   
   time which nist1, etc, are filled
   
   :Default: 
   :Dimension: (nsteps)
   :Group: Timary
   :Type: double
   :Unit: s
.. py:attribute:: yldnmx
   
   max. rate-of-change, yldot/yl
   
   :Default: 
   :Dimension: (nsteps)
   :Group: Timary
   :Type: double
   :Unit: 1/s
.. py:attribute:: iyldnmx
   
   index of vector yldnmx
   
   :Default: 
   :Dimension: (nsteps)
   :Group: Timary
   :Type: integer
   :Unit: 
.. py:attribute:: ni_stor
   
   density for rundt output times
   
   :Default: 
   :Dimension: (n_stor,0:nx+1,0:ny+1,nisp)
   :Group: Timary
   :Type: double
   :Unit: 1/m**3
.. py:attribute:: up_stor
   
   par vel for rundt output times
   
   :Default: 
   :Dimension: (n_stor,0:nx+1,0:ny+1,nisp)
   :Group: Timary
   :Type: double
   :Unit: m/s
.. py:attribute:: te_stor
   
   Te for rundt output times
   
   :Default: 
   :Dimension: (n_stor,0:nx+1,0:ny+1)
   :Group: Timary
   :Type: double
   :Unit: eV
.. py:attribute:: ti_stor
   
   Ti for rundt output times
   
   :Default: 
   :Dimension: (n_stor,0:nx+1,0:ny+1)
   :Group: Timary
   :Type: double
   :Unit: eV
.. py:attribute:: ng_stor
   
   ng for rundt output times
   
   :Default: 
   :Dimension: (n_stor,0:nx+1,0:ny+1,ngsp)
   :Group: Timary
   :Type: double
   :Unit: 1/m**3
.. py:attribute:: phi_stor
   
   phi for rundt times
   
   :Default: 
   :Dimension: (n_stor,0:nx+1,0:ny+1)
   :Group: Timary
   :Type: double
   :Unit: V
.. py:attribute:: tim_stor
   
   output times for rundt
   
   :Default: 
   :Dimension: (n_stor)
   :Group: Timary
   :Type: double
   :Unit: t
.. py:attribute:: nfe_stor
   
   num func evals in rundt interval
   
   :Default: 
   :Dimension: (n_stor)
   :Group: Timary
   :Type: double
   :Unit: 
.. py:attribute:: dtreal_stor
   
   dtreal for rundt output times
   
   :Default: 
   :Dimension: (n_stor)
   :Group: Timary
   :Type: double
   :Unit: t
.. py:attribute:: savefname
   
   name of pfb save file pfdt_'savefname'
   
   :Default: 
   :Dimension: 
   :Group: Timary
   :Type: character(5)
   :Unit: 
.. py:attribute:: mi
   
   ion mass in kg, calculated from minu
   
   :Default: 
   :Dimension: (1:nisp)
   :Group: Compla
   :Type: double
   :Unit: kg
.. py:attribute:: zi
   
   ion charge number, calc. from ziin
   
   :Default: 
   :Dimension: (1:nisp)
   :Group: Compla
   :Type: double
   :Unit: 
.. py:attribute:: mg
   
   gas species mass, calc. fr minu
   
   :Default: 
   :Dimension: (1:ngsp)
   :Group: Compla
   :Type: double
   :Unit: kg
.. py:attribute:: facmg
   
   scale factor for mg to recov old case
   
   :Default: 
   :Dimension: (1:31)
   :Group: Compla
   :Type: double
   :Unit: 
.. py:attribute:: znucl
   
   tot. nucl. charge, calc. from znuclin
   
   :Default: 
   :Dimension: (1:nisp)
   :Group: Compla
   :Type: integer
   :Unit: 
.. py:attribute:: ni
   
   ion density in primary cell (ix,iy)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m^-3
.. py:attribute:: lni
   
   log(ion dens) in prim. cell (ix,iy)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m^-3
.. py:attribute:: nm
   
   mass density [nm(,,1) is sum, exclud.
   gas, if nusp=1, isimpon=5] in cell
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: kg*m^-3
.. py:attribute:: nz2
   
   sum of ni*zi**2 over all ion species
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: m^-3
.. py:attribute:: uu
   
   ratio ion-flux/density at x-face;
   if orthog mesh, poloidal ion velocity
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: uup
   
   poloidal ion vel (|| flow contrib)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: up
   
   par ion vel if full mom eqn on
   (mass-dens. avg if isimpon = 5)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: upi
   
   inter. par ion vel even if force bal
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: upifmb
   
   par ion vel fmombal if isimpon=5
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: uz
   
   toroidal ion vel in pol X rad direct
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: v2
   
   vel normal to parallel & rad. direc.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: v2xgp
   
   v2 ion vel for v2x_gradx_P eng terms
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: v2ce
   
   portion of v2 from ExB
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: v2cb
   
   portion of v2 from grad_B
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: ve2cb
   
   electron v2 from grad_B
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: v2cd
   
   portion of ion v2 from grad_PxB
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: ve2cd
   
   portion of elec v2 from grad_PxB
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: q2cd
   
   ion heat flux from grad_PxB
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: v2rd
   
   portion of v2 from resistive drift
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: v2dd
   
   portion of v2 from anomalous drift
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: vy
   
   radial ion velocity
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: vygp
   
   radial ion vel for vy_grady_P eng terms
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: vytan
   
   radial ion vel.*tan(vtag) on x-face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: vygtan
   
   radial gas grad-T vel.*tan(vtag) on
   x-face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: vyce
   
   portion of vy from ExB
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: vycb
   
   portion of vy from grad_B
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: veycb
   
   electron vy from grad_B
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: vycp
   
   ion vy from grad_PixB
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: veycp
   
   electron vy from grad_PeXB
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: vyrd
   
   portion of vy from resistive drift
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: vydd
   
   portion of vy from anomalous drift
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: vyavis
   
   rad vel from anom perp vis (ExB,P)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: vex
   
   Poloidal electron velocity
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: upe
   
   parallel electron velocity
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: vep
   
   old parallel electron velocity-remove
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: ve2
   
   old '2' electron velocity-remove
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: vey
   
   Radial electron velocity
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: vycf
   
   radial vel from class. viscosity
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: vycr
   
   radial vel from class. thermal force
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: te
   
   electron temperature in primary cell
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: J
.. py:attribute:: ti
   
   ion temperature in primary cell
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: J
.. py:attribute:: ng
   
   gas density in primary cell (ix,iy)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Compla
   :Type: double
   :Unit: m^-3
.. py:attribute:: lng
   
   log(gas dens) in prim. cell (ix,iy)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Compla
   :Type: double
   :Unit: m^-3
.. py:attribute:: uug
   
   ratio gas-flux/density at x-face;
   if orthog mesh, poloidal gas velocity
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: vyg
   
   radial gas velocity
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: tg
   
   gas temperature in primary cell
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Compla
   :Type: double
   :Unit: J
.. py:attribute:: istgcon
   
   =0, set tg(,,i)=rtg2ti*ti; if >0, set
   tg=(1-istgcon)*rtg2ti*ti+istgcon*tgas*ev
   
   :Default: 
   :Dimension: (6)
   :Group: Compla
   :Type: double
   :Unit: 
.. py:attribute:: tev
   
   ion temperature at vertex of cell
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: J
.. py:attribute:: niv
   
   ion dens up-right vert[rm,zm(,,4)]
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m^-3
.. py:attribute:: upv
   
   ion par vel up-right vert[rm,zm(,,4)]
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: ngv
   
   gas dens up-right vert[rm,zm(,,4)]
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Compla
   :Type: double
   :Unit: m^-3
.. py:attribute:: tiv
   
   ion temperature at vertex of cell
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: J
.. py:attribute:: niy0
   
   ion density below y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m^-3
.. py:attribute:: niy1
   
   ion density above y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m^-3
.. py:attribute:: niy0s
   
   old ion density below y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m^-3
.. py:attribute:: niy1s
   
   old ion density above y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m^-3
.. py:attribute:: ney0
   
   elec density below y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: m^-3
.. py:attribute:: ney1
   
   elec density above y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: m^-3
.. py:attribute:: nity0
   
   total ion density below y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: m^-3
.. py:attribute:: nity1
   
   total ion density above y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: m^-3
.. py:attribute:: tey0
   
   elec temp below y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: eV
.. py:attribute:: tey1
   
   elec temp above y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: eV
.. py:attribute:: tiy0
   
   ion temp below y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: eV
.. py:attribute:: tiy1
   
   ion temp above y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: eV
.. py:attribute:: tiy0s
   
   old ion temp below y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: eV
.. py:attribute:: tiy1s
   
   old ion temp above y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: eV
.. py:attribute:: tgy0
   
   atom temp below y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Compla
   :Type: double
   :Unit: eV
.. py:attribute:: tgy1
   
   atom temp above y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Compla
   :Type: double
   :Unit: eV
.. py:attribute:: ngy0
   
   gas density below y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Compla
   :Type: double
   :Unit: m^-3
.. py:attribute:: ngy1
   
   gas density above y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Compla
   :Type: double
   :Unit: m^-3
.. py:attribute:: pgy0
   
   gas pressure below y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Compla
   :Type: double
   :Unit: J/m^3
.. py:attribute:: pgy1
   
   gas pressure above y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Compla
   :Type: double
   :Unit: J/m^3
.. py:attribute:: pg
   
   gas pressure at cell center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Compla
   :Type: double
   :Unit: J/m^3
.. py:attribute:: phiy0
   
   potential below y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: V
.. py:attribute:: phiy1
   
   potential above y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: V
.. py:attribute:: phiy0s
   
   old potential below y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: V
.. py:attribute:: phiy1s
   
   old potential above y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: V
.. py:attribute:: pr
   
   total pressure at center of cell
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: J/m^3
.. py:attribute:: prev
   
   elec pressure at vertex of cell
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: J/m^3
.. py:attribute:: prtv
   
   total pressure at vertex of cell
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: J/m^3
.. py:attribute:: pri
   
   ion plasma pressure
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: J/m^3
.. py:attribute:: priv
   
   ion pressure at vertex of cells
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: J/m^3
.. py:attribute:: priy0
   
   ion pressure below y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: J/m^3
.. py:attribute:: priy1
   
   ion pressure above y-face center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: J/m^3
.. py:attribute:: pre
   
   el. plasma pressure
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: J/m^3
.. py:attribute:: ne
   
   electron dens in primary cell (ix,iy)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: m^-3
.. py:attribute:: nit
   
   tot ion dens in primary cell (ix,iy)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: m^-3
.. py:attribute:: nginit
   
   init gas dens in primary cell (ix,iy)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: m^-3
.. py:attribute:: phi
   
   potential in primary cell (ix,iy)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: V
.. py:attribute:: phiv
   
   potential at vertex of cell
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: V
.. py:attribute:: zeff
   
   Z_effective charge in cell (ix,iy)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: 
.. py:attribute:: netap
   
   ne*parallel resistivity
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: 
.. py:attribute:: znot
   
   =Sum(n_z * Z^2)/n_i in cell
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: 
.. py:attribute:: zimpc
   
   Zimp (avg-ion model) in cell (ix,iy)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: 
.. py:attribute:: nil
   
   ion density at last output
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m^-3
.. py:attribute:: upl
   
   parallel ion velocity at last output
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: tel
   
   electron temperature at last output
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: J
.. py:attribute:: til
   
   ion temperature at last output
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: J
.. py:attribute:: ngl
   
   gas density at last output
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Compla
   :Type: double
   :Unit: m^-3
.. py:attribute:: phil
   
   potential at last output
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: V
.. py:attribute:: upxpt
   
   parallel velocity at x-point
   
   :Default: 
   :Dimension: (1:nusp,1:nxpt)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: nixpt
   
   ion density at x-point
   
   :Default: 
   :Dimension: (1:nusp,1:nxpt)
   :Group: Compla
   :Type: double
   :Unit: m^-3
.. py:attribute:: visyxpt
   
   ion viscosity at x-point
   
   :Default: 
   :Dimension: (1:nusp,1:nxpt)
   :Group: Compla
   :Type: double
   :Unit: 
.. py:attribute:: vyhxpt
   
   horiz. ion drift vel. at x-point
   
   :Default: 
   :Dimension: (1:nusp,1:nxpt)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: vyvxpt
   
   vert. ion drift vel. at x-point
   
   :Default: 
   :Dimension: (1:nusp,1:nxpt)
   :Group: Compla
   :Type: double
   :Unit: m/s
.. py:attribute:: fmihxpt
   
   horiz. mom. flux at x-point
   
   :Default: 
   :Dimension: (1:nusp,1:nxpt)
   :Group: Compla
   :Type: double
   :Unit: Nwt
.. py:attribute:: fmivxpt
   
   vert. mom. flux at x-point
   
   :Default: 
   :Dimension: (1:nusp,1:nxpt)
   :Group: Compla
   :Type: double
   :Unit: Nwt
.. py:attribute:: rtaux
   
   Norm. poloidal neutral line-dens.,
   Ly-a opacity to plates
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: 1e-16 m^-2
.. py:attribute:: rtauy
   
   Norm. radial neutral line-dens.,
   norm. Ly-a opacity to radial wall
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: 1e-16 m^-2
.. py:attribute:: rtau
   
   Min. norm neutral line-dens.,
   min. Ly-a opacity; min(rtaux,rtauy)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: 1e-16 m^-2
.. py:attribute:: betap
   
   poloidal plasma beta
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Compla
   :Type: double
   :Unit: 
.. py:attribute:: fetx
   
   total energy flow through a poloidal cell face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Postproc
   :Type: double
   :Unit: W
.. py:attribute:: fety
   
   total energy flow through a radial cell face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Postproc
   :Type: double
   :Unit: W
.. py:attribute:: pdrift
   
   power in bringing new ion to flow velocity
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Postproc
   :Type: double
   :Unit: W/m^3
.. py:attribute:: peirad
   
   tot. power lost by electrons and ions in
   rad., ion. and dissoc.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Postproc
   :Type: double
   :Unit: W
.. py:attribute:: png2ni
   
   power exchange bwt. neutral and ion
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Postproc
   :Type: double
   :Unit: W
.. py:attribute:: pmomv
   
   power exchange bwt. neutal & ion from flow
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Postproc
   :Type: double
   :Unit: W/m^3
.. py:attribute:: jdote
   
   power from J.E heating
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Postproc
   :Type: double
   :Unit: W
.. py:attribute:: engerr
   
   local error in power balance
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Postproc
   :Type: double
   :Unit: 
.. py:attribute:: iion
   
   net ionization current per gas isotope
   
   :Default: 
   :Dimension: (ngsp)
   :Group: Postproc
   :Type: double
   :Unit: A
.. py:attribute:: irecomb
   
   net recombination current of gas isotope
   
   :Default: 
   :Dimension: (ngsp)
   :Group: Postproc
   :Type: double
   :Unit: A
.. py:attribute:: icxgas
   
   net charge exchange current of gas isotope
   
   :Default: 
   :Dimension: (ngsp)
   :Group: Postproc
   :Type: double
   :Unit: A
.. py:attribute:: pradimpt
   
   net impurity photon rad. loss
   
   :Default: 
   :Dimension: (ngsp)
   :Group: Postproc
   :Type: double
   :Unit: W
.. py:attribute:: pradimp
   
   rad. loss for each impurity charge state
   
   :Default: 
   :Dimension: (0:nzspmx,ngsp-1)
   :Group: Postproc
   :Type: double
   :Unit: W
.. py:attribute:: pwr_plth
   
   hydrog rad pwr flux on divertor plate
   
   :Default: 
   :Dimension: (0:ny+1,2*nxpt)
   :Group: Postproc
   :Type: double
   :Unit: W/m**2
.. py:attribute:: pwr_pltz
   
   impur rad pwr flux on divertor plate
   
   :Default: 
   :Dimension: (0:ny+1,2*nxpt)
   :Group: Postproc
   :Type: double
   :Unit: W/m**2
.. py:attribute:: pwr_wallh
   
   hydrog rad pwr flux on outer wall
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Postproc
   :Type: double
   :Unit: W/m**2
.. py:attribute:: pwr_wallz
   
   impur rad pwr flux on outer wall
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Postproc
   :Type: double
   :Unit: W/m**2
.. py:attribute:: pwr_pfwallh
   
   hydrog rad pwr flux on PF wall
   
   :Default: 
   :Dimension: (0:nx+1,nxpt)
   :Group: Postproc
   :Type: double
   :Unit: W/m**2
.. py:attribute:: pwr_pfwallz
   
   impur rad pwr flux on PF wall
   
   :Default: 
   :Dimension: (0:nx+1,nxpt)
   :Group: Postproc
   :Type: double
   :Unit: W/m**2
.. py:attribute:: sdelb
   
   elec pwr flux to left div
   
   :Default: 
   :Dimension: (0:ny+1,1:nxpt)
   :Group: Postproc
   :Type: double
   :Unit: W/m**2
.. py:attribute:: sderb
   
   elec ion pwr flux to right div
   
   :Default: 
   :Dimension: (0:ny+1,1:nxpt)
   :Group: Postproc
   :Type: double
   :Unit: W/m**2
.. py:attribute:: sdilb
   
   tot ion pwr flux to left div
   
   :Default: 
   :Dimension: (0:ny+1,1:nxpt)
   :Group: Postproc
   :Type: double
   :Unit: W/m**2
.. py:attribute:: sdirb
   
   tot ion pwr flux to right div
   
   :Default: 
   :Dimension: (0:ny+1,1:nxpt)
   :Group: Postproc
   :Type: double
   :Unit: W/m**2
.. py:attribute:: sbindlb
   
   tot bind eng pwr flux to left div
   
   :Default: 
   :Dimension: (0:ny+1,1:nxpt)
   :Group: Postproc
   :Type: double
   :Unit: W/m**2
.. py:attribute:: sbindrb
   
   tot bind eng pwr flux to right div
   
   :Default: 
   :Dimension: (0:ny+1,1:nxpt)
   :Group: Postproc
   :Type: double
   :Unit: W/m**2
.. py:attribute:: sdrlb
   
   tot rad pwr flux to left div
   
   :Default: 
   :Dimension: (0:ny+1,1:nxpt)
   :Group: Postproc
   :Type: double
   :Unit: W/m**2
.. py:attribute:: sdrrb
   
   tot rad flux to right div
   
   :Default: 
   :Dimension: (0:ny+1,1:nxpt)
   :Group: Postproc
   :Type: double
   :Unit: W/m**2
.. py:attribute:: sdtlb
   
   tot pwr flux to left div
   
   :Default: 
   :Dimension: (0:ny+1,1:nxpt)
   :Group: Postproc
   :Type: double
   :Unit: W/m**2
.. py:attribute:: sdtrb
   
   tot pwr flux to right div
   
   :Default: 
   :Dimension: (0:ny+1,1:nxpt)
   :Group: Postproc
   :Type: double
   :Unit: W/m**2
.. py:attribute:: gdilb
   
   particle flux to left div
   
   :Default: 
   :Dimension: (0:ny+1,nisp,nxpt)
   :Group: Postproc
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: gdirb
   
   particle flux to right div
   
   :Default: 
   :Dimension: (0:ny+1,nisp,nxpt)
   :Group: Postproc
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: engilb
   
   ave ion energy to left div
   
   :Default: 
   :Dimension: (0:ny+1,nisp,nxpt)
   :Group: Postproc
   :Type: double
   :Unit: Volts
.. py:attribute:: engirb
   
   ave ion energy to right div
   
   :Default: 
   :Dimension: (0:ny+1,nisp,nxpt)
   :Group: Postproc
   :Type: double
   :Unit: Volts
.. py:attribute:: gwalli
   
   particle flux to left div
   
   :Default: 
   :Dimension: (0:nx+1,nisp)
   :Group: Postproc
   :Type: double
   :Unit: 1/m**2s
.. py:attribute:: engwalli
   
   ave ion energy to outer wall
   
   :Default: 
   :Dimension: (0:nx+1,nisp)
   :Group: Postproc
   :Type: double
   :Unit: Volts
.. py:attribute:: swallr
   
   radiation pwr flux to outer wall
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Postproc
   :Type: double
   :Unit: W/m**2
.. py:attribute:: swalli
   
   ion pwr flux to outer wall
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Postproc
   :Type: double
   :Unit: W/m**2
.. py:attribute:: swalle
   
   elec pwr flux to outer wall
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Postproc
   :Type: double
   :Unit: W/m**2
.. py:attribute:: swbind
   
   binding energy flux to outer wall
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Postproc
   :Type: double
   :Unit: W/m**2
.. py:attribute:: swallt
   
   total pwr flux to outer wall
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Postproc
   :Type: double
   :Unit: W/m**2
.. py:attribute:: spfwallr
   
   radiation pwr flux to PF wall
   
   :Default: 
   :Dimension: (0:nx+1,nxpt)
   :Group: Postproc
   :Type: double
   :Unit: W/m**2
.. py:attribute:: pwrsore
   
   power src into electrons in cell ix,iy
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Volsrc
   :Type: double
   :Unit: W
.. py:attribute:: pwrsori
   
   power src into ions in cell ix,iy
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Volsrc
   :Type: double
   :Unit: W
.. py:attribute:: volpsor
   
   current src into ions in cell ix,iy
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Volsrc
   :Type: double
   :Unit: 1/s
.. py:attribute:: volmsor
   
   up mom src in cell ix,iy
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Volsrc
   :Type: double
   :Unit: kg m/s**2
.. py:attribute:: voljcsor
   
   uniform core-region curr sor. in ix,iy
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Volsrc
   :Type: double
   :Unit: A
.. py:attribute:: volpsorg
   
   curr source for gas in cell ix,iy
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Volsrc
   :Type: double
   :Unit: 1/s
.. py:attribute:: psgov_use
   
   user-specified gas source
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Volsrc
   :Type: double
   :Unit: 1/m**3 s
.. py:attribute:: ivolcur
   
   total volume current
   
   :Default: 
   :Dimension: (1:nisp)
   :Group: Volsrc
   :Type: double
   :Unit: A
.. py:attribute:: mvolcur
   
   total volume parallel mom. curr.
   
   :Default: 
   :Dimension: (1:nisp)
   :Group: Volsrc
   :Type: double
   :Unit: kgA m/s
.. py:attribute:: effvng
   
   normalizing factor of gas source; calc
   
   :Default: 
   :Dimension: (1:ngsp)
   :Group: Volsrc
   :Type: double
   :Unit: m
.. py:attribute:: ivolcurg
   
   tot. volumn gas source strength
   
   :Default: 
   :Dimension: (1:ngsp)
   :Group: Volsrc
   :Type: double
   :Unit: A
.. py:attribute:: z0ng
   
   axial or x loc. of gas particle profile
   
   :Default: 
   :Dimension: (1:ngsp)
   :Group: Volsrc
   :Type: double
   :Unit: m
.. py:attribute:: r0ng
   
   rad. or y loc. of gas particle profile
   
   :Default: 
   :Dimension: (1:ngsp)
   :Group: Volsrc
   :Type: double
   :Unit: m
.. py:attribute:: zwng
   
   axial or y Gaussian 1/2 width gas prtcl
   
   :Default: 
   :Dimension: (1:ngsp)
   :Group: Volsrc
   :Type: double
   :Unit: m
.. py:attribute:: rwng
   
   rad. or y Gaussian 1/2 width gas prtcl
   
   :Default: 
   :Dimension: (1:ngsp)
   :Group: Volsrc
   :Type: double
   :Unit: m
.. py:attribute:: b02d
   
   net B-field scale fac. =b0+b0_use present iter.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Bfield
   :Type: double
   :Unit: 
.. py:attribute:: b0old
   
   net B-field scale factor at last iter.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Bfield
   :Type: double
   :Unit: 
.. py:attribute:: b0_use
   
   spatial B-field scale factor; user input
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Bfield
   :Type: double
   :Unit: 
.. py:attribute:: rbpol
   
   major radius*poloidal magnetic field
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Bfield
   :Type: double
   :Unit: mT
.. py:attribute:: btot
   
   total magnetic field strength
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Bfield
   :Type: double
   :Unit: T
.. py:attribute:: rbfbt
   
   ratio bphi/btot at density cell center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Bfield
   :Type: double
   :Unit: 
.. py:attribute:: rbfbt2
   
   ratio bphi/btot**2 at density cell center
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Bfield
   :Type: double
   :Unit: 1/T
.. py:attribute:: curvrby
   
   curvature drift factor on y-face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Bfield
   :Type: double
   :Unit: 1/mT
.. py:attribute:: curvrb2
   
   curvature drift factor on x-face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Bfield
   :Type: double
   :Unit: 1/mT
.. py:attribute:: gradby
   
   grad_B drift (p_perp) drift factor, y-face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Bfield
   :Type: double
   :Unit: 1/mT
.. py:attribute:: gradb2
   
   grad_B drift (p_perp) drift factor, x-face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Bfield
   :Type: double
   :Unit: 1/mT
.. py:attribute:: dbm2dx
   
   pol deriv of 1/B**2 on x-face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Bfield
   :Type: double
   :Unit: 1/T**2m
.. py:attribute:: dbm2dy
   
   rad deriv of 1/B**2 on y-face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Bfield
   :Type: double
   :Unit: 1/T**2m
.. py:attribute:: bfacxrozh
   
   [1-B**2/B**2_ave], x-face; Rozhansky
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Bfield
   :Type: double
   :Unit: 
.. py:attribute:: bfacyrozh
   
   [1-B**2/B**2_ave], y-face; Rozhansky
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Bfield
   :Type: double
   :Unit: 
.. py:attribute:: ni0
   
   old ion density
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Oldpla
   :Type: double
   :Unit: m^-3
.. py:attribute:: ng0
   
   old neutral density
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Oldpla
   :Type: double
   :Unit: m^-3
.. py:attribute:: te0
   
   old electron temperature
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Oldpla
   :Type: double
   :Unit: J
.. py:attribute:: ti0
   
   old ion temperature
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Oldpla
   :Type: double
   :Unit: J
.. py:attribute:: phi0
   
   old electrostatic potential
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Oldpla
   :Type: double
   :Unit: V
.. py:attribute:: up0
   
   old parallel velocity
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Oldpla
   :Type: double
   :Unit: m/s
.. py:attribute:: vy0
   
   old radial velocity
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Oldpla
   :Type: double
   :Unit: m/s
.. py:attribute:: fqp
   
   pol proj of par cur, east face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: Amp
.. py:attribute:: fq2
   
   pol proj of 2 cur, east face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: Amp
.. py:attribute:: fqx
   
   net poloidal current, east face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: Amp
.. py:attribute:: fqxb
   
   poloidal cur from grad_B, east face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: Amp
.. py:attribute:: fdiaxlb
   
   left boundary Dia current for bc
   
   :Default: 
   :Dimension: (0:ny+1,1:nxpt)
   :Group: Comflo
   :Type: double
   :Unit: Amp
.. py:attribute:: fdiaxrb
   
   right boundary Dia current for bc
   
   :Default: 
   :Dimension: (0:ny+1,1:nxpt)
   :Group: Comflo
   :Type: double
   :Unit: Amp
.. py:attribute:: floxebgt
   
   BxgradTe diamag part floxe (-> feex)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: W
.. py:attribute:: floxibgt
   
   BxgradTi diamag part floxi (-> feex)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Comflo
   :Type: double
   :Unit: W
.. py:attribute:: fqy
   
   net radial current, north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: Amp
.. py:attribute:: fqyb
   
   radial current from grad_B, north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: Amp
.. py:attribute:: fqyn
   
   radial cur from cx coll, north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: Amp
.. py:attribute:: fqym
   
   radial cur from inertia, north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: Amp
.. py:attribute:: fqymi
   
   spec rad cur from inertia, north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Comflo
   :Type: double
   :Unit: Amp
.. py:attribute:: fqya
   
   anomalous visc rad cur, north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: Amp
.. py:attribute:: fqydt
   
   time-dep inertial rad cur, north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: Amp
.. py:attribute:: fqydti
   
   spec time-dep inert rad cur, north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Comflo
   :Type: double
   :Unit: Amp
.. py:attribute:: fqyao
   
   old anom mobil rad current, north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: Amp
.. py:attribute:: fqyae
   
   anom mobil rad current for electrons, north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: Amp
.. py:attribute:: fqyai
   
   anom mobil rad current for ions, north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: Amp
.. py:attribute:: fqyd
   
   diamag radial current; north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: Amp
.. py:attribute:: fqygp
   
   net radial curr. uses grad_P, north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: Amp
.. py:attribute:: fq2d
   
   diamag 2-current; east face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: Amp
.. py:attribute:: fqypneo
   
   rad-cur from neo particle flux
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: Amp
.. py:attribute:: fq2pneo
   
   2-cur from neo particle flux
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: Amp
.. py:attribute:: fqyqneo
   
   rad-cur from neo heat flux
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: Amp
.. py:attribute:: fq2qneo
   
   2-cur from neo heat flux
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: Amp
.. py:attribute:: fnix
   
   ion poloidal current, east face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Comflo
   :Type: double
   :Unit: 1/s
.. py:attribute:: fnixcb
   
   ion grad-B pol. current, east face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Comflo
   :Type: double
   :Unit: 1/s
.. py:attribute:: fniy
   
   ion radial current, north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Comflo
   :Type: double
   :Unit: 1/s
.. py:attribute:: fniy4ord
   
   4th ord ion radial current, north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Comflo
   :Type: double
   :Unit: 1/s
.. py:attribute:: fniycb
   
   ion grad-B rad. current, north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Comflo
   :Type: double
   :Unit: 1/s
.. py:attribute:: flnix
   
   ion poloidal log-current, east face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Comflo
   :Type: double
   :Unit: 1/s
.. py:attribute:: flniy
   
   ion radial log-current, north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Comflo
   :Type: double
   :Unit: 1/s
.. py:attribute:: fmix
   
   ion poloidal momentum current,east face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nusp)
   :Group: Comflo
   :Type: double
   :Unit: Nwt
.. py:attribute:: fmiy
   
   ion radial momentum current, north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nusp)
   :Group: Comflo
   :Type: double
   :Unit: Nwt
.. py:attribute:: fmixy
   
   nonorthog ion pol. mom. curr., east f.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nusp)
   :Group: Comflo
   :Type: double
   :Unit: Nwt
.. py:attribute:: fmity
   
   rad flux of cross-field tor. mom*R/Bp;
   nisp dimen, not nusp as for pot eqn
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Comflo
   :Type: double
   :Unit: 
.. py:attribute:: fmgx
   
   pol. neutral mom. current, east face ### IJ 2016/10/11
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,ngsp)
   :Group: Comflo
   :Type: double
   :Unit: Nwt
.. py:attribute:: fmgy
   
   rad. neutral mom. current, north face ### IJ 2016/10/11
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,ngsp)
   :Group: Comflo
   :Type: double
   :Unit: Nwt
.. py:attribute:: feex
   
   poloidal electron thermal current,
   east face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: J/s
.. py:attribute:: feey
   
   radial electron thermal current,
   north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: J/s
.. py:attribute:: feexy
   
   nonorthog elec. pol. therm cur, east f.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: J/s
.. py:attribute:: feey4ord
   
   elec. pol. kye4order therm cur, east f.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: J/s
.. py:attribute:: feix
   
   poloidal ion thermal current, east face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: J/s
.. py:attribute:: feiy
   
   radial ion thermal current, north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: J/s
.. py:attribute:: fegx
   
   poloidal neutral thermal current, east face ### IJ 2016/09/22
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,ngsp)
   :Group: Comflo
   :Type: double
   :Unit: J/s
.. py:attribute:: fegy
   
   radial neutral thermal current, north face ### IJ 2016/09/22
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,ngsp)
   :Group: Comflo
   :Type: double
   :Unit: J/s
.. py:attribute:: qipar
   
   parallel conductive ion heat flux
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nisp)
   :Group: Comflo
   :Type: double
   :Unit: J/m**2s
.. py:attribute:: qgpar
   
   parallel conductive gas heat flux
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,ngsp)
   :Group: Comflo
   :Type: double
   :Unit: J/m**2s
.. py:attribute:: fniycbo
   
   fniy cor. iy=0 bdry for grad_B, grad_P
   
   :Default: 
   :Dimension: (0:nx+1,1:nisp)
   :Group: Comflo
   :Type: double
   :Unit: 1/s
.. py:attribute:: feiycbo
   
   feiy cor. iy=0 bdry for grad_B, grad_P
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Comflo
   :Type: double
   :Unit: J/s
.. py:attribute:: feeycbo
   
   feey cor. iy=0 bdry for grad_B, grad_P
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Comflo
   :Type: double
   :Unit: J/s
.. py:attribute:: feixy
   
   nonorthog ion pol. thermal cur, east f.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: J/s
.. py:attribute:: feiy4ord
   
   ion pol. kyi4order therm cur, east f.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: J/s
.. py:attribute:: fngx
   
   neutral polodial current, east face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Comflo
   :Type: double
   :Unit: 1/s
.. py:attribute:: fngx4ord
   
   4th ord gas radial current, north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Comflo
   :Type: double
   :Unit: 1/s
.. py:attribute:: flngx
   
   neutral pol. log-current, east face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Comflo
   :Type: double
   :Unit: 1/s
.. py:attribute:: fngxs
   
   neutral pol cur w/o fngxy, east face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Comflo
   :Type: double
   :Unit: 1/s
.. py:attribute:: fngy
   
   neutral radial current, north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Comflo
   :Type: double
   :Unit: 1/s
.. py:attribute:: fngy4ord
   
   4th ord gas radial current, north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Comflo
   :Type: double
   :Unit: 1/s
.. py:attribute:: flngy
   
   neutral radial log-current, north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Comflo
   :Type: double
   :Unit: 1/s
.. py:attribute:: fngxy
   
   nonorthog gas pol. cur., east face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Comflo
   :Type: double
   :Unit: 1/s
.. py:attribute:: flngxy
   
   nonorthog gas pol.log-cur., east face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Comflo
   :Type: double
   :Unit: 1/s
.. py:attribute:: fngyx
   
   nonorthog gas rad. cur., north face
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Comflo
   :Type: double
   :Unit: 1/s
.. py:attribute:: fnixtot
   
   total poloidal ion cur.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: 1/s
.. py:attribute:: fniytot
   
   total radial ion cur.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Comflo
   :Type: double
   :Unit: 1/s
.. py:attribute:: idxn
   
   index of yl vector for ni(ix,iy)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Indexes
   :Type: integer
   :Unit: 
.. py:attribute:: idxg
   
   index of yl vector for ng(ix,iy,ig)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Indexes
   :Type: integer
   :Unit: 
.. py:attribute:: idxtg
   
   index of yl vector for tg(ix,iy,ig)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Indexes
   :Type: integer
   :Unit: 
.. py:attribute:: idxu
   
   index of yl vector for up(ix,iy)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nusp)
   :Group: Indexes
   :Type: integer
   :Unit: 
.. py:attribute:: idxti
   
   index of yl vector for ti(ix,iy)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Indexes
   :Type: integer
   :Unit: 
.. py:attribute:: idxte
   
   index of yl vector for te(ix,iy)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Indexes
   :Type: integer
   :Unit: 
.. py:attribute:: idxphi
   
   index of yl vector for phi(ix,iy)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Indexes
   :Type: integer
   :Unit: 
.. py:attribute:: ivfirst
   
   first eqn number iv at (ix,iy)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Indexes
   :Type: integer
   :Unit: 
.. py:attribute:: igyl
   
   ix,iy indices for vector yl(iv)
   
   :Default: 
   :Dimension: (neqmx,2)
   :Group: Indexes
   :Type: integer
   :Unit: 
.. py:attribute:: iseqalg
   
   flag(=1) for eqn being algebraic
   
   :Default: 
   :Dimension: (neqmx)
   :Group: Indexes
   :Type: integer
   :Unit: 
.. py:attribute:: isvarup
   
   flag(=1) for variable being up
   
   :Default: 
   :Dimension: (numvar)
   :Group: Indexes
   :Type: integer
   :Unit: 
.. py:attribute:: isvarphi
   
   flag(=1) for variable being phi
   
   :Default: 
   :Dimension: (numvar)
   :Group: Indexes
   :Type: integer
   :Unit: 
.. py:attribute:: hu
   
   present timestep
   
   :Default: 
   :Dimension: (nsteps,ngrid)
   :Group: Stat
   :Type: double
   :Unit: 1/s
.. py:attribute:: gpe
   
   ratio, linear to nonlinear iter, nli/nni
   
   :Default: 
   :Dimension: (nsteps,ngrid)
   :Group: Stat
   :Type: double
   :Unit: 
.. py:attribute:: npe
   
   cumulative number of precond. eval.
   
   :Default: 
   :Dimension: (nsteps,ngrid)
   :Group: Stat
   :Type: integer
   :Unit: 
.. py:attribute:: nps
   
   cumulative number of precond. solves
   
   :Default: 
   :Dimension: (nsteps,ngrid)
   :Group: Stat
   :Type: integer
   :Unit: 
.. py:attribute:: nfe
   
   cumulative number of RHS eval.
   
   :Default: 
   :Dimension: (nsteps,ngrid)
   :Group: Stat
   :Type: integer
   :Unit: 
.. py:attribute:: nst
   
   cum. number of steps taken
   
   :Default: 
   :Dimension: (nsteps,ngrid)
   :Group: Stat
   :Type: integer
   :Unit: 
.. py:attribute:: nni
   
   cumulative number of nonlinear iter.
   
   :Default: 
   :Dimension: (nsteps,ngrid)
   :Group: Stat
   :Type: integer
   :Unit: 
.. py:attribute:: nli
   
   cumulative number of linear iter.
   
   :Default: 
   :Dimension: (nsteps,ngrid)
   :Group: Stat
   :Type: integer
   :Unit: 
.. py:attribute:: nje
   
   cumulative number of Jacobian eval.
   
   :Default: 
   :Dimension: (nsteps,ngrid)
   :Group: Stat
   :Type: integer
   :Unit: 
.. py:attribute:: ncfn
   
   number of nonlinear converg. failures
   
   :Default: 
   :Dimension: (nsteps,ngrid)
   :Group: Stat
   :Type: integer
   :Unit: 
.. py:attribute:: ncfl
   
   number of linear converg. failures
   
   :Default: 
   :Dimension: (nsteps,ngrid)
   :Group: Stat
   :Type: integer
   :Unit: 
.. py:attribute:: nqu
   
   method order last used
   
   :Default: 
   :Dimension: (nsteps,ngrid)
   :Group: Stat
   :Type: integer
   :Unit: 
.. py:attribute:: iddas
   
   idid for daspk
   
   :Default: 
   :Dimension: (nsteps,ngrid)
   :Group: Stat
   :Type: integer
   :Unit: 
.. py:attribute:: eqmxer
   
   eqn. number giving maximum error in lsode
   
   :Default: 
   :Dimension: (nsteps,ngrid)
   :Group: Stat
   :Type: integer
   :Unit: 
.. py:attribute:: npsn
   
   cum. number of Jacobian solves for Newton
   
   :Default: 
   :Dimension: (ngrid)
   :Group: Stat
   :Type: integer
   :Unit: 
.. py:attribute:: njen
   
   cum. Jacobian evals. for Newton iter.
   
   :Default: 
   :Dimension: (ngrid)
   :Group: Stat
   :Type: integer
   :Unit: 
.. py:attribute:: newbcl
   
   switch on new sheath model
   
   :Default: 
   :Dimension: (2)
   :Group: Poten
   :Type: integer
   :Unit: 
.. py:attribute:: newbcr
   
   switch on new sheath model
   
   :Default: 
   :Dimension: (2)
   :Group: Poten
   :Type: integer
   :Unit: 
.. py:attribute:: bcel
   
   electron sheath energy transmission factor
   on the left boundary
   
   :Default: 
   :Dimension: (0:ny+1,nxpt)
   :Group: Poten
   :Type: double
   :Unit: 
.. py:attribute:: bcer
   
   electron sheath energy transmission factor
   on the right boundary
   
   :Default: 
   :Dimension: (0:ny+1,nxpt)
   :Group: Poten
   :Type: double
   :Unit: 
.. py:attribute:: bcil
   
   ion sheath energy transmission factor
   on the left boundary
   
   :Default: 
   :Dimension: (0:ny+1,nxpt)
   :Group: Poten
   :Type: double
   :Unit: 
.. py:attribute:: bcir
   
   ion sheath energy transmission factor
   on the right boundary
   
   :Default: 
   :Dimension: (0:ny+1,nxpt)
   :Group: Poten
   :Type: double
   :Unit: 
.. py:attribute:: kappal
   
   sheath pot'l drop on left boundary, phi/Te
   
   :Default: 
   :Dimension: (0:ny+1,nxpt)
   :Group: Poten
   :Type: double
   :Unit: 
.. py:attribute:: kappar
   
   sheath pot'l drop on right boundary, phi/Te
   
   :Default: 
   :Dimension: (0:ny+1,nxpt)
   :Group: Poten
   :Type: double
   :Unit: 
.. py:attribute:: bctype
   
   /0,ny*0,0/
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Poten
   :Type: integer
   :Unit: 
.. py:attribute:: phi0r
   
   plate pot'l at right poloidal boundary
   
   :Default: 
   :Dimension: (0:ny+1,nxpt)
   :Group: Poten
   :Type: double
   :Unit: V
.. py:attribute:: phi0l
   
   plate pot'l at left poloidal boundary
   
   :Default: 
   :Dimension: (0:ny+1,nxpt)
   :Group: Poten
   :Type: double
   :Unit: V
.. py:attribute:: capx
   
   /ny*0.0/
   
   :Default: 
   :Dimension: (1:ny)
   :Group: Poten
   :Type: double
   :Unit: 
.. py:attribute:: dphi_iy1
   
   /(nx+2)*0./ #incremental phi at iy=1 to have
   Te=constant for second phi BC
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Poten
   :Type: double
   :Unit: V
.. py:attribute:: kincorlb
   
   kinetic corr. factor for elec part. loss, left b
   
   :Default: 
   :Dimension: (0:ny+1,nxpt)
   :Group: Poten
   :Type: double
   :Unit: 
.. py:attribute:: kincorrb
   
   kinetic corr. factor for elec part. loss, right b
   
   :Default: 
   :Dimension: (0:ny+1,nxpt)
   :Group: Poten
   :Type: double
   :Unit: 
.. py:attribute:: ex
   
   poloidal electric field
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Gradients
   :Type: double
   :Unit: V/m
.. py:attribute:: ey
   
   radial electric field
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Gradients
   :Type: double
   :Unit: V/m
.. py:attribute:: gpix
   
   X-gradient of ion pressure
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Gradients
   :Type: double
   :Unit: Pa/m
.. py:attribute:: gpiy
   
   Y-gradient of ion pressure
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Gradients
   :Type: double
   :Unit: Pa/m
.. py:attribute:: gpex
   
   X-gradient of el. pressure
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Gradients
   :Type: double
   :Unit: Pa/m
.. py:attribute:: gpey
   
   Y-gradient of el. pressure
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Gradients
   :Type: double
   :Unit: Pa/m
.. py:attribute:: gprx
   
   X-gradient of total pressure
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Gradients
   :Type: double
   :Unit: Pa/m
.. py:attribute:: gpry
   
   Y-gradient of total pressure
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Gradients
   :Type: double
   :Unit: Pa/m
.. py:attribute:: gtex
   
   X-gradient of el. temperature
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Gradients
   :Type: double
   :Unit: J/m
.. py:attribute:: gtey
   
   Y-gradient of el. temperature
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Gradients
   :Type: double
   :Unit: J/m
.. py:attribute:: gtix
   
   X-gradient of ion temperature
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Gradients
   :Type: double
   :Unit: J/m
.. py:attribute:: gtiy
   
   Y-gradient of ion temperature
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Gradients
   :Type: double
   :Unit: J/m
.. py:attribute:: frice
   
   Electron parallel Coulomb friction
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Cfric
   :Type: double
   :Unit: J/m**4
.. py:attribute:: frici
   
   Ion parallel Coulomb friction
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nisp)
   :Group: Cfric
   :Type: double
   :Unit: J/m**4
.. py:attribute:: fricnrl
   
   NRL ion par fric ni*mi*nu*(up1-up2)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nusp)
   :Group: Cfric
   :Type: double
   :Unit: J/m**4
.. py:attribute:: isalfecalc
   
   =1 for internal calc of alfe
   
   :Default: 
   :Dimension: (1:nisp)
   :Group: Cfric
   :Type: integer
   :Unit: 
.. py:attribute:: isbetaicalc
   
   =1 for internal calc of betai
   
   :Default: 
   :Dimension: (1:nisp)
   :Group: Cfric
   :Type: integer
   :Unit: 
.. py:attribute:: alfe
   
   grad_Te thm force coeff isalfecalc=0
   
   :Default: 
   :Dimension: (1:nisp)
   :Group: Cfric
   :Type: double
   :Unit: 
.. py:attribute:: betai
   
   grad_Ti thm force coeff isbetaicalc=0
   
   :Default: 
   :Dimension: (1:nisp)
   :Group: Cfric
   :Type: double
   :Unit: 
.. py:attribute:: inewton
   
   =1 for Newton iter., =0 for time-dependent
   reset=1 internally if svrpkg=nksol or newton
   
   :Default: 
   :Dimension: (30)
   :Group: Grid
   :Type: integer
   :Unit: 
.. py:attribute:: ijac
   
   
   
   :Default: 
   :Dimension: (ngrid)
   :Group: Grid
   :Type: integer
   :Unit: 
.. py:attribute:: w
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Wkspace
   :Type: double
   :Unit: 
.. py:attribute:: w0
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Wkspace
   :Type: double
   :Unit: 
.. py:attribute:: w1
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Wkspace
   :Type: double
   :Unit: 
.. py:attribute:: w2
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Wkspace
   :Type: double
   :Unit: 
.. py:attribute:: w3
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Wkspace
   :Type: double
   :Unit: 
.. py:attribute:: flox
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Locflux
   :Type: double
   :Unit: 
.. py:attribute:: floy
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Locflux
   :Type: double
   :Unit: 
.. py:attribute:: conx
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Locflux
   :Type: double
   :Unit: 
.. py:attribute:: cony
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Locflux
   :Type: double
   :Unit: 
.. py:attribute:: floxe
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Locflux
   :Type: double
   :Unit: 
.. py:attribute:: floye
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Locflux
   :Type: double
   :Unit: 
.. py:attribute:: floxi
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Locflux
   :Type: double
   :Unit: 
.. py:attribute:: floyi
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Locflux
   :Type: double
   :Unit: 
.. py:attribute:: floxg
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Locflux
   :Type: double
   :Unit: 
.. py:attribute:: floyg
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Locflux
   :Type: double
   :Unit: 
.. py:attribute:: fgtdx
   
   scale factor for gas grad-x T vel
   
   :Default: 
   :Dimension: (0:nx+1)
   :Group: Locflux
   :Type: double
   :Unit: 
.. py:attribute:: fgtdy
   
   scale factor for gas grad-x T vel
   
   :Default: 
   :Dimension: (0:ny+1)
   :Group: Locflux
   :Type: double
   :Unit: 
.. py:attribute:: conxe
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Locflux
   :Type: double
   :Unit: 
.. py:attribute:: conye
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Locflux
   :Type: double
   :Unit: 
.. py:attribute:: conxi
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Locflux
   :Type: double
   :Unit: 
.. py:attribute:: conyi
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Locflux
   :Type: double
   :Unit: 
.. py:attribute:: conxg
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Locflux
   :Type: double
   :Unit: 
.. py:attribute:: conyg
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Locflux
   :Type: double
   :Unit: 
.. py:attribute:: floxge
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Locflux
   :Type: double
   :Unit: 
.. py:attribute:: floyge
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Locflux
   :Type: double
   :Unit: 
.. py:attribute:: conxge
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Locflux
   :Type: double
   :Unit: 
.. py:attribute:: conyge
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Locflux
   :Type: double
   :Unit: 
.. py:attribute:: visx
   
   poloidal viscosity coeff.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Conduc
   :Type: double
   :Unit: kg/m s
.. py:attribute:: visy
   
   radial viscosity coeff.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Conduc
   :Type: double
   :Unit: kg/m s
.. py:attribute:: hcxe
   
   poloidal elec. therm. conduct.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: 1/m s
.. py:attribute:: hcye
   
   radial elec. therm. conduct.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: 1/m s
.. py:attribute:: hcxij
   
   j-species pol. ion therm. conduct.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Conduc
   :Type: double
   :Unit: 1/m s
.. py:attribute:: hcyij
   
   j-species rad. ion therm. conduct.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Conduc
   :Type: double
   :Unit: 1/m s
.. py:attribute:: hcxg
   
   j-species pol. gas therm. conduct.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Conduc
   :Type: double
   :Unit: 1/m s
.. py:attribute:: hcyg
   
   j-species rad. gas therm. conduct.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Conduc
   :Type: double
   :Unit: 1/m s
.. py:attribute:: hcxi
   
   summed pol. ion+neut therm. conduct.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: 1/m s
.. py:attribute:: hcxineo
   
   neocl. pol. ion+neut therm. conduct.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: 1/m s
.. py:attribute:: hcyi
   
   summed rad. ion+neut therm. conduct.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: 1/m s
.. py:attribute:: hcxn
   
   poloidal neutral therm. conduct.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: 1/m s
.. py:attribute:: hcyn
   
   radial neutral therm. conduct.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: 1/m s
.. py:attribute:: kxbohm
   
   spatially depend. diff. on x-face
   set by user; Bohm if isbohmcalc=1
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kybohm
   
   spatially depend. diff. on y-face
   set by user; Bohm if isbohmcalc=1
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: m**2/s
.. py:attribute:: vybohm
   
   spatially depend. convect. y-vel
   set user if isbohmcalc=0; else =0
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: m/s
.. py:attribute:: dif_use
   
   spatially depend. diff; if
   isbohmcalc=1, user input if all
   facbni+facbup+facbee+facbei =0,
   or kybohm if facbni, etc. > 0;
   if isbohmcalc=2, then
   D = difni*kybohm/(difni+kybohm)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Conduc
   :Type: double
   :Unit: m**2/s
.. py:attribute:: difp_use
   
   for gen pr diff; see dif_use comment
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Conduc
   :Type: double
   :Unit: m**2/s
.. py:attribute:: dif2_use
   
   for dif2; see dif_use comment
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Conduc
   :Type: double
   :Unit: m**2/s
.. py:attribute:: tray_use
   
   for travis; see dif_use comment
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Conduc
   :Type: double
   :Unit: m**2/s
.. py:attribute:: trax_use
   
   pol. analog to tra_use
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Conduc
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kye_use
   
   for kye; see dif_use comment
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kyi_use
   
   for kyi; see dif_use comment
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kxe_use
   
   user elec pol. heat cond
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kxi_use
   
   user ion pol. heat cond.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kxg_use
   
   user gas pol. heat cond.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Conduc
   :Type: double
   :Unit: m**2/s
.. py:attribute:: kyg_use
   
   user gas rad. heat cond.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Conduc
   :Type: double
   :Unit: m**2/s
.. py:attribute:: dutm_use
   
   for difutm; see dif_use comment
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Conduc
   :Type: double
   :Unit: m**2/s
.. py:attribute:: vy_use
   
   user-set rad vel;for isbohmcalc=0
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Conduc
   :Type: double
   :Unit: m/s
.. py:attribute:: vyup_use
   
   user-set conv vel of ion || vel, up
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: m/s
.. py:attribute:: vyte_use
   
   user-set rad elec eng vel
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: m/s
.. py:attribute:: vyti_use
   
   user-set rad ion eng vel
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: m/s
.. py:attribute:: fniyos_use
   
   user-set particle flux
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Conduc
   :Type: double
   :Unit: 1/s m**2
.. py:attribute:: feeyosn_use
   
   user-set Te energy flux
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: J/s m**2
.. py:attribute:: feiyosn_use
   
   user-set Ti energy flux
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: J/s m**2
.. py:attribute:: vy_cft
   
   calc vy from fniyos_use (fix flux)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Conduc
   :Type: double
   :Unit: m/s
.. py:attribute:: vyte_cft
   
   calc vyte from feeyos_use (fix flux)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: m/s
.. py:attribute:: vyti_cft
   
   calc vyte from feiyos_use (fix flux)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: m/s
.. py:attribute:: nuiz
   
   ionization rate (=ne*sigma*v)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,ngsp)
   :Group: Conduc
   :Type: double
   :Unit: 1/s
.. py:attribute:: nucx
   
   charge-exchg rate for neut(sigv*ni)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,ngsp)
   :Group: Conduc
   :Type: double
   :Unit: 1/s
.. py:attribute:: nucxi
   
   charge-exchg rate for ion (sigv*ng)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nisp)
   :Group: Conduc
   :Type: double
   :Unit: 1/s
.. py:attribute:: nueli
   
   elast scatt rate for ion (sigv*ng)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nisp)
   :Group: Conduc
   :Type: double
   :Unit: 1/s
.. py:attribute:: nuelg
   
   elast scatt rate for gas (sigv*nimp)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,ngsp)
   :Group: Conduc
   :Type: double
   :Unit: 1/s
.. py:attribute:: nuix
   
   fnuizx*nuiz+fnucxx*nucx
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,ngsp)
   :Group: Conduc
   :Type: double
   :Unit: 1/s
.. py:attribute:: nurc
   
   recombination rate
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,ngsp)
   :Group: Conduc
   :Type: double
   :Unit: 1/s
.. py:attribute:: nuvl
   
   vol loss rate, ~cs/l_parloss for 1-D
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nisp)
   :Group: Conduc
   :Type: double
   :Unit: 1/s
.. py:attribute:: cfvli
   
   /nisp*0./#scal fac for individ ion rate nuvl
   
   :Default: 
   :Dimension: (nisp)
   :Group: Conduc
   :Type: double
   :Unit: 
.. py:attribute:: eqp
   
   Te,i equipart. fact; needs *(Te-Ti)*vol
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: 1/m**3
.. py:attribute:: eqpg
   
   Tg,i equipart. fact; needs *(Tg-Ti)*vol
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,ngsp)
   :Group: Conduc
   :Type: double
   :Unit: 1/m**3
.. py:attribute:: engcoolm
   
   cool rate ion/atoms by mols if ishymol=1
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: J/s
.. py:attribute:: eeli
   
   electron energy loss per ionization
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: J
.. py:attribute:: eta1
   
   Braginskii ion visc coeff eta_1
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: J-s/m**3
.. py:attribute:: rtaue
   
   Brag. R coeff (t_e/me)/(w_ce*t_e)**2
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: s/kg
.. py:attribute:: dclass_e
   
   classical elec perp heat conduc.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: m**2/s
.. py:attribute:: dclass_i
   
   classical ion perp heat conduc.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Conduc
   :Type: double
   :Unit: m**2/s
.. py:attribute:: visxneo
   
   Braginskii eta_0 neo-modified
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Conduc
   :Type: double
   :Unit: kg/m s
.. py:attribute:: visvol_v
   
   vel-based viscosity in (n*m*up)^dot eqn
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Conduc
   :Type: double
   :Unit: 
.. py:attribute:: visvol_q
   
   heat-flux-based viscosity (n*m*up)^dot eqn
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Conduc
   :Type: double
   :Unit: 
.. py:attribute:: nuii
   
   Braginski nuii coll freq.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Conduc
   :Type: double
   :Unit: 
.. py:attribute:: nuiistar
   
   neoclassical nuii coll freq.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Conduc
   :Type: double
   :Unit: 
.. py:attribute:: alfneo
   
   neoclassical factor for q-based visc.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Conduc
   :Type: double
   :Unit: 
.. py:attribute:: k2neo
   
   neoclassical coeff reducing therm cond
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Conduc
   :Type: double
   :Unit: 
.. py:attribute:: ktneo
   
   neoclassical coeff of grad Ti
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Conduc
   :Type: double
   :Unit: 
.. py:attribute:: snic
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Rhsides
   :Type: double
   :Unit: 
.. py:attribute:: sniv
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Rhsides
   :Type: double
   :Unit: 
.. py:attribute:: psorc
   
   cell ctr ioniz. sor plasma (>0)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Rhsides
   :Type: double
   :Unit: part/s
.. py:attribute:: psor
   
   cell ave ioniz. sor plasma (>0)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Rhsides
   :Type: double
   :Unit: part/s
.. py:attribute:: psort
   
   ioniz. source for plasma (>0)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Rhsides
   :Type: double
   :Unit: part/s
.. py:attribute:: psorxrc
   
   cell ctr cx &recomb. for ions (<0)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Rhsides
   :Type: double
   :Unit: part/s
.. py:attribute:: psorxr
   
   cell ave cx &recomb. for ions (<0)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Rhsides
   :Type: double
   :Unit: part/s
.. py:attribute:: psor_tmpov
   
   work array for psor,etc for ave
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Rhsides
   :Type: double
   :Unit: part/s
.. py:attribute:: psorgc
   
   cell ctr ioniz. sor neutral (<0)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Rhsides
   :Type: double
   :Unit: part/s
.. py:attribute:: psorg
   
   cell ave ioniz. sor neutral (<0)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Rhsides
   :Type: double
   :Unit: part/s
.. py:attribute:: psorrgc
   
   cell ctr recomb. source for neutrals
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Rhsides
   :Type: double
   :Unit: part/s
.. py:attribute:: psorrg
   
   cell ave recomb. source for neutrals
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Rhsides
   :Type: double
   :Unit: part/s
.. py:attribute:: psorcxgc
   
   cell ctr cx source for neutrals
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Rhsides
   :Type: double
   :Unit: part/s
.. py:attribute:: psorcxg
   
   cell ave cx source for neutrals
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Rhsides
   :Type: double
   :Unit: part/s
.. py:attribute:: psori
   
   impurity gas source
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Rhsides
   :Type: double
   :Unit: part/s
.. py:attribute:: psordis
   
   diss. source of hydrogen
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Rhsides
   :Type: double
   :Unit: part/s
.. py:attribute:: psorbgg
   
   diag artific neut backg source
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Rhsides
   :Type: double
   :Unit: part/s
.. py:attribute:: psorbgz
   
   diag artific impur backg source
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Rhsides
   :Type: double
   :Unit: part/s
.. py:attribute:: erliz
   
   H rad'n loss for ioniz'n
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Rhsides
   :Type: double
   :Unit: J/s
.. py:attribute:: erlrc
   
   H rad'n loss for recom'n
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Rhsides
   :Type: double
   :Unit: J/s
.. py:attribute:: vsoreec
   
   cell ctr tot elec vol eng source
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Rhsides
   :Type: double
   :Unit: J/s
.. py:attribute:: vsoree
   
   cell ave tot elec vol eng source
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Rhsides
   :Type: double
   :Unit: J/s
.. py:attribute:: pwrebkg
   
   elec energy backgrd source; limits te~tebg
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Rhsides
   :Type: double
   :Unit: W/m**3
.. py:attribute:: pwribkg
   
   ion energy backgrd source; limits ti~tibg
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Rhsides
   :Type: double
   :Unit: W/m**3
.. py:attribute:: wjdote
   
   Joule heating rate
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Rhsides
   :Type: double
   :Unit: J/s
.. py:attribute:: smoc
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nusp)
   :Group: Rhsides
   :Type: double
   :Unit: 
.. py:attribute:: smov
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nusp)
   :Group: Rhsides
   :Type: double
   :Unit: 
.. py:attribute:: msor
   
   ioniz. mom. source for ions
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Rhsides
   :Type: double
   :Unit: kg-m/s**2
.. py:attribute:: msorxr
   
   cx&recomb. mom. sink for ions
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Rhsides
   :Type: double
   :Unit: kg-m/s**2
.. py:attribute:: seec
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Rhsides
   :Type: double
   :Unit: 
.. py:attribute:: seev
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Rhsides
   :Type: double
   :Unit: 
.. py:attribute:: seic
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Rhsides
   :Type: double
   :Unit: 
.. py:attribute:: seiv
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Rhsides
   :Type: double
   :Unit: 
.. py:attribute:: resco
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: Rhsides
   :Type: double
   :Unit: 
.. py:attribute:: resng
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Rhsides
   :Type: double
   :Unit: 
.. py:attribute:: reseg
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:ngsp)
   :Group: Rhsides
   :Type: double
   :Unit: 
.. py:attribute:: resmo
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nusp)
   :Group: Rhsides
   :Type: double
   :Unit: 
.. py:attribute:: resee
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Rhsides
   :Type: double
   :Unit: 
.. py:attribute:: resei
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Rhsides
   :Type: double
   :Unit: 
.. py:attribute:: resphi
   
   
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Rhsides
   :Type: double
   :Unit: 
.. py:attribute:: ng_mc
   
   neutral gas density from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: part/m**3
.. py:attribute:: ng_mc_rsd
   
   neutral gas density rsd from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: #
.. py:attribute:: pg_mc
   
   neutral gas pressure from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: Pa
.. py:attribute:: pg_mc_rsd
   
   neutral gas pressure rsd from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: #
.. py:attribute:: ng_ue
   
   neutral gas density from Monte-Carlo-Neutrals model, blended with fluid result
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: part/m**3
.. py:attribute:: ng_ue_rsd
   
   neutral gas density rsd from Monte-Carlo-Neutrals model, blended with fluid result
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: #
.. py:attribute:: upg_ue
   
   neutral gas parallel velocity from Monte-Carlo-Neutrals model, blended with fluid result
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: m/s
.. py:attribute:: upg_ue_rsd
   
   neutral gas parallel velocity rsd from Monte-Carlo-Neutrals model, blended with fluid result
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: #
.. py:attribute:: pg_ue
   
   neutral gas pressure from Monte-Carlo-Neutrals model, blended with fluid result
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: Pa
.. py:attribute:: pg_ue_rsd
   
   neutral gas pressure rsd from Monte-Carlo-Neutrals model, blended with fluid result
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: #
.. py:attribute:: tg_ue
   
   neutral gas temperature from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: J
.. py:attribute:: tg_ue_rsd
   
   neutral gas temperature rsd from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: #
.. py:attribute:: sng_ue
   
   neutral particle source density (convective only)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: part/m**3-s
.. py:attribute:: smg_ue
   
   neutral parallel momentum source density
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: N/m**3
.. py:attribute:: seg_ue
   
   neutral energy source density (convective only)
   ## Vectors ###
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: W/m**3
.. py:attribute:: jng_mc
   
   neutral gas particle flux density from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl,3)
   :Group: MCN_sources
   :Type: double
   :Unit: part/s*m**2
.. py:attribute:: jng_mc_rsd
   
   relative standard deviation
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl,3)
   :Group: MCN_sources
   :Type: double
   :Unit: part/s*m**2
.. py:attribute:: jng_ue
   
   neutral gas particle flux density interpolated to UEDGE grid
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl,3)
   :Group: MCN_sources
   :Type: double
   :Unit: part/s*m**2
.. py:attribute:: jng_ue_rsd
   
   relative standard deviation
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl,3)
   :Group: MCN_sources
   :Type: double
   :Unit: #
.. py:attribute:: vg_mc
   
   neutral gas velocity from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl,3)
   :Group: MCN_sources
   :Type: double
   :Unit: m/s
.. py:attribute:: vg_mc_rsd
   
   relative standard deviation
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl,3)
   :Group: MCN_sources
   :Type: double
   :Unit: #
.. py:attribute:: vg_ue
   
   neutral gas velocity interpolated to UEDGE grid
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl,3)
   :Group: MCN_sources
   :Type: double
   :Unit: m/s
.. py:attribute:: vg_ue_rsd
   
   relative standard deviation
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl,3)
   :Group: MCN_sources
   :Type: double
   :Unit: #
.. py:attribute:: fngx_mc
   
   blended poloidal neutral gas particle flux from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: part/s
.. py:attribute:: fngx_mc_rsd
   
   relative standard deviation
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: #
.. py:attribute:: fngy_mc
   
   blended poloidal neutral gas particle flux from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: part/s
.. py:attribute:: fngy_mc_rsd
   
   relative standard deviation
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: #
.. py:attribute:: fngx_ue
   
   blended poloidal neutral gas particle flux on UEDGE grid
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: part/s
.. py:attribute:: fngx_ue_rsd
   
   relative standard deviation
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: #
.. py:attribute:: fngy_ue
   
   blended poloidal neutral gas particle flux on UEDGE grid
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: part/s
.. py:attribute:: fngy_ue_rsd
   
   relative standard deviation
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: #
.. py:attribute:: fmgx_mc
   
   blended poloidal neutral gas momentum flux from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: N
.. py:attribute:: fmgx_mc_rsd
   
   relative standard deviation
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: #
.. py:attribute:: fmgy_mc
   
   blended radial neutral gas momentum flux from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: N
.. py:attribute:: fmgy_mc_rsd
   
   relative standard deviation
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: #
.. py:attribute:: fmgx_ue
   
   blended poloidal neutral gas momentum flux on UEDGE grid
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: N
.. py:attribute:: fmgx_ue_rsd
   
   relative standard deviation
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: #
.. py:attribute:: fmgy_ue
   
   blended radial neutral gas momentum flux on UEDGE grid
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: N
.. py:attribute:: fmgy_ue_rsd
   
   relative standard deviation
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: #
.. py:attribute:: fmgxy_ue
   
   blended poloidal neutral gas momentum flux on nonorthog. UEDGE grid
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: N
.. py:attribute:: fmgxy_ue_rsd
   
   relative standard deviation
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: #
.. py:attribute:: jeg_mc
   
   neutral gas energy flux density vector from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl,3)
   :Group: MCN_sources
   :Type: double
   :Unit: W/m**2
.. py:attribute:: jeg_mc_rsd
   
   relative standard deviation
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl,3)
   :Group: MCN_sources
   :Type: double
   :Unit: #
.. py:attribute:: jeg_ue
   
   neutral gas energy flux density vector interpolated to UEDGE grid
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl,3)
   :Group: MCN_sources
   :Type: double
   :Unit: W/m**2
.. py:attribute:: jeg_ue_rsd
   
   relative standard deviation
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl,3)
   :Group: MCN_sources
   :Type: double
   :Unit: #
.. py:attribute:: fegx_mc
   
   blended poloidal neutral gas heat flux from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: W
.. py:attribute:: fegx_mc_rsd
   
   relative standard deviation
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: #
.. py:attribute:: fegy_mc
   
   blended radial neutral gas heat flux from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: W
.. py:attribute:: fegy_mc_rsd
   
   relative standard deviation
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: #
.. py:attribute:: fegx_ue
   
   blended poloidal neutral gas heat flux on UEDGE grid
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: W
.. py:attribute:: fegx_ue_rsd
   
   relative standard deviation
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: #
.. py:attribute:: fegy_ue
   
   blended radial neutral gas heat flux on UEDGE grid
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: W
.. py:attribute:: fegy_ue_rsd
   
   relative standard deviation
   ## Tensors ###
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: #
.. py:attribute:: stressg_mc
   
   neutral gas stress tensor from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl,1:3,1:3)
   :Group: MCN_sources
   :Type: double
   :Unit: Pa
.. py:attribute:: stressg_mc_rsd
   
   neutral gas stress tensor from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl,1:3,1:3)
   :Group: MCN_sources
   :Type: double
   :Unit: Pa
.. py:attribute:: stressg_ue
   
   neutral gas stress tensor interpolated to UEDGE grid
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl,1:3,1:3)
   :Group: MCN_sources
   :Type: double
   :Unit: Pa
.. py:attribute:: stressg_ue_rsd
   
   neutral gas stress tensor interpolated to UEDGE grid
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl,1:3,1:3)
   :Group: MCN_sources
   :Type: double
   :Unit: Pa
.. py:attribute:: pxz_mc
   
   neutral gas pressure from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: Pa
.. py:attribute:: pxz_mc_rsd
   
   neutral gas pressure rsd from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,nfl)
   :Group: MCN_sources
   :Type: double
   :Unit: Pa
.. py:attribute:: mcnsor_ni
   
   ion particle source from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp,1:nstra)
   :Group: MCN_sources
   :Type: double
   :Unit: part/s
.. py:attribute:: mcnsor_up
   
   ion parallel momentum source from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp,1:nstra)
   :Group: MCN_sources
   :Type: double
   :Unit: kg-m/s**2
.. py:attribute:: mcnsor_te
   
   electron thermal energy source from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nstra)
   :Group: MCN_sources
   :Type: double
   :Unit: J/s
.. py:attribute:: mcnsor_ti
   
   ion thermal energy source from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nstra)
   :Group: MCN_sources
   :Type: double
   :Unit: J/s
.. py:attribute:: mcncurr
   
   neutral source current from each strata in Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (1:nstra)
   :Group: MCN_sources
   :Type: double
   :Unit: part/s
.. py:attribute:: uesor_ni
   
   scaled ion particle source from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: MCN_sources
   :Type: double
   :Unit: part/s
.. py:attribute:: uesor_up
   
   scaled ion parallel momentum source from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: MCN_sources
   :Type: double
   :Unit: kg-m/s**2
.. py:attribute:: uesor_te
   
   scaled electron thermal energy source from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: MCN_sources
   :Type: double
   :Unit: J/s
.. py:attribute:: uesor_ti
   
   scaled ion thermal energy source from Monte-Carlo-Neutrals model
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: MCN_sources
   :Type: double
   :Unit: J/s
.. py:attribute:: uecurr
   
   neutral source current from each strata according to UEDGE plasma model
   
   :Default: 
   :Dimension: (1:nstra)
   :Group: MCN_sources
   :Type: double
   :Unit: part/s
.. py:attribute:: olduecurr
   
   neutral source current from each strata according to UEDGE plasma model
   
   :Default: 
   :Dimension: (1:nstra)
   :Group: MCN_sources
   :Type: double
   :Unit: part/s
.. py:attribute:: strascal
   
   scaling factor for plasma source terms due to each strata
   
   :Default: 
   :Dimension: (1:nstra)
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: wsor
   
   normalization constant for plasma source terms from EIRENE file 'fort.32'
   
   :Default: 
   :Dimension: (1:nstra)
   :Group: MCN_sources
   :Type: double
   :Unit: part/s
.. py:attribute:: esor
   
   unused constant from EIRENE file 'fort.32'
   
   :Default: 
   :Dimension: (1:nstra)
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: sni
   
   normalized ion particle sources from EIRENE file 'fort.32'
   or absolute ion particle source from DEGAS2
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nfl,1:nstra)
   :Group: MCN_sources
   :Type: double
   :Unit: part/s
.. py:attribute:: smo
   
   normalized ion parallel momentum sources from EIRENE file 'fort.32',
   or absolute ion parallel momentum source from DEGAS2
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nfl,1:nstra)
   :Group: MCN_sources
   :Type: double
   :Unit: kg-m/s**2
.. py:attribute:: smor
   
   'radial' component of ion momentum source from DEGAS2
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nfl,1:nstra)
   :Group: MCN_sources
   :Type: double
   :Unit: kg-m/s**2
.. py:attribute:: smophi
   
   'toroidal' component of ion momentum source from DEGAS2
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nfl,1:nstra)
   :Group: MCN_sources
   :Type: double
   :Unit: kg-m/s**2
.. py:attribute:: smoz
   
   'vertical' component of ion momentum source from DEGAS2
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nfl,1:nstra)
   :Group: MCN_sources
   :Type: double
   :Unit: kg-m/s**2
.. py:attribute:: see
   
   normalized electron energy source from EIRENE file 'fort.32'
   or absolute electron energy source from DEGAS2
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nstra)
   :Group: MCN_sources
   :Type: double
   :Unit: J/s
.. py:attribute:: sei
   
   normalized ion energy source from EIRENE file 'fort.32'
   or absolute ion energy source from DEGAS2
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nstra)
   :Group: MCN_sources
   :Type: double
   :Unit: J/s
.. py:attribute:: labela
   
   data from Monte Carlo neutrals code:
   
   :Default: 
   :Dimension: (1:12)
   :Group: MCN_sources
   :Type: character(8)
   :Unit: 
.. py:attribute:: labelm
   
   data from Monte Carlo neutrals code:
   
   :Default: 
   :Dimension: (1:12)
   :Group: MCN_sources
   :Type: character(8)
   :Unit: 
.. py:attribute:: labeli
   
   data from Monte Carlo neutrals code:
   
   :Default: 
   :Dimension: (1:12)
   :Group: MCN_sources
   :Type: character(8)
   :Unit: 
.. py:attribute:: naf
   
   data from Monte Carlo neutrals code:
   atomic neutral density
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:natmi)
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: taf
   
   data from Monte Carlo neutrals code:
   atomic neutral temperature
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:natmi)
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: nmf
   
   data from Monte Carlo neutrals code:
   molecular neutral density
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:nmoli)
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: tmf
   
   data from Monte Carlo neutrals code:
   molecular neutral temperature
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:nmoli)
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: ntf
   
   data from Monte Carlo neutrals code:
   molecular ion density
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:nioni)
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: ttf
   
   data from Monte Carlo neutrals code:
   molecular ion temperature
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:nioni)
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: fnax
   
   data from Monte Carlo neutrals code:
   x-particle flux of atomic neutrals
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:natmi)
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: fnmx
   
   data from Monte Carlo neutrals code:
   x-particle flux of molecular neutrals
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:nmoli)
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: fntx
   
   data from Monte Carlo neutrals code:
   x-particle flux of molecular ions
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:nioni)
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: fnay
   
   data from Monte Carlo neutrals code:
   y-particle flux of atomic neutrals
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:natmi)
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: fnmy
   
   data from Monte Carlo neutrals code:
   y-particle flux of molecular neutrals
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:nmoli)
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: fnty
   
   data from Monte Carlo neutrals code:
   y-particle flux of molecular ions
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:nioni)
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: fnaz
   
   data from Monte Carlo neutrals code:
   z-particle flux of atomic neutrals
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:natmi)
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: fnmz
   
   data from Monte Carlo neutrals code:
   z-particle flux of molecular neutrals
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:nmoli)
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: fntz
   
   data from Monte Carlo neutrals code:
   z-particle flux of molecular ions
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:nioni)
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: feay
   
   data from EIRENE file fort.44: y-energy flux of atomic neutrals
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:natmi)
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: femy
   
   data from EIRENE file fort.44: y-energy flux of molecular neutrals
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:nmoli)
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: feax
   
   data from EIRENE file fort.44: x-energy flux of atomic neutrals
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:natmi)
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: femx
   
   data from EIRENE file fort.44: x-energy flux of molecular neutrals
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:nmoli)
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: hatm
   
   data from EIRENE file fort.44: h-alpha radiation from atomic neutrals
   
   :Default: 
   :Dimension: (1:nxf,1:nyf)
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: hmol
   
   data from EIRENE file fort.44: h-alpha radiation from molecular neutrals
   
   :Default: 
   :Dimension: (1:nxf,1:nyf)
   :Group: MCN_sources
   :Type: double
   :Unit: 
.. py:attribute:: labelmc
   
   labels for Monte Carlo species
   
   :Default: 
   :Dimension: (1:12)
   :Group: MCN_test
   :Type: character(8)
   :Unit: 
.. py:attribute:: nmc
   
   density from Monte Carlo neutrals code
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:nmcsp)
   :Group: MCN_test
   :Type: double
   :Unit: 
.. py:attribute:: tmc
   
   temperature from Monte Carlo neutrals code
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:nmcsp)
   :Group: MCN_test
   :Type: double
   :Unit: 
.. py:attribute:: fnmcx
   
   x-component of particle flux from Monte Carlo neutrals code
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:nmcsp)
   :Group: MCN_test
   :Type: double
   :Unit: 
.. py:attribute:: fnmcy
   
   y-component of particle flux from Monte Carlo neutrals code
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:nmcsp)
   :Group: MCN_test
   :Type: double
   :Unit: 
.. py:attribute:: fnmcz
   
   z-component of particle flux from Monte Carlo neutrals code
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:nmcsp)
   :Group: MCN_test
   :Type: double
   :Unit: 
.. py:attribute:: femcx
   
   x-component of energy flux from Monte Carlo neutrals code
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:nmcsp)
   :Group: MCN_test
   :Type: double
   :Unit: 
.. py:attribute:: femcy
   
   y-component of energy flux from Monte Carlo neutrals code
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:nmcsp)
   :Group: MCN_test
   :Type: double
   :Unit: 
.. py:attribute:: femcz
   
   z-component of energy flux from Monte Carlo neutrals code
   
   :Default: 
   :Dimension: (1:nxf,1:nyf,1:nmcsp)
   :Group: MCN_test
   :Type: double
   :Unit: 
.. py:attribute:: v2c
   
   v2 velocity component at cell centers
   
   :Default: 
   :Dimension: (1:nx,1:ny,1:nisp)
   :Group: MCN_bkgd
   :Type: double
   :Unit: m/s
.. py:attribute:: vyc
   
   vy velocity component at cell centers
   
   :Default: 
   :Dimension: (1:nx,1:ny,1:nisp)
   :Group: MCN_bkgd
   :Type: double
   :Unit: m/s
.. py:attribute:: upc
   
   up velocity component at cell centers
   
   :Default: 
   :Dimension: (1:nx,1:ny,1:nisp)
   :Group: MCN_bkgd
   :Type: double
   :Unit: m/s
.. py:attribute:: uuc
   
   uu velocity component at cell centers
   
   :Default: 
   :Dimension: (1:nx,1:ny,1:nisp)
   :Group: MCN_bkgd
   :Type: double
   :Unit: m/s
.. py:attribute:: utc
   
   ut velocity component at cell centers
   
   :Default: 
   :Dimension: (1:nx,1:ny,1:nisp)
   :Group: MCN_bkgd
   :Type: double
   :Unit: m/s
.. py:attribute:: vr
   
   vr velocity component at cell centers
   
   :Default: 
   :Dimension: (1:nx,1:ny,1:nisp)
   :Group: MCN_bkgd
   :Type: double
   :Unit: m/s
.. py:attribute:: vphi
   
   vphi velocity component at cell centers
   
   :Default: 
   :Dimension: (1:nx,1:ny,1:nisp)
   :Group: MCN_bkgd
   :Type: double
   :Unit: m/s
.. py:attribute:: vz
   
   vz velocity component at cell centers
   
   :Default: 
   :Dimension: (1:nx,1:ny,1:nisp)
   :Group: MCN_bkgd
   :Type: double
   :Unit: m/s
.. py:attribute:: v2tg1
   
   v2 velocity component at target plate number 1 (ix=0)
   
   :Default: 
   :Dimension: (1:ny,1:nisp)
   :Group: MCN_bkgd
   :Type: double
   :Unit: m/s
.. py:attribute:: vytg1
   
   vy velocity component at target plate number 1 (ix=0)
   
   :Default: 
   :Dimension: (1:ny,1:nisp)
   :Group: MCN_bkgd
   :Type: double
   :Unit: m/s
.. py:attribute:: uptg1
   
   up velocity component at target plate number 1 (ix=0)
   
   :Default: 
   :Dimension: (1:ny,1:nisp)
   :Group: MCN_bkgd
   :Type: double
   :Unit: m/s
.. py:attribute:: uutg1
   
   uu velocity component at target plate number 1 (ix=0)
   
   :Default: 
   :Dimension: (1:ny,1:nisp)
   :Group: MCN_bkgd
   :Type: double
   :Unit: m/s
.. py:attribute:: uttg1
   
   ut velocity component at target plate number 1 (ix=0)
   
   :Default: 
   :Dimension: (1:ny,1:nisp)
   :Group: MCN_bkgd
   :Type: double
   :Unit: m/s
.. py:attribute:: vrtg1
   
   vr velocity component at target plate number 1 (ix=0)
   
   :Default: 
   :Dimension: (1:ny,1:nisp)
   :Group: MCN_bkgd
   :Type: double
   :Unit: m/s
.. py:attribute:: vphitg1
   
   vphi velocity component at target plate number 1 (ix=0)
   
   :Default: 
   :Dimension: (1:ny,1:nisp)
   :Group: MCN_bkgd
   :Type: double
   :Unit: m/s
.. py:attribute:: vztg1
   
   vz velocity component at target plate number 1 (ix=0)
   
   :Default: 
   :Dimension: (1:ny,1:nisp)
   :Group: MCN_bkgd
   :Type: double
   :Unit: m/s
.. py:attribute:: v2tg2
   
   v2 velocity component at target plate number 2 (ix=nx)
   
   :Default: 
   :Dimension: (1:ny,1:nisp)
   :Group: MCN_bkgd
   :Type: double
   :Unit: m/s
.. py:attribute:: vytg2
   
   vy velocity component at target plate number 2 (ix=nx)
   
   :Default: 
   :Dimension: (1:ny,1:nisp)
   :Group: MCN_bkgd
   :Type: double
   :Unit: m/s
.. py:attribute:: uptg2
   
   up velocity component at target plate number 2 (ix=nx)
   
   :Default: 
   :Dimension: (1:ny,1:nisp)
   :Group: MCN_bkgd
   :Type: double
   :Unit: m/s
.. py:attribute:: uutg2
   
   uu velocity component at target plate number 2 (ix=nx)
   
   :Default: 
   :Dimension: (1:ny,1:nisp)
   :Group: MCN_bkgd
   :Type: double
   :Unit: m/s
.. py:attribute:: uttg2
   
   ut velocity component at target plate number 2 (ix=nx)
   
   :Default: 
   :Dimension: (1:ny,1:nisp)
   :Group: MCN_bkgd
   :Type: double
   :Unit: m/s
.. py:attribute:: vrtg2
   
   vr velocity component at target plate number 2 (ix=nx)
   
   :Default: 
   :Dimension: (1:ny,1:nisp)
   :Group: MCN_bkgd
   :Type: double
   :Unit: m/s
.. py:attribute:: vphitg2
   
   vphi velocity component at target plate number 2 (ix=nx)
   
   :Default: 
   :Dimension: (1:ny,1:nisp)
   :Group: MCN_bkgd
   :Type: double
   :Unit: m/s
.. py:attribute:: vztg2
   
   vz velocity component at target plate number 2 (ix=nx)
   
   :Default: 
   :Dimension: (1:ny,1:nisp)
   :Group: MCN_bkgd
   :Type: double
   :Unit: m/s
.. py:attribute:: uedgecmd
   
   uedge command
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(16)
   :Unit: 
.. py:attribute:: uedgescript
   
   uedge script to run
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(16)
   :Unit: 
.. py:attribute:: uedgefile
   
   uedge output file
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(16)
   :Unit: 
.. py:attribute:: uedgesave
   
   uedge save file
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(16)
   :Unit: 
.. py:attribute:: eirenecmd
   
   eirene command
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(16)
   :Unit: 
.. py:attribute:: eirenefile
   
   eirene output file
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(16)
   :Unit: 
.. py:attribute:: degas2cmd
   
   degas2 MC executable
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(32)
   :Unit: 
.. py:attribute:: degas2mpi
   
   degas2 MC executable for use with MPI
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(32)
   :Unit: 
.. py:attribute:: degas2file
   
   degas2 output file
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(32)
   :Unit: 
.. py:attribute:: gecmd
   
   degas2 readgeometry executable
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(32)
   :Unit: 
.. py:attribute:: geufile
   
   geometry input file for degas2
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(32)
   :Unit: 
.. py:attribute:: gedfile
   
   degas2 readgeometry output file
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(32)
   :Unit: 
.. py:attribute:: bkcmd
   
   degas2 readbackground executable
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(32)
   :Unit: 
.. py:attribute:: bkufile
   
   uedge output file for degas2 readbackground
   NOTE: the same file(grid) must also be used for readgeometry, as specified in the geufile
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(32)
   :Unit: 
.. py:attribute:: bkdfile
   
   degas2 readbackground output file
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(32)
   :Unit: 
.. py:attribute:: degas2outcmd
   
   degas2 outputbrowser executable
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(32)
   :Unit: 
.. py:attribute:: degas2outscript
   
   degas2 outputbrowser input file
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(32)
   :Unit: 
.. py:attribute:: degas2outsh
   
   sed script to clean up output files
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(32)
   :Unit: 
.. py:attribute:: mcnflights
   
   number of mc pseudo-particle trajectories
   
   :Default: 
   :Dimension: (1:10)
   :Group: Ext_neutrals
   :Type: integer
   :Unit: 
.. py:attribute:: ncsetcmd
   
   netcdf file editor command
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(32)
   :Unit: 
.. py:attribute:: ncsetvar
   
   variable to edit in bkdfile netcdf file
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(32)
   :Unit: 
.. py:attribute:: mpicmd
   
   MPI command
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(64)
   :Unit: 
.. py:attribute:: npopt
   
   option to specify # procs
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(8)
   :Unit: 
.. py:attribute:: runid
   
   description of run
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(80)
   :Unit: 
.. py:attribute:: neut_output_dir
   
   output directory
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(64)
   :Unit: 
.. py:attribute:: neut_ng_file
   
   neutral density file
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(64)
   :Unit: 
.. py:attribute:: neut_pg_file
   
   neutral pressure file
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(64)
   :Unit: 
.. py:attribute:: neut_jng1_file
   
   neutral particle flux: R
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(64)
   :Unit: 
.. py:attribute:: neut_jng2_file
   
   neutral particle flux: T
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(64)
   :Unit: 
.. py:attribute:: neut_jng3_file
   
   neutral particle flux: Z
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(64)
   :Unit: 
.. py:attribute:: neut_pg11_file
   
   neutral stress: RR
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(64)
   :Unit: 
.. py:attribute:: neut_pg22_file
   
   neutral stress: TT
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(64)
   :Unit: 
.. py:attribute:: neut_pg33_file
   
   neutral stress: ZZ
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(64)
   :Unit: 
.. py:attribute:: neut_pg23_file
   
   neutral stress: TZ
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(64)
   :Unit: 
.. py:attribute:: neut_pg31_file
   
   neutral stress: ZR
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(64)
   :Unit: 
.. py:attribute:: neut_pg12_file
   
   neutral stress: RT
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(64)
   :Unit: 
.. py:attribute:: neut_jeg1_file
   
   neutral heat flux: R
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(64)
   :Unit: 
.. py:attribute:: neut_jeg2_file
   
   neutral heat flux: T
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(64)
   :Unit: 
.. py:attribute:: neut_jeg3_file
   
   neutral heat flux: Z
   
   :Default: 
   :Dimension: 
   :Group: Ext_neutrals
   :Type: character(64)
   :Unit: 
.. py:attribute:: pnc_histfile
   
   file name
   
   :Default: 
   :Dimension: 
   :Group: PNC_params
   :Type: character(64)
   :Unit: 
.. py:attribute:: pnc_savefile
   
   default pdb filename for saving pnc data
   
   :Default: 
   :Dimension: 
   :Group: PNC_params
   :Type: character(64)
   :Unit: 
.. py:attribute:: pnc_balancefile
   
   file to store diagnositic info for each step
   # test preconditioner with alternate neutrals model
   
   :Default: 
   :Dimension: 
   :Group: PNC_params
   :Type: character(64)
   :Unit: 
.. py:attribute:: pnc_cfparvis
   
   factor for parallel visc. in preconditioner
   
   :Default: 
   :Dimension: (1:31)
   :Group: PNC_params
   :Type: double
   :Unit: 
.. py:attribute:: pnc_cftravis
   
   factor for perp. visc. in preconditioner
   
   :Default: 
   :Dimension: (1:31)
   :Group: PNC_params
   :Type: double
   :Unit: 
.. py:attribute:: pnc_cfni
   
   factor for ni in preconditioner
   
   :Default: 
   :Dimension: (1:31)
   :Group: PNC_params
   :Type: double
   :Unit: 
.. py:attribute:: pnc_cfup
   
   factor for up in preconditioner
   
   :Default: 
   :Dimension: (1:31)
   :Group: PNC_params
   :Type: double
   :Unit: 
.. py:attribute:: ni_pnc
   
   ion density in primary cell (ix,iy) at last pnc step
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: PNC_data
   :Type: double
   :Unit: 1/m**3
.. py:attribute:: up_pnc
   
   parallel velocity in primary cell (ix,iy) at last pnc step
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nisp)
   :Group: PNC_data
   :Type: double
   :Unit: m/s
.. py:attribute:: ti_pnc
   
   ion temperature in primary cell (ix,iy) at last pnc step
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: PNC_data
   :Type: double
   :Unit: J
.. py:attribute:: te_pnc
   
   electron temperature in primary cell (ix,iy) at last pnc step
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: PNC_data
   :Type: double
   :Unit: J
.. py:attribute:: phi_pnc
   
   potential in primary cell (ix,iy) at last pnc step
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: PNC_data
   :Type: double
   :Unit: V
.. py:attribute:: ng_pnc
   
   neutral density in primary cell (ix,iy) at last pnc step
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nfl)
   :Group: PNC_data
   :Type: double
   :Unit: 1/m**3
.. py:attribute:: upg_pnc
   
   parallel neutral velocity in primary cell (ix,iy) at last pnc step
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nfl)
   :Group: PNC_data
   :Type: double
   :Unit: m/s
.. py:attribute:: tg_pnc
   
   neutral temperature in primary cell (ix,iy) at last pnc step
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nfl)
   :Group: PNC_data
   :Type: double
   :Unit: J
.. py:attribute:: sng_pnc
   
   neutral particle source at last pnc step
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nfl)
   :Group: PNC_data
   :Type: double
   :Unit: 1/s
.. py:attribute:: smg_pnc
   
   neutral momentum source at last pnc step
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nfl)
   :Group: PNC_data
   :Type: double
   :Unit: N
.. py:attribute:: seg_pnc
   
   neutral energy source at last pnc step
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nfl)
   :Group: PNC_data
   :Type: double
   :Unit: W
.. py:attribute:: sni_pnc
   
   density source in primary cell (ix,iy) at last pnc step
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nfl,1:nstra)
   :Group: PNC_data
   :Type: double
   :Unit: 1/s
.. py:attribute:: smor_pnc
   
   radial momentum source in primary cell (ix,iy) at last pnc step
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nfl,1:nstra)
   :Group: PNC_data
   :Type: double
   :Unit: N
.. py:attribute:: smophi_pnc
   
   toroidal momentum source in primary cell (ix,iy) at last pnc step
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nfl,1:nstra)
   :Group: PNC_data
   :Type: double
   :Unit: N
.. py:attribute:: smoz_pnc
   
   vertical momentum source in primary cell (ix,iy) at last pnc step
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nfl,1:nstra)
   :Group: PNC_data
   :Type: double
   :Unit: N
.. py:attribute:: sei_pnc
   
   ion energy source in primary cell (ix,iy) at last pnc step
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nstra)
   :Group: PNC_data
   :Type: double
   :Unit: W
.. py:attribute:: see_pnc
   
   electron energy source in primary cell (ix,iy) at last pnc step
   fngx_pnc(0:nx+1,0:ny+1,1:nfl)			_real	[1/s]		#neutral particle flux
   fngy_pnc(0:nx+1,0:ny+1,1:nfl)			_real	[1/s]		#neutral particle flux
   feg_pnc(0:nx+1,0:ny+1,1:3,1:nfl)		_real	[W]			#total neutral heat flux
   vg_pnc(0:nx+1,0:ny+1,1:3,1:nfl)		_real	[m/s]		#neutral velocity
   qg_pnc(0:nx+1,0:ny+1,1:3,1:nfl)		_real	[W/m**2]	#neutral heat flux
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,1:nstra)
   :Group: PNC_data
   :Type: double
   :Unit: W
.. py:attribute:: psorold
   
   unpert. ioniz. sources
   
   :Default: 
   :Dimension: (1:nisp)
   :Group: Save_terms
   :Type: double
   :Unit: part/s
.. py:attribute:: psorxrold
   
   unpert. recom. & cx sources
   
   :Default: 
   :Dimension: (1:nisp)
   :Group: Save_terms
   :Type: double
   :Unit: part/s
.. py:attribute:: msorold
   
   unpert. ioniz. mom. sources
   
   :Default: 
   :Dimension: (1:nisp)
   :Group: Save_terms
   :Type: double
   :Unit: kg-m/s**2
.. py:attribute:: msorxrold
   
   unpert. recom. & cx mom. sources
   
   :Default: 
   :Dimension: (1:nisp)
   :Group: Save_terms
   :Type: double
   :Unit: kg-m/s**2
.. py:attribute:: ylodt
   
   primary variables for ODE's at last output
   
   :Default: 
   :Dimension: (neqmx)
   :Group: Time_dep_nwt
   :Type: double
   :Unit: 
.. py:attribute:: dtoptx
   
   spatial-depend. time step, min. in a cell
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Time_dep_nwt
   :Type: double
   :Unit: 
.. py:attribute:: dtoptv
   
   variable-dependent time step; each var. diff
   
   :Default: 
   :Dimension: (neqmx)
   :Group: Time_dep_nwt
   :Type: double
   :Unit: 
.. py:attribute:: dtuse
   
   time step used based on model_dt value
   
   :Default: 
   :Dimension: (neqmx)
   :Group: Time_dep_nwt
   :Type: double
   :Unit: 
.. py:attribute:: rconds
   
   condition numbers of Jacobians
   
   :Default: 
   :Dimension: (300,ngrid)
   :Group: Condition_number
   :Type: double
   :Unit: 
.. py:attribute:: jac
   
   Nonzero entries of the Jacobian matrix.
   This array, together with jacj and jaci,
   contain the Jacobian in compressed sparse
   row format.
   
   :Default: 
   :Dimension: (nnzmx)
   :Group: Jacobian
   :Type: double
   :Unit: 
.. py:attribute:: jaci
   
   Nonzero structure of Jacobian matrix jac.
   jaci(i+1) - jaci(i) = no. of nonzeros
   in row i of jac.
   
   :Default: 
   :Dimension: (neqp1)
   :Group: Jacobian
   :Type: integer
   :Unit: 
.. py:attribute:: jacj
   
   Column indices of nonzero entries in jac.
   
   :Default: 
   :Dimension: (nnzmx)
   :Group: Jacobian
   :Type: integer
   :Unit: 
.. py:attribute:: rcsc
   
   Nonzero entries of the Jacobian matrix.
   This array, together with jcsc and icsc,
   contain the Jacobian in compressed sparse
   column format.
   
   :Default: 
   :Dimension: (nnzmx)
   :Group: Jacobian_csc
   :Type: double
   :Unit: 
.. py:attribute:: jcsc
   
   Nonzero structure of Jacobian matrix rcsc.
   jcsc(j+1) - jcsc(j) = no. of nonzeros
   in column j of rcsc.
   
   :Default: 
   :Dimension: (neq+1)
   :Group: Jacobian_csc
   :Type: integer
   :Unit: 
.. py:attribute:: icsc
   
   Row indices of nonzero entries in rcsc.
   
   :Default: 
   :Dimension: (nnzmx)
   :Group: Jacobian_csc
   :Type: integer
   :Unit: 
.. py:attribute:: yldot_pert
   
   Perturbed yldot within Jac_calc (diagnostic)
   
   :Default: 
   :Dimension: (neqmx)
   :Group: Jacobian_csc
   :Type: double
   :Unit: 
.. py:attribute:: yldot_unpt
   
   Initial yldot with Jac_calc (diagnostic)
   
   :Default: 
   :Dimension: (neqmx)
   :Group: Jacobian_csc
   :Type: double
   :Unit: 
.. py:attribute:: jac1
   
   Nonzero elements of Jacobian
   
   :Default: 
   :Dimension: (nnz1mx)
   :Group: Jacobian_part
   :Type: double
   :Unit: 
.. py:attribute:: ia1
   
   Row indices of elements in jac1, or
   nonzero structure of jac1 in csr format
   
   :Default: 
   :Dimension: (nnz1mx)
   :Group: Jacobian_part
   :Type: integer
   :Unit: 
.. py:attribute:: ja1
   
   Column indices of elements in jac1
   
   :Default: 
   :Dimension: (nnz1mx)
   :Group: Jacobian_part
   :Type: integer
   :Unit: 
.. py:attribute:: perm
   
   Integer array containing the permutation
   used in reordering the rows and columns of
   the Jacobian matrix.
   
   :Default: 
   :Dimension: (neq)
   :Group: Jacreorder
   :Type: integer
   :Unit: 
.. py:attribute:: qperm
   
   Integer array holding the inverse of the
   permutation in array perm.
   
   :Default: 
   :Dimension: (neq)
   :Group: Jacreorder
   :Type: integer
   :Unit: 
.. py:attribute:: levels
   
   Work array used by the bfs reordering
   subroutine. See subroutine bfs for
   more details.
   
   :Default: 
   :Dimension: (neq)
   :Group: Jacreorder
   :Type: integer
   :Unit: 
.. py:attribute:: mask
   
   Work array used by the bfs reordering
   subroutine. See bfs subroutine.
   
   :Default: 
   :Dimension: (neq)
   :Group: Jacreorder
   :Type: integer
   :Unit: 
.. py:attribute:: jacfull
   
   
   
   :Default: 
   :Dimension: (neq,neq)
   :Group: Jacobian_full
   :Type: double
   :Unit: 
.. py:attribute:: premeth
   
   type of preconditioning used in the
   linear iteration:
   ='banded' means use full banded jacobian as
   preconditioner.
   ='ilut' means use ilut preconditioning.
   ='inel' means use INEL ILU preconditioning
   
   :Default: 
   :Dimension: 
   :Group: Preconditioning
   :Type: character(8)
   :Unit: 
.. py:attribute:: adiag
   
   diagonals of the Jacobian matrix
   
   :Default: 
   :Dimension: (neq,ndiagmx)
   :Group: Nonzero_diagonals
   :Type: double
   :Unit: 
.. py:attribute:: siginel
   
   work array used by INEL precond5
   
   :Default: 
   :Dimension: (neq)
   :Group: Nonzero_diagonals
   :Type: double
   :Unit: 
.. py:attribute:: fmuinel
   
   work array used by INEL precond5
   
   :Default: 
   :Dimension: (neq)
   :Group: Nonzero_diagonals
   :Type: double
   :Unit: 
.. py:attribute:: rwkd
   
   work array used by cdiagsrt
   
   :Default: 
   :Dimension: (ndiagmx)
   :Group: Nonzero_diagonals
   :Type: double
   :Unit: 
.. py:attribute:: iwkd1
   
   number of nonzeros in each diagonal
   
   :Default: 
   :Dimension: (2*neq-1)
   :Group: Nonzero_diagonals
   :Type: integer
   :Unit: 
.. py:attribute:: iwkd2
   
   work array used by cdiagsrt
   
   :Default: 
   :Dimension: (ndiagmx)
   :Group: Nonzero_diagonals
   :Type: integer
   :Unit: 
.. py:attribute:: ngscal
   
   ratio of initial gas density to ion dens
   
   :Default: 
   :Dimension: (6)
   :Group: UEint
   :Type: double
   :Unit: 
.. py:attribute:: nibeg
   
   initial ion density
   
   :Default: 
   :Dimension: (1:31)
   :Group: UEint
   :Type: double
   :Unit: 
.. py:attribute:: minu
   
   ion mass in units of proton mass (AMU)
   
   :Default: 
   :Dimension: (1:31)
   :Group: UEint
   :Type: double
   :Unit: 
.. py:attribute:: ziin
   
   ion charge read in, used to reset zi in
   group Compla which gets erased on gallot
   
   :Default: 
   :Dimension: (1:31)
   :Group: UEint
   :Type: double
   :Unit: 
.. py:attribute:: znuclin
   
   total nuclear charge of ion (i.d. isotope)
   
   :Default: 
   :Dimension: (1:31)
   :Group: UEint
   :Type: integer
   :Unit: 
.. py:attribute:: pyrestart_file
   
   Python file that can also be used to restart
   
   :Default: 
   :Dimension: 
   :Group: UEint
   :Type: character(80)
   :Unit: 
.. py:attribute:: uedge_savefile
   
   default pdb filename for saving uedge data
   
   :Default: 
   :Dimension: 
   :Group: Interp
   :Type: character(64)
   :Unit: 
.. py:attribute:: ixlbo
   
   prev. grid value for ixlb
   
   :Default: 
   :Dimension: (1:nxpt)
   :Group: Interp
   :Type: integer
   :Unit: 
.. py:attribute:: ixpt1o
   
   prev. grid value for ixpt1
   
   :Default: 
   :Dimension: (1:nxpt)
   :Group: Interp
   :Type: integer
   :Unit: 
.. py:attribute:: ixpt2o
   
   prev. grid value for ixpt2
   
   :Default: 
   :Dimension: (1:nxpt)
   :Group: Interp
   :Type: integer
   :Unit: 
.. py:attribute:: ixrbo
   
   prev. grid value for ixrb
   
   :Default: 
   :Dimension: (1:nxpt)
   :Group: Interp
   :Type: integer
   :Unit: 
.. py:attribute:: ixst
   
   starting ix for 6 poloid region interp
   
   :Default: 
   :Dimension: (1:6)
   :Group: Interp
   :Type: integer
   :Unit: 
.. py:attribute:: ixsto
   
   value of ixst on previous grid
   
   :Default: 
   :Dimension: (1:6)
   :Group: Interp
   :Type: integer
   :Unit: 
.. py:attribute:: ixend
   
   end ix for 6 poloid region interp
   
   :Default: 
   :Dimension: (1:6)
   :Group: Interp
   :Type: integer
   :Unit: 
.. py:attribute:: ixendo
   
   value of ixend on previous grid
   
   :Default: 
   :Dimension: (1:6)
   :Group: Interp
   :Type: integer
   :Unit: 
.. py:attribute:: xnrmo
   
   norm. x-grd; old x-grid, old y-grid
   
   :Default: 
   :Dimension: (0:nxold+1,0:nyold+1)
   :Group: Interp
   :Type: double
   :Unit: 
.. py:attribute:: xvnrmo
   
   norm. xv-grd; old x-grid, old y-grid
   
   :Default: 
   :Dimension: (0:nxold+1,0:nyold+1)
   :Group: Interp
   :Type: double
   :Unit: 
.. py:attribute:: xnrmox
   
   norm. x-grd;nxold grd interp. to new ny
   
   :Default: 
   :Dimension: (0:nxold+1,0:ny+1)
   :Group: Interp
   :Type: double
   :Unit: 
.. py:attribute:: xvnrmox
   
   norm. xv-grd;nxold grd interp.to new ny
   
   :Default: 
   :Dimension: (0:nxold+1,0:ny+1)
   :Group: Interp
   :Type: double
   :Unit: 
.. py:attribute:: xnrmnx
   
   norm. x-grd; second intermed. grid
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Interp
   :Type: double
   :Unit: 
.. py:attribute:: xvnrmnx
   
   norm. xv-grd; second intermed. grid
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Interp
   :Type: double
   :Unit: 
.. py:attribute:: ynrmo
   
   norm. y-grd; old x-grid, old y-grid
   
   :Default: 
   :Dimension: (0:nxold+1,0:nyold+1)
   :Group: Interp
   :Type: double
   :Unit: 
.. py:attribute:: yvnrmo
   
   norm. yv-grd; old x-grid, old y-grid
   
   :Default: 
   :Dimension: (0:nxold+1,0:nyold+1)
   :Group: Interp
   :Type: double
   :Unit: 
.. py:attribute:: ynrmox
   
   norm. y-grd; old x-grid, new y-grid
   
   :Default: 
   :Dimension: (0:nxold+1,0:ny+1)
   :Group: Interp
   :Type: double
   :Unit: 
.. py:attribute:: yvnrmox
   
   norm. yv-grd; old xv-grid, new y-grid
   
   :Default: 
   :Dimension: (0:nxold+1,0:ny+1)
   :Group: Interp
   :Type: double
   :Unit: 
.. py:attribute:: ynrmnx
   
   norm. y-grd; second intermed. grid
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Interp
   :Type: double
   :Unit: 
.. py:attribute:: yvnrmnx
   
   norm. yv-grd; second intermed. grid
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Interp
   :Type: double
   :Unit: 
.. py:attribute:: wrkint
   
   wrk array; vars on old x, new y grid
   
   :Default: 
   :Dimension: (0:nxold+1,0:ny+1)
   :Group: Interp
   :Type: double
   :Unit: 
.. py:attribute:: wrkint2
   
   wrk array; vars on second interm. grid
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Interp
   :Type: double
   :Unit: 
.. py:attribute:: ixmg
   
   ix index used for (ixo,iy) pt.
   
   :Default: 
   :Dimension: (0:nxold+1,0:ny+1)
   :Group: Interp
   :Type: integer
   :Unit: 
.. py:attribute:: iyomg
   
   iyo index used for (ixo,iy) pt.
   
   :Default: 
   :Dimension: (0:nxold+1,0:ny+1)
   :Group: Interp
   :Type: integer
   :Unit: 
.. py:attribute:: ixvmg
   
   ixv index used for (ixvo,iy) pt.
   
   :Default: 
   :Dimension: (0:nxold+1,0:ny+1)
   :Group: Interp
   :Type: integer
   :Unit: 
.. py:attribute:: iyvomg
   
   iyvo index used for (ixvo,iy) pt.
   
   :Default: 
   :Dimension: (0:nxold+1,0:ny+1)
   :Group: Interp
   :Type: integer
   :Unit: 
.. py:attribute:: ix2g
   
   ix index for sec. interm. (ix,iy) pt.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Interp
   :Type: integer
   :Unit: 
.. py:attribute:: iy2g
   
   iy index for sec. interm. (ixo,iy) pt.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Interp
   :Type: integer
   :Unit: 
.. py:attribute:: ixv2g
   
   ixv index for sec. interm.(ixvo,iy) pt.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Interp
   :Type: integer
   :Unit: 
.. py:attribute:: iyv2g
   
   iyv index for sec.interm.(ixvo,iy) pt.
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Interp
   :Type: integer
   :Unit: 
.. py:attribute:: nis
   
   ion dens at last success. calc
   
   :Default: 
   :Dimension: (0:nxold+1,0:nyold+1,1:nisp)
   :Group: Interp
   :Type: double
   :Unit: m^-3
.. py:attribute:: tes
   
   elec. temp at last success. calc
   
   :Default: 
   :Dimension: (0:nxold+1,0:nyold+1)
   :Group: Interp
   :Type: double
   :Unit: J
.. py:attribute:: tis
   
   ion temp at last success. calc
   
   :Default: 
   :Dimension: (0:nxold+1,0:nyold+1)
   :Group: Interp
   :Type: double
   :Unit: J
.. py:attribute:: tgs
   
   gas temp at last success. calc
   
   :Default: 
   :Dimension: (0:nxold+1,0:nyold+1,1:ngsp)
   :Group: Interp
   :Type: double
   :Unit: J
.. py:attribute:: phis
   
   potential at last success. calc
   
   :Default: 
   :Dimension: (0:nxold+1,0:nyold+1)
   :Group: Interp
   :Type: double
   :Unit: V
.. py:attribute:: ups
   
   parall. vel at last success. calc
   
   :Default: 
   :Dimension: (0:nxold+1,0:nyold+1,1:nisp)
   :Group: Interp
   :Type: double
   :Unit: m/s
.. py:attribute:: ngs
   
   gas dens at last success. calc.
   
   :Default: 
   :Dimension: (0:nxold+1,0:nyold+1,1:ngsp)
   :Group: Interp
   :Type: double
   :Unit: m^-3
.. py:attribute:: afracs
   
   rel. imp. frac at last succ. calc
   
   :Default: 
   :Dimension: (0:nxold+1,0:nyold+1)
   :Group: Interp
   :Type: double
   :Unit: 
.. py:attribute:: nisg
   
   global array for nis
   
   :Default: 
   :Dimension: (0:nxoldg+1,0:nyoldg+1,1:nisp)
   :Group: Global_vars
   :Type: double
   :Unit: m^-3
.. py:attribute:: tesg
   
   global array for tes
   
   :Default: 
   :Dimension: (0:nxoldg+1,0:nyoldg+1)
   :Group: Global_vars
   :Type: double
   :Unit: J
.. py:attribute:: tisg
   
   global array for tis
   
   :Default: 
   :Dimension: (0:nxoldg+1,0:nyoldg+1)
   :Group: Global_vars
   :Type: double
   :Unit: J
.. py:attribute:: tgsg
   
   global array for tgs
   
   :Default: 
   :Dimension: (0:nxoldg+1,0:nyoldg+1)
   :Group: Global_vars
   :Type: double
   :Unit: J
.. py:attribute:: phisg
   
   global array for phis
   
   :Default: 
   :Dimension: (0:nxoldg+1,0:nyoldg+1)
   :Group: Global_vars
   :Type: double
   :Unit: V
.. py:attribute:: upsg
   
   global array for ups
   
   :Default: 
   :Dimension: (0:nxoldg+1,0:nyoldg+1,1:nisp)
   :Group: Global_vars
   :Type: double
   :Unit: m/s
.. py:attribute:: ngsg
   
   global array for ngs
   
   :Default: 
   :Dimension: (0:nxoldg+1,0:nyoldg+1,1:ngsp)
   :Group: Global_vars
   :Type: double
   :Unit: m^-3
.. py:attribute:: afracsg
   
   global array for afracs
   
   :Default: 
   :Dimension: (0:nxoldg+1,0:nyoldg+1)
   :Group: Global_vars
   :Type: double
   :Unit: 
.. py:attribute:: ipassin
   
   integer input variables to be passed
   
   :Default: 
   :Dimension: (1:100)
   :Group: Global_input
   :Type: integer
   :Unit: 
.. py:attribute:: rpassin
   
   real input variables to be passed
   
   :Default: 
   :Dimension: (1:100)
   :Group: Global_input
   :Type: double
   :Unit: 
.. py:attribute:: cpassin
   
   character input variables to be passed
   
   :Default: 
   :Dimension: (1:30)
   :Group: Global_input
   :Type: character(8)
   :Unit: 
.. py:attribute:: ndleg
   
   number of x-domains in nxleg regions
   
   :Default: 
   :Dimension: (1:10,1:2)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: ndxcore
   
   number of x-domains in nxcore(,1:2) regions
   
   :Default: 
   :Dimension: (1:10)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: ndycore
   
   number of y-domains in core
   
   :Default: 
   :Dimension: (1:10)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: ndysol
   
   number of y-domains in sol
   
   :Default: 
   :Dimension: (1:10)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: idxpt
   
   PF/core domains with up touching X-point
   
   :Default: 
   :Dimension: (1:2)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: ixmin
   
   min global ix for given domain
   
   :Default: 
   :Dimension: (32)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: ixmax
   
   max global ix for given domain
   
   :Default: 
   :Dimension: (32)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: iymin
   
   min global iy for given domain
   
   :Default: 
   :Dimension: (32)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: iymax
   
   max global iy for given domain
   
   :Default: 
   :Dimension: (32)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: ixmnbcg
   
   B.C. type at ix=ixmin bdry;=0 inter.,=1 ex.
   
   :Default: 
   :Dimension: (32)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: ixmxbcg
   
   B.C. type at ix=ixmax bdry;=0 inter.,=1 ex.
   
   :Default: 
   :Dimension: (32)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: iymnbcg
   
   B.C. type at iy=iymin bdry;=0 inter.,=1 ex.
   
   :Default: 
   :Dimension: (32)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: iymxbcg
   
   B.C. type at iy=iymax bdry;=0 inter.,=1 ex.
   
   :Default: 
   :Dimension: (32)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: ncell
   
   number of cells for given domain
   
   :Default: 
   :Dimension: (32)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: idxp1g
   
   domain to the right of given domain (ix+1)
   
   :Default: 
   :Dimension: (32)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: idxm1g
   
   domain to the left of given domain (ix-1)
   
   :Default: 
   :Dimension: (32)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: idyp1g
   
   domain above given domain (iy+1)
   
   :Default: 
   :Dimension: (32)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: idym1g
   
   domain below given domain (iy-1)
   
   :Default: 
   :Dimension: (32)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: idcorng
   
   domains touching corners; from lower left,
   numbering as in rm,zm: (l,r bot=1,2; top=3,4)
   
   :Default: 
   :Dimension: (32,1:4)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: ixpt1g
   
   ixpt1 for a given domain
   
   :Default: 
   :Dimension: (32)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: ixpt2g
   
   ixpt2 for a given domain
   
   :Default: 
   :Dimension: (32)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: iysptrxg
   
   iysptrx for a given domain
   
   :Default: 
   :Dimension: (32)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: vrsend
   
   real array used for passing global data via MPI
   
   :Default: 
   :Dimension: (nvrsend)
   :Group: Indices_domain_dcg
   :Type: double
   :Unit: 
.. py:attribute:: visend
   
   int array used for passing global data via MPI
   
   :Default: 
   :Dimension: (nvisend)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: neq_locg
   
   number of vars per domain
   
   :Default: 
   :Dimension: (32)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: ispwrbc
   
   =1 for core pwr flux BC if corresp to ixpt2g
   
   :Default: 
   :Dimension: (32)
   :Group: Indices_domain_dcg
   :Type: integer
   :Unit: 
.. py:attribute:: ivcum
   
   counter to build yl-local to yl-global map
   
   :Default: 
   :Dimension: (32)
   :Group: Indices_loc_glob_map
   :Type: integer
   :Unit: 
.. py:attribute:: ivloc2sdg
   
   map loc-var to glob-var, single domain
   
   :Default: 
   :Dimension: (neqmx,32)
   :Group: Indices_loc_glob_map
   :Type: integer
   :Unit: 
.. py:attribute:: ivloc2mdg
   
   map loc-var to glob-var, mult domain
   
   :Default: 
   :Dimension: (neqmx,32)
   :Group: Indices_loc_glob_map
   :Type: integer
   :Unit: 
.. py:attribute:: ivl2gstnl
   
   1st arg loc-eqn number;
   2nd arg poss Jac vars - global-mp; 3rd arg domain
   
   :Default: 
   :Dimension: (neq_locgmx,9*numvar,32)
   :Group: Indices_loc_glob_map
   :Type: integer
   :Unit: 
.. py:attribute:: iellast
   
   last meaningful entry into ivl2gstnl
   
   :Default: 
   :Dimension: (neqmx,32)
   :Group: Indices_loc_glob_map
   :Type: integer
   :Unit: 
.. py:attribute:: idcorn
   
   domains touching corners; from lower left,
   numbering as in rm,zm: (l,r bot=1,2; top=3,4)
   
   :Default: 
   :Dimension: (1:4)
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: iv_totbdy
   
   number of elems. in bdry messages vrsendl
   
   :Default: 
   :Dimension: (1:8)
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: typebdyi
   
   mpi tags for bdry iv_totbdy along edges
   
   :Default: 
   :Dimension: (1:4)
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: typecni
   
   mpi tags for bdry iv_totbdy at corners
   
   :Default: 
   :Dimension: (1:4)
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: typebdy
   
   mpi tags for bdry vrsendl along edges
   
   :Default: 
   :Dimension: (1:4)
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: typecn
   
   mpi tags for bdry vrsendl at corners
   
   :Default: 
   :Dimension: (1:4)
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: vrsendl
   
   real array used for passing local data via MPI
   
   :Default: 
   :Dimension: (nvrsendl)
   :Group: Indices_domain_dcl
   :Type: double
   :Unit: 
.. py:attribute:: visendl
   
   int array used for passing local data via MPI
   
   :Default: 
   :Dimension: (nvisendl)
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: ivloc2sdgl
   
   maps loc-var to glob-var, single domain
   
   :Default: 
   :Dimension: (nvisendl)
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: ivloc2mdgl
   
   maps loc-var to glob-var, mult domain
   
   :Default: 
   :Dimension: (nvisendl)
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: ivl2gstnll
   
   1st arg loc-eqn number;
   2nd arg poss Jac vars-global-mp
   
   :Default: 
   :Dimension: (neq_locl,9*numvarl)
   :Group: Indices_domain_dcl
   :Type: integer
   :Unit: 
.. py:attribute:: ylold
   
   
   
   :Default: 
   :Dimension: (neqmx)
   :Group: Jacaux
   :Type: double
   :Unit: 
.. py:attribute:: yldot1
   
   
   
   :Default: 
   :Dimension: (neqmx)
   :Group: Jacaux
   :Type: double
   :Unit: 
.. py:attribute:: yldot0
   
   
   
   :Default: 
   :Dimension: (neqmx)
   :Group: Jacaux
   :Type: double
   :Unit: 
.. py:attribute:: fnormnw
   
   
   
   :Default: 
   :Dimension: (neqmx)
   :Group: Jacaux
   :Type: double
   :Unit: 
.. py:attribute:: ysave
   
   last two yl's in Newton (1,) most recent
   
   :Default: 
   :Dimension: (2,neqmx)
   :Group: Newtaux
   :Type: double
   :Unit: 
.. py:attribute:: ycor
   
   
   
   :Default: 
   :Dimension: (neqmx)
   :Group: Newtaux
   :Type: double
   :Unit: 
.. py:attribute:: saux2
   
   maximum update allowed in Newton
   
   :Default: 
   :Dimension: (nmaxnewt)
   :Group: Newtaux
   :Type: double
   :Unit: 
.. py:attribute:: sumf
   
   ave value of right-hand-sides after Newton
   
   :Default: 
   :Dimension: (0:nmaxnewt)
   :Group: Newtaux
   :Type: double
   :Unit: 
.. py:attribute:: irwd
   
   
   
   :Default: 
   :Dimension: (nmaxnewt,2)
   :Group: Newtaux
   :Type: integer
   :Unit: 
.. py:attribute:: rwdmax
   
   
   
   :Default: 
   :Dimension: (nmaxnewt,2)
   :Group: Newtaux
   :Type: double
   :Unit: 
.. py:attribute:: rwdmin
   
   
   
   :Default: 
   :Dimension: (nmaxnewt,2)
   :Group: Newtaux
   :Type: double
   :Unit: 
.. py:attribute:: sumnew1
   
   average change in variables for Newton iter.
   
   :Default: 
   :Dimension: (nmaxnewt)
   :Group: Newtaux
   :Type: double
   :Unit: 
.. py:attribute:: sumr1dy
   
   
   
   :Default: 
   :Dimension: (nmaxnewt)
   :Group: Newtaux
   :Type: double
   :Unit: 
.. py:attribute:: rcn
   
   radial position of density cell
   
   :Default: 
   :Dimension: (0:nxm+1,0:nym+1)
   :Group: RZ_cell_info
   :Type: double
   :Unit: m
.. py:attribute:: zcn
   
   vertical position of density cell
   
   :Default: 
   :Dimension: (0:nxm+1,0:nym+1)
   :Group: RZ_cell_info
   :Type: double
   :Unit: m
.. py:attribute:: rfn
   
   radial position of density face
   
   :Default: 
   :Dimension: (-1:nxm+1,-1:nym+1)
   :Group: RZ_cell_info
   :Type: double
   :Unit: m
.. py:attribute:: zfn
   
   vertical position of density face
   
   :Default: 
   :Dimension: (-1:nxm+1,-1:nym+1)
   :Group: RZ_cell_info
   :Type: double
   :Unit: m
.. py:attribute:: rcv
   
   radial position of velocity cell
   
   :Default: 
   :Dimension: (0:nxm+1,0:nym+1)
   :Group: RZ_cell_info
   :Type: double
   :Unit: m
.. py:attribute:: zcv
   
   vertical position of velocity cell
   
   :Default: 
   :Dimension: (0:nxm+1,0:nym+1)
   :Group: RZ_cell_info
   :Type: double
   :Unit: m
.. py:attribute:: rfv
   
   radial position of velocity face
   
   :Default: 
   :Dimension: (0:nxm+2,0:nym+2)
   :Group: RZ_cell_info
   :Type: double
   :Unit: m
.. py:attribute:: zfv
   
   vertical position of velocity face
   
   :Default: 
   :Dimension: (0:nxm+2,0:nym+2)
   :Group: RZ_cell_info
   :Type: double
   :Unit: m
.. py:attribute:: nzloc
   
   imp. dens. for each Z at one grid cell
   
   :Default: 
   :Dimension: (0:nzspmx)
   :Group: Imprad
   :Type: double
   :Unit: /m**3
.. py:attribute:: impradloc
   
   rad. power loss density for each Z at one grid cell
   
   :Default: 
   :Dimension: (0:nzspmx)
   :Group: Imprad
   :Type: double
   :Unit: Watts/m**3
.. py:attribute:: pwrzec
   
   elec energy loss via impurities at cell-cntr
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Imprad
   :Type: double
   :Unit: Watts/m**3
.. py:attribute:: pwrze
   
   elec energy loss via impurities; cell-ave
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Imprad
   :Type: double
   :Unit: Watts/m**3
.. py:attribute:: pradc
   
   cell ctr total impurity radiation
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Imprad
   :Type: double
   :Unit: Watts/m**3
.. py:attribute:: pradcff
   
   cell ctr impurity radiation (fixed-fraction)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Imprad
   :Type: double
   :Unit: Watts/m**3
.. py:attribute:: prad
   
   cell ave total impurity radiation
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Imprad
   :Type: double
   :Unit: Watts/m**3
.. py:attribute:: pradzc
   
   cell ctr imp rad due to each imp. ch. state
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,0:nzspmx,1:ngsp-1)
   :Group: Imprad
   :Type: double
   :Unit: Watts/m**3
.. py:attribute:: pradz
   
   cell ave imp rad due to each imp. ch. state
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1,0:nzspmx,1:ngsp-1)
   :Group: Imprad
   :Type: double
   :Unit: Watts/m**3
.. py:attribute:: na
   
   atomic density of impurity (=afrac*ne)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Imprad
   :Type: double
   :Unit: /m**3
.. py:attribute:: ntau
   
   confinement parameter for impurity (=atau*ne)
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Imprad
   :Type: double
   :Unit: sec/m**3
.. py:attribute:: nratio
   
   ratio of neutrals to electrons
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Imprad
   :Type: double
   :Unit: 
.. py:attribute:: afrac
   
   atomic impur conc; set internally to afracs
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Imprad
   :Type: double
   :Unit: 
.. py:attribute:: atau
   
   lifetime of impurity
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Imprad
   :Type: double
   :Unit: sec
.. py:attribute:: tau1
   
   time to escape to inboard divertor plate
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Imprad
   :Type: double
   :Unit: 
.. py:attribute:: tau2
   
   time to escape to outboard divertor plate
   
   :Default: 
   :Dimension: (0:nx+1,0:ny+1)
   :Group: Imprad
   :Type: double
   :Unit: 
.. py:attribute:: fnzysi
   
   profiles along inner wall
   
   :Default: 
   :Dimension: (0:nx+1,nzspt)
   :Group: Impurity_source_flux
   :Type: double
   :Unit: 
.. py:attribute:: fnzyso
   
   profiles along outer wall
   
   :Default: 
   :Dimension: (0:nx+1,nzspt)
   :Group: Impurity_source_flux
   :Type: double
   :Unit: 
.. py:attribute:: natomic
   
   maximum charge state of each isotope
   
   :Default: 
   :Dimension: (1:5)
   :Group: Reduced_ion_interface
   :Type: integer
   :Unit: 
.. py:attribute:: amu
   
   atomic mass, relative to proton
   
   :Default: 
   :Dimension: (1:misotope)
   :Group: Reduced_ion_interface
   :Type: double
   :Unit: none
.. py:attribute:: tempa
   
   temperature
   
   :Default: 
   :Dimension: (1:misotope)
   :Group: Reduced_ion_interface
   :Type: double
   :Unit: J
.. py:attribute:: qneut
   
   parallel heat flux of neutral
   
   :Default: 
   :Dimension: (1:misotope)
   :Group: Reduced_ion_interface
   :Type: double
   :Unit: J/m**2-s
.. py:attribute:: uneut
   
   parallel flow speed of neutral
   
   :Default: 
   :Dimension: (1:misotope)
   :Group: Reduced_ion_interface
   :Type: double
   :Unit: m/s
.. py:attribute:: den
   
   density
   
   :Default: 
   :Dimension: (1:misotope,0:nchstate)
   :Group: Reduced_ion_interface
   :Type: double
   :Unit: 1/m**3
.. py:attribute:: gradp
   
   parallel pressure grad
   
   :Default: 
   :Dimension: (1:misotope,1:nchstate)
   :Group: Reduced_ion_interface
   :Type: double
   :Unit: J/m**4
.. py:attribute:: gradt
   
   parallel temp gradient
   
   :Default: 
   :Dimension: (1:misotope,1:nchstate)
   :Group: Reduced_ion_interface
   :Type: double
   :Unit: J/m**4
.. py:attribute:: friction
   
   parallel friction force
   
   :Default: 
   :Dimension: (1:misotope,1:nchstate)
   :Group: Reduced_ion_interface
   :Type: double
   :Unit: J/m**4
.. py:attribute:: friccomp
   
   par friction components
   friccomp(,,1)~ upi-upj
   friccomp(,,2)~ qcond
   friccomp(,,3)~ h;higher mom
   friccomp(,,4)~ caplam;elec?
   friccomp(,,5)~ ioniz/recomb
   
   :Default: 
   :Dimension: (1:misotope,1:nchstate,1:5)
   :Group: Reduced_ion_interface
   :Type: double
   :Unit: J/m**4
.. py:attribute:: nuion
   
   ionization rate
   
   :Default: 
   :Dimension: (1:misotope,0:nchstate)
   :Group: Reduced_ion_interface
   :Type: double
   :Unit: 1/s
.. py:attribute:: nurec
   
   recombination rate
   
   :Default: 
   :Dimension: (1:misotope,1:nchstate)
   :Group: Reduced_ion_interface
   :Type: double
   :Unit: 1/s
.. py:attribute:: qcond
   
   parallel heat flux
   
   :Default: 
   :Dimension: (1:misotope,1:nchstate)
   :Group: Reduced_ion_interface
   :Type: double
   :Unit: J/m**2-s
.. py:attribute:: ucond
   
   parallel flow speed
   
   :Default: 
   :Dimension: (1:misotope,1:nchstate)
   :Group: Reduced_ion_interface
   :Type: double
   :Unit: m/s
.. py:attribute:: dztot
   
   total local isotope density
   
   :Default: 
   :Dimension: (1:misotope)
   :Group: Reduced_ion_interface
   :Type: double
   :Unit: 1/m**3
.. py:attribute:: iwork
   
   integer work array
   
   :Default: 
   :Dimension: (liw)
   :Group: Solver_work_arrays
   :Type: integer
   :Unit: 
.. py:attribute:: rwork
   
   real work array
   
   :Default: 
   :Dimension: (lrw)
   :Group: Solver_work_arrays
   :Type: double
   :Unit: 
.. py:attribute:: iwwp
   
   integer work array
   
   :Default: 
   :Dimension: (liwp)
   :Group: Jac_work_arrays
   :Type: integer
   :Unit: 
.. py:attribute:: wwp
   
   real work array
   
   :Default: 
   :Dimension: (lwp)
   :Group: Jac_work_arrays
   :Type: double
   :Unit: 
.. py:attribute:: rwk1
   
   
   
   :Default: 
   :Dimension: (neq+1)
   :Group: Temporary_work_arrays
   :Type: double
   :Unit: 
.. py:attribute:: rwk2
   
   
   
   :Default: 
   :Dimension: (neq+1)
   :Group: Temporary_work_arrays
   :Type: double
   :Unit: 
.. py:attribute:: iwk1
   
   
   
   :Default: 
   :Dimension: (neq+1)
   :Group: Temporary_work_arrays
   :Type: integer
   :Unit: 
.. py:attribute:: iwk2
   
   
   
   :Default: 
   :Dimension: (neq+1)
   :Group: Temporary_work_arrays
   :Type: integer
   :Unit: 
.. py:attribute:: iwk3
   
   
   
   :Default: 
   :Dimension: (neq+1)
   :Group: Temporary_work_arrays
   :Type: integer
   :Unit: 
.. py:attribute:: nezag
   
   Zagorski's electron density
   
   :Default: 
   :Dimension: (imx+1,imy+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: nizag
   
   Zagorski's hydrogen ion density
   
   :Default: 
   :Dimension: (imx+1,imy+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: tezag
   
   Zagorski's electron temperature
   
   :Default: 
   :Dimension: (imx+1,imy+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: tizag
   
   Zagorski's electron temperature
   
   :Default: 
   :Dimension: (imx+1,imy+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: vizag
   
   Zagorski's velocity 1
   
   :Default: 
   :Dimension: (imx+1,imy+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: uizag
   
   Zagorski's velocity 2
   
   :Default: 
   :Dimension: (imx+1,imy+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: v0zag
   
   Zagorski's velocity 3
   
   :Default: 
   :Dimension: (imx+1,imy+1,5)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: u0zag
   
   Zagorski's velocity 4
   
   :Default: 
   :Dimension: (imx+1,imy+1,5)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: zeffzag
   
   Zagorski's Zeff
   
   :Default: 
   :Dimension: (imx+1,imy+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: elfzag
   
   Zagorski's elf
   
   :Default: 
   :Dimension: (imx+1,imy+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: vezag
   
   Zagorski's ve
   
   :Default: 
   :Dimension: (imx+1,imy+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: uezag
   
   Zagorski's ue
   
   :Default: 
   :Dimension: (imx+1,imy+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: curxzag
   
   Zagorski's curx
   
   :Default: 
   :Dimension: (imx+1,imy+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: curyzag
   
   Zagorski's cury
   
   :Default: 
   :Dimension: (imx+1,imy+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: t0zag
   
   Zagorski's t0
   
   :Default: 
   :Dimension: (imx+1,imy+1,5)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: n0zag
   
   Zagorski's n0
   
   :Default: 
   :Dimension: (imx+1,imy+1,5)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: tz0zag
   
   Zagorski's tz0
   IMPUR common block
   
   :Default: 
   :Dimension: (imx+1,imy+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: nzzag
   
   Zagorski's nz
   
   :Default: 
   :Dimension: (imx+1,imy+1,lnst+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: vzzag
   
   Zagorski's vz
   
   :Default: 
   :Dimension: (imx+1,imy+1,lnst+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: uzzag
   
   Zagorski's uz
   GEOMETRY common block
   
   :Default: 
   :Dimension: (imx+1,imy+1,lnst+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: rxzag
   
   Zagorski's rx
   
   :Default: 
   :Dimension: (imx+1,imy+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: ryzag
   
   Zagorski's ry
   
   :Default: 
   :Dimension: (imx+1,imy+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: ggzag
   
   Zagorski's gg
   
   :Default: 
   :Dimension: (imx+1,imy+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: bratiozag
   
   Zagorski's bratio
   
   :Default: 
   :Dimension: (imx+1,imy+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: viparzag
   
   Zagorski's vipar
   
   :Default: 
   :Dimension: (imx+1,imy+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: veparzag
   
   Zagorski's vepar
   
   :Default: 
   :Dimension: (imx+1,imy+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: vzparzag
   
   Zagorski's vzpar
   NET common block
   
   :Default: 
   :Dimension: (imx+1,imy+1,lnst+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: xzag
   
   Zagorski's x
   
   :Default: 
   :Dimension: (imx+2)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: yzag
   
   Zagorski's y
   WARIANT common block
   
   :Default: 
   :Dimension: (imy+2)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: fe
   
   
   
   :Default: 
   :Dimension: (imx+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: fi
   
   
   
   :Default: 
   :Dimension: (imx+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: sdod
   
   
   
   :Default: 
   :Dimension: (imx+1,imy+1,lnst+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: yielh
   
   
   
   :Default: 
   :Dimension: (imx+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: yielz
   
   
   
   :Default: 
   :Dimension: (imx+1,lnst+1)
   :Group: Zag_output
   :Type: double
   :Unit: 
.. py:attribute:: uedge_ver
   
   
   
   :Default: 
   :Dimension: 
   :Group: Ident_vars
   :Type: character(80)
   :Unit: 
.. py:attribute:: uedge_date
   
   
   
   :Default: 
   :Dimension: 
   :Group: Ident_vars
   :Type: character(80)
   :Unit: 
.. py:attribute:: logfname
   
   name of the log file to which to write
   
   :Default: 
   :Dimension: 
   :Group: Logging
   :Type: character(64)
   :Unit: 
.. py:attribute:: xcz
   
   
   
   :Default: 
   :Dimension: (1:nxx)
   :Group: Convdiffeqn
   :Type: double
   :Unit: 
.. py:attribute:: xfz
   
   
   
   :Default: 
   :Dimension: (1:nxx)
   :Group: Convdiffeqn
   :Type: double
   :Unit: 
.. py:attribute:: vrz
   
   
   
   :Default: 
   :Dimension: (1:nxx)
   :Group: Convdiffeqn
   :Type: double
   :Unit: 
.. py:attribute:: drz
   
   
   
   :Default: 
   :Dimension: (1:nxx)
   :Group: Convdiffeqn
   :Type: double
   :Unit: 
.. py:attribute:: dens
   
   
   
   :Default: 
   :Dimension: (1:nxx)
   :Group: Convdiffeqn
   :Type: double
   :Unit: 
.. py:attribute:: vrhs
   
   
   
   :Default: 
   :Dimension: (1:nxx)
   :Group: Convdiffeqn
   :Type: double
   :Unit: 
.. py:attribute:: drhs
   
   
   
   :Default: 
   :Dimension: (1:nxx)
   :Group: Convdiffeqn
   :Type: double
   :Unit: 
.. py:attribute:: gampz
   
   
   
   :Default: 
   :Dimension: (1:nxx)
   :Group: Convdiffeqn
   :Type: double
   :Unit: 
.. py:attribute:: gampzt
   
   
   
   :Default: 
   :Dimension: (1:nxx,1:ntim)
   :Group: Convdiffeqn
   :Type: double
   :Unit: 
.. py:attribute:: nnt
   
   solution at output times
   
   :Default: 
   :Dimension: (1:nxx,1:ntim)
   :Group: Convdiffeqn
   :Type: double
   :Unit: 
.. py:attribute:: timo
   
   output times
   
   :Default: 
   :Dimension: (1:ntim)
   :Group: Convdiffeqn
   :Type: double
   :Unit: 
