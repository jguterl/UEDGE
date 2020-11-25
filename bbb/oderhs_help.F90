subroutine capsfscal(sfscal,neq,yl)
integer ix,iy,ifld,igsp
real,intent(inout) sfscal(neq)
integer,intent(in):: neq
real,intent(in) yl(*)
Use Jacaux, only: yldot0
Use Dim, only: nx,ny,ngsp,nisp
USe Indexes, only: idxn,idxg
use capfloor,only: icapsfscal,acapsfscal
use compla, only: ni,ng
use UEpar,only: nzbackg,ngbackg

if (icapsfscal.gt.0) then
call rhsnk (neq, yl, yldot0)

do ix=1,nx
    do iy=1,ny
        do ifld=1,nisp
    ii=idxn(ix,iy,ifld)
    sfscal(ii)=sfscal(ii)*(acapsfscal*nzbackg(ifld)/ni(ix,iy,ifld))**icapsfscal
    enddo
       do igsp=1,ngsp
          ii=idxg(ix,iy,igsp)
          sfscal(ii)=sfscal(ii)*(acapsfscal*ngbackg(igsp)/ng(ix,iy,igsp))**icapsfscal
       enddo
    enddo

    enddo
endif
end subroutine capsfscal
