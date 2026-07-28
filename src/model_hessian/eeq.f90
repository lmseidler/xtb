! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
! Copyright (C) 2026 Leopold M. Seidler
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

!> Shared utilities used by model Hessian implementations.
!>
!> The EEQ contribution is an optional additive term, not a standalone
!> model Hessian.

module xtb_modelhessian_eeq
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants
   use xtb_mctc_convert, only : aatoau
   use xtb_type_param, only : chrg_parameter
   implicit none

   public :: c6
   public :: vander
   public :: min_fk
   public :: ixyz
   public :: jnd
   public :: ind
   public :: rcutoff
   public :: itabrow
   public :: outofp2
   public :: trsn2
   public :: strtch2
   public :: bend2
   public :: getvdw_hess
   public :: fk_lindh
   public :: fk_swart
   public :: fk_vdw
   public :: add_eeq_hessian
   private

!  van-der-Waals radii used in the D2 model
   real(wp),parameter :: vander(86) = aatoau * (/ &
     0.91_wp,0.92_wp, & ! H, He
     0.75_wp,1.28_wp,1.35_wp,1.32_wp,1.27_wp,1.22_wp,1.17_wp,1.13_wp, & ! Li-Ne
     1.04_wp,1.24_wp,1.49_wp,1.56_wp,1.55_wp,1.53_wp,1.49_wp,1.45_wp, & ! Na-Ar
     1.35_wp,1.34_wp, & ! K, Ca
     1.42_wp,1.42_wp,1.42_wp,1.42_wp,1.42_wp, & ! Sc-Zn
     1.42_wp,1.42_wp,1.42_wp,1.42_wp,1.42_wp, &
     1.50_wp,1.57_wp,1.60_wp,1.61_wp,1.59_wp,1.57_wp, & ! Ga-Kr
     1.48_wp,1.46_wp, & ! Rb, Sr
     1.49_wp,1.49_wp,1.49_wp,1.49_wp,1.49_wp, & ! Y-Cd
     1.49_wp,1.49_wp,1.49_wp,1.49_wp,1.49_wp, &
     1.52_wp,1.64_wp,1.71_wp,1.72_wp,1.72_wp,1.71_wp, & ! In-Xe
     2.00_wp,2.00_wp, &
     2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp, & ! La-Yb
     2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp, &
     2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp, & ! Lu-Hg
     2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp, &
     2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp /) ! Tl-Rn
!  C6 coefficients used in the D2 model
   real(wp),parameter :: c6(86) = (/&
      0.14_wp, 0.08_wp, & ! H,He
      1.61_wp, 1.61_wp, 3.13_wp, 1.75_wp, 1.23_wp, 0.70_wp, 0.75_wp, 0.63_wp, &
      5.71_wp, 5.71_wp,10.79_wp, 9.23_wp, 7.84_wp, 5.57_wp, 5.07_wp, 4.61_wp, &
     10.80_wp,10.80_wp, & ! K,Ca
     10.80_wp,10.80_wp,10.80_wp,10.80_wp,10.80_wp, & ! Sc-Zn
     10.80_wp,10.80_wp,10.80_wp,10.80_wp,10.80_wp, &
     16.99_wp,17.10_wp,16.37_wp,12.64_wp,12.47_wp,12.01_wp, & ! Ga-Kr
     24.67_wp,24.67_wp, & ! Rb,Sr
     24.67_wp,24.67_wp,24.67_wp,24.67_wp,24.67_wp, & ! Y-Cd
     24.67_wp,24.67_wp,24.67_wp,24.67_wp,24.67_wp, &
     37.32_wp,38.71_wp,38.44_wp,31.74_wp,31.50_wp,29.99_wp, & ! In-Xe
     50.00_wp,50.00_wp, & ! Cs,Ba
     50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp, & ! La-Yb
     50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp, &
     50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp, & ! Lu-Hg
     50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp, &
     50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp /) ! Tl-Rn

   real(wp),parameter :: min_fk = 1.0e-3_wp

contains

pure elemental function ixyz(i,iatom)
   integer :: ixyz
   integer,intent(in) :: i,iatom
   ixyz = (iatom-1)*3 + i
end function ixyz

pure elemental function jnd(i,j)
   integer :: jnd
   integer,intent(in) :: i,j
   jnd = i*(i-1)/2 +j
end function jnd

pure elemental function ind(i,iatom,j,jatom)
   integer :: ind
   integer,intent(in) :: i,iatom,j,jatom
   ind=jnd(max(ixyz(i,iatom),ixyz(j,jatom)),min(ixyz(i,iatom),ixyz(j,jatom)))
end function ind

pure function rcutoff(xyz,katom,latom,rcut)
   logical  :: rcutoff
   real(wp),intent(in) :: xyz(3,*)
   real(wp),intent(in) :: rcut
   real(wp) :: rkl(3),rkl2
   integer, intent(in) :: katom,latom
   rcutoff=.false.
   rkl=xyz(:,kAtom)-xyz(:,lAtom)
   rkl2 = sum(rkl**2)
   if(rkl2.gt.rcut) rcutoff=.true.
end function rcutoff

pure elemental function itabrow(i)
   integer :: itabrow
   integer,intent(in) :: i

   ! NOTE: rows 4-7 all return 3 — every element past Ar (Z>18) collapses
   ! to "row 3". Intentional? Periods 4,5,6,7 (K-Rn, Z=19-86) should be
   ! rows 4,5,6,7. Carried over verbatim from legacy src/lindh.f90 iTabRow.
   ! Affects Lindh/Swart/GFF model Hessian force-constant scaling via
   ! fk_lindh/fk_swart which branch on itabrow return value.
   itabrow=0
   if (i.gt. 0 .and. i.le. 2) then
      itabrow=1
   else if (i.gt. 2 .and. i.le.10) then
      itabrow=2
   else if (i.gt.10 .and. i.le.18) then
      itabrow=3
   else if (i.gt.18 .and. i.le.36) then
      itabrow=3
   else if (i.gt.36 .and. i.le.54) then
      itabrow=3
   else if (i.gt.54 .and. i.le.86) then
      itabrow=3
   else if (i.gt.86) then
      itabrow=3
   end if

   return
end function itabrow

pure subroutine outofp2(xyz,teta,bt)
   real(wp),intent(out) :: teta
   real(wp),intent(out) :: bt(3,4)
   real(wp),intent(in)  :: xyz(3,4)
   real(wp) :: r1(3),r2(3),r3(3)
   real(wp) :: q41,q42,q43,e41(3),e42(3),e43(3)
   real(wp) :: cosfi1,fi1,dfi1,cosfi2,fi2,dfi2,cosfi3,fi3,dfi3
   real(wp) :: c14(3,3),br14(3,3)
   real(wp) :: r42(3),r43(3)
   integer  :: ix,iy,iz
!  4 -> 1 (bond)
   r1=xyz(:,1)-xyz(:,4)
   q41=norm2(r1)
   e41 = r1 / q41
!  4 -> 2 (bond in plane)
   r2=xyz(:,2)-xyz(:,4)
   q42=norm2(r2)
   e42 = r2 / q42
!  4 -> 3 (bond in plane)
   r3=xyz(:,3)-xyz(:,4)
   q43=norm2(r3)
   e43 = r3 / q43
!
!  get the angle between e43 and e42
!
   cosfi1 = dot_product(e43,e42)

   fi1=acos(cosfi1)
   dfi1 = 180.d0 * fi1 / pi
!
!  dirty exit! this happens when an earlier structure is ill defined.
!
   if (abs(fi1-pi).lt.1.0d-13) then
      teta=0.0_wp
      bt = 0.0_wp
      return
   end if
!
!  get the angle between e41 and e43
!
   cosfi2 = dot_product(e41,e43)

   fi2=acos(cosfi2)
   dfi2 = 180.d0 * fi2 / pi
!
!  get the angle between e41 and e42
!
   cosfi3 = dot_product(e41,e42)

   fi3=acos(cosfi3)
   dfi3 = 180.d0 * fi3 / pi
!
!  the first two centers are trivially
!
   c14(:,1) = xyz(:,1)
   c14(:,2) = xyz(:,4)
!
!  the 3rd is
!
   r42=xyz(:,2)-xyz(:,4)
   r43=xyz(:,3)-xyz(:,4)
   c14(1,3)=r42(2)*r43(3)-r42(3)*r43(2)
   c14(2,3)=r42(3)*r43(1)-r42(1)*r43(3)
   c14(3,3)=r42(1)*r43(2)-r42(2)*r43(1)
!
!  exit if 2-3-4 are collinear
!  (equivalent to the above check, but this is more concrete)
!
   if ((c14(1,3)**2+c14(2,3)**2+c14(3,3)**2).lt.1.0d-10) then
      teta=0.0d0
      bt = 0.0_wp
      return
   end if
   c14(1,3)=c14(1,3)+xyz(1,4)
   c14(2,3)=c14(2,3)+xyz(2,4)
   c14(3,3)=c14(3,3)+xyz(3,4)

   call bend2(c14,teta,br14)

   teta = teta - 0.5_wp*pi
!
!--compute the wdc matrix
!
   do ix = 1, 3
      iy = mod(ix+1, 4)+(ix+1)/4
      iz = mod(iy+1, 4)+(iy+1)/4

      bt(ix,1) = - br14(ix,1)
      bt(ix,2) =   r43(iz)*br14(iy,3) - r43(iy)*br14(iz,3)
      bt(ix,3) = - r42(iz)*br14(iy,3) + r42(iy)*br14(iz,3)

      bt(ix,4) = - (bt(ix,1)+bt(ix,2)+bt(ix,3))

   end do

   bt = -bt
end subroutine outofp2

pure subroutine trsn2(xyz,tau,bt)
   real(wp),intent(out) :: bt(3,4)
   real(wp),intent(out) :: tau
   real(wp),intent(in)  :: xyz(3,4)
   real(wp) :: rij(3),rij1,brij(3,2)
   real(wp) :: rjk(3),rjk1,brjk(3,2)
   real(wp) :: rkl(3),rkl1,brkl(3,2)
   real(wp) :: bf2(3,3),fi2,sinfi2,cosfi2
   real(wp) :: bf3(3,3),fi3,sinfi3,cosfi3
   real(wp) :: costau,sintau
   integer  :: ix,iy,iz
   call strtch2(xyz(1,1),rij1,brij)
   call strtch2(xyz(1,2),rjk1,brjk)
   call strtch2(xyz(1,3),rkl1,brkl)
   call bend2(xyz(1,1),fi2,bf2)
   sinfi2=sin(fi2)
   cosfi2=cos(fi2)
   call bend2(xyz(1,2),fi3,bf3)
   sinfi3=sin(fi3)
   cosfi3=cos(fi3)
   costau = ( ( brij(2,1)*brjk(3,2) - brij(3,1)*brjk(2,2) ) * &
               ( brjk(2,1)*brkl(3,2) - brjk(3,1)*brkl(2,2) ) + &
               ( brij(3,1)*brjk(1,2) - brij(1,1)*brjk(3,2) ) * &
               ( brjk(3,1)*brkl(1,2) - brjk(1,1)*brkl(3,2) ) + &
               ( brij(1,1)*brjk(2,2) - brij(2,1)*brjk(1,2) ) * &
               ( brjk(1,1)*brkl(2,2) - brjk(2,1)*brkl(1,2) ) ) &
             / (sinfi2*sinfi3)
   sintau = ( brij(1,2) * (brjk(2,1)*brkl(3,2)-brjk(3,1)*brkl(2,2)) &
             + brij(2,2) * (brjk(3,1)*brkl(1,2)-brjk(1,1)*brkl(3,2)) &
             + brij(3,2) * (brjk(1,1)*brkl(2,2)-brjk(2,1)*brkl(1,2)) ) &
             / (sinfi2*sinfi3)
   tau = atan2(sintau,costau)
   if (abs(tau).eq.pi) tau=pi
   do ix = 1, 3
      iy=ix+1
      if (iy.gt.3) iy=iy-3
      iz=iy+1
      if (iz.gt.3) iz=iz-3
      bt(ix,1) = (brij(iy,2)*brjk(iz,2)-brij(iz,2)*brjk(iy,2)) &
         &           / (rij1*sinfi2**2)
      bt(ix,4) = (brkl(iy,1)*brjk(iz,1)-brkl(iz,1)*brjk(iy,1)) &
         &           / (rkl1*sinfi3**2)
      bt(ix,2) = -( (rjk1-rij1*cosfi2) * bt(ix,1) &
         &             +         rkl1*cosfi3  * bt(ix,4))/rjk1
      bt(ix,3) = - ( bt(ix,1)+bt(ix,2)+bt(ix,4))
   end do
end subroutine trsn2

pure subroutine strtch2(xyz,avst,b)
   real(wp),intent(out) :: b(3,2)
   real(wp),intent(in)  :: xyz(3,2)
   real(wp) :: r(3)
   real(wp) :: rr
   real(wp),intent(out) :: avst
   r=xyz(:,2)-xyz(:,1)
   rr=norm2(r)
   avst=rr
   b(:,1)=-r/rr
   b(:,2)=-b(:,1)
end subroutine strtch2

pure subroutine bend2(xyz,fir,bf)
   real(wp),intent(out) :: bf(3,3)
   real(wp),intent(in)  :: xyz(3,3)
   real(wp) :: brij(3,2)
   real(wp) :: brjk(3,2)
   real(wp) :: co,crap
   real(wp),intent(out) :: fir
   real(wp) :: si
   real(wp) :: rij1,rjk1
   integer  :: i
   call strtch2(xyz(1,1),rij1,brij)
   call strtch2(xyz(1,2),rjk1,brjk)
   co=0.0_wp
   crap=0.0_wp
   do i = 1, 3
      co=co+brij(i,1)*brjk(i,2)
      crap=crap+(brjk(i,2)+brij(i,1))**2
   end do
   if (sqrt(crap).lt.1.0d-6) then
      fir=pi-asin(sqrt(crap))
      si=sqrt(crap)
   else
      fir=acos(co)
      si=sqrt(1.0_wp-co**2)
   end if
   if (abs(fir-pi).lt.1.0d-13) then
      fir=pi
      return
   end if
   do i = 1, 3
      bf(i,1)= (co*brij(i,1)-brjk(i,2))/(si*rij1)
      bf(i,3)= (co*brjk(i,2)-brij(i,1))/(si*rjk1)
      bf(i,2)=-(bf(i,1)+bf(i,3))
   end do
end subroutine bend2

pure elemental subroutine getvdwxy(rx,ry,rz, c66, s6,r0, vdw)
   !cc Ableitung nach rx und ry
   real(wp),intent(in)  :: rx,ry,rz,c66,s6,r0
   real(wp),intent(out) :: vdw
   real(wp) :: t1,t2,t3,t4,t5,t6,t7,t11,t12,t16,t17,t25,t26,t35
   real(wp) :: t40,t41,t43,t44,t56,avdw
   ! write(*,*) 's6:', s6
   avdw=20.0
   t1 = s6 * C66
   t2 = rx ** 2
   t3 = ry ** 2
   t4 = rz ** 2
   t5 = t2 + t3 + t4
   t6 = t5 ** 2
   t7 = t6 ** 2
   t11 = sqrt(t5)
   t12 = 0.1D1 / r0
   t16 = exp(-avdw * (t11 * t12 - 0.1D1))
   t17 = 0.1D1 + t16
   t25 = t17 ** 2
   t26 = 0.1D1 / t25
   t35 = 0.1D1 / t7
   t40 = avdw ** 2
   t41 = r0 ** 2
   t43 = t40 / t41
   t44 = t16 ** 2
   t56 = -0.48D2 * t1 / t7 / t5 / t17 * rx * ry + 0.13D2 * t1 / t11 / &
      t7 * t26 * rx * avdw * t12 * ry * t16 - 0.2D1 * t1 * t35 / t25 / &
      t17 * t43 * rx * t44 * ry + t1 * t35 * t26 * t43 * rx * ry * t16
   vdw=t56
end subroutine getvdwxy

pure elemental subroutine getvdwxx(rx, ry, rz, c66, s6, r0, vdw)
   !cc Ableitung nach rx und rx
   real(wp),intent(in)  :: rx,ry,rz,c66,s6,r0
   real(wp),intent(out) :: vdw
   real(wp) :: t1,t2,t3,t4,t5,t6,t7,t10,t11,t15,t16,t17,t24,t25,t29
   real(wp) :: t33,t41,t42,t44,t45,t62,avdw
   avdw=20.0
   ! write(*,*) 's6:', s6
   t1 = s6 * C66
   t2 = rx ** 2
   t3 = ry ** 2
   t4 = rz ** 2
   t5 = t2 + t3 + t4
   t6 = t5 ** 2
   t7 = t6 ** 2
   t10 = sqrt(t5)
   t11 = 0.1D1 / r0
   t15 = exp(-avdw * (t10 * t11 - 0.1D1))
   t16 = 0.1D1 + t15
   t17 = 0.1D1 / t16
   t24 = t16 ** 2
   t25 = 0.1D1 / t24
   t29 = t11 * t15
   t33 = 0.1D1 / t7
   t41 = avdw ** 2
   t42 = r0 ** 2
   t44 = t41 / t42
   t45 = t15 ** 2
   t62 = -0.48D2 * t1 / t7 / t5 * t17 * t2 + 0.13D2 * t1 / t10 / t7 * &
      t25 * t2 * avdw * t29 + 0.6D1 * t1 * t33 * t17 - 0.2D1 * t1 * t33 &
      / t24 / t16 * t44 * t2 * t45 - t1 / t10 / t6 / t5 * t25 * avdw * &
      t29 + t1 * t33 * t25 * t44 * t2 * t15
   vdw=t62
end subroutine getvdwxx

!> Unified D2 dispersion Hessian block for a pair.
pure subroutine getvdw_hess(vec, c66, s6, r0, vdw)
   real(wp), intent(in) :: vec(3)
   real(wp), intent(in) :: c66, s6, r0
   real(wp), intent(out) :: vdw(3,3)

   call getvdwxx(vec(1), vec(2), vec(3), c66, s6, r0, vdw(1,1))
   call getvdwxy(vec(1), vec(2), vec(3), c66, s6, r0, vdw(1,2))
   call getvdwxy(vec(1), vec(3), vec(2), c66, s6, r0, vdw(1,3))
   call getvdwxx(vec(2), vec(1), vec(3), c66, s6, r0, vdw(2,2))
   call getvdwxy(vec(2), vec(3), vec(1), c66, s6, r0, vdw(2,3))
   call getvdwxx(vec(3), vec(1), vec(2), c66, s6, r0, vdw(3,3))
   vdw(2,1) = vdw(1,2)
   vdw(3,1) = vdw(1,3)
   vdw(3,2) = vdw(2,3)
end subroutine getvdw_hess

pure elemental function fk_lindh(alpha,r0,r2) result(gmm)
   real(wp),intent(in) :: alpha,r0,r2
   real(wp) :: gmm
   gmm = exp(alpha*(r0**2 - r2))
end function fk_lindh

pure elemental function fk_swart(alpha,r0,r2) result(gmm)
   real(wp),intent(in) :: alpha,r0,r2
   real(wp) :: gmm
   gmm = exp(-alpha*(sqrt(r2)/r0 - 1.0_wp))
end function fk_swart

pure elemental function fk_vdw(alpha,r0,r2) result(gmm)
   real(wp),intent(in) :: alpha,r0,r2
   real(wp) :: gmm
   gmm = exp(-alpha*(r0 - sqrt(r2))**2)
end function fk_vdw

!> Add EEQ contribution to an existing packed Hessian.
subroutine add_eeq_hessian(n,at,xyz,chrg,chrgeq,kq,hess)

!! ------------------------------------------------------------------------
!  Input
!! ------------------------------------------------------------------------
   integer, intent(in)    :: n                ! number of atoms
   integer, intent(in)    :: at(n)            ! ordinal numbers
   real(wp),intent(in)    :: xyz(3,n)         ! geometry
   real(wp),intent(in)    :: chrg             ! total charge
   real(wp),intent(in)    :: kq               ! scaling parameter
   type(chrg_parameter),intent(in) :: chrgeq  ! charge model
!! ------------------------------------------------------------------------
!  Output
!! ------------------------------------------------------------------------
   real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
   real(wp),allocatable   :: hessian(:,:,:,:) ! molecular hessian of IES

!  π itself
   real(wp),parameter :: pi = 3.1415926535897932384626433832795029_wp
!  √π
   real(wp),parameter :: sqrtpi = sqrt(pi)
!  √(2/π)
   real(wp),parameter :: sqrt2pi = sqrt(2.0_wp/pi)
!
!! ------------------------------------------------------------------------
!  charge model
!! ------------------------------------------------------------------------
   integer  :: m ! dimension of the Lagrangian
   real(wp),allocatable :: Amat(:,:)
   real(wp),allocatable :: Xvec(:)
   real(wp),allocatable :: Ainv(:,:)
   real(wp),allocatable :: dAmat(:,:,:)
   real(wp),allocatable :: dqdr(:,:,:)

!! ------------------------------------------------------------------------
!  local variables
!! ------------------------------------------------------------------------
   integer  :: i,j,k,l
   real(wp) :: r,rij(3),r2
   real(wp) :: gamij,gamij2
   real(wp) :: arg,arg2,tmp,dtmp
   real(wp) :: lambda
   real(wp) :: es,expterm,erfterm
   real(wp) :: htmp,rxr(3,3)
   real(wp) :: rcovij,rr

!! ------------------------------------------------------------------------
!  scratch variables
!! ------------------------------------------------------------------------
   real(wp),allocatable :: alpha(:)
   real(wp),allocatable :: xtmp(:)
   real(wp),allocatable :: atmp(:,:)

!! ------------------------------------------------------------------------
!  Lapack work variables
!! ------------------------------------------------------------------------
   integer, allocatable :: ipiv(:)
   real(wp),allocatable :: temp(:)
   real(wp),allocatable :: work(:)
   integer  :: lwork
   integer  :: info
   real(wp) :: test(1)

!! ------------------------------------------------------------------------
!  initizialization
!! ------------------------------------------------------------------------
   m    = n+1
   allocate( ipiv(m), source = 0 )
   allocate( Amat(m,m), Xvec(m), alpha(n), dqdr(3,n,m), source = 0.0_wp )

!! ------------------------------------------------------------------------
!  set up the A matrix and X vector
!! ------------------------------------------------------------------------
!  αi -> alpha(i), ENi -> xi(i), κi -> kappa(i), Jii -> gam(i)
!  γij = 1/√(αi+αj)
!  Xi  = -ENi + κi·√CNi
!  Aii = Jii + 2/√π·γii
!  Aij = erf(γij·Rij)/Rij = 2/√π·F0(γ²ij·R²ij)
!! ------------------------------------------------------------------------
!  prepare some arrays
!$omp parallel default(none) &
!$omp shared(n,at,chrgeq) &
!$omp private(i) &
!$omp shared(Xvec,alpha)
!$omp do schedule(dynamic)
   do i = 1, n
      Xvec(i) = -chrgeq%en(i)
      alpha(i) = chrgeq%alpha(i)**2
   enddo
!$omp enddo
!$omp endparallel

!$omp parallel default(none) &
!$omp shared(n,at,xyz,chrgeq,alpha) &
!$omp private(i,j,r,gamij) &
!$omp shared(Amat)
!$omp do schedule(dynamic)
   ! prepare A matrix
   do i = 1, n
      ! EN of atom i
      do j = 1, i-1
         r = sqrt(sum((xyz(:,j) - xyz(:,i))**2))
         gamij = 1.0_wp/sqrt(alpha(i)+alpha(j))
         Amat(j,i) = erf(gamij*r)/r
         Amat(i,j) = Amat(j,i)
      enddo
      Amat(i,i) = chrgeq%gam(i) + sqrt2pi/sqrt(alpha(i))
   enddo
!$omp enddo
!$omp endparallel

!! ------------------------------------------------------------------------
!  solve the linear equations to obtain partial charges
!! ------------------------------------------------------------------------
   Amat(m,1:m) = 1.0_wp
   Amat(1:m,m) = 1.0_wp
   Amat(m,m  ) = 0.0_wp
   Xvec(m)     = chrg
   ! generate temporary copy
   allocate( Atmp(m,m), source = Amat )
   allocate( Xtmp(m),   source = Xvec )

   ! assume work space query, set best value to test after first dsysv call
   call dsysv('u', m, 1, Atmp, m, ipiv, Xtmp, m, test, -1, info)
   lwork = int(test(1))
   allocate( work(lwork), source = 0.0_wp )

   call dsysv('u',m,1,Atmp,m,ipiv,Xtmp,m,work,lwork,info)
   if(info > 0) call raise('E','(goedecker_solve) DSYSV failed')

   if(abs(sum(Xtmp(:n))-chrg) > 1.e-6_wp) &
      call raise('E','(goedecker_solve) charge constrain error')
   !print'(3f20.14)',Xtmp

!! ------------------------------------------------------------------------
!  calculate isotropic electrostatic (IES) energy
!! ------------------------------------------------------------------------
!  E = ∑i (ENi - κi·√CNi)·qi + ∑i (Jii + 2/√π·γii)·q²i
!      + ½ ∑i ∑j,j≠i qi·qj·2/√π·F0(γ²ij·R²ij)
!    = q·(½A·q - X)
!! ------------------------------------------------------------------------
!   work(:m) = Xvec
!   call dsymv('u',m,0.5_wp,Amat,m,Xtmp,1,-1.0_wp,work,1)
!   es = dot_product(Xtmp,work(:m))
!   energy = es + energy

!! ------------------------------------------------------------------------
!  calculate molecular gradient of the IES energy
!! ------------------------------------------------------------------------
!  dE/dRj -> g(:,j), ∂Xi/∂Rj -> -dcn(:,i,j), ½∂Aij/∂Rj -> dAmat(:,j,i)
!  dE/dR = (½∂A/∂R·q - ∂X/∂R)·q
!  ∂Aij/∂Rj = ∂Aij/∂Ri
!! ------------------------------------------------------------------------
   allocate( dAmat(3,n,m), source = 0.0_wp )
!$omp parallel default(none) &
!$omp shared(n,xyz,alpha,Amat,Xtmp) &
!$omp private(i,j,rij,r2,gamij,arg,dtmp) &
!$omp reduction(+:dAmat)
!$omp do schedule(dynamic)
   do i = 1, n
      do j = 1, i-1
         rij = xyz(:,i) - xyz(:,j)
         r2 = sum(rij**2)
         gamij = 1.0_wp/sqrt(alpha(i) + alpha(j))
         arg = gamij**2*r2
         dtmp = 2.0_wp*gamij*exp(-arg)/(sqrtpi*r2)-Amat(j,i)/r2
         dAmat(:,i,i) = +dtmp*rij*Xtmp(j) + dAmat(:,i,i)
         dAmat(:,j,j) = -dtmp*rij*Xtmp(i) + dAmat(:,j,j)
         dAmat(:,i,j) = +dtmp*rij*Xtmp(i)
         dAmat(:,j,i) = -dtmp*rij*Xtmp(j)
      enddo
   enddo
!$omp enddo
!$omp endparallel

!! ------------------------------------------------------------------------
!  invert the A matrix using a Bunch-Kaufman factorization
!  A⁻¹ = (L·D·L^T)⁻¹ = L^T·D⁻¹·L
!! ------------------------------------------------------------------------
   allocate( Ainv(m,m), source = Amat )

   ! assume work space query, set best value to test after first dsytrf call
   call dsytrf('L',m,Ainv,m,ipiv,test,-1,info)
   if (int(test(1)) > lwork) then
      deallocate(work)
      lwork=int(test(1))
      allocate( work(lwork), source = 0.0_wp )
   endif

   ! Bunch-Kaufman factorization A = L*D*L**T
   call dsytrf('L',m,Ainv,m,ipiv,work,lwork,info)
   if(info > 0)then
      call raise('E', '(goedecker_inversion) DSYTRF failed')
   endif

   ! A⁻¹ from factorized L matrix, save lower part of A⁻¹ in Ainv matrix
   ! Ainv matrix is overwritten with lower triangular part of A⁻¹
   call dsytri('L',m,Ainv,m,ipiv,work,info)
   if (info > 0) then
      call raise('E', '(goedecker_inversion) DSYTRI failed')
   endif

   ! symmetrizes A⁻¹ matrix from lower triangular part of inverse matrix
   do i = 1, m
      do j = i+1, m
         Ainv(i,j)=Ainv(j,i)
      enddo
   enddo

!! ------------------------------------------------------------------------
!  calculate gradient of the partial charge w.r.t. the nuclear coordinates
!! ------------------------------------------------------------------------
   !call dsymm('r','l',3*n,m,-1.0_wp,Ainv,m,dAmat,3*n,1.0_wp,dqdr,3*n)
   call dgemm('n','n',3*n,m,m,-1.0_wp,dAmat,3*n,Ainv,m,1.0_wp,dqdr,3*n)
   !print'(/,"analytical gradient")'
   !print'(3f20.14)',dqdr(:,:,:n)

!! ------------------------------------------------------------------------
!  molecular Hessian calculation
!! ------------------------------------------------------------------------
   do i = 1, n
      do j = 1, i-1
         rij = xyz(:,j) - xyz(:,i)
         r2 = sum(rij**2)
         r = sqrt(r2)
         gamij = 1.0_wp/sqrt(alpha(i)+alpha(j))
         gamij2 = gamij**2
         arg2 = gamij2 * r2
         arg = sqrt(arg2)
         erfterm = Xtmp(i)*Xtmp(j)*erf(arg)/r
         expterm = Xtmp(i)*Xtmp(j)*2*gamij*exp(-arg2)/sqrtpi
         ! ∂²(qAq)/(∂Ri∂Rj):
         ! ∂²(qAq)/(∂Xi∂Xi) = (1-3X²ij/R²ij-2γ²ijX²ij) 2γij/√π exp[-γ²ij·R²ij]/R²ij
         !                  - (R²ij-3X²ij) erf[γij·Rij]/R⁵ij
         ! ∂²(qAq)/(∂Xi∂Xj) = (R²ij-3X²ij) erf[γij·Rij]/R⁵ij
         !                  - (1-3X²ij/R²ij-2γ²ijX²ij) 2γij/√π exp[-γ²ij·R²ij]/R²ij
         ! ∂²(qAq)/(∂Xi∂Yi) = 3X²ij erf[γij·Rij]/R⁵ij
         !                  - (3X²ij/R²ij+2γ²ijX²ij) 2γij/√π exp[-γ²ij·R²ij]/R²ij
         ! ∂²(qAq)/(∂Xi∂Yj) = (3X²ij/R²ij+2γ²ijX²ij) 2γij/√π exp[-γ²ij·R²ij]/R²ij
         !                  - 3X²ij erf[γij·Rij]/R⁵ij
         rxr(1,1) = erfterm * ( 3*rij(1)**2/r2**2 - 1.0_wp/r2 ) &
                  - expterm * ( 3*rij(1)**2/r2**2 + 2*gamij2*rij(1)**2/r2 - 1/r2 )
         rxr(2,2) = erfterm * ( 3*rij(2)**2/r2**2 - 1.0_wp/r2 ) &
                  - expterm * ( 3*rij(2)**2/r2**2 + 2*gamij2*rij(2)**2/r2 - 1/r2 )
         rxr(3,3) = erfterm * ( 3*rij(3)**2/r2**2 - 1.0_wp/r2 ) &
                  - expterm * ( 3*rij(3)**2/r2**2 + 2*gamij2*rij(3)**2/r2 - 1/r2 )
         rxr(2,1) = erfterm * 3*rij(2)*rij(1)/r2**2 &
                  - expterm * ( 3*rij(2)*rij(1)/r2**2 + 2*gamij2*rij(2)*rij(1)/r2 )
         rxr(3,1) = erfterm * 3*rij(3)*rij(1)/r2**2 &
                  - expterm * ( 3*rij(3)*rij(1)/r2**2 + 2*gamij2*rij(3)*rij(1)/r2 )
         rxr(3,2) = erfterm * 3*rij(3)*rij(2)/r2**2 &
                  - expterm * ( 3*rij(3)*rij(2)/r2**2 + 2*gamij2*rij(3)*rij(2)/r2 )

         do k = 1, m
            rxr(1,1) = rxr(1,1) + 0.5_wp*dqdr(1,i,k)*dAmat(1,j,k) &
                                + 0.5_wp*dqdr(1,j,k)*dAmat(1,i,k)
            rxr(2,1) = rxr(2,1) + 0.5_wp*dqdr(2,i,k)*dAmat(1,j,k) &
                                + 0.5_wp*dqdr(2,j,k)*dAmat(1,i,k)
            rxr(3,1) = rxr(3,1) + 0.5_wp*dqdr(3,i,k)*dAmat(1,j,k) &
                                + 0.5_wp*dqdr(3,j,k)*dAmat(1,i,k)
            rxr(2,2) = rxr(2,2) + 0.5_wp*dqdr(2,i,k)*dAmat(2,j,k) &
                                + 0.5_wp*dqdr(2,j,k)*dAmat(2,i,k)
            rxr(3,2) = rxr(3,2) + 0.5_wp*dqdr(3,i,k)*dAmat(2,j,k) &
                                + 0.5_wp*dqdr(3,j,k)*dAmat(2,i,k)
            rxr(3,3) = rxr(3,3) + 0.5_wp*dqdr(3,i,k)*dAmat(3,j,k) &
                                + 0.5_wp*dqdr(3,j,k)*dAmat(3,i,k)
         enddo
         ! symmetrize
         rxr(1,2) = rxr(2,1)
         rxr(1,3) = rxr(3,1)
         rxr(2,3) = rxr(3,2)

         ! save diagonal elements for atom i
         hess(ind(1,i,1,i))=hess(ind(1,i,1,i)) + kq*rxr(1,1)
         hess(ind(2,i,1,i))=hess(ind(2,i,1,i)) + kq*rxr(2,1)
         hess(ind(2,i,2,i))=hess(ind(2,i,2,i)) + kq*rxr(2,2)
         hess(ind(3,i,1,i))=hess(ind(3,i,1,i)) + kq*rxr(3,1)
         hess(ind(3,i,2,i))=hess(ind(3,i,2,i)) + kq*rxr(3,2)
         hess(ind(3,i,3,i))=hess(ind(3,i,3,i)) + kq*rxr(3,3)
         ! save elements between atom i and atom j
         hess(ind(1,i,1,j))=hess(ind(1,i,1,j)) - kq*rxr(1,1)
         hess(ind(1,i,2,j))=hess(ind(1,i,2,j)) - kq*rxr(2,1)
         hess(ind(1,i,3,j))=hess(ind(1,i,3,j)) - kq*rxr(3,1)
         hess(ind(2,i,1,j))=hess(ind(2,i,1,j)) - kq*rxr(2,1)
         hess(ind(2,i,2,j))=hess(ind(2,i,2,j)) - kq*rxr(2,2)
         hess(ind(2,i,3,j))=hess(ind(2,i,3,j)) - kq*rxr(3,2)
         hess(ind(3,i,1,j))=hess(ind(3,i,1,j)) - kq*rxr(3,1)
         hess(ind(3,i,2,j))=hess(ind(3,i,2,j)) - kq*rxr(3,2)
         hess(ind(3,i,3,j))=hess(ind(3,i,3,j)) - kq*rxr(3,3)
         ! save diagonal elements for atom j
         hess(ind(1,j,1,j))=hess(ind(1,j,1,j)) + kq*rxr(1,1)
         hess(ind(2,j,1,j))=hess(ind(2,j,1,j)) + kq*rxr(2,1)
         hess(ind(2,j,2,j))=hess(ind(2,j,2,j)) + kq*rxr(2,2)
         hess(ind(3,j,1,j))=hess(ind(3,j,1,j)) + kq*rxr(3,1)
         hess(ind(3,j,2,j))=hess(ind(3,j,2,j)) + kq*rxr(3,2)
         hess(ind(3,j,3,j))=hess(ind(3,j,3,j)) + kq*rxr(3,3)
      enddo
   enddo

   ! ∂²(qA)/(∂Ri∂q)·∂q/∂Rj
  ! hessian = hessian + reshape(matmul(reshape(dqdr,(/3*n,m/)),&
  !    transpose(reshape(dAmat,(/3*n,m/)))),(/3,n,3,n/))
   !call dgemm('n','t',3*n,m,3*n,+1.0_wp,dqdr,3*n,dAmat,3*n,1.0_wp,hessian,3*n)
   !call dgemm('n','t',3*n,m,3*n,+1.0_wp,dAmat,3*n,dqdr,3*n,1.0_wp,hessian,3*n)

end subroutine add_eeq_hessian

end module xtb_modelhessian_eeq
