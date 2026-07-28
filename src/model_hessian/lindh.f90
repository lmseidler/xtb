! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
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

!! ========================================================================
!  This module implements various model Hessians, the actual definition
!  of the model Hessian is given in the description of each header, the
!  implementation following tries to do its very best to transform it back
!  from redundant internal coordinates (where the model is defined) to
!  Cartesian coordinates which are then used in any imaginable fashion
!  later on in the actual optimization.                        - SAW190131
!! ========================================================================

module xtb_modelhessian_lindh
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi
   use xtb_chargemodel, only : new_charge_model_2019
   use xtb_bmatrix, only : bmat_bond, bmat_torsion, bmat_accum_packed, &
      & bmat_accum_pairblock_packed
   use xtb_modelhessian_type, only : TModelHessian
   use xtb_modelhessian_eeq, only : c6, vander, itabrow, rcutoff, outofp2, &
      & getvdw_hess, fk_lindh, fk_vdw, add_eeq_hessian
   use xtb_type_param, only : chrg_parameter
   implicit none (type, external)

   private

   type :: TLindhParameters
      real(wp) :: rav(3, 3)
      real(wp) :: aav(3, 3)
      real(wp) :: dav(3, 3)
      real(wp) :: outofplane_dispersion
   end type TLindhParameters

   type(TLindhParameters), parameter :: lindh_d2_parameters = TLindhParameters( &
      rav=reshape([ &
         1.3500_wp, 2.1000_wp, 2.5300_wp, &
         2.1000_wp, 2.8700_wp, 3.4000_wp, &
         2.5300_wp, 3.4000_wp, 3.4000_wp], [3, 3]), &
      aav=reshape([ &
         1.0000_wp, 0.3949_wp, 0.3949_wp, &
         0.3949_wp, 0.2800_wp, 0.2800_wp, &
         0.3949_wp, 0.2800_wp, 0.2800_wp], [3, 3]), &
      dav=reshape([ &
         0.0000_wp, 0.0000_wp, 0.0000_wp, &
         0.0000_wp, 0.0000_wp, 0.0000_wp, &
         0.0000_wp, 0.0000_wp, 0.0000_wp], [3, 3]), &
      outofplane_dispersion=1.0_wp)

   type(TLindhParameters), parameter :: lindh_parameters = TLindhParameters( &
      rav=reshape([ &
         1.3500_wp, 2.1000_wp, 2.5300_wp, &
         2.1000_wp, 2.8700_wp, 3.8000_wp, &
         2.5300_wp, 3.8000_wp, 4.5000_wp], [3, 3]), &
      aav=reshape([ &
         1.0000_wp, 0.3949_wp, 0.3949_wp, &
         0.3949_wp, 0.2800_wp, 0.1200_wp, &
         0.3949_wp, 0.1200_wp, 0.0600_wp], [3, 3]), &
      dav=reshape([ &
         0.0000_wp, 3.6000_wp, 3.6000_wp, &
         3.6000_wp, 5.3000_wp, 5.3000_wp, &
         3.6000_wp, 5.3000_wp, 5.3000_wp], [3, 3]), &
      outofplane_dispersion=0.0_wp)

   type, abstract, extends(TModelHessian) :: TLindhModelHessianBase
   contains
      procedure :: stretch
      procedure :: bend
      procedure :: torsion
      procedure :: outofplane
      procedure :: add_charge
      procedure(lindh_parameter_provider), deferred, private :: get_parameters
   end type TLindhModelHessianBase

   type, public, extends(TLindhModelHessianBase) :: TLindhModelHessian
   contains
      procedure, private :: get_parameters => get_lindh_parameters
   end type TLindhModelHessian

   type, public, extends(TLindhModelHessianBase) :: TLindhD2ModelHessian
   contains
      procedure, private :: get_parameters => get_lindh_d2_parameters
   end type TLindhD2ModelHessian

    abstract interface
       pure function lindh_parameter_provider(self) result(parameters)
          import :: TLindhModelHessianBase, TLindhParameters
          implicit none (type, external)
         class(TLindhModelHessianBase), intent(in) :: self
         type(TLindhParameters) :: parameters
      end function lindh_parameter_provider
   end interface

contains

pure function get_lindh_parameters(self) result(parameters)
   class(TLindhModelHessian), intent(in) :: self
   type(TLindhParameters) :: parameters

   parameters = lindh_parameters
end function get_lindh_parameters

pure function get_lindh_d2_parameters(self) result(parameters)
   class(TLindhD2ModelHessian), intent(in) :: self
   type(TLindhParameters) :: parameters

   parameters = lindh_d2_parameters
end function get_lindh_d2_parameters

subroutine add_charge(self, xyz, n, hess, at, kq)
   class(TLindhModelHessianBase), intent(in) :: self
   integer, intent(in) :: n
   real(wp), intent(in) :: xyz(3, n)
   real(wp), intent(inout) :: hess((3*n)*(3*n+1)/2)
   integer, intent(in) :: at(n)
   real(wp), intent(in) :: kq

   type(chrg_parameter) :: chrgeq

   call new_charge_model_2019(chrgeq, n, at)
   call add_eeq_hessian(n, at, xyz, 0.0_wp, chrgeq, kq, hess)
end subroutine add_charge

!! ========================================================================
!  Lindh's Model Hessian augmented with D2
!! ------------------------------------------------------------------------
!  Implemented after:
!  Lindh, R., Bernhardsson, A., Karlström, G., & Malmqvist, P.-Å. (1995).
!  On the use of a Hessian model function in molecular geometry optimizations.
!  Chem. Phys. Lett., 241(4), 423–428. doi:10.1016/0009-2614(95)00646-l
!
!  gij = exp[αij(R²ref - R²ij)]
!  kij   = rkr·gij
!  kijk  = rkf·gij·gjk
!  kijkl = rkt·gij·gjk·gkl
!
!  Originally Lindh proposed (we tweaked those a little bit):
!  rkr = 0.45, rkf = 0.15, rkt = 0.005
!
!  the reference distances are divided by rows in the PSE:
!  rAv:        1        2        3       aAv:        1        2        3
!    1    1.3500   2.1000   2.5300         1    1.0000   0.3949   0.3949
!    2    2.1000   2.8700   3.4000         2    0.3949   0.2800   0.2800
!    3    2.5300   3.4000   3.4000         3    0.3949   0.2800   0.2800
!
!  This Hessian is additionally augmented with D2, please note that D2
!  is not implemented in atomic units and requires some magical conversion
!  factor somewhere hidden in the implementation below.
!! ------------------------------------------------------------------------
!! ========================================================================
!  Lindh's Model Hessian updated around 2007
!! ------------------------------------------------------------------------
!  R. Lindh, personal communication.
!
!  gij = exp[αij(R²ref - R²ij)]
!  dij = exp[-4·(Rvdw - Rij)²]
!  kij   = rkr·gij + rkd·dij
!  kijk  = rkf·(gij+½·rkd/rkr·dij)·(gjk+½·rkd/rkr·djk)
!  kijkl = rkt·(gij+½·rkd/rkr·dij)·(gjk+½·rkd/rkr·djk)·(gkl+½·rkd/rkr·dkl)
!
!  parameters tweaked by R. Lindh in 2007:
!  rkr = 0.45, rkf = 0.10, rkt = 0.0025, rko = 0.16, rkd = 0.05
!
!  the reference distances are divided by rows in the PSE:
!  rAv:        1        2        3       aAv:        1        2        3
!    1    1.3500   2.1000   2.5300         1    1.0000   0.3949   0.3949
!    2    2.1000   2.8700   3.8000         2    0.3949   0.2800   0.1200
!    3    2.5300   3.8000   4.5000         3    0.3949   0.1200   0.0600
!
!  dAv:        1        2        3
!    1    0.0000   3.6000   3.6000
!    2    3.6000   5.3000   5.3000
!    3    3.6000   5.3000   5.3000
!
!! ------------------------------------------------------------------------
pure subroutine stretch(self, xyz, n, hess, at, kr, kd, s6, lcutoff, rcut)
   class(TLindhModelHessianBase), intent(in) :: self
   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
   real(wp),intent(in)    :: kr
   real(wp),intent(in)    :: kd
   real(wp),intent(in)    :: s6
   logical, intent(inout) :: lcutoff(n,n)
   real(wp),intent(in)    :: rcut

   integer  :: i,ir,j,jr
   type(TLindhParameters) :: parameters
   real(wp) :: aav(3,3), rav(3,3), dav(3,3)
   real(wp) :: vec(3), rij2, r0, d0
   real(wp) :: alpha,gmm
   real(wp) :: c6i,c6j,c6ij,rv
   real(wp) :: vdw(3,3)
   real(wp) :: bmat6(6)

   parameters = self%get_parameters()
   aav = parameters%aav
   rav = parameters%rav
   dav = parameters%dav

!! ------------------------------------------------------------------------
!  Hessian for stretch
!! ------------------------------------------------------------------------
   stretch_iAt: do i = 1, n
      ir = itabrow(at(i))

      stretch_jAt: do j = 1, i-1
         jr=itabrow(at(j))

         ! save for later
         lcutoff(i,j) = rcutoff(xyz,i,j,rcut)
         lcutoff(j,i) = lcutoff(i,j)

         vec = xyz(:,i) - xyz(:,j)
         rij2 = dot_product(vec, vec)
         r0 = rav(ir,jr)
         d0 = dav(ir,jr)
         alpha=aav(ir,jr)

         !cccccc vdwx ccccccccccccccccccccccccccccccccc
         c6i=c6(at(i))
         c6j=c6(at(j))
         c6ij=sqrt(c6i*c6j)
         rv=(vander(at(i))+vander(at(j)))

         call getvdw_hess(vec, c6ij, s6, rv, vdw)
         !cccccc ende vdwx ccccccccccccccccccccccccccccccc

         gmm = kr*fk_lindh(alpha,r0,rij2) &
            + kr*kd * fk_vdw(4.0_wp,d0,rij2)

         !gmm = max(gmm,min_fk)

         ! pure stretch: gmm * B^T B
         bmat6 = bmat_bond(vec)
         call bmat_accum_packed(n, hess, [i, j], bmat6, gmm)

         ! D2 dispersion Cartesian second derivative
         call bmat_accum_pairblock_packed(n, hess, i, j, vdw)

      end do stretch_jAt
   end do stretch_iAt

end subroutine stretch

pure subroutine bend(self, xyz, n, hess, at, force_constant, kd, lcutoff)
   class(TLindhModelHessianBase), intent(in) :: self
   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
   real(wp),intent(in)    :: force_constant
   real(wp),intent(in)    :: kd
   logical, intent(in)    :: lcutoff(n,n)

   integer  :: i,ir,j,jr,m,mr,ii
   type(TLindhParameters) :: parameters
   real(wp) :: aav(3,3), rav(3,3), dav(3,3)
   real(wp),parameter :: rzero = 1.0e-10_wp
   real(wp) :: xij,yij,zij,rij2,rrij,r1
   real(wp) :: xmi,ymi,zmi,rmi2,rmi,r0mi,ami,d0mj,gmi
   real(wp) :: xmj,ymj,zmj,rmj2,rmj,r0mj,amj,d0mi,gmj
   real(wp) :: test,gij,rl2,rl,rmidotrmj
   real(wp) :: sinphi,cosphi,costhetax,costhetay,costhetaz
   real(wp) :: alpha
   real(wp) :: si(3),sj(3),sm(3),x(2),y(2),z(2)
   real(wp) :: bmat9(9)

   parameters = self%get_parameters()
   aav = parameters%aav
   rav = parameters%rav
   dav = parameters%dav

!! ------------------------------------------------------------------------
!  Hessian for bending
!! ------------------------------------------------------------------------
   bend_mAt: do m = 1, n
      mr=itabrow(at(m))
      bend_iAt: do i = 1, n
          if (i == m) cycle bend_iAt
         ir=itabrow(at(i))
         if(lcutoff(i,m)) cycle bend_iAt

         xmi=(xyz(1,i)-xyz(1,m))
         ymi=(xyz(2,i)-xyz(2,m))
         zmi=(xyz(3,i)-xyz(3,m))
         rmi2 = xmi**2 + ymi**2 + zmi**2
         rmi=sqrt(rmi2)
         r0mi=rav(mr,ir)
         d0mi=dav(mr,ir)
         ami=aav(mr,ir)

         bend_jAt: do j = 1, i-1
             if (j == m) cycle bend_jAt
            jr=itabrow(at(j))
            if(lcutoff(j,i)) cycle bend_jAt
            if(lcutoff(j,m)) cycle bend_jAt

            xmj=(xyz(1,j)-xyz(1,m))
            ymj=(xyz(2,j)-xyz(2,m))
            zmj=(xyz(3,j)-xyz(3,m))
            rmj2 = xmj**2 + ymj**2 + zmj**2
            rmj=sqrt(rmj2)
            r0mj=rav(mr,jr)
            d0mj=dav(mr,jr)
            amj=aav(mr,jr)

            ! test if zero angle
            test=xmi*xmj+ymi*ymj+zmi*zmj
            test=test/(rmi*rmj)
            if (abs(test-1.0_wp) < 1.0e-12_wp) cycle bend_jAt

            xij=(xyz(1,j)-xyz(1,i))
            yij=(xyz(2,j)-xyz(2,i))
            zij=(xyz(3,j)-xyz(3,i))
            rij2 = xij**2 + yij**2 + zij**2
            rrij=sqrt(rij2)

            gmi = fk_lindh(ami,r0mi,rmi2) &
                + 0.5_wp*kd * fk_vdw(4.0_wp,d0mi,rmi2)
            gmj = fk_lindh(amj,r0mj,rmj2) &
                + 0.5_wp*kd * fk_vdw(4.0_wp,d0mj,rmj2)

            gij = force_constant*gmi*gmj

            rl2=(ymi*zmj-zmi*ymj)**2+(zmi*xmj-xmi*zmj)**2+(xmi*ymj-ymi*xmj)**2

            if(rl2 < 1.e-14_wp) then
               rl=0.0_wp
            else
               rl=sqrt(rl2)
            end if

            !gij = max(gij,min_fk)

            if ((rmj > rzero) .and. (rmi > rzero) .and. (rrij > rzero)) then
               sinphi=rl/(rmj*rmi)
               rmidotrmj=xmi*xmj+ymi*ymj+zmi*zmj
               cosphi=rmidotrmj/(rmj*rmi)
               ! none linear case
               if (sinphi > rzero) then
                  si(1)=(xmi/rmi*cosphi-xmj/rmj)/(rmi*sinphi)
                  si(2)=(ymi/rmi*cosphi-ymj/rmj)/(rmi*sinphi)
                  si(3)=(zmi/rmi*cosphi-zmj/rmj)/(rmi*sinphi)
                  sj(1)=(cosphi*xmj/rmj-xmi/rmi)/(rmj*sinphi)
                  sj(2)=(cosphi*ymj/rmj-ymi/rmi)/(rmj*sinphi)
                  sj(3)=(cosphi*zmj/rmj-zmi/rmi)/(rmj*sinphi)
                   sm(1)=-si(1)-sj(1)
                   sm(2)=-si(2)-sj(2)
                   sm(3)=-si(3)-sj(3)
                   bmat9 = [si, sm, sj]
                   call bmat_accum_packed(n, hess, [i, m, j], bmat9, gij)
                else
                  ! linear case
                  if ((abs(ymi) > rzero) .or. (abs(xmi) > rzero)) then
                     x(1)=-ymi
                     y(1)=xmi
                     z(1)=0.0_wp
                     x(2)=-xmi*zmi
                     y(2)=-ymi*zmi
                     z(2)=xmi*xmi+ymi*ymi
                  else
                     x(1)=1.0_wp
                     y(1)=0.0_wp
                     z(1)=0.0_wp
                     x(2)=0.0_wp
                     y(2)=1.0_wp
                     z(2)=0.0_wp
                  end if
                  do ii=1,2
                     r1=sqrt(x(ii)**2+y(ii)**2+z(ii)**2)
                     costhetax=x(ii)/r1
                     costhetay=y(ii)/r1
                     costhetaz=z(ii)/r1
                     si(1)=-costhetax/rmi
                     si(2)=-costhetay/rmi
                     si(3)=-costhetaz/rmi
                     sj(1)=-costhetax/rmj
                     sj(2)=-costhetay/rmj
                     sj(3)=-costhetaz/rmj
                     sm(1)=-(si(1)+sj(1))
                     sm(2)=-(si(2)+sj(2))
                     sm(3)=-(si(3)+sj(3))
                     bmat9 = [si, sm, sj]
                     call bmat_accum_packed(n, hess, [i, m, j], bmat9, gij)
                   end do

               end if
            end if

         end do bend_jAt
      end do bend_iAt
   end do bend_mAt

end subroutine bend

subroutine torsion(self, xyz, n, hess, at, force_constant, kd, lcutoff)
   class(TLindhModelHessianBase), intent(in) :: self
   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
   real(wp),intent(in)    :: force_constant
   real(wp),intent(in)    :: kd
   logical, intent(in)    :: lcutoff(n,n)

   integer  :: i,ir,j,jr,k,kr,l,lr,ij,kl
   type(TLindhParameters) :: parameters
   real(wp) :: aav(3,3), rav(3,3), dav(3,3)
!  allow only angles in the range of 35-145
   real(wp),parameter :: a35 = (35.0_wp/180.0_wp)* pi
   real(wp),parameter :: cosfi_max=cos(a35)
   real(wp) :: txyz(3,4),c(3,4)
   real(wp) :: rij(3),rij0,aij,rij2,d0ij,gij
   real(wp) :: rjk(3),rjk0,ajk,rjk2,d0jk,gjk
   real(wp) :: rkl(3),rkl0,akl,rkl2,d0kl,gkl
   real(wp) :: cosfi2,cosfi3,cosfi4
   real(wp) :: beta,tij
   real(wp) :: si(3),sj(3),sk(3),sl(3)
   real(wp) :: brow12(12)

   parameters = self%get_parameters()
   aav = parameters%aav
   rav = parameters%rav
   dav = parameters%dav

!! ------------------------------------------------------------------------
!  Hessian for torsion
!! ------------------------------------------------------------------------
   torsion_jAt: do j = 1,n
      jr=itabrow(at(j))
      txyz(:,2)=xyz(:,j)
      torsion_kAt: do k = 1, n
         if (k == j) cycle torsion_kAt
         kr=itabrow(at(k))
         if(lcutoff(k,j)) cycle torsion_kAt
         txyz(:,3) = xyz(:,k)
         torsion_iAt: do i = 1, n
            ij=n*(j-1)+i
            if (i == j) cycle torsion_iAt
            if (i == k) cycle torsion_iAt
            ir=itabrow(at(i))
            if(lcutoff(i,k)) cycle torsion_iAt
            if(lcutoff(i,j)) cycle torsion_iAt

            txyz(:,1)=xyz(:,i)
            torsion_lAt: do l = 1, n
               kl=n*(k-1)+l
               if (ij <= kl) cycle torsion_lAt
               if (l == i)   cycle torsion_lAt
               if (l == j)   cycle torsion_lAt
               if (l == k)   cycle torsion_lAt
               lr=itabrow(at(l))
!
               if(lcutoff(l,i)) cycle torsion_lAt
               if(lcutoff(l,k)) cycle torsion_lAt
               if(lcutoff(l,j)) cycle torsion_lAt

               txyz(:,4)=xyz(:,l)

               rij=xyz(:,i)-xyz(:,j)
               d0ij=dav(ir,jr)
               rij0=rav(ir,jr)
               aij =aav(ir,jr)

               rjk=xyz(:,j)-xyz(:,k)
               d0jk=dav(jr,kr)
               rjk0=rav(jr,kr)
               ajk =aav(jr,kr)

               rkl=xyz(:,k)-xyz(:,l)
               d0kl=dav(kr,lr)
               rkl0=rav(kr,lr)
               akl =aav(kr,lr)

               rij2=dot_product(rij,rij)
               rjk2=dot_product(rjk,rjk)
               rkl2=dot_product(rkl,rkl)

                cosfi2=dot_product(rij,rjk)/sqrt(rij2*rjk2)
                if (abs(cosfi2) > cosfi_max) cycle torsion_lAt
                cosfi3=dot_product(rkl,rjk)/sqrt(rkl2*rjk2)
                if (abs(cosfi3) > cosfi_max) cycle torsion_lAt

               gij = fk_lindh(aij,rij0,rij2) &
                  + 0.5_wp*kd * fk_vdw(4.0_wp,d0ij,rij2)
               gjk = fk_lindh(ajk,rjk0,rjk2) &
                  + 0.5_wp*kd * fk_vdw(4.0_wp,d0jk,rjk2)
               gkl = fk_lindh(akl,rkl0,rkl2) &
                  + 0.5_wp*kd * fk_vdw(4.0_wp,d0kl,rkl2)

               tij = force_constant * gij*gjk*gkl

               !tij = max(tij,10*min_fk)

               c = bmat_torsion(txyz)
               si = c(:,1)
               sj = c(:,2)
               sk = c(:,3)
               sl = c(:,4)
               brow12 = [si, sj, sk, sl]
               call bmat_accum_packed(n, hess, [i, j, k, l], brow12, tij)

            end do torsion_lAt
         end do torsion_iAt
      end do torsion_kAt
   end do torsion_jAt

end subroutine torsion

pure subroutine outofplane(self, xyz, n, hess, at, force_constant, kd, lcutoff)
   class(TLindhModelHessianBase), intent(in) :: self
   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
   real(wp),intent(in)    :: force_constant
   real(wp),intent(in)    :: kd
   logical, intent(in)    :: lcutoff(n,n)

   integer  :: i,ir,j,jr,k,kr,l,lr
   type(TLindhParameters) :: parameters
   real(wp) :: aav(3,3), rav(3,3), dav(3,3), outofplane_kd
   real(wp) :: txyz(3,4),c(3,4)
   real(wp) :: rij(3),rij0,aij,rij2,gij,d0ij
   real(wp) :: rik(3),rik0,aik,rik2,gik,d0ik
   real(wp) :: ril(3),ril0,ail,ril2,gil,d0il
   real(wp) :: cosfi2,cosfi3,cosfi4
   real(wp) :: beta,tij,tau
   real(wp) :: si(3),sj(3),sk(3),sl(3)
   real(wp) :: brow12(12)

   parameters = self%get_parameters()
   aav = parameters%aav
   rav = parameters%rav
   dav = parameters%dav
   outofplane_kd = parameters%outofplane_dispersion*kd

!! ------------------------------------------------------------------------
!  Hessian for out-of-plane
!! ------------------------------------------------------------------------
   outofplane_iAt: do i = 1, n
      ir = itabrow(at(i))
      txyz(:,4) = xyz(:,i)
      outofplane_jAt: do j = 1, n
         if (j == i) cycle outofplane_jAt
         if(lcutoff(j,i)) cycle outofplane_jAt
         jr = itabrow(at(j))
         txyz(:,1) = xyz(:,j)
         outofplane_kAt: do k = 1, n
            if (i == k) cycle outofplane_kAt
            if (j == k) cycle outofplane_kAt
            if(lcutoff(k,i)) cycle outofplane_kAt
            if(lcutoff(k,j)) cycle outofplane_kAt
            kr = itabrow(at(k))
            txyz(:,2) = xyz(:,k)
            outofplane_lAt: do l = 1, n
               lr = itabrow(at(l))
               txyz(:,3) = xyz(:,l)
               if (l == i) cycle outofplane_lAt
               if (l == j) cycle outofplane_lAt
               if (l == k) cycle outofplane_lAt
               if(lcutoff(l,i)) cycle outofplane_lAt
               if(lcutoff(l,k)) cycle outofplane_lAt
               if(lcutoff(l,j)) cycle outofplane_lAt

               rij=xyz(:,i)-xyz(:,j)
               d0ij=dav(ir,jr)
               rij0=rav(ir,jr)
               aij =aav(ir,jr)

               rik=xyz(:,i)-xyz(:,k)
               d0ik=dav(ir,kr)
               rik0=rav(ir,kr)
               aik =aav(ir,kr)

               ril=xyz(:,i)-xyz(:,l)
               d0il=dav(ir,lr)
               ril0=rav(ir,lr)
               ail =aav(ir,lr)

               rij2=dot_product(rij,rij)
               rik2=dot_product(rik,rik)
               ril2=dot_product(ril,ril)

                cosfi2=dot_product(rij,rik)/sqrt(rij2*rik2)
                if (abs(abs(cosfi2)-1.0_wp) < 1.0e-1_wp) cycle outofplane_lAt
                cosfi3=dot_product(rij,ril)/sqrt(rij2*ril2)
                if (abs(abs(cosfi3)-1.0_wp) < 1.0e-1_wp) cycle outofplane_lAt
                cosfi4=dot_product(rik,ril)/sqrt(rik2*ril2)
                if (abs(abs(cosfi4)-1.0_wp) < 1.0e-1_wp) cycle outofplane_lAt

               gij = fk_lindh(aij,rij0,rij2) &
                  + 0.5_wp*outofplane_kd * fk_vdw(4.0_wp,d0ij,rij2)
               gik = fk_lindh(aik,rik0,rik2) &
                  + 0.5_wp*outofplane_kd * fk_vdw(4.0_wp,d0ik,rik2)
               gil = fk_lindh(ail,ril0,ril2) &
                  + 0.5_wp*outofplane_kd * fk_vdw(4.0_wp,d0il,ril2)

               tij = force_constant * gij*gik*gil

               !tij = max(tij,10*min_fk)

                call outofp2(txyz,tau,c)
                if (abs(tau) > 45.0_wp*(pi/180.0_wp)) cycle outofplane_lAt

               si = c(:,4)
               sj = c(:,1)
               sk = c(:,2)
               sl = c(:,3)
               brow12 = [si, sj, sk, sl]
               call bmat_accum_packed(n, hess, [i, j, k, l], brow12, tij)

            end do outofplane_lAt
         end do outofplane_kAt
      end do outofplane_jAt
   end do outofplane_iAt

end subroutine outofplane

end module xtb_modelhessian_lindh
