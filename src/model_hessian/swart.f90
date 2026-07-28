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

module xtb_modelhessian_swart
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi
   use xtb_mctc_param, only : covalent_radius_2009
   use xtb_chargemodel, only : new_charge_model_2019
   use xtb_bmatrix, only : bmat_bond, bmat_angle, bmat_linbend, bmat_torsion, &
      & bmat_outofplane, oop_angle, bmat_accum_packed, bmat_accum_pairblock_packed
   use xtb_modelhessian_type, only : TModelHessian
   use xtb_modelhessian_eeq, only : c6, vander, rcutoff, getvdw_hess, &
      & fk_swart, fk_vdw, add_eeq_hessian
   use xtb_type_param, only : chrg_parameter
   implicit none (type, external)

   private

   real(wp), parameter :: swart_rcov(size(covalent_radius_2009)) = &
      covalent_radius_2009
   ! Current parameterization uses equal values, but rvdw remains a distinct input.
   real(wp), parameter :: swart_rvdw(size(covalent_radius_2009)) = &
      covalent_radius_2009

   type, public, extends(TModelHessian) :: TSwartModelHessian
   contains
      procedure :: stretch
      procedure :: bend
      procedure :: torsion
      procedure :: outofplane
      procedure :: add_charge
   end type TSwartModelHessian

contains

!! ========================================================================
!  Swart's Model Hessian augmented with D2
!! ------------------------------------------------------------------------
!  Implemented after:
!  M. Swart, F. M. Bickelhaupt, Int. J. Quantum Chem., 2006, 106, 2536–2544.
!  DOI:10.1002/qua.21049
!
!  gij = exp[-(Rij/Cij-1)]
!  kij   = rkr·gij
!  kijk  = rkf·gij·gjk
!  kijkl = rkt·gij·gjk·gkl
!
!  The proposed force constants by Swart are:
!  rkr = 0.35, rkf = 0.15, rkt = 0.005
!
!  This Hessian is additionally augmented with D2, please note that D2
!  is not implemented in atomic units and requires some magical conversion
!  factor somewhere hidden in the implementation below.
!! ------------------------------------------------------------------------
subroutine add_charge(self, xyz, n, hess, at, kq)
   class(TSwartModelHessian), intent(in) :: self
   integer, intent(in) :: n
   real(wp), intent(in) :: xyz(3, n)
   real(wp), intent(inout) :: hess((3*n)*(3*n+1)/2)
   integer, intent(in) :: at(n)
   real(wp), intent(in) :: kq

   type(chrg_parameter) :: chrgeq

   call new_charge_model_2019(chrgeq, n, at)
   call add_eeq_hessian(n, at, xyz, 0.0_wp, chrgeq, kq, hess)
end subroutine add_charge

pure subroutine stretch(self, xyz, n, hess, at, kr, kd, s6, lcutoff, rcut)
   class(TSwartModelHessian), intent(in) :: self
   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
   real(wp),intent(in)    :: kr
   real(wp),intent(in)    :: kd
   real(wp),intent(in)    :: s6
   logical, intent(inout) :: lcutoff(n,n)
   real(wp),intent(in)    :: rcut

   integer  :: i,j
   real(wp) :: rij(3), rij2, r0, d0
   real(wp) :: gmm
   real(wp) :: c6i,c6j,c6ij,rv
   real(wp) :: vdw(3,3)
   real(wp) :: bmat6(6)

!! ------------------------------------------------------------------------
!  Hessian for stretch
!! ------------------------------------------------------------------------
   stretch_iAt: do i = 1, n

      stretch_jAt: do j = 1, i-1

         ! save for later
         lcutoff(i,j) = rcutoff(xyz,i,j,rcut)
         lcutoff(j,i) = lcutoff(i,j)

         rij = xyz(:,i) - xyz(:,j)
         rij2 = dot_product(rij, rij)
         r0 = swart_rcov(at(i))+swart_rcov(at(j))
         d0 = swart_rvdw(at(i))+swart_rvdw(at(j))

         !cccccc vdwx ccccccccccccccccccccccccccccccccc
         c6i = c6(at(i))
         c6j = c6(at(j))
         c6ij = sqrt(c6i*c6j)
         rv = vander(at(i)) + vander(at(j))

         call getvdw_hess(rij, c6ij, s6, rv, vdw)
         !cccccc ende vdwx ccccccccccccccccccccccccccccccc

         gmm = kr*fk_swart(1.0_wp,  r0,rij2) &
            + kr*kd * fk_vdw(5.0_wp,d0,rij2)

         !gmm = max(gmm,min_fk)

         ! pure stretch: gmm * B^T B
         bmat6 = bmat_bond(rij)
         call bmat_accum_packed(n, hess, [i, j], bmat6, gmm)

         ! D2 dispersion Cartesian second derivative
         call bmat_accum_pairblock_packed(n, hess, i, j, vdw)

      end do stretch_jAt
   end do stretch_iAt

end subroutine stretch

pure subroutine bend(self, xyz, n, hess, at, force_constant, kd, lcutoff)
   class(TSwartModelHessian), intent(in) :: self
   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
   real(wp),intent(in)    :: force_constant
   real(wp),intent(in)    :: kd
   logical, intent(in)    :: lcutoff(n,n)

   integer  :: i,j,m,ii
   real(wp),parameter :: rzero = 1.0e-10_wp
   real(wp) :: vec_mi(3), vec_mj(3), rmi2, rmi, r0mi, d0mi, gmi
   real(wp) :: rmj2, rmj, r0mj, d0mj, gmj
   real(wp) :: test, gij, rl2, rl, crv(3)
   real(wp) :: sinphi
   real(wp) :: bmat9(9), bmat29(2,9)

!! ------------------------------------------------------------------------
!  Hessian for bending
!! ------------------------------------------------------------------------
   bend_mAt: do m = 1, n
      bend_iAt: do i = 1, n
         if (i==m) cycle bend_iAt
         if (lcutoff(i,m)) cycle bend_iAt

         vec_mi = xyz(:,i) - xyz(:,m)
         rmi2 = dot_product(vec_mi, vec_mi)
         rmi = sqrt(rmi2)
         r0mi = swart_rcov(at(m)) + swart_rcov(at(i))
         d0mi = swart_rvdw(at(m)) + swart_rvdw(at(i))

         bend_jAt: do j = 1, i-1
            if (j==m) cycle bend_jAt
            if (lcutoff(j,i)) cycle bend_jAt
            if (lcutoff(j,m)) cycle bend_jAt

            vec_mj = xyz(:,j) - xyz(:,m)
            rmj2 = dot_product(vec_mj, vec_mj)
            rmj = sqrt(rmj2)
            r0mj = swart_rcov(at(m)) + swart_rcov(at(j))
            d0mj = swart_rvdw(at(m)) + swart_rvdw(at(j))

            ! test if zero angle
            test = dot_product(vec_mi, vec_mj) / (rmi * rmj)
            if (abs(test - 1.0_wp) < 1.0e-12_wp) cycle bend_jAt

            gmi = fk_swart(1.0_wp, r0mi, rmi2) &
               + 0.5_wp*kd * fk_vdw(5.0_wp, d0mi, rmi2)
            gmj = fk_swart(1.0_wp, r0mj, rmj2) &
               + 0.5_wp*kd * fk_vdw(5.0_wp, d0mj, rmj2)

            gij = force_constant * gmi * gmj

            crv(1) = vec_mi(2)*vec_mj(3) - vec_mi(3)*vec_mj(2)
            crv(2) = vec_mi(3)*vec_mj(1) - vec_mi(1)*vec_mj(3)
            crv(3) = vec_mi(1)*vec_mj(2) - vec_mi(2)*vec_mj(1)
            rl2 = dot_product(crv, crv)

            if (rl2 < 1.0e-14_wp) then
               rl = 0.0_wp
            else
               rl = sqrt(rl2)
            end if

            !gij = max(gij,min_fk)

            if ((rmj > rzero) .and. (rmi > rzero)) then
               sinphi = rl / (rmj * rmi)
               ! none linear case
               if (sinphi > rzero) then
                  ! shared Wilson B row for non-linear angle (i-m-j)
                  bmat9 = bmat_angle(vec_mi, vec_mj)
                  call bmat_accum_packed(n, hess, [i, m, j], bmat9, gij)
               else
                  ! linear case: shared Wilson B rows for linear bend (two rows)
                  bmat29 = bmat_linbend(vec_mi, vec_mj)
                  do ii = 1, 2
                     call bmat_accum_packed(n, hess, [i, m, j], bmat29(ii, :), gij)
                  end do
               end if
            end if

         end do bend_jAt
      end do bend_iAt
   end do bend_mAt

end subroutine bend

pure subroutine torsion(self, xyz, n, hess, at, force_constant, kd, lcutoff)
   class(TSwartModelHessian), intent(in) :: self
   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
   real(wp),intent(in)    :: force_constant
   real(wp),intent(in)    :: kd
   logical, intent(in)    :: lcutoff(n,n)

   integer  :: i,j,k,l,ij,kl
!  allow only angles in the range of 35-145
   real(wp),parameter :: a35 = (35.0_wp/180.0_wp)* pi
   real(wp),parameter :: cosfi_max=cos(a35)
   real(wp) :: txyz(3,4),c(3,4)
   real(wp) :: rij(3),rij0,aij,rij2,d0ij,gij
   real(wp) :: rjk(3),rjk0,ajk,rjk2,d0jk,gjk
   real(wp) :: rkl(3),rkl0,akl,rkl2,d0kl,gkl
   real(wp) :: cosfi2,cosfi3,cosfi4
   real(wp) :: beta,tij
   real(wp) :: brow12(12)

!! ------------------------------------------------------------------------
!  Hessian for torsion
!! ------------------------------------------------------------------------
   torsion_jAt: do j = 1,n
      txyz(:,2)=xyz(:,j)
      torsion_kAt: do k = 1, n
         if (k==j) cycle torsion_kAt
         if(lcutoff(k,j)) cycle torsion_kAt
         txyz(:,3) = xyz(:,k)
         torsion_iAt: do i = 1, n
            ij=n*(j-1)+i
            if (i==j) cycle torsion_iAt
            if (i==k) cycle torsion_iAt
            if(lcutoff(i,k)) cycle torsion_iAt
            if(lcutoff(i,j)) cycle torsion_iAt

            txyz(:,1)=xyz(:,i)
            torsion_lAt: do l = 1, n
               kl=n*(l-1)+k
               if (ij<=kl) cycle torsion_lAt
               if (l==i)   cycle torsion_lAt
               if (l==j)   cycle torsion_lAt
               if (l==k)   cycle torsion_lAt
!
               if(lcutoff(l,i)) cycle torsion_lAt
               if(lcutoff(l,k)) cycle torsion_lAt
               if(lcutoff(l,j)) cycle torsion_lAt

               txyz(:,4)=xyz(:,l)

               rij=xyz(:,i)-xyz(:,j)
               d0ij=swart_rvdw(at(i))+swart_rvdw(at(j))
               rij0=swart_rcov(at(i))+swart_rcov(at(j))

               rjk=xyz(:,j)-xyz(:,k)
               d0jk=swart_rvdw(at(j))+swart_rvdw(at(k))
               rjk0=swart_rcov(at(j))+swart_rcov(at(k))

               rkl=xyz(:,k)-xyz(:,l)
               d0kl=swart_rvdw(at(k))+swart_rvdw(at(l))
               rkl0=swart_rcov(at(k))+swart_rcov(at(l))

               rij2=dot_product(rij,rij)
               rjk2=dot_product(rjk,rjk)
               rkl2=dot_product(rkl,rkl)

                cosfi2=dot_product(rij,rjk)/sqrt(rij2*rjk2)
                if (abs(cosfi2) > cosfi_max) cycle torsion_lAt
                cosfi3=dot_product(rkl,rjk)/sqrt(rkl2*rjk2)
                if (abs(cosfi3) > cosfi_max) cycle torsion_lAt

               gij = fk_swart(1.0_wp,rij0,rij2) &
                  + 0.5_wp*kd * fk_vdw(5.0_wp,d0ij,rij2)
               gjk = fk_swart(1.0_wp,rjk0,rjk2) &
                  + 0.5_wp*kd * fk_vdw(5.0_wp,d0jk,rjk2)
               gkl = fk_swart(1.0_wp,rkl0,rkl2) &
                  + 0.5_wp*kd * fk_vdw(5.0_wp,d0kl,rkl2)

               tij = force_constant * gij*gjk*gkl

                !tij = max(tij,10*min_fk)

                c = bmat_torsion(txyz)
                brow12 = [c(:,1), c(:,2), c(:,3), c(:,4)]
                call bmat_accum_packed(n, hess, [i, j, k, l], brow12, tij)

             end do torsion_lAt
          end do torsion_iAt
       end do torsion_kAt
    end do torsion_jAt

end subroutine torsion

pure subroutine outofplane(self, xyz, n, hess, at, force_constant, kd, lcutoff)
   class(TSwartModelHessian), intent(in) :: self
   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
   real(wp),intent(in)    :: force_constant
   real(wp),intent(in)    :: kd
   logical, intent(in)    :: lcutoff(n,n)

   integer  :: i,ir,j,jr,k,kr,l,lr
   real(wp) :: txyz(3,4),c(3,4)
   real(wp) :: rij(3),rij0,d0ij,rij2,gij
   real(wp) :: rik(3),rik0,d0ik,rik2,gik
   real(wp) :: ril(3),ril0,d0il,ril2,gil
   real(wp) :: cosfi2,cosfi3,cosfi4
   real(wp) :: beta,tij,tau
   real(wp) :: brow12(12)

!! ------------------------------------------------------------------------
!  Hessian for out-of-plane
!! ------------------------------------------------------------------------
   outofplane_iAt: do i = 1, n
      txyz(:,4) = xyz(:,i)
      outofplane_jAt: do j = 1, n
         if (j==i) cycle outofplane_jAt
         if(lcutoff(j,i)) cycle outofplane_jAt
         txyz(:,1) = xyz(:,j)
         outofplane_kAt: do k = 1, n
            if (i==k) cycle outofplane_kAt
            if (j==k) cycle outofplane_kat
            if(lcutoff(k,i)) cycle outofplane_kAt
            if(lcutoff(k,j)) cycle outofplane_kAt
            txyz(:,2) = xyz(:,k)
            outofplane_lAt: do l = 1, n
               txyz(:,3) = xyz(:,l)
               if (l==i)   cycle outofplane_lAt
               if (l==j)   cycle outofplane_lAt
               if (l==k)   cycle outofplane_lAt
               if(lcutoff(l,i)) cycle outofplane_lAt
               if(lcutoff(l,k)) cycle outofplane_lAt
               if(lcutoff(l,j)) cycle outofplane_lAt

               rij=xyz(:,i)-xyz(:,j)
               rij0=swart_rcov(at(i))+swart_rcov(at(j))
               d0ij=swart_rvdw(at(i))+swart_rvdw(at(j))

               rik=xyz(:,i)-xyz(:,k)
               rik0=swart_rcov(at(i))+swart_rcov(at(k))
               d0ik=swart_rvdw(at(i))+swart_rvdw(at(k))

               ril=xyz(:,i)-xyz(:,l)
               ril0=swart_rcov(at(i))+swart_rcov(at(l))
               d0il=swart_rvdw(at(i))+swart_rvdw(at(l))

               rij2=dot_product(rij,rij)
               rik2=dot_product(rik,rik)
               ril2=dot_product(ril,ril)

                cosfi2=dot_product(rij,rik)/sqrt(rij2*rik2)
                if (abs(abs(cosfi2)-1.0_wp) < 1.0e-1_wp) cycle outofplane_lAt
                cosfi3=dot_product(rij,ril)/sqrt(rij2*ril2)
                if (abs(abs(cosfi3)-1.0_wp) < 1.0e-1_wp) cycle outofplane_lAt
                cosfi4=dot_product(rik,ril)/sqrt(rik2*ril2)
                if (abs(abs(cosfi4)-1.0_wp) < 1.0e-1_wp) cycle outofplane_lAt

               gij = fk_swart(1.0_wp,rij0,rij2) &
                  + 0.5_wp*kd * fk_vdw(5.0_wp,d0ij,rij2)
               gik = fk_swart(1.0_wp,rik0,rik2) &
                  + 0.5_wp*kd * fk_vdw(5.0_wp,d0ik,rik2)
               gil = fk_swart(1.0_wp,ril0,ril2) &
                  + 0.5_wp*kd * fk_vdw(5.0_wp,d0il,ril2)

               tij = force_constant * gij*gik*gil

               !tij = max(tij,10*min_fk)

                 tau = oop_angle(txyz)
                 if (abs(tau) > 45.0_wp*(pi/180.0_wp)) cycle outofplane_lAt

                c = bmat_outofplane(txyz)
                brow12 = [c(:,4), c(:,1), c(:,2), c(:,3)]
                call bmat_accum_packed(n, hess, [i, j, k, l], brow12, tij)

            end do outofplane_lAt
          end do outofplane_kAt
       end do outofplane_jAt
    end do outofplane_iAt

end subroutine outofplane

end module xtb_modelhessian_swart
