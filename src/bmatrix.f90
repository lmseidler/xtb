! This file is part of xtb.
!
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
module xtb_bmatrix
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi
   use xtb_mctc_blas, only : mctc_dot
   implicit none

   private
   public :: bmat_bond, bmat_angle, bmat_linbend, &
      & bmat_torsion, bmat_outofplane

   character(len=*), parameter :: source = 'xtb_bmatrix'

contains

!> Wilson B matrix for bonds
pure function bmat_bond(vec) result(bmat)
   real(wp), intent(in) :: vec(3)

   real(wp) :: l, bmat(6)

   bmat = 0.0_wp
   l = norm2(vec)

   bmat(1:3) = vec(:) / l
   bmat(4:6) = -vec(:) / l
end function bmat_bond

!> Wilson B matrix for non-linear angles
function bmat_angle(vec1, vec2) result(bmat)
   real(wp), intent(in) :: vec1(3), vec2(3)
   real(wp) :: bmat(9)
   real(wp) :: l1, l2, nvec1(3), nvec2(3)
   real(wp) :: dl(2, 6), dnvec(2, 3, 6), dinprod(9)
   real(wp) :: dot_n1n2
   integer :: ii

   l1 = norm2(vec1)
   l2 = norm2(vec2)
   nvec1 = vec1 / l1
   nvec2 = vec2 / l2

   dl = 0.0_wp
   dl(1, 1:3) = nvec1
   dl(1, 4:6) = -nvec1
   dl(2, 1:3) = nvec2
   dl(2, 4:6) = -nvec2

   dnvec = 0.0_wp
   do ii = 1, 6
      dnvec(1, 1:3, ii) = -nvec1 * dl(1, ii) / l1
      dnvec(2, 1:3, ii) = -nvec2 * dl(2, ii) / l2
   end do
   do ii = 1, 3
      dnvec(1, ii, ii) = dnvec(1, ii, ii) + 1.0_wp / l1
      dnvec(2, ii, ii) = dnvec(2, ii, ii) + 1.0_wp / l2
      dnvec(1, ii, ii+3) = dnvec(1, ii, ii+3) - 1.0_wp / l1
      dnvec(2, ii, ii+3) = dnvec(2, ii, ii+3) - 1.0_wp / l2
   end do

   dinprod = 0.0_wp
   do ii = 1, 3
      dinprod(ii) = mctc_dot(dnvec(1, :, ii), nvec2)
      dinprod(ii+3) = mctc_dot(dnvec(1, :, ii+3), nvec2) + mctc_dot(dnvec(2, :, ii+3), nvec1)
      dinprod(ii+6) = mctc_dot(dnvec(2, :, ii), nvec1)
   end do

   dot_n1n2 = mctc_dot(nvec1, nvec2)
   bmat = -dinprod / sqrt(max(1.0e-15_wp, 1.0_wp-dot_n1n2**2))
end function bmat_angle

!> Wilson B matrix for linear angles
function bmat_linbend(vec1, vec2) result(bmat)
   real(wp), intent(in) :: vec1(3), vec2(3)
   real(wp) :: bmat(2, 9)
   real(wp) :: l1, l2, nvec1(3), nvec2(3)
   real(wp) :: vn(3), vn2(3), nvn
   real(wp), parameter :: xaxis(3) = [1.0_wp, 0.0_wp, 0.0_wp], yaxis(3) = [0.0_wp, 1.0_wp, 0.0_wp]

   l1 = norm2(vec1)
   l2 = norm2(vec2)
   nvec1 = vec1 / l1
   nvec2 = vec2 / l2

   vn(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
   vn(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
   vn(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)
   nvn = norm2(vn)

   if (nvn < 1.0e-15_wp) then
      vn = xaxis - mctc_dot(xaxis, vec1) / l1**2 * vec1
      nvn = norm2(vn)
      if (nvn < 1.0e-15_wp) then
         vn = yaxis - mctc_dot(yaxis, vec1) / l1**2 * vec1
         nvn = norm2(vn)
      end if
   end if
   vn = vn / nvn

   vn2(1) = (vec1(2)-vec2(2)) * vn(3) - (vec1(3)-vec2(3)) * vn(2)
   vn2(2) = (vec1(3)-vec2(3)) * vn(1) - (vec1(1)-vec2(1)) * vn(3)
   vn2(3) = (vec1(1)-vec2(1)) * vn(2) - (vec1(2)-vec2(2)) * vn(1)
   vn2 = vn2 / norm2(vn2)

   bmat = 0.0_wp
   bmat(2, 1:3) = vn / l1
   bmat(2, 7:9) = vn / l2
   bmat(2, 4:6) = -bmat(2, 1:3) - bmat(2, 7:9)
   bmat(1, 1:3) = vn2 / l1
   bmat(1, 7:9) = vn2 / l2
   bmat(1, 4:6) = -bmat(1, 1:3) - bmat(1, 7:9)
end function bmat_linbend

!> Wilson B matrix for torsions (4 atoms)
pure subroutine bmat_torsion(xyz, tau, bt)
   real(wp), intent(in) :: xyz(3, 4)
   real(wp), intent(out) :: tau
   real(wp), intent(out) :: bt(3, 4)

   real(wp) :: rij1, rjk1, rkl1
   real(wp) :: brij(3, 2), brjk(3, 2), brkl(3, 2)
   real(wp) :: bf2(3, 3), fi2, sinfi2, cosfi2
   real(wp) :: bf3(3, 3), fi3, sinfi3, cosfi3
   real(wp) :: costau, sintau
   integer :: ix, iy, iz

    call bmat_strtch(xyz(:, 1:2), rij1, brij)
    call bmat_strtch(xyz(:, 2:3), rjk1, brjk)
    call bmat_strtch(xyz(:, 3:4), rkl1, brkl)
    call bmat_bend(xyz(:, 1:3), fi2, bf2)
    sinfi2 = sin(fi2)
    cosfi2 = cos(fi2)
    call bmat_bend(xyz(:, 2:4), fi3, bf3)
   sinfi3 = sin(fi3)
   cosfi3 = cos(fi3)

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
   tau = atan2(sintau, costau)
   if (abs(tau) == pi) tau = pi

   do ix = 1, 3
      iy = ix + 1
      if (iy > 3) iy = iy - 3
      iz = iy + 1
      if (iz > 3) iz = iz - 3
      bt(ix, 1) = (brij(iy, 2)*brjk(iz, 2) - brij(iz, 2)*brjk(iy, 2)) &
         & / (rij1 * sinfi2**2)
      bt(ix, 4) = (brkl(iy, 1)*brjk(iz, 1) - brkl(iz, 1)*brjk(iy, 1)) &
         & / (rkl1 * sinfi3**2)
      bt(ix, 2) = -((rjk1 - rij1*cosfi2) * bt(ix, 1) &
         & + rkl1*cosfi3 * bt(ix, 4)) / rjk1
      bt(ix, 3) = -(bt(ix, 1) + bt(ix, 2) + bt(ix, 4))
   end do
end subroutine bmat_torsion

!> Wilson B matrix for out-of-plane angles (4 atoms, atom 4 central)
pure subroutine bmat_outofplane(xyz, teta, bt)
   real(wp), intent(in) :: xyz(3, 4)
   real(wp), intent(out) :: teta
   real(wp), intent(out) :: bt(3, 4)

   real(wp) :: r1(3), r2(3), r3(3)
   real(wp) :: q41, q42, q43, e41(3), e42(3), e43(3)
   real(wp) :: cosfi1, fi1, cosfi2, fi2, cosfi3, fi3
   real(wp) :: c14(3, 3), br14(3, 3)
   real(wp) :: r42(3), r43(3)
   integer :: ix, iy, iz

!  4 -> 1 (bond out of plane)
   r1 = xyz(:, 1) - xyz(:, 4)
   q41 = norm2(r1)
   e41 = r1 / q41
!  4 -> 2 (bond in plane)
   r2 = xyz(:, 2) - xyz(:, 4)
   q42 = norm2(r2)
   e42 = r2 / q42
!  4 -> 3 (bond in plane)
   r3 = xyz(:, 3) - xyz(:, 4)
   q43 = norm2(r3)
   e43 = r3 / q43

!  angle between e43 and e42
   cosfi1 = dot_product(e43, e42)
   fi1 = acos(cosfi1)

!  dirty exit when 2-3-4 collinear
   if (abs(fi1 - pi) < 1.0e-13_wp) then
      teta = 0.0_wp
      bt = 0.0_wp
      return
   end if

!  angle between e41 and e43
   cosfi2 = dot_product(e41, e43)
   fi2 = acos(cosfi2)

!  angle between e41 and e42
   cosfi3 = dot_product(e41, e42)
   fi3 = acos(cosfi3)

!  first two centers
   c14(:, 1) = xyz(:, 1)
   c14(:, 2) = xyz(:, 4)

!  3rd center: cross product of r42 x r43
   r42 = xyz(:, 2) - xyz(:, 4)
   r43 = xyz(:, 3) - xyz(:, 4)
   c14(1, 3) = r42(2)*r43(3) - r42(3)*r43(2)
   c14(2, 3) = r42(3)*r43(1) - r42(1)*r43(3)
   c14(3, 3) = r42(1)*r43(2) - r42(2)*r43(1)

!  exit if 2-3-4 collinear
   if (sum(c14(:, 3)**2) < 1.0e-10_wp) then
      teta = 0.0_wp
      bt = 0.0_wp
      return
   end if
   c14(:, 3) = c14(:, 3) + xyz(:, 4)

   call bmat_bend(c14, teta, br14)
   teta = teta - 0.5_wp * pi

!  compute the WDC matrix
   do ix = 1, 3
      iy = mod(ix + 1, 4) + (ix + 1) / 4
      iz = mod(iy + 1, 4) + (iy + 1) / 4

      bt(ix, 1) = -br14(ix, 1)
      bt(ix, 2) = r43(iz)*br14(iy, 3) - r43(iy)*br14(iz, 3)
      bt(ix, 3) = -r42(iz)*br14(iy, 3) + r42(iy)*br14(iz, 3)
      bt(ix, 4) = -(bt(ix, 1) + bt(ix, 2) + bt(ix, 3))
   end do

   bt = -bt
end subroutine bmat_outofplane

!> Stretch B row for 2-atom fragment (private helper)
pure subroutine bmat_strtch(xyz, avst, b)
   real(wp), intent(in) :: xyz(3, 2)
   real(wp), intent(out) :: avst
   real(wp), intent(out) :: b(3, 2)

   real(wp) :: r(3), rr

   r = xyz(:, 2) - xyz(:, 1)
   rr = norm2(r)
   avst = rr
   b(:, 1) = -r / rr
   b(:, 2) = -b(:, 1)
end subroutine bmat_strtch

!> Bend B row for 3-atom fragment (private helper)
pure subroutine bmat_bend(xyz, fir, bf)
   real(wp), intent(in) :: xyz(3, 3)
   real(wp), intent(out) :: fir
   real(wp), intent(out) :: bf(3, 3)

   real(wp) :: brij(3, 2), brjk(3, 2)
   real(wp) :: co, crap, si
   real(wp) :: rij1, rjk1
   integer :: i

    call bmat_strtch(xyz(:, 1:2), rij1, brij)
    call bmat_strtch(xyz(:, 2:3), rjk1, brjk)
   co = 0.0_wp
   crap = 0.0_wp
   do i = 1, 3
      co = co + brij(i, 1)*brjk(i, 2)
      crap = crap + (brjk(i, 2) + brij(i, 1))**2
   end do
   if (sqrt(crap) < 1.0e-6_wp) then
      fir = pi - asin(sqrt(crap))
      si = sqrt(crap)
   else
      fir = acos(co)
      si = sqrt(1.0_wp - co**2)
   end if
   if (abs(fir - pi) < 1.0e-13_wp) then
      fir = pi
      return
   end if
   do i = 1, 3
      bf(i, 1) = (co*brij(i, 1) - brjk(i, 2)) / (si*rij1)
      bf(i, 3) = (co*brjk(i, 2) - brij(i, 1)) / (si*rjk1)
      bf(i, 2) = -(bf(i, 1) + bf(i, 3))
   end do
end subroutine bmat_bend

end module xtb_bmatrix