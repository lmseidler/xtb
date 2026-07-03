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
   use xtb_mctc_math, only : crossProd
   implicit none

   private
   public :: bmat_bond, bmat_angle, bmat_linbend, &
      & bmat_torsion, bmat_outofplane

   character(len=*), parameter :: source = 'xtb_bmatrix'

   !> Tolerance for near-linear/collinear detection (radians)
   real(wp), parameter :: tol_collinear = 1.0e-13_wp
   !> Tolerance for near-zero vector norms
   real(wp), parameter :: tol_norm = 1.0e-15_wp
   !> Tolerance for cross-product magnitude (collinearity guard)
   real(wp), parameter :: tol_cross = 1.0e-10_wp
   !> Tolerance for near-linear bend detection
   real(wp), parameter :: tol_bend = 1.0e-6_wp
   !> Smallest positive argument avoiding division by zero
   real(wp), parameter :: eps_sqrt = 1.0e-15_wp

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
pure function bmat_angle(vec1, vec2) result(bmat)
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
      dinprod(ii) = dot_product(dnvec(1, :, ii), nvec2)
      dinprod(ii+3) = dot_product(dnvec(1, :, ii+3), nvec2) &
         & + dot_product(dnvec(2, :, ii+3), nvec1)
      dinprod(ii+6) = dot_product(dnvec(2, :, ii), nvec1)
   end do

   dot_n1n2 = dot_product(nvec1, nvec2)
   bmat = -dinprod / sqrt(max(eps_sqrt, 1.0_wp-dot_n1n2**2))
end function bmat_angle

!> Wilson B matrix for linear angles
pure function bmat_linbend(vec1, vec2) result(bmat)
   real(wp), intent(in) :: vec1(3), vec2(3)

   real(wp) :: bmat(2, 9)
   real(wp) :: l1, l2, nvec1(3), nvec2(3)
   real(wp) :: vn(3), vn2(3), nvn
   real(wp), parameter :: xaxis(3) = [1.0_wp, 0.0_wp, 0.0_wp], yaxis(3) = [0.0_wp, 1.0_wp, 0.0_wp]

   l1 = norm2(vec1)
   l2 = norm2(vec2)
   nvec1 = vec1 / l1
   nvec2 = vec2 / l2

   vn = crossProd(vec1, vec2)
   nvn = norm2(vn)

   if (nvn < tol_norm) then
      vn = xaxis - dot_product(xaxis, vec1) / l1**2 * vec1
      nvn = norm2(vn)
      if (nvn < tol_norm) then
         vn = yaxis - dot_product(yaxis, vec1) / l1**2 * vec1
         nvn = norm2(vn)
      end if
   end if
   vn = vn / nvn

   vn2 = crossProd(vec1-vec2, vn)
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
pure function bmat_torsion(xyz) result(bt)
   real(wp), intent(in) :: xyz(3, 4)

   real(wp) :: bt(3, 4)
   real(wp) :: rij1, rjk1, rkl1
   real(wp) :: brij(3, 2), brjk(3, 2), brkl(3, 2)
   real(wp) :: fi2, sinfi2, cosfi2
   real(wp) :: fi3, sinfi3, cosfi3
   integer :: ix, iy, iz

   rij1 = norm2(xyz(:, 2)-xyz(:, 1))
   rjk1 = norm2(xyz(:, 3)-xyz(:, 2))
   rkl1 = norm2(xyz(:, 4)-xyz(:, 3))
   brij = bmat_strtch(xyz(:, 1:2))
   brjk = bmat_strtch(xyz(:, 2:3))
   brkl = bmat_strtch(xyz(:, 3:4))

   fi2 = bend_angle(xyz(:, 1:3))
   sinfi2 = sin(fi2)
   cosfi2 = cos(fi2)
   fi3 = bend_angle(xyz(:, 2:4))
   sinfi3 = sin(fi3)
   cosfi3 = cos(fi3)

   do ix = 1, 3
      iy = ix + 1
      if (iy > 3) iy = iy - 3
      iz = iy + 1
      if (iz > 3) iz = iz - 3
      bt(ix, 1) = (brij(iy, 2)*brjk(iz, 2)-brij(iz, 2)*brjk(iy, 2)) &
         & / (rij1*sinfi2**2)
      bt(ix, 4) = (brkl(iy, 1)*brjk(iz, 1)-brkl(iz, 1)*brjk(iy, 1)) &
         & / (rkl1*sinfi3**2)
      bt(ix, 2) = -((rjk1-rij1*cosfi2)*bt(ix, 1) &
         & + rkl1 * cosfi3 * bt(ix, 4)) / rjk1
      bt(ix, 3) = -(bt(ix, 1)+bt(ix, 2)+bt(ix, 4))
   end do
end function bmat_torsion

!> Wilson B matrix for out-of-plane angles (4 atoms, atom 4 central)
pure function bmat_outofplane(xyz) result(bt)
   real(wp), intent(in) :: xyz(3, 4)

   real(wp) :: bt(3, 4)
   real(wp) :: r2(3), r3(3)
   real(wp) :: q42, q43, e42(3), e43(3)
   real(wp) :: cosfi1, fi1
   real(wp) :: c14(3, 3), br14(3, 3)
   real(wp) :: r42(3), r43(3)
   integer :: ix, iy, iz

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
   if (abs(fi1-pi) < tol_collinear) then
      bt = 0.0_wp
      return
   end if

!  first two centers
   c14(:, 1) = xyz(:, 1)
   c14(:, 2) = xyz(:, 4)

!  3rd center: cross product of r42 x r43
   r42 = xyz(:, 2) - xyz(:, 4)
   r43 = xyz(:, 3) - xyz(:, 4)
   c14(:, 3) = crossProd(r42, r43)

!  exit if 2-3-4 collinear
   if (sum(c14(:, 3)**2) < tol_cross) then
      bt = 0.0_wp
      return
   end if
   c14(:, 3) = c14(:, 3) + xyz(:, 4)

   br14 = bmat_bend(c14)

!  compute the WDC matrix
   do ix = 1, 3
      iy = mod(ix+1, 4) + (ix+1) / 4
      iz = mod(iy+1, 4) + (iy+1) / 4

      bt(ix, 1) = -br14(ix, 1)
      bt(ix, 2) = r43(iz) * br14(iy, 3) - r43(iy) * br14(iz, 3)
      bt(ix, 3) = -r42(iz) * br14(iy, 3) + r42(iy) * br14(iz, 3)
      bt(ix, 4) = -(bt(ix, 1)+bt(ix, 2)+bt(ix, 3))
   end do

   bt = -bt
end function bmat_outofplane

!> Stretch B row for 2-atom fragment (private helper)
pure function bmat_strtch(xyz) result(b)
   real(wp), intent(in) :: xyz(3, 2)

   real(wp) :: b(3, 2)
   real(wp) :: r(3), rr

   r = xyz(:, 2) - xyz(:, 1)
   rr = norm2(r)
   b(:, 1) = -r / rr
   b(:, 2) = -b(:, 1)
end function bmat_strtch

!> Bend angle for 3-atom fragment (private helper)
pure function bend_angle(xyz) result(fir)
   real(wp), intent(in) :: xyz(3, 3)

   real(wp) :: fir
   real(wp) :: brij(3, 2), brjk(3, 2)
   real(wp) :: co, crap
   integer :: i

   brij = bmat_strtch(xyz(:, 1:2))
   brjk = bmat_strtch(xyz(:, 2:3))
   co = 0.0_wp
   crap = 0.0_wp
   do i = 1, 3
      co = co + brij(i, 1) * brjk(i, 2)
      crap = crap + (brjk(i, 2)+brij(i, 1))**2
   end do
   if (sqrt(crap) < tol_bend) then
      fir = pi - asin(sqrt(crap))
   else
      fir = acos(co)
   end if
   if (abs(fir-pi) < tol_collinear) then
      fir = pi
   end if
end function bend_angle

!> Bend B row for 3-atom fragment (private helper)
pure function bmat_bend(xyz) result(bf)
   real(wp), intent(in) :: xyz(3, 3)

   real(wp) :: bf(3, 3)
   real(wp) :: brij(3, 2), brjk(3, 2)
   real(wp) :: co, crap, fir, si
   real(wp) :: rij1, rjk1
   integer :: i

   brij = bmat_strtch(xyz(:, 1:2))
   brjk = bmat_strtch(xyz(:, 2:3))
   rij1 = norm2(xyz(:, 2)-xyz(:, 1))
   rjk1 = norm2(xyz(:, 3)-xyz(:, 2))
   co = 0.0_wp
   crap = 0.0_wp
   do i = 1, 3
      co = co + brij(i, 1) * brjk(i, 2)
      crap = crap + (brjk(i, 2)+brij(i, 1))**2
   end do
   if (sqrt(crap) < tol_bend) then
      fir = pi - asin(sqrt(crap))
      si = sqrt(crap)
   else
      fir = acos(co)
      si = sqrt(1.0_wp-co**2)
   end if
   bf = 0.0_wp
   if (abs(fir-pi) < tol_collinear) then
      return
   end if
   do i = 1, 3
      bf(i, 1) = (co*brij(i, 1)-brjk(i, 2)) / (si*rij1)
      bf(i, 3) = (co*brjk(i, 2)-brij(i, 1)) / (si*rjk1)
      bf(i, 2) = -(bf(i, 1)+bf(i, 3))
   end do
end function bmat_bend

end module xtb_bmatrix
